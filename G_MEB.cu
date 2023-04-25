#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unordered_set>
#include <string>
#include <fstream>
#include <iostream>

/** Version using one block per vertex in V, 2d arrays, and shared memory**/

#define PAGELOCKED
#define WRITE_TO_FILE

// #define MARVEL
#define ARXIV
// #define HUMAN
// #define MOVIE


#ifdef MOVIE
#define N_v 943
#define N_u 1682
#define N_E 99906
#endif

#ifdef HUMAN
#define N_v 4025
#define N_u 96
#define N_E 40197
#endif


/* Marvel Dataset*/
#ifdef MARVEL
#define N_v 6486
#define N_u 12942
#define N_E 96662
#endif

#ifdef ARXIV
#define N_v 16726
#define N_u 22015
#define N_E 58595
#endif



// #define BLOCKS_PER_BATCH 64
#define BICLIQUE_ROW_LENGTH (((N_v+N_u)/32)+1)

#define ROW_OUTPUT_SIZE (N_E)*10
#define COL_OUTPUT_SIZE (N_E)*10

struct edge {
    int u;
    int v;
};

/** Set k-th bit **/
__host__ __device__ void set_bit(uint32_t *map, int k){
    map[k/32] = map[k/32] | (1 << (k % 32));
}


/** Set k-th bit to a value **/
__host__ __device__ void set_val(uint32_t *map, int k, int bit){
    map[k/32] = map[k/32] | (bit << (k % 32));
}

/** Clear k-th bit **/
__host__ __device__ void clear_bit(uint32_t *map, int k){
    map[k/32] = map[k/32] & ~(1 << (k%32));
}

/** Return k-th bit **/
__host__ __device__ uint32_t get_bit(uint32_t *map , int k){
    return (map[k/32] & ((uint32_t)1 << k%32)) >> k%32 ;
}

/** Generate the indices of where each vertex in set V begins its edge list
 *  relative to E_V (edge list for vertices in V)
 *  **/
__global__ void generate_EV_indices(int d_E_V[N_E], int d_v_indices[N_v]){
    int current_vertex = -1;
    for (int i = 0; i < N_E; i++)
    {
        if (d_E_V[i] > current_vertex)
        {
            current_vertex++;
            d_v_indices[current_vertex] = i;
        }
        
    }
    
}

/** generate_EV_indices analog for EU **/

__global__ void generate_EU_indices(int d_E_U[N_E], int d_u_indices[N_u]){
    int current_vertex = -1;
    for (int i = 0; i < N_E; i++)
    {
        if (d_E_U[i] > current_vertex)
        {
            current_vertex++;
            d_u_indices[current_vertex] = i;
        }
        
    }
    
}


/** Create a bitmap of V-U connections (N_v rows by N_u columns)
 *  **/
__global__ void make_bitmap(uint32_t *d_bitmap, int *d_E_V, int *d_E_U, int *d_v_indices, size_t bitmap_pitch, uint32_t *d_bitmap_T, size_t bitmap_T_pitch, int *freq_u, int *freq_v){

    for(int i = 0; i < N_E; i++ ){
        // printf("i: %d d_E_U: %d\n", i, d_E_U[i]);
        set_bit((uint32_t *)((char *)d_bitmap+d_E_V[i]*bitmap_pitch), d_E_U[i]);
        set_bit((uint32_t *)((char *)d_bitmap_T+d_E_U[i]*bitmap_T_pitch), d_E_V[i]);
        freq_u[d_E_U[i]]++;
        freq_v[d_E_V[i]]++;
    }    

    
}

/*** Row-wise traversal of bitmap to find maximum edge biclique ***/

__global__ void row_traversal(int vertices_start, uint32_t *d_bitmap, size_t bitmap_pitch,  int *d_E_U, int *d_v_indices, uint32_t *d_output, int stream_num, int *output_row, int *maximal_edge, int *freq_u){
    int vertex_V = blockIdx.x + vertices_start;
    int thread_num = threadIdx.x;
    int threadsPerBlock = blockDim.x;

    if (vertex_V >= N_v){return;}
    
    
    __shared__ uint32_t main_row[(N_u/32) + 1]; // Shared memory row for comparison
    __shared__ char check_rows[N_v];

    for(int i = thread_num; i < (N_u/32) + 1; i+=threadsPerBlock){ // Populate main row for comparison
        main_row[i] = ((uint32_t *)((char *)d_bitmap+vertex_V*bitmap_pitch))[i];
    }
        
    __syncthreads();  // Barrier to advance
    int start = d_v_indices[vertex_V];
    int end = N_E;
    if(vertex_V < N_v - 1){
        end = d_v_indices[vertex_V+1];
    }

    // Determine if node should be visited (only visit 2-neighbors)
    for (int i = thread_num; i < N_v; i+=threadsPerBlock)
    {
        int sum = 0;
        for (int j = start; j < end; j++)
        {
            if((get_bit(main_row, d_E_U[j]) & get_bit((uint32_t *)((char *)d_bitmap+i*bitmap_pitch), d_E_U[j])) != 0){
                sum++; 
            }
        }
        if (sum > 0){
            check_rows[i] = (char) 1;
        }else{
            check_rows[i] = (char) 0;
        }
    }

    __syncthreads();

    // Main processing loop
    int frequency;
    for (int i = thread_num; i < N_v; i+=threadsPerBlock)
    {
        if (check_rows[i] == (char) 0 || i < vertex_V){continue;}

        uint32_t resultant[BICLIQUE_ROW_LENGTH];  // Resultant bit vector for biclique to be generated
        frequency = 0;
        for (int j = 0; j < BICLIQUE_ROW_LENGTH; j++)
        {
            resultant[j] = 0;
        }
        

        set_bit(resultant, vertex_V);
        set_bit(resultant, i);


        
        // Set resultant bit for all shared neighbors between node j and block node
        int num_edges = end - start;
        for (int j = start; j < end; j++)
        {
            if (get_bit(main_row, d_E_U[j]) & get_bit((uint32_t *)((char *)d_bitmap+i*bitmap_pitch), d_E_U[j]))
            {
                set_bit(resultant, N_v+d_E_U[j]);
            }
            
        }
        

        int equals_prior = 0;
        // Find all 2-neighbors that have same shared 1-neighbors
        for (int j = 0; j < i; j++)
        {
            if (check_rows[j] == (char) 0){continue;}
            int num_same = 0;
            int set_j = 1;
            for (int k = start; k < end; k++)
            {
                int current_row = (get_bit(main_row, d_E_U[k]) & get_bit((uint32_t *)((char *)d_bitmap+i*bitmap_pitch), d_E_U[k]));
                int comparison = get_bit(main_row, d_E_U[k]) & get_bit((uint32_t *)((char *)d_bitmap+j*bitmap_pitch), d_E_U[k]);
                if ((current_row == comparison))
                {
                    num_same++;
                }
                if (current_row != comparison && current_row == 1)
                {
                    set_j = 0;
                }
            }
            if (num_same == num_edges)
            {
                equals_prior = 1;
                break;
            }
            set_val(resultant, j, set_j);
            frequency+=set_j;
        }
        if (equals_prior == 1)
        {
            continue;
        }
        

        for (int j = i; j < N_v; j++)
        {
            if (check_rows[j] == (char) 0){continue;}
            int set_j = 1;
            for (int k = start; k < end; k++)
            {
                int current_row = (get_bit(main_row, d_E_U[k]) & get_bit((uint32_t *)((char *)d_bitmap+i*bitmap_pitch), d_E_U[k]));
                int comparison = get_bit(main_row, d_E_U[k]) & get_bit((uint32_t *)((char *)d_bitmap+j*bitmap_pitch), d_E_U[k]);
                if (current_row != comparison && current_row == 1)
                {
                    set_j = 0;
                }
                
            }
            set_val(resultant, j, set_j);
            frequency+=set_j;
        }
        
        if (!equals_prior)
        {

            
            int insert_row = atomicAdd(&output_row[stream_num], 1);
            int num_v = 0;
            int num_u = 0;
            for (int j = 0; j < N_v+N_u; j++)
            {
                int value = get_bit(resultant,j);
                set_val(&d_output[insert_row*BICLIQUE_ROW_LENGTH], j, value);
                if (j < N_v && value == 1)
                {
                    num_v++;
                }else if(j>=N_v && value == 1){
                    num_u++;
                    if (frequency == freq_u[j-N_v])
                    {
                        atomicCAS(&freq_u[j-N_v], frequency, 0);
                    }
                    
                }
                
            }
            int old = atomicMax(&maximal_edge[2], num_u*num_v);
            // Maximum edge biclique check
            if (old < num_u*num_v)
            {
                atomicExch(&maximal_edge[0], num_v);
                atomicExch(&maximal_edge[1], num_u);
                atomicExch(&maximal_edge[3], insert_row);
            }   
        }   
    }
}

/*** Column-wise traversal. Analog of row_traversal. Check row_traversal comments for details. ***/
__global__ void col_traversal(int vertices_start, uint32_t *d_bitmap, size_t bitmap_pitch,  int *d_E_U, int *d_v_indices, uint32_t *d_output, int stream_num, int *output_row, int *maximal_edge){
    int vertex_V = blockIdx.x + vertices_start;
    int thread_num = threadIdx.x;
    int threadsPerBlock = blockDim.x;

    if (vertex_V >= N_u){return;}
    
    
    __shared__ uint32_t main_row[(N_v/32) + 1]; // Shared memory row for comparison
    __shared__ char check_rows[N_u];

    for(int i = thread_num; i < (N_v/32) + 1; i+=threadsPerBlock){ // Populate main row for comparison
        main_row[i] = ((uint32_t *)((char *)d_bitmap+vertex_V*bitmap_pitch))[i];
    }
        
    __syncthreads();  // Barrier to advance
    int start = d_v_indices[vertex_V];
    int end = N_E;
    if(vertex_V < N_u - 1){
        end = d_v_indices[vertex_V+1];
    }
    for (int i = thread_num; i < N_u; i+=threadsPerBlock)
    {
        int sum = 0;
        for (int j = start; j < end; j++)
        {
            if((get_bit(main_row, d_E_U[j]) & get_bit((uint32_t *)((char *)d_bitmap+i*bitmap_pitch), d_E_U[j])) != 0){
                sum++; 
            }
        }
        if (sum > 0){
            check_rows[i] = (char) 1;
        }else{
            check_rows[i] = (char) 0;
        }        
    }

    __syncthreads();
    
    for (int i = thread_num; i < N_u; i+=threadsPerBlock)
    {
        if (check_rows[i] == (char) 0 || i < vertex_V){continue;}

        uint32_t resultant[BICLIQUE_ROW_LENGTH];
        
        for (int j = 0; j < BICLIQUE_ROW_LENGTH; j++)
        {
            resultant[j] = 0;
        }
        

        set_bit(resultant, N_v+vertex_V);
        set_bit(resultant, N_v+i);

        int num_edges = end - start;
        int overlap_count = 0;
        for (int j = start; j < end; j++)
        {
            if (get_bit(main_row, d_E_U[j]) & get_bit((uint32_t *)((char *)d_bitmap+i*bitmap_pitch), d_E_U[j]))
            {
                set_bit(resultant, d_E_U[j]);
                overlap_count++;
            }
            
        }

        if (overlap_count <= 2)
        {
            continue;
        }
        
        

        int equals_prior = 0;
        for (int j = 0; j < i; j++)
        {
            if (check_rows[j] == (char) 0){continue;}
            int num_same = 0;
            int set_j = 1;
            for (int k = start; k < end; k++)
            {
                int current_row = (get_bit(main_row, d_E_U[k]) & get_bit((uint32_t *)((char *)d_bitmap+i*bitmap_pitch), d_E_U[k]));
                int comparison = get_bit(main_row, d_E_U[k]) & get_bit((uint32_t *)((char *)d_bitmap+j*bitmap_pitch), d_E_U[k]);
                if ((current_row == comparison))
                {
                    num_same++;
                }
                if (current_row != comparison && current_row == 1)
                {
                    set_j = 0;
                }
            }
            if (num_same == num_edges)
            {
                equals_prior = 1;
                break;
            }
            set_val(resultant, N_v+j, set_j);
            
        }

        if (equals_prior == 1)
        {
            continue;
        }
        

        for (int j = i; j < N_u; j++)
        {
            if (check_rows[j] == (char) 0){continue;}
            int set_j = 1;
            for (int k = start; k < end; k++)
            {
                int current_row = (get_bit(main_row, d_E_U[k]) & get_bit((uint32_t *)((char *)d_bitmap+i*bitmap_pitch), d_E_U[k]));
                int comparison = get_bit(main_row, d_E_U[k]) & get_bit((uint32_t *)((char *)d_bitmap+j*bitmap_pitch), d_E_U[k]);
                if (current_row != comparison && current_row == 1)
                {
                    set_j = 0;
                }
                
            }
            set_val(resultant, N_v+j, set_j);
        }
        
        if (!equals_prior)
        {
            
            
            int insert_row = atomicAdd(&output_row[stream_num], 1);
            
            int num_v = 0;
            int num_u = 0;
            for (int j = 0; j < N_v+N_u; j++)
            {
                int value = get_bit(resultant,j);
                set_val(&d_output[insert_row*BICLIQUE_ROW_LENGTH], j, value);
                if (j < N_v && value == 1)
                {
                    num_v++;
                }else if(j>=N_v && value == 1){
                    num_u++;
                }
                
            }
            int old = atomicMax(&maximal_edge[2], num_u*num_v);
            if (old < num_u*num_v)
            {
                atomicExch(&maximal_edge[0], num_v);
                atomicExch(&maximal_edge[1], num_u);
                atomicExch(&maximal_edge[3], insert_row);
            } 
        }   
    }
}


int sort_comp(const void * a, const void * b){
    int x = *((int*) a);
    int y = *((int*) b);
    return x - y ;
}

int edge_comp_u(const void * a, const void * b){
    const struct edge *E_A = (edge *)a;
    const struct edge *E_B = (edge *)b;
    
    
    if (E_A->u - E_B->u != 0)
    {
        return E_A->u - E_B->u;
    }else{
        return E_A->v - E_B->v;
    }
}

int edge_comp_v(const void * a, const void * b){
    const struct edge *E_A = (edge *)a;
    const struct edge *E_B = (edge *)b;
    
    if (E_A->v - E_B->v != 0)
    {
        return E_A->v - E_B->v;
    }else{
        return E_A->u - E_B->u;
    }
}

int main(int argc, char**  argv){

    int i, j;

    /** CUDA Configurations **/

    if (argc != 3 ){    
        printf("Please enter arguments as follows: ./MBE_fci_reduce <nThreads> <nStreams>");
    }
    
    
    int threads_per_block = atoi(argv[1]);

    int threads_density = threads_per_block;
    // printf("Threads from density: %d\n", threads_density);
    dim3 gridDim_v(N_v, 1, 1);
    dim3 gridDim_u(N_u, 1, 1);
    dim3 blockDim(threads_per_block, 1, 1);
    dim3 blockDim_comp(threads_density, 1, 1);

    int nStreams = atoi(argv[2]);
    // int nStreams = STREAMS;

    /** Host/Device variables **/

    // Edge lists
    #ifdef DEBUG
    int temp_E_V[N_E] = {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 4, 4};
    int temp_E_U[N_E] = {0, 1, 5, 0, 1, 3, 0, 1, 4, 2, 4, 2, 3, 4, 5, 6};
    #endif

    #ifdef PAGELOCKED
    int *E_V;
    int *E_U;
    #endif

    #ifndef PAGELOCKED
    int E_V[N_E];
    int E_U[N_E];
    #endif

    int *d_E_V;
    int *d_E_U;

    int *E_V_u;
    int *E_U_u;
    int *dE_V_u;
    int *dE_U_u;

    struct edge edges_by_u[N_E];
    struct edge edges_by_v[N_E];

    int *freq_u;
    int *freq_v;

    // Bitmaps
    uint32_t *d_bitmap;
    uint32_t *d_bitmap_T;

    size_t bitmap_pitch;
    size_t bitmap_T_pitch;


    // Output arrays
    uint32_t *h_row_output;
    uint32_t *h_col_output;
    uint32_t *d_row_output;
    uint32_t *d_col_output;

    // CUDA Streams for compare and find
    cudaStream_t *stream = (cudaStream_t*)malloc(2*nStreams*sizeof(cudaStream_t));

    // E_V index arrays for faster accesses
    int *d_v_indices;

    int *d_u_indices;

    int *output_row;
    int *output_col;

    int *maximal_edge;

    /** Input Dataset **/

    #ifdef MARVEL
    char in_file_string[]= "formatted_marvel.txt";
    #endif
    #ifdef ARXIV
    char in_file_string[]= "formatted_arxiv.txt";
    #endif
    #ifdef HUMAN
    char in_file_string[] = "formatted_human.txt";
    #endif
    #ifdef MOVIE
    char in_file_string[] = "formatted_movie.txt";
    #endif


    /** Allocate host variables **/

    // Using page_locked
    #ifdef PAGELOCKED
    cudaMallocHost(&E_V, sizeof(int)*N_E);
    cudaMallocHost(&E_U, sizeof(int)*N_E);
    #endif


    /** Read input file **/
    #ifdef DEBUG
    for ( i = 0; i < N_E; i++)
    {
        E_V[i] = temp_E_V[i];
        E_U[i] = temp_E_U[i];
    }
    #endif

    

    #ifndef DEBUG
    printf("Reading input file\n");

    FILE *in_file = fopen(in_file_string, "r");
    int fill = 0;
    while (fscanf(in_file, "%d,%d", &E_V[fill], &E_U[fill])!=EOF)
    {
        fill++;
    }
    fclose(in_file);
    printf("Read input file\n");
    #endif

    for ( i = 0; i < N_E; i++)
    {
        
        edges_by_u[i].u = E_U[i];
        edges_by_u[i].v = E_V[i];
        edges_by_v[i].u = E_U[i];
        edges_by_v[i].v = E_V[i];
    }

    qsort(edges_by_u, N_E, sizeof(edge), edge_comp_u);
    qsort(edges_by_v, N_E, sizeof(edge), edge_comp_v);

    cudaMallocHost(&E_U_u, sizeof(int)*N_E);
    cudaMallocHost(&E_V_u, sizeof(int)*N_E);
    for ( i = 0; i < N_E; i++)
    {
        E_U_u[i] = edges_by_u[i].u;
        E_V_u[i] = edges_by_u[i].v;
        E_U[i] = edges_by_v[i].u;
        E_V[i] = edges_by_v[i].v;
        // printf("(%d,%d) ", E_U[i], E_V[i]);
    }
    // printf("\n");
    cudaEvent_t start,stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);

    int er = 0;
    /** Allocate Device Variables **/
    er+=cudaMallocPitch(&d_bitmap, &bitmap_pitch, ((N_u/32)+1)*sizeof(uint32_t), N_v);
    er+=cudaMallocPitch(&d_bitmap_T, &bitmap_T_pitch, ((N_v/32)+1)*sizeof(uint32_t), N_u);

    er+=cudaMalloc(&d_E_V, sizeof(int)*N_E);
    er+=cudaMalloc(&d_E_U, sizeof(int)*N_E);
    er+=cudaMalloc(&d_v_indices, sizeof(int)*N_v);

    er+=cudaMalloc(&d_u_indices, sizeof(int)*N_u);
    er+=cudaMalloc(&dE_V_u, sizeof(int)*N_E);
    er+=cudaMalloc(&dE_U_u, sizeof(int)*N_E);
    
    er+=cudaMemcpy(d_E_V, E_V, sizeof(int) * N_E, cudaMemcpyHostToDevice);
    if (er == 1){printf("error with E_V memcpy!\n");}
    er+=cudaMemcpy(d_E_U, E_U, sizeof(int) * N_E, cudaMemcpyHostToDevice);
    if (er == 1){printf("error with E_U memcpy!\n");}
    er+=cudaMemset2D(d_bitmap_T, bitmap_T_pitch, 0, ((N_v/32)+1)*sizeof(uint32_t), N_u);
    er+=cudaMemset2D(d_bitmap, bitmap_pitch, 0, ((N_u/32)+1)*sizeof(uint32_t), N_v);
    er+=cudaMemset(d_v_indices, 0, sizeof(int)*N_v);
    er+=cudaMemcpy(dE_V_u, E_V_u, sizeof(int) * N_E, cudaMemcpyHostToDevice);
    er+=cudaMemcpy(dE_U_u, E_U_u, sizeof(int) * N_E, cudaMemcpyHostToDevice);

    er+=cudaMallocManaged(&output_row, sizeof(int)*nStreams);
    er+=cudaMallocManaged(&output_col, sizeof(int)*nStreams);
    er+=cudaMallocManaged(&maximal_edge, sizeof(int)*4);

    er+=cudaMallocManaged(&freq_u, sizeof(int)*N_u);
    er+=cudaMallocManaged(&freq_v, sizeof(int)*N_v);

    if (er > 0)                 // Error in device memory copy/setting
    {
        printf("Exiting due to memory copy/setting error\n");
        cudaFree(d_bitmap_T);
        cudaFreeHost(E_V);
        cudaFreeHost(E_U);
        cudaFree(d_bitmap);
        cudaFree(d_E_U);
        cudaFree(d_E_V);
        cudaFree(d_v_indices);
        exit(1);
    }

    for ( i = 0; i < nStreams*2; i++)
    {
        cudaStreamCreate(&stream[i]);
    }
    
    
    /** Generate indices for faster accesses **/

    generate_EV_indices<<<1,1>>>(d_E_V, d_v_indices);
    generate_EU_indices<<<1,1>>>(dE_U_u, d_u_indices);

    cudaDeviceSynchronize();

    /** Make bitmap **/

    make_bitmap<<<1, 1>>>(d_bitmap, d_E_V, d_E_U, d_v_indices, bitmap_pitch, d_bitmap_T, bitmap_T_pitch, freq_u, freq_v);


    if (er > 0)                 // Error in device memory copy/setting
    {
        printf("Exiting due to memory copy/setting error for rows/cols store\n");

        cudaFreeHost(E_V);
        cudaFreeHost(E_U);

        cudaFree(d_bitmap_T);
        cudaFree(d_bitmap);
        cudaFree(d_E_U);
        cudaFree(d_E_V);
        cudaFree(d_v_indices);
        exit(1);
    }


    #ifdef PAGELOCKED
    cudaFreeHost(E_V);
    cudaFreeHost(E_U);
    #endif
    cudaDeviceSynchronize();    

    er+=cudaMallocHost(&h_row_output, sizeof(uint32_t)*ROW_OUTPUT_SIZE*BICLIQUE_ROW_LENGTH);
    er+=cudaMallocHost(&h_col_output, sizeof(uint32_t)*COL_OUTPUT_SIZE*BICLIQUE_ROW_LENGTH);
    er+=cudaMalloc(&d_row_output, sizeof(uint32_t)*ROW_OUTPUT_SIZE*BICLIQUE_ROW_LENGTH);
    er+=cudaMalloc(&d_col_output, sizeof(uint32_t)*COL_OUTPUT_SIZE*BICLIQUE_ROW_LENGTH);

    er+=cudaMemset(d_row_output, 0, sizeof(uint32_t)*ROW_OUTPUT_SIZE*BICLIQUE_ROW_LENGTH);
    er+=cudaMemset(d_col_output, 0, sizeof(uint32_t)*COL_OUTPUT_SIZE*BICLIQUE_ROW_LENGTH);

    // printf("Starting row/col traversal\n");
    int vblocks_per_stream = (int)ceil((double)N_v/(double)nStreams);
    int ublocks_per_stream = (int)ceil((double)N_u/(double)nStreams);

    int row_size = (int)ceil((double)ROW_OUTPUT_SIZE/(double)nStreams);
    int col_size = (int)ceil((double)COL_OUTPUT_SIZE/(double)nStreams);

    for (int i = 0; i < nStreams; i++)
    {
        output_row[i] = i*row_size;
        output_col[i] = i*col_size;
    }

    for (int i = 0; i < 4; i++)
    {
        maximal_edge[i] = 0;
    }
    
    
    int *r_starts = (int *) malloc (sizeof(int)*nStreams);
    int *c_starts = (int *) malloc (sizeof(int)*nStreams);
    for (int i = 0; i < nStreams; i++)
    {
        int row_start = output_row[i];
        int col_start = output_col[i];
        r_starts[i] = row_start;
        c_starts[i] = col_start;

        int output_end = (i+1)*row_size;
        int col_end = (i+1)*col_size;
        if (output_end > ROW_OUTPUT_SIZE)
        {
            output_end = ROW_OUTPUT_SIZE;
        }
        if (col_end > COL_OUTPUT_SIZE)
        {
            col_end = COL_OUTPUT_SIZE;
        }
        
        
        int v_start = i *vblocks_per_stream;
        int v_end = (i+1)*vblocks_per_stream;
        if (v_end >= N_v)
        {
            v_end = N_v;
        }
        if (v_start < N_v )
        {
            dim3 rowGrid(v_end-v_start, 1, 1);
            row_traversal<<<rowGrid, blockDim_comp, 0, stream[i]>>>(v_start, d_bitmap, bitmap_pitch, d_E_U, d_v_indices, d_row_output, i, output_row, maximal_edge, freq_u);
            cudaMemcpyAsync(&h_row_output[output_row[i]*BICLIQUE_ROW_LENGTH], &d_row_output[output_row[i]*BICLIQUE_ROW_LENGTH], (output_end-output_row[i])*sizeof(uint32_t)*BICLIQUE_ROW_LENGTH, cudaMemcpyDeviceToHost, stream[i]);
        }

        int u_start = i *ublocks_per_stream;
        int u_end = (i+1)*ublocks_per_stream;
        if (u_end >= N_u)
        {
            u_end = N_u;
        }
        if (u_start < N_u)
        {
            dim3 colGrid(u_end-u_start, 1, 1);
            col_traversal<<<colGrid, blockDim_comp, 0, stream[2*i]>>>(u_start, d_bitmap_T, bitmap_T_pitch, dE_V_u, d_u_indices, d_col_output, i, output_col, maximal_edge);
            cudaMemcpyAsync(&h_col_output[output_col[i]*BICLIQUE_ROW_LENGTH], &d_col_output[output_col[i]*BICLIQUE_ROW_LENGTH], (col_end-output_col[i])*sizeof(uint32_t)*BICLIQUE_ROW_LENGTH, cudaMemcpyDeviceToHost, stream[2*i]);
        }
        
        
        
    }
    
    
    cudaDeviceSynchronize();


    cudaEventRecord(stop);

    cudaEventSynchronize(stop);
    float ms = 0;
    cudaEventElapsedTime(&ms, start, stop);

    printf("Finished GPU portion\n");
    printf("GPU Running Time: %.6f ms\n", ms);

    int amount_left = 0;
    for (int i = 0; i < N_u; i++)
    {
        if(freq_u[i]>0){amount_left++;}
    }
    // printf("Amount left: %d\n", amount_left);

    printf("===Maximal Edge Biclique===\n");
    printf("%d %d %d\n", maximal_edge[0], maximal_edge[1], maximal_edge[2]);

    for (int i = 0; i < 2*nStreams; ++i)
    {
        cudaStreamDestroy(stream[i]);
    }

    #ifdef WRITE_TO_FILE
    std::ofstream out_file("bi_output.txt");
    #endif
    // FILE *out_file;
    // out_file = fopen("bi_output.txt", "w");
    int num_bicliques = 0;
    using namespace std;
    unordered_set <string> biclique_hash;


    for ( i = 0; i < ROW_OUTPUT_SIZE; i++)
    {
        string bi_string ="";
        int num_v = 0;
        int num_u = 0;
        for ( j = 0; j < N_v+N_u; j++)
        {
            if (get_bit(&h_row_output[i*BICLIQUE_ROW_LENGTH], j) == 1)
            {
                
                // char prefix;
                if (j < N_v)
                {
                    num_v++;
                    bi_string+='v';
                    bi_string+=to_string(j);
                    // fprintf(out_file, "v%d ", j);
                } else{
                    num_u++;
                    bi_string+='u';
                    bi_string+=to_string(j-N_v);
                    // fprintf(out_file, "u%d ", j-N_v);
                }
                bi_string+=' ';
            }
            
        }
        if (bi_string =="")
        {
            continue;
        }
        
        if (biclique_hash.find(bi_string) == biclique_hash.end() && num_v > 0 && num_u > 0)
        {
            biclique_hash.insert(bi_string);
            num_bicliques++;
            #ifdef WRITE_TO_FILE
            out_file << bi_string;
            out_file << '\n';
            #endif
        }
        

        // fprintf(out_file,"\n");
    }
    int from_col = 0;
    for ( i = 0; i < COL_OUTPUT_SIZE; i++)
    {
        string bi_string = "";
        int num_v = 0;
        int num_u = 0;
        for ( j = 0; j < N_v+N_u; j++)
        {
            if (get_bit(&h_col_output[i*BICLIQUE_ROW_LENGTH], j) == 1)
            {
                // char prefix;
                if (j < N_v)
                {
                    num_v++;
                    bi_string+='v';
                    bi_string+=to_string(j);
                    // fprintf(out_file, "v%d ", j);
                } else{
                    num_u++;
                    bi_string+='u';
                    bi_string+=to_string(j-N_v);
                    // fprintf(out_file, "u%d ", j-N_v);
                }
                bi_string+=' ';
            }
            
        }
        if (biclique_hash.find(bi_string) == biclique_hash.end())
        {
            biclique_hash.insert(bi_string);
            num_bicliques++;
            from_col++;
            #ifdef WRITE_TO_FILE
            out_file << bi_string;
            out_file << '\n';
            #endif
        }
        

    //     // fprintf(out_file,"\n");
    }
    printf("Number of bicliques found: %d\n", num_bicliques);
    #ifdef WRITE_TO_FILE
    out_file.close();
    #endif


    cudaFree(d_bitmap);
    cudaFree(d_E_U);
    cudaFree(d_E_V);
    cudaFree(d_v_indices);

    cudaFree(d_bitmap_T);
}