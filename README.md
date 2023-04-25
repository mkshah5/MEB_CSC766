# MEB_CSC766
Maximum Edge Biclique problem on GPU

## Dependencies
- CUDA 11+
- nvcc

## Building
- Run `make` to build
- To test different datasets, define dataset name at top of G_MEB.cu (see file for example)

## Running
- Run `./G_MEB <nThreads> <nStreams>`
- nThreads specifies number of threads per CUDA block
- nStreams specifies number of CUDA streams
- If `WRITE_TO_FILE` is defined, all enumerated maximal bicliques will be output to the file `bi_output.txt`. 

## Expected Output
G_MEB will output the size of the maximum edge biclique as follows:
```
===Maximal Edge Biclique===

left_size right_size num_edges
```

`left_size` is the number of nodes in the left partition. `right_size` is the number of nodes in the right partition. `num_edges` is the number of edges in the maximum edge biclique.

