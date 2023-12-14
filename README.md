# NMP
Algorithm to solve the Nestedness Maximization Problem



### HOW TO COMPILE nmp.c ###

open terminal and type: gcc -o nmp nmp.c -lm -O3

------------------------------------------------------------------------------------

### HOW TO EXECUTE nmp ###

open a terminal and type: ./nmp <FILE CONTAINING THE BIPARTITE NETWORK>

where <FILE CONTAINING THE BIPARTITE NETWORK> must be formatted as in the following example:

0,0,0,14,1
0,0,0,0,2
0,1,0,0,0

Each row must contain exactly the same number of entries (5 in the example above). 
The commas "," are important. 
Entries must be non-negative. 

EXAMPLE: ./nmp sampleNetwork.csv


***NOTE: for better results do the following:

1) open the source code "nmp.c"  
2) go to the main() and modify the variable beta1. The default value is: beta1 = 20.0. Choose bigger value for better result,  although convergence may slow down.
3) save and exit nmp.c and recompile

------------------------------------------------------------------------------------

### OUTPUT OF THE PROGRAM ###
1) The program prints at screen the size of the matrix and the ground state energy

2) The program creates two files:

2a) One with the results of the algorithm, named "result_sampleNetwork.csv.txt" 
2b) One with the nested adjacency matrix, named "Anest_sampleNetwork.csv.txt" in the format of an edge list:    

Node_i Node_j 1
Node_i Node_k 1
...

------------------------------------------------------------------------------------

### HOW PLOT THE NESTED ADJACENCY MATRIX ###

open gnuplot and type: plot 'Anest_sampleNetwork.csv.txt' u 1:2:3 w p pt 5 palette



