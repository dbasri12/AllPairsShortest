# Distributed-Linear-Equation-Solver
CMPT 431 Project: Parallel Algorithm Implementation Using C++ Threads and MPI

4. All-Pairs Shortest Path: Described in the textbook (DC Chapter 5).

## Dependencies

All code can be run on Head Node (cs-cloud.cs.surrey.sfu.ca) except

```
test_scripts/plot_results.py
```
which needs seaborn to run

### Local Dependencies
1. Distributed Algorithms:  
      [MPI][1]  
2. Result plotting:  
   [seaborn and its dependencies][2]  
   ```
   sudo apt-get install python3 python3-dev python3-pip
   pip3 install seaborn
   ```

## Input File Generation
Already supplied in [input_graphs][3]

## Steps to Run Algorithms (serial, parallel, distributed)

### Slurm with Main Test Scripts

Starting in Main Directory on Head Node
1. Running Serial and Parallel SSSP
   ```
   sbatch correctness_test.sh
   ```
2. Running Distributed SSSP
   ```
   sbatch MPI_test.sh
   ```

### Slurm Using Script Generation

Starting in main directory
```
make
cd test_scripts/
python3 cluster_script_generator.py
cd ../
```

### Local
Starting in Main Directory on Local Machine
1. Running Serial and Parallel SSSP
   ```
   sbatch correctness_test.sh
   ```
2. Running Distributed SSSP
   ```
   sbatch MPI_test.sh
   ```

[1]: https://mpitutorial.com/tutorials/installing-mpich2/ "MPICH2 Local Installation Guide"
[2]: https://seaborn.pydata.org/installing.html "seaborn Installation Guide"
[3]: ../../tree/main/input_graphs