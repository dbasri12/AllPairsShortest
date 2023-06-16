#!/bin/bash
#
#SBATCH --cpus-per-task=4
#SBATCH --time=10:00
#SBATCH --mem=1G
#SBATCH --partition=slow

make DEBUG=1

echo " "
srun ./SSSP --sourceNode 0 --strategy 1 --inputFile "input_graphs/GG_Test"
echo " "
srun ./SSSP --sourceNode 0 --strategy 2 --inputFile "input_graphs/GG_Test"
echo " "
srun ./SSSP_parallel --nThreads 4 --sourceNode 0  --strategy 1 --inputFile "input_graphs/GG_Test"
echo " "
srun ./SSSP_parallel --nThreads 4 --sourceNode 0  --strategy 2 --inputFile "input_graphs/GG_Test"
echo " "
srun ./SSSP_parallel --nThreads 4 --sourceNode 0  --strategy 3 --granularity 3 --inputFile "input_graphs/GG_Test"
echo " "
srun ./SSSP --sourceNode 2 --strategy 1 --inputFile "input_graphs/GG_Test"
echo " "
srun ./SSSP --sourceNode 2 --strategy 2 --inputFile "input_graphs/GG_Test"
echo " "
srun ./SSSP_parallel --nThreads 4 --sourceNode 2  --strategy 1 --inputFile "input_graphs/GG_Test"
echo " "
srun ./SSSP_parallel --nThreads 4 --sourceNode 2  --strategy 2 --inputFile "input_graphs/GG_Test"
echo " "
srun ./SSSP_parallel --nThreads 4 --sourceNode 2  --strategy 3 --granularity 3 --inputFile "input_graphs/GG_Test"

make clean