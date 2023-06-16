#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=5G
#SBATCH --time=07:00
#SBATCH --partition=slow

make DEBUG=1

srun ./SSSP_MPI --sourceNode 0 --strategy 1 --inputFile "input_graphs/GG_Test"
echo " "
srun ./SSSP_MPI --sourceNode 0 --strategy 2 --inputFile "input_graphs/GG_Test"
echo " "
srun ./SSSP_MPI --sourceNode 2 --strategy 1 --inputFile "input_graphs/GG_Test"
echo " "
srun ./SSSP_MPI --sourceNode 2 --strategy 2 --inputFile "input_graphs/GG_Test"

make clean