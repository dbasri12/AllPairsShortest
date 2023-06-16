#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=10:00
#SBATCH --partition=slow

srun ./SSSP_MPI --sourceNode 3466 --strategy 2 --inputFile "input_graphs/CA-GrQc"
srun ./SSSP_MPI --sourceNode 3466 --strategy 1 --inputFile "input_graphs/CA-GrQc"
