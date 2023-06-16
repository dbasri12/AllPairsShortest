#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00
#SBATCH --mem=4G
#SBATCH --partition=slow

srun ./SSSP --sourceNode 3466 --strategy 1 --inputFile "input_graphs/CA-GrQc"
srun ./SSSP --sourceNode 3466 --strategy 2 --inputFile "input_graphs/CA-GrQc"
