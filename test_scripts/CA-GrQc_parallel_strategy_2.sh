#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00
#SBATCH --mem=4G
#SBATCH --partition=slow

srun ./SSSP_parallel --nThreads 1 --sourceNode 3466 --strategy 2 --granularity 200 --inputFile "input_graphs/CA-GrQc"
srun ./SSSP_parallel --nThreads 2 --sourceNode 3466 --strategy 2 --granularity 200 --inputFile "input_graphs/CA-GrQc"
srun ./SSSP_parallel --nThreads 4 --sourceNode 3466 --strategy 2 --granularity 200 --inputFile "input_graphs/CA-GrQc"
srun ./SSSP_parallel --nThreads 8 --sourceNode 3466 --strategy 2 --granularity 200 --inputFile "input_graphs/CA-GrQc"
