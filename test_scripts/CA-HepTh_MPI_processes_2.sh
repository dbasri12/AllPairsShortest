#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=4G
#SBATCH --time=10:00
#SBATCH --partition=slow

srun ./SSSP_MPI --sourceNode 24325 --strategy 2 --inputFile "input_graphs/CA-HepTh"
srun ./SSSP_MPI --sourceNode 24325 --strategy 1 --inputFile "input_graphs/CA-HepTh"
