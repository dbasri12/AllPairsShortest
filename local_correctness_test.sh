#!/bin/bash

make DEBUG=1

echo " "
./SSSP --sourceNode 0 --strategy 1 --inputFile "input_graphs/GG_Test"
echo " "
./SSSP --sourceNode 0 --strategy 2 --inputFile "input_graphs/GG_Test"
echo " "
./SSSP_parallel --nThreads 4 --sourceNode 0  --strategy 1 --inputFile "input_graphs/GG_Test"
echo " "
./SSSP_parallel --nThreads 4 --sourceNode 0  --strategy 2 --inputFile "input_graphs/GG_Test"
echo " "
./SSSP_parallel --nThreads 4 --sourceNode 0  --strategy 3 --granularity 3 --inputFile "input_graphs/GG_Test"
echo " "
./SSSP --sourceNode 2 --strategy 1 --inputFile "input_graphs/GG_Test"
echo " "
./SSSP --sourceNode 2 --strategy 2 --inputFile "input_graphs/GG_Test"
echo " "
./SSSP_parallel --nThreads 4 --sourceNode 2  --strategy 1 --inputFile "input_graphs/GG_Test"
echo " "
./SSSP_parallel --nThreads 4 --sourceNode 2  --strategy 2 --inputFile "input_graphs/GG_Test"
echo " "
./SSSP_parallel --nThreads 4 --sourceNode 2  --strategy 3 --granularity 3 --inputFile "input_graphs/GG_Test"

make clean