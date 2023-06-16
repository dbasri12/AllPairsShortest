import sys
import os
from defaults import *


# Generates all parallel SSSP scripts for execution on cluster for an input graph with a given source node
def main(input_graph_name, source_node=0):
    # SERIAL SCRIPT
    f = open(f"{input_graph_name}_serial.sh", "w")
    f.write(
f'''#!/bin/bash
#
#SBATCH --cpus-per-task={1}
#SBATCH --time={time[Computation.SERIAL]}
#SBATCH --mem={mem[Computation.SERIAL]}G
#SBATCH --partition={partition[Computation.SERIAL]}

''')
    for strategy in range(1, strategies[Computation.SERIAL] + 1):
        f.write(f'srun ./SSSP --sourceNode {source_node} --strategy {strategy} --inputFile "input_graphs/{input_graph_name}"\n')
    f.close()

    # ALL PARALLEL SCRIPTS
    for strategy in range(1, strategies[Computation.PARALLEL] + 1):
        f = open(f"{input_graph_name}_parallel_strategy_{strategy}.sh", "w")
        f.write(
f'''#!/bin/bash
#
#SBATCH --cpus-per-task={workers[-1]}
#SBATCH --time={time[Computation.PARALLEL]}
#SBATCH --mem={mem[Computation.PARALLEL]}G
#SBATCH --partition={partition[Computation.PARALLEL]}

''')
        for threads in workers:
            f.write(f'srun ./SSSP_parallel --nThreads {threads} --sourceNode {source_node} --strategy {strategy} --granularity {granularity} --inputFile "input_graphs/{input_graph_name}"\n')
        f.close()
    
    # ALL DISTRIBUTED SCRIPTS
    for processes in workers:
        f = open(f"{input_graph_name}_MPI_processes_{processes}.sh", "w")
        f.write(
f'''#!/bin/bash
#
#SBATCH --cpus-per-task={MPI_cpus}
#SBATCH --nodes={MPI_nodes}
#SBATCH --ntasks={processes}
#SBATCH --mem={mem[Computation.DISTRIBUTED]}G
#SBATCH --time={time[Computation.DISTRIBUTED]}
#SBATCH --partition={partition[Computation.DISTRIBUTED]}

''')
        for strategy in range(strategies[Computation.DISTRIBUTED], 0, -1):
            f.write(f'srun ./SSSP_MPI --sourceNode {source_node} --strategy {strategy} --inputFile "input_graphs/{input_graph_name}"\n')
        f.close()


# Give either:
# No Arguments  :   Default list of input graphs on source node 0
# One Argument  :   Default list of input graphs on given source node
# More Arguments:   Source node using rest of arguments as input graphs
if __name__ == '__main__':
    if len(sys.argv) == 1:
        for i in default_input_graph_paths:
            main(i)
    elif len(sys.argv) == 2:
        for i in default_input_graph_paths:
            main(i, sys.argv[1])
    else:
        for i in range(2, len(sys.argv)):
            main(sys.argv[i], sys.argv[1])
