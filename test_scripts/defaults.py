from enum import IntEnum

default_input_graph_paths = ["CA-GrQc", "CA-HepTh"]

class Computation(IntEnum):
    SERIAL = 0
    PARALLEL = 1
    DISTRIBUTED = 2

granularity = 100

strategies = [2, 3, 2]

repeat = len(strategies)

workers = [1, 2, 4, 8]

time = ["10:00"] * repeat
mem = [4] * repeat
partition = ["slow"] * repeat

MPI_cpus = 1
MPI_nodes = 1