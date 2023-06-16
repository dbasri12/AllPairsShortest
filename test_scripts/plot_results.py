import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn
seaborn.set()

from defaults import *

CAGrQc_serial = [7.02592, 0.000930071]

CAGrQc_parallel_strategy_1 = [8.299242,4.452193,2.541873,2.373522]

CAGrQc_parallel_strategy_2 = [0.001560,0.007331,0.005717,0.021456]

CAGrQc_parallel_strategy_3 = [8.832310,13.419844,9.555574,21.476067]

CAGrQc_MPI_processes_1 = [110.035]
CAGrQc_MPI_processes_2 = []
CAGrQc_MPI_processes_4 = []
CAGrQc_MPI_processes_8 = []


CAHepTh_serial = [43.7203, 0.00216007]

CAHepTh_parallel_strategy_1 = [49.278850, 25.301551, 13.552267, 12.832467]

CAHepTh_parallel_strategy_2 = [0.004112, 0.010091, 0.016396, 0.037774]

CAHepTh_parallel_strategy_3 = [54.190907, 73.597726, 96.827684, 67.410748]

CAHepTh_MPI_strategy_1 = [178.053, 163.012, 150.210, 148.2]
CAHepTh_MPI_strategy_2 = [160.460, 155.223, 147.334, 147.2]

def main():
    # np.linspace(1, 8, CAGrQc_serial[0], endpoint=True)
    plt.plot(np.linspace(1, 8, 4, endpoint=True), np.zeros(4), ':', label='Bellman-Ford')
    plt.plot(CAGrQc_serial[1], '--', label='Dijkstra')
    plt.plot(workers, CAHepTh_parallel_strategy_1, label='Parallel Bellman-Ford')
    plt.plot(workers, CAHepTh_parallel_strategy_2, '-.', label='Parallel Dijkstra')
    plt.title('CA-HepTh Benchmarks')
    plt.xlabel('Cores/Threads')
    plt.ylabel('Execution Time (s)')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
