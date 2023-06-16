#include <stdlib.h>
#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <vector>
#include <queue>
#include <functional>
#include <algorithm>
#include <utility>

#include "core/graph.h"
#include "core/utils.h"

#define DEFAULT_STRATEGY    "1"

/******************************************************************************
 **              Bellman-Ford Algorithm
 ******************************************************************************/

void calcDistance(Graph &g, int *length, int start,int end) {
       // Update the shortest path to the source for each vertex
       for (uintV u = start; u < end; u++) {
           uintE neighbors = g.vertices_[u].getOutDegree();
           for (uintE i = 0; i < neighbors; i++) {
                uintV v = g.vertices_[u].getOutNeighbor(i);

                // length != INT_MAX to avoid bit roll-over when add 1
                if (length[u] != INT_MAX && length[u] + 1 < length[v]) {
                    length[v] = length[u] + 1;
                }
           }
       }
}

void SSSP_BellmanFord(Graph &g, const int &source_node, const int &world_rank, const int &world_size) {
    uintV n = g.n_;
    int *length = new int[n];

    // Initialization
    for (int i = 0; i < n; i++) {
        length[i] = INT_MAX;
    }
    length[source_node] = 0;

    uint startx;
    uint endx;
    int ver_per_t = n / world_size;
    if(world_rank == world_size - 1) {
        int excess = n % world_size;
        startx = world_rank * ver_per_t;
        endx = startx + ver_per_t + excess;
    } else {
        startx = world_rank * ver_per_t;
        endx = startx + ver_per_t;
    }
    
    timer t , t2;
    t.start();
    t2.start();
    double time_taken = 0.0;
    double time_taken2 = 0.0;

    for (uintV round = 1; round < n-1; round++) {
        calcDistance(g, length, startx, endx);
        for(int i=0; i<n; i++) {
            int local_temp[2];
            local_temp[0] = length[i];
            local_temp[1] = world_rank;

            int global_temp[2];
            MPI_Allreduce(&local_temp,&global_temp,1,MPI_2INT,MPI_MINLOC,MPI_COMM_WORLD);
            if(global_temp[0] < local_temp[0]) {
                length[i] = global_temp[0];
            }
        }   
    }
    time_taken=t.stop();

#ifdef PRINT_DEBUG
    for (int j = startx; j < endx; j++) {
            printf("%d, %d, %d, %d, %.3f\n", world_rank, source_node, j, length[j], time_taken);
    }
#endif

    MPI_Barrier(MPI_COMM_WORLD);
    time_taken2=t2.stop();

    if(world_rank==0){
#ifdef PRINT_DEBUG
        printf("Total time taken (s) : %.3f\n", time_taken2);
#else
        printf("%.3f\n", time_taken2);
#endif
    }

    delete[] length;

    return;
}

/******************************************************************************
 **              Dijkstra Algorithm
 ******************************************************************************/

typedef std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> PriorityQueue;

uintV findMin(PriorityQueue &pq, const int &rank, const int &world_size) {
    std::pair<int, int> ver_pair;
    if (!pq.empty()) {
        ver_pair = pq.top();
    } else {
        ver_pair = std::make_pair(INT_MAX, 0);
    }

    // printf("Before gather %d\n", rank);
    std::vector<std::pair<int, int>> recv_buf(world_size); 
    MPI_Gather(&ver_pair, sizeof(std::pair<int, int>), MPI_BYTE, 
                (void*)recv_buf.data(), sizeof(std::pair<int, int>), 
                MPI_BYTE, 0, MPI_COMM_WORLD);
    
    // Find the min pair in the root
    if (rank == 0) {
        ver_pair = *std::min_element(recv_buf.begin(), recv_buf.end(), [](const std::pair<int,int> &a, const std::pair<int,int> &b) {
            return a.first < b.first;
        });
    }
    MPI_Bcast(&ver_pair, sizeof(std::pair<int, int>), MPI_BYTE, 0, MPI_COMM_WORLD);

    if (ver_pair.second == pq.top().second) {
        pq.pop();
    }

    return ver_pair.second;
}

void relax(PriorityQueue &pq, const uintV &v, const uintV &u, int *length) {
    uint temp = length[v] + 1;
    if (length[u] > temp) {
        length[u] = temp;
        
        PriorityQueue temp;
        while (!pq.empty()) {
            std::pair<int, int> q = pq.top();
            pq.pop();
            if (q.second == u) {
                temp.push({length[u], u});
            } else {
                temp.push(q);
            }
        }
        // std::swap(temp, pq);
        pq = temp;
    }
}

void SSSP_Dijkstra(Graph &g, const int &source, const int &world_size, const int &world_rank) {
    uintV n = g.n_;
    
    int *length = new int [n];
    bool *visited = new bool[n];
    for (int i = 0; i < n; i++) {
        length[i] = INT_MAX;
        visited[i] = false;
    }

    int numOfParts = world_size;
    if (world_size > n) {
        numOfParts = n;
    }

    if (world_rank >= numOfParts) {
        return;
    }

    uintV *start = new uintV[numOfParts];
    uintV *end = new uintV[numOfParts];

    // Decomposes the vertices to each process
    uintV ver_per_p = n/numOfParts;
    for (int i = 0; i < numOfParts; i++) {
        start[i] = i*ver_per_p;
        end[i] = start[i] + ver_per_p;
    }
    end[numOfParts-1] = n;

    length[source] = 0;

    PriorityQueue pq;
    for (uintV v = start[world_rank]; v < end[world_rank]; v++) {
        pq.push({length[v], v});
    }

    int global_pq_size;
    int pq_size = pq.size();
    MPI_Allreduce(&pq_size, &global_pq_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    timer t1;
    timer global_time;
    double time_taken = 0.0;
    double time_taken2 = 0.0;

    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    global_time.start();
    t1.start();
    // Update the shortest path to the source for each vertex
    while (global_pq_size > 0) {
        uintV v = findMin(pq, world_rank, world_size);

        // If the vertex has already been visited, don't do it again
        if (visited[v]) continue;
        visited[v] = true;

        uint neighbor = g.vertices_[v].getOutDegree();
        for (int i = 0; i < neighbor; i++) {
            uintV u = g.vertices_[v].getOutNeighbor(i);
            relax(pq, v, u, length);
        }

        pq_size = pq.size();
        MPI_Allreduce(&pq_size, &global_pq_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }

    time_taken = t1.stop();
    // -------------------------------------------------------------------

#ifdef PRINT_DEBUG
    // Print the statistics for each thread
    // Example output for 2 threads:
    // process_id, source, dest, length, time_taken
    // 0, 0, 0, 0, 0.12
    // 1, 0, 1, 1, 0.12
      for (int i = 0; i < n; i++) {
        if (i >= start[world_rank] && i < end[world_rank]) {
            printf("%d, %d, %d, %d, %.3f\n", world_rank, source, i, length[i], time_taken);
        }
    }
#endif

    MPI_Barrier(MPI_COMM_WORLD);
    time_taken2 = global_time.stop();

    if(world_rank == 0){
#ifdef PRINT_DEBUG
        printf("Total time taken (s) : %.3f\n", time_taken2);
#else
        printf("%.3f\n", time_taken2);
#endif
    }


    delete[] length;
    delete[] visited;
    delete[] start;
    delete[] end;

    return;
}

/******************************************************************************
 **              Main
 ******************************************************************************/
int main(int argc, char *argv[]) {
  // Initialize command line arguments

  cxxopts::Options options(
        "page_rank_push",
        "Calculate page_rank using serial and parallel execution");
    options.add_options(
        "",
        {
            {"sourceNode", "The source node used to calculate distance to all other nodes",
                cxxopts::value<uint>()->default_value("0")},
            {"inputFile", "Input graph file path",
                cxxopts::value<std::string>()->default_value("/scratch/input_graphs/roadNet-CA")},
            {"strategy", "Strategy for decomposition",
                cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
        });

    auto cl_options = options.parse(argc, argv);
    uint source_node = cl_options["sourceNode"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();
    uint strategy = cl_options["strategy"].as<uint>();

    MPI_Init(NULL,NULL);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);

    if(world_rank==0){
#ifdef PRINT_DEBUG
        std::cout << std::fixed;
        //std::cout << "Number of Threads : " << n_threads << std::endl;
        std::cout << "Source node : " << source_node << std::endl;
        std::cout << "Strategy: " << strategy << std::endl;
#endif
    }

    Graph g;
    if(world_rank==0){
#ifdef PRINT_DEBUG
        std::cout << "Reading graph\n";
#endif
    }
    g.readGraphFromBinary<int>(input_file_path);

    if(world_rank==0){
#ifdef PRINT_DEBUG
        std::cout << "Created graph\n";
        std::cout << "MPI Bellman-Ford"<< std::endl;
#endif
#ifdef PRINT_DEBUG
        std::cout<<"processor_id, source, dest, length, time_taken"<<std::endl;
#endif
    }

    switch(strategy) {
        case 1:
#ifdef PRINT_DEBUG
            if(world_rank==0)
                printf("Bellman-Ford algorithm\n");
#endif
            SSSP_BellmanFord(g, source_node, world_rank, world_size);
            break;
        case 2:
#ifdef PRINT_DEBUG
            if(world_rank==0)
                printf("Dijkstra algorithm\n");
#endif
            SSSP_Dijkstra(g, source_node, world_size, world_rank);
            break;
        default:
            if(world_rank==0)
                printf("Please choose a valid strategy\n");
            break;
    }

    MPI_Finalize();
    return 0;

}




