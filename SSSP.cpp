#include <stdlib.h>

#include <iomanip>
#include <iostream>
#include <vector>
#include <queue>
#include <functional>

#include "core/graph.h"
#include "core/utils.h"

#define DEFAULT_STRATEGY    "1"

void SSSP_BellmanFord(Graph &g, int source) {
    uintV n = g.n_;
    uintE m = g.m_;
    int *length = new int[n];

    // Initialization
    for (int i = 0; i < n; i++) {
        length[i] = INT_MAX;
    }
    length[source] = 0;

    timer t1;
    double time_taken = 0.0;

    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    t1.start();
    for (int round = 1; round < n-1; round++) {
       // Update the shortest path to the source for each vertex
       for (int u = 0; u < n; u++) {
           uintE neighbors = g.vertices_[u].getOutDegree();
           for (int i = 0; i < neighbors; i++) {
                uintV v = g.vertices_[u].getOutNeighbor(i);

                // length != INT_MAX to avoid bit roll-over when add 1
                if (length[u] != INT_MAX && length[u] + 1 < length[v]) {
                    length[v] = length[u] + 1;
                }
           }
       }
    }

    time_taken = t1.stop();
    // -------------------------------------------------------------------
#ifdef PRINT_DEBUG
    // Print the statistics for each thread
    // Example output for 2 threads:
    // thread_id, source, dest, length, time_taken
    // 0, 0, 0, 0, 0.12
    // 1, 0, 1, 1, 0.12
    std::cout << "thread_id, source, dest, length, time_taken" << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << "0, " << source << ", " << i  <<", " << length[i] << ", " << time_taken << std::endl;
    }
#endif

    std::cout << "Time taken (in seconds) : " << time_taken << "\n";

    delete[] length;

    return;
}

void SSSP_Dijkstra(Graph &g, int source) {
    uintV n = g.n_;
    uintE m = g.m_;
    
    std::vector<int> length(n, INT_MAX);
    std::vector<int> visited(n, false);
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;

    length[source] = 0;
    pq.push({0, source});

    timer t1;
    double time_taken = 0.0;

    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    t1.start();
    // Update the shortest path to the source for each vertex
    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();
        if (visited[u]) continue;
        visited[u] = true;
        uintE neighbors = g.vertices_[u].getOutDegree();
        for (int i = 0; i < neighbors; i++) {
            uintV v = g.vertices_[u].getOutNeighbor(i);

            // length != INT_MAX to avoid bit roll-over when add 1
            if (length[u] != INT_MAX && length[u] + 1 < length[v]) {
                length[v] = length[u] + 1;
                pq.push({length[v], v});
            }
        }
    }

    time_taken = t1.stop();
    // -------------------------------------------------------------------
#ifdef PRINT_DEBUG
    // Print the statistics for each thread
    // Example output for 2 threads:
    // thread_id, source, dest, length, time_taken
    // 0, 0, 0, 0, 0.12
    // 1, 0, 1, 1, 0.12
    std::cout << "thread_id, source, dest, length, time_taken" << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << "0, " << source << ", " << i  <<", " << length[i] << ", " << time_taken << std::endl;
    }
#endif
   
    std::cout << "Time taken (in seconds) : " << time_taken << "\n";

    return;
}

int main(int argc, char *argv[]) {
    cxxopts::Options options(
        "page_rank_push",
        "Calculate page_rank using serial and parallel execution");
    options.add_options(
        "",
        {
            {"sourceNode", "The source node used to calculate distance to all other nodes",
                cxxopts::value<uint>()->default_value("0")},
            {"inputFile", "Input graph file path",
            cxxopts::value<std::string>()->default_value(
                "/scratch/input_graphs/roadNet-CA")},
                {"strategy", "Strategy for decomposition",
                cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
        });

    auto cl_options = options.parse(argc, argv);
    uint source_node = cl_options["sourceNode"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();
    uint strategy = cl_options["strategy"].as<uint>();

    std::cout << "Source node : " << source_node << std::endl;
    std::cout << "Strategy: " << strategy << std::endl;

    Graph g;
    std::cout << "Reading graph\n";
    g.readGraphFromBinary<int>(input_file_path);
    std::cout << "Created graph\n";

    switch (strategy) {
    case 1:
        std::cout << "Serial Bellman-Ford SSSP" << std::endl;
        SSSP_BellmanFord(g, source_node);
        break;
    case 2:
        std::cout << "Serial Dijkstra's SSSP" << std::endl;
        SSSP_Dijkstra(g, source_node);
        break;
    }

    return 0;
}
