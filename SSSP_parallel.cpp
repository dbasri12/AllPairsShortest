#include <stdlib.h>

#include <iomanip>
#include <iostream>
#include <thread>
#include <atomic>
#include <vector>
#include <queue>
#include <functional>
#include <mutex>

#include "core/graph.h"
#include "core/utils.h"

#define DEFAULT_GRANULARITY "0"
#define DEFAULT_STRATEGY    "1"

typedef struct ThreadData {
    uintV start;
    uintV end;
    double time_taken;
} Data_T;

/******************************************************************************
 **              Bellman-Ford Vertex-based decomposition
 ******************************************************************************/

void BellmanFordCalcDistance(Graph &g, const uintV n, std::atomic<int> *length, Data_T &data, CustomBarrier &barrier) {
    timer t;
    t.start();
    for (uintV round = 1; round < n-1; round++) {
       // Update the shortest path to the source for each vertex
       for (uintV u = data.start; u < data.end; u++) {
           uintE neighbors = g.vertices_[u].getOutDegree();
           for (uintE i = 0; i < neighbors; i++) {
                uintV v = g.vertices_[u].getOutNeighbor(i);

                // length != INT_MAX to avoid bit roll-over when add 1
                int weight_v = length[v].load();
                int weight_u = length[u].load();
                if (weight_u != INT_MAX && weight_u + 1 < weight_v) {
                    while(!length[v].compare_exchange_weak(weight_v, weight_u + 1,
                                        std::memory_order_release,
                                        std::memory_order_relaxed))
                    ; // the body of the loop is empty
                }
           }
       }
       barrier.wait();
    }
    data.time_taken = t.stop();
}

void SSSP_BellmanFord_vertexDecomp(Graph &g, const int source, const int n_threads) {
    uintV n = g.n_;
    std::atomic<int> *length = new std::atomic<int>[n];

    // Initialization
    for (int i = 0; i < n; i++) {
        if (i != source) {
            atomic_init(&length[i], INT_MAX);
        }
    }
    atomic_init(&length[source], 0);

    // If there are more threads than the number of vertices,
    // then limit the number of threads used to the number
    // of edges
    // Reason: We can't process a fraction of a vertex, can we?
    int thread_num = n_threads;
    if (n_threads > n) {
        std::cout << "There number of threads requested exceeded the number of vertices. Default the number of threads to: " << thread_num << std::endl;
        thread_num = n;
    }

    // Thread data
    Data_T *data = new Data_T[thread_num];
    std::thread *ts = new std::thread[thread_num];
    CustomBarrier barrier(thread_num);

    // Set the starting vertex and calculate
    // the number of edges per thread
    int ver_per_t = n/thread_num;
    for (int i = 0; i < thread_num; i++) {
        data[i].start = i*ver_per_t;
        data[i].end = data[i].start + ver_per_t;
    }
    data[thread_num].end = n;

    timer t1;
    double time_taken = 0.0;

    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    t1.start();
    for (int i = 0; i < thread_num; i++) {
        ts[i] = std::thread(BellmanFordCalcDistance, std::ref(g), n, std::ref(length),
                                std::ref(data[i]), std::ref(barrier));
    }

    for (int i = 0; i < thread_num; i++) {
        ts[i].join();
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
    for (int i = 0; i < n_threads; i++) {
        for (int j = data[i].start; j < data[i].end; j++) {
            std::cout << i << ", " << source << ", " << j <<", " << length[j] << ", " << data[i].time_taken << std::endl;
        }
    }
    std::cout << "Time taken (in seconds) : ";
#endif

    std::cout << time_taken << "\n";

    delete[] length;
    delete[] data;
    delete[] ts;
    return;
}

/******************************************************************************
 **              Dijkstra's Parallel Algorithm
 ******************************************************************************/

void DijkstraCalcDistance(Graph &g, const uintV n, int n_threads, std::atomic<int> *&length, std::atomic<int> &done_counter, std::mutex &pq_mutex, std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> &pq, double &time_taken) {
    timer t_time;
    t_time.start();

    while(true) {
        if (pq.empty()) {
            done_counter++;
            while (pq.empty()) {
                if (done_counter.load() == n_threads) break;
            }
            if (done_counter.load() == n_threads) break;
            done_counter--;
            continue;
        }
        pq_mutex.lock();
        if (pq.empty()) {
            pq_mutex.unlock();
            continue;
        }
        int u = pq.top().second;
        pq.pop();
        pq_mutex.unlock();
        uintE neighbors = g.vertices_[u].getOutDegree();
        for (uintE i = 0; i < neighbors; i++) {
            uintV v = g.vertices_[u].getOutNeighbor(i);

            // uVal != INT_MAX to avoid bit roll-over when add 1
            int uVal = length[u].load();
            int vVal = length[v].load();
            if (uVal != INT_MAX && uVal + 1 < vVal && length[v].compare_exchange_weak(vVal, uVal + 1)) {
                pq_mutex.lock();
                pq.push({uVal + 1, v});
                pq_mutex.unlock();
            }
        }
    }

    time_taken = t_time.stop();
}

void SSSP_Dijkstra_vertexDecomp(Graph &g, const int source, const int n_threads) {
    uintV n = g.n_;

    std::atomic<int> *length = new std::atomic<int>[n];
    std::atomic<int> done_counter(0);
    std::mutex pq_mutex;
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;

    // Initialization
    for (int i = 0; i < n; i++) {
        if (i != source) {
            atomic_init(&length[i], INT_MAX);
        }
    }
    atomic_init(&length[source], 0);
    pq.push({0, source});

    // If there are more threads than the number of vertices,
    // then limit the number of threads used to the number
    // of edges
    // Reason: We can't process a fraction of a vertex, can we?
    int thread_num = n_threads;
    if (n_threads > n) {
        std::cout << "There number of threads requested exceeded the number of vertices. Default the number of threads to: " << thread_num << std::endl;
        thread_num = n;
    }

    // Thread data
    std::thread *ts = new std::thread[thread_num];
    double *t_time = new double[thread_num];

    timer t1;
    double time_taken = 0.0;

    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    t1.start();
    for (int i = 0; i < n_threads; i++) {
        ts[i] = std::thread(DijkstraCalcDistance, std::ref(g), n, thread_num, std::ref(length),
                                std::ref(done_counter), std::ref(pq_mutex), std::ref(pq), std::ref(t_time[i]));
    }

    for (int i = 0; i < thread_num; i++) {
        ts[i].join();
    }
    time_taken = t1.stop();
    // -------------------------------------------------------------------
#ifdef PRINT_DEBUG
    // Print the statistics for each thread
    // Example output for 2 threads:
    // thread_id, source, dest, length, time_taken
    // 0, 0, 0, 0, 0.12
    // 1, 0, 1, 1, 0.12
    std::cout << "source, dest, length" << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << source << ", " << i <<", " << length[i] << std::endl;
    }
    std::cout << "Time taken (in seconds) : ";
#endif

    std::cout << time_taken << "\n";

    delete[] length;
    delete[] ts;
    delete[] t_time;
    return;
}

/******************************************************************************
 **              Vertex-based coarse grained decomposition
 ******************************************************************************/

uintV getNextVertexToBeProcessed(const uint &granularity, const uintV &n, std::atomic<uintV> &counter) {
    uintV temp = counter.fetch_add(granularity, std::memory_order_relaxed);

    if(temp >= n) {
      return -1;
    }

    return temp;
}

void calcDistanceDynamic(Graph &g, const uintV n, std::atomic<int> *&length, std::atomic<uintV> &vertex_counter, const uint granularity, double &time_taken) {
    timer t_time;
    t_time.start();

    const uint k = granularity;
    for (uintV round = 1; round < n-1; round++) {
        while(true) {
            uintV u = getNextVertexToBeProcessed(granularity, n, vertex_counter);

            if(u == -1) 
                break;

            for (uint j = 0; j < k && u < n; j++, u++) {
               uintE neighbors = g.vertices_[u].getOutDegree();
                for (uintE i = 0; i < neighbors; i++) {
                    uintV v = g.vertices_[u].getOutNeighbor(i);

                    // uVal != INT_MAX to avoid bit roll-over when add 1
                    for(int uVal = length[u].load(), vVal = length[v].load(); uVal != INT_MAX && uVal + 1 < vVal && !length[v].compare_exchange_weak(vVal, uVal + 1);) {;}
                }
            }
        }

        // Reset counter
        if(vertex_counter.load() != 0) {
            vertex_counter = 0;
        }
    }

    time_taken = t_time.stop();
}

void SSSP_dynamicDecomp(Graph &g, const int source, const int n_threads, uint granularity) {
    uintV n = g.n_;
    std::atomic<int> *length = new std::atomic<int>[n];
    std::atomic<uintV> vertex_counter(0);

    // Initialization
    for (int i = 0; i < n; i++) {
        if (i != source) {
            atomic_init(&length[i], INT_MAX);
        }
    }
    atomic_init(&length[source], 0);

    // If there are more threads than the number of vertices,
    // then limit the number of threads used to the number
    // of edges
    // Reason: We can't process a fraction of a vertex, can we?
    int thread_num = n_threads;
    if (n_threads > n) {
        std::cout << "There number of threads requested exceeded the number of vertices. Default the number of threads to: " << thread_num << std::endl;
        thread_num = n;
    }

    // Thread data
    std::thread *ts = new std::thread[thread_num];
    double *t_time = new double[thread_num];

    timer t1;
    double time_taken = 0.0;

    // Create threads and distribute the work across T threads
    // -------------------------------------------------------------------
    t1.start();
    for (int i = 0; i < n_threads; i++) {
        ts[i] = std::thread(calcDistanceDynamic, std::ref(g), n, std::ref(length), std::ref(vertex_counter), 
                                granularity, std::ref(t_time[i]));
    }

    for (int i = 0; i < n_threads; i++) {
        ts[i].join();
    }
    time_taken = t1.stop();
    // -------------------------------------------------------------------
#ifdef PRINT_DEBUG
    // Print the statistics for each thread
    // Example output for 2 threads:
    // thread_id, source, dest, length, time_taken
    // 0, 0, 0, 0, 0.12
    // 1, 0, 1, 1, 0.12
    std::cout << "source, dest, length" << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << source << ", " << i <<", " << length[i] << std::endl;
    }
    std::cout << "Time taken (in seconds) : ";
#endif

    std::cout << time_taken << "\n";

    delete[] length;
    delete[] ts;
    delete[] t_time;
    return;
}

/******************************************************************************
 **              Main
 ******************************************************************************/

int main(int argc, char *argv[]) {
    cxxopts::Options options(
        "page_rank_push",
        "Calculate page_rank using serial and parallel execution");
    options.add_options(
        "",
        {
            {"nThreads", "Number of Threads",
                cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_THREADS)},
            {"sourceNode", "The source node used to calculate distance to all other nodes",
                cxxopts::value<uint>()->default_value("0")},
            {"inputFile", "Input graph file path",
                cxxopts::value<std::string>()->default_value("/scratch/input_graphs/roadNet-CA")},
            {"strategy", "Strategy for decomposition",
                cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
            {"granularity", "The granularity of vertex decomposition",
                cxxopts::value<uint>()->default_value(DEFAULT_GRANULARITY)},
        });

    auto cl_options = options.parse(argc, argv);
    uint n_threads = cl_options["nThreads"].as<uint>();
    uint source_node = cl_options["sourceNode"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();
    uint strategy = cl_options["strategy"].as<uint>();
    uint granularity = cl_options["granularity"].as<uint>();
    if (strategy >= 1 && strategy <= 2)
        granularity = 1;

    std::cout << std::fixed;
#ifdef PRINT_DEBUG
    std::cout << "Number of Threads : " << n_threads << std::endl;
    std::cout << "Source node : " << source_node << std::endl;
    std::cout << "Strategy: " << strategy << std::endl;
    std::cout << "Granularity: " << granularity << std::endl;
#endif

    Graph g;
#ifdef PRINT_DEBUG
    std::cout << "Reading graph\n";
#endif
    g.readGraphFromBinary<int>(input_file_path);
#ifdef PRINT_DEBUG
    std::cout << "Created graph\n";
#endif

    switch (strategy) {
    case 1:
#ifdef PRINT_DEBUG
        std::cout << "Bellman-Ford Vertex-based decomposition" << std::endl;
#endif
        SSSP_BellmanFord_vertexDecomp(g, source_node, n_threads);
        break;
    case 2:
#ifdef PRINT_DEBUG
        std::cout << "Dijkstra's Parallel Algorithm" << std::endl;
#endif
        SSSP_Dijkstra_vertexDecomp(g, source_node, n_threads);
        break;
    case 3:
#ifdef PRINT_DEBUG
        std::cout << "Bellman-Ford Vertex-based coarse grained decomposition" << std::endl;
#endif
        SSSP_dynamicDecomp(g, source_node, n_threads, granularity);
        break;
    default:
        std::cout << "Please choose a valid strategy for decomposition\n";
        break;
    }

    return 0;
}
