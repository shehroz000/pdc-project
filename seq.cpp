#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <map>
#include <set>
#include <string> // For std::to_string
#include <algorithm> // For std::remove_if, std::find_if
#include <cmath>     // For std::abs

// --- METIS Header ---
// You MUST have METIS installed and this header accessible.
// Ensure your compiler can find it (e.g., using -I/path/to/metis/include)
#include <metis.h>

// Define infinity - choose a value larger than any possible path sum
const double INF = std::numeric_limits<double>::infinity();
const int NO_PREDECESSOR = -1;

// Represents an edge: target vertex and weight (for SSSP)
struct Edge {
    int to;
    double weight;
    // For METIS graph construction, we might need to track original edge index or if it's a reverse edge
};

// Graph representation: Adjacency list (for SSSP)
using Graph = std::vector<std::vector<Edge>>;

// Structure for Priority Queue elements: {distance, vertex}
using PrioQueueEntry = std::pair<double, int>;

// --- SSSP Functions (dijkstra, update_sssp, print_paths) ---
// (These functions are copied from the end of Step 2's response)
// (Scroll down to see them in their original place if needed; for brevity, not repeated inline here,
//  but assume they are present in your .cpp file)
void dijkstra(int start_node, const Graph& graph, std::vector<double>& dist, std::vector<int>& pred) {
    int num_vertices = graph.size();
    dist.assign(num_vertices, INF);
    pred.assign(num_vertices, NO_PREDECESSOR);

    dist[start_node] = 0.0;
    std::priority_queue<PrioQueueEntry, std::vector<PrioQueueEntry>, std::greater<PrioQueueEntry>> pq;
    pq.push({0.0, start_node});

    while (!pq.empty()) {
        double d = pq.top().first;
        int u = pq.top().second;
        pq.pop();

        if (d > dist[u]) {
            continue;
        }

        for (const Edge& edge : graph[u]) {
            int v = edge.to;
            double weight = edge.weight;
            if (dist[u] != INF && dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                pred[v] = u;
                pq.push({dist[v], v});
            }
        }
    }
}

void update_sssp(
    int u_changed, int v_changed, double new_weight,
    int start_node, Graph& graph, std::vector<double>& dist, std::vector<int>& pred
) {
    std::cout << "\nUpdating edge (" << u_changed << " -> " << v_changed
              << ") with new weight " << (new_weight == INF ? "INF (deleted)" : std::to_string(new_weight))
              << std::endl;

    double old_weight = INF;
    bool edge_existed = false;

    auto& u_edges = graph[u_changed];
    auto it_edge = std::find_if(u_edges.begin(), u_edges.end(),
                                [v_changed](const Edge& e){ return e.to == v_changed; });

    if (it_edge != u_edges.end()) {
        old_weight = it_edge->weight;
        edge_existed = true;
        if (new_weight == INF) { // Deletion
            u_edges.erase(it_edge);
        } else { // Weight update
            it_edge->weight = new_weight;
        }
    }

    if (!edge_existed && new_weight != INF) {
        std::cout << "  Edge did not exist. Inserting edge (" << u_changed << " -> " << v_changed
                  << ") with weight " << new_weight << std::endl;
        graph[u_changed].push_back({v_changed, new_weight});
        old_weight = INF;
    } else if (!edge_existed && new_weight == INF) {
         std::cout << "  Attempting to delete non-existent edge. No change." << std::endl;
         return;
    }

    std::priority_queue<PrioQueueEntry, std::vector<PrioQueueEntry>, std::greater<PrioQueueEntry>> pq;

    if (new_weight < old_weight) {
        std::cout << "  Detected Weight Decrease/Insertion." << std::endl;
        if (dist[u_changed] != INF && dist[u_changed] + new_weight < dist[v_changed]) {
            std::cout << "  Potential improvement for vertex " << v_changed << " (" << dist[v_changed]
                      << " -> " << dist[u_changed] + new_weight << "). Adding to update queue." << std::endl;
            dist[v_changed] = dist[u_changed] + new_weight;
            pred[v_changed] = u_changed;
            pq.push({dist[v_changed], v_changed});
        } else {
             std::cout << "  No immediate improvement for vertex " << v_changed << ". No update triggered by this edge directly." << std::endl;
        }
    }
    else if (new_weight > old_weight || (new_weight == INF && edge_existed) ) { // Check edge_existed for true deletion case
         std::cout << "  Detected Weight Increase/Deletion." << std::endl;
         double tolerance = 1e-9;
         bool was_on_shortest_path = (pred[v_changed] == u_changed) &&
                                    (old_weight != INF && std::abs((dist[u_changed] + old_weight) - dist[v_changed]) < tolerance);

        if (was_on_shortest_path) {
             std::cout << "  Edge (" << u_changed << "->" << v_changed << ") was on the shortest path to "
                       << v_changed << ". Re-evaluation needed." << std::endl;
             dist[v_changed] = INF;
             pred[v_changed] = NO_PREDECESSOR;
             pq.push({INF, v_changed});
             std::cout << "  Invalidated distance for " << v_changed << " and added to PQ for re-evaluation." << std::endl;
        } else {
             std::cout << "  Edge was not on the shortest path to " << v_changed << " or change is not impactful. No update triggered by this edge directly." << std::endl;
        }
    } else {
         std::cout << "  Weight unchanged or no effective change. No update needed." << std::endl;
         return;
    }

    int relaxation_count = 0;
    while (!pq.empty()) {
        double d = pq.top().first;
        int u = pq.top().second;
        pq.pop();

        if (d > dist[u] && d != INF) { // d != INF allows INF to be processed for re-evaluation
             continue;
        }
        
        // If d is INF, it means u needs re-evaluation. Attempt to find its best current path.
        // This is a simplification; paper's method might be more direct.
        if (d == INF) {
            double current_min_dist_for_u = INF;
            int current_pred_for_u = NO_PREDECESSOR;
            // This requires iterating through ALL nodes to find incoming edges to u
            // This part is inefficient and highlights the need for reverse graph or specific
            // algorithm for handling nodes that become unreachable from their previous parent.
            // The paper's specific algorithm for handling increases is CRITICAL here.
            // For this placeholder, we'll assume something magically updates dist[u]
            // or it remains INF. The propagation below will then work from this state.
            // A more robust placeholder would be to re-run Dijkstra from 'start_node' if many
            // nodes get INF here, or use the paper's method for finding alternative parents.
            // For now, we just ensure dist[u] is used as is.
            if(dist[u] == INF) { // if after external check it's still INF
                std::cout << "  Skipping propagation from vertex " << u << " as it remains disconnected (INF)." << std::endl;
                continue;
            }
            std::cout << "  Re-evaluating from vertex " << u << " (current dist=" << dist[u] << ")" << std::endl;
        } else {
           std::cout << "  Processing updates from vertex " << u << " (dist=" << d << ")" << std::endl;
        }

        for (const Edge& edge : graph[u]) {
            int v = edge.to;
            double weight = edge.weight;
            if (dist[u] != INF && dist[u] + weight < dist[v]) {
                relaxation_count++;
                 std::cout << "    Updating " << v << ": " << dist[v] << " -> " << dist[u] + weight << std::endl;
                dist[v] = dist[u] + weight;
                pred[v] = u;
                pq.push({dist[v], v});
            }
        }
    }
     std::cout << "Update propagation finished. Relaxations performed: " << relaxation_count << std::endl;
}

void print_paths(int start_node, const std::vector<double>& dist, const std::vector<int>& pred) {
    int n = dist.size();
    std::cout << "Distances and paths from node " << start_node << ":" << std::endl;
    for (int i = 0; i < n; ++i) {
        std::cout << "  Vertex " << i << ": Dist = ";
        if (dist[i] == INF) {
            std::cout << "INF";
        } else {
            std::cout << dist[i];
        }
        std::cout << ", Pred = " << pred[i] << std::endl;
    }
}


// --- METIS Integration Functions ---

// Converts our Graph to CSR format for METIS.
// For METIS, we treat the graph as undirected and use unit edge weights.
void convert_to_metis_csr(const Graph& directed_graph,
                          idx_t& nvtxs,
                          std::vector<idx_t>& xadj,
                          std::vector<idx_t>& adjncy,
                          std::vector<idx_t>& adjwgt) {
    nvtxs = directed_graph.size();
    xadj.clear();
    adjncy.clear();
    adjwgt.clear();

    xadj.resize(nvtxs + 1);

    // Create a temporary adjacency list for an undirected version
    std::vector<std::vector<std::pair<int, int>>> undirected_adj(nvtxs);
    for (idx_t u = 0; u < nvtxs; ++u) {
        for (const auto& edge : directed_graph[u]) {
            idx_t v = edge.to;
            // Add edge u-v (and v-u implicitly for METIS structure by adding symmetrically)
            // Ensure no duplicates if original graph could have u->v and v->u
            bool found_uv = false;
            for(const auto& p : undirected_adj[u]) if(p.first == v) found_uv = true;
            if(!found_uv) undirected_adj[u].push_back({v, 1}); // Unit weight

            bool found_vu = false;
            for(const auto& p : undirected_adj[v]) if(p.first == u) found_vu = true;
            if(!found_vu) undirected_adj[v].push_back({u, 1}); // Unit weight
        }
    }
    
    // Sort adjacency lists for consistent CSR if METIS prefers/requires
    for (idx_t i = 0; i < nvtxs; ++i) {
        std::sort(undirected_adj[i].begin(), undirected_adj[i].end());
    }


    idx_t current_edge_idx = 0;
    xadj[0] = 0;
    for (idx_t i = 0; i < nvtxs; ++i) {
        for (const auto& edge_pair : undirected_adj[i]) {
            adjncy.push_back(edge_pair.first); // Neighbor
            adjwgt.push_back(edge_pair.second); // Unit weight (1)
            current_edge_idx++;
        }
        xadj[i + 1] = current_edge_idx;
    }
}


// --- Main Function: Example Usage ---
int main() {
    // Example Graph (5 vertices)
    int num_vertices = 5;
    Graph graph(num_vertices);

    // Add edges: {to, weight}
    graph[0].push_back({1, 10.0});
    graph[0].push_back({3, 5.0});
    graph[1].push_back({2, 1.0});
    graph[1].push_back({3, 2.0});
    graph[2].push_back({4, 4.0});
    graph[3].push_back({1, 3.0}); // Note: 0->3 then 3->1 (original example)
    graph[3].push_back({2, 9.0});
    graph[3].push_back({4, 2.0});
    // graph[4].push_back({0, 7.0}); // Cycle back, removing for simpler acyclic example first
    // graph[4].push_back({2, 6.0});

    int start_node = 0;
    std::vector<double> dist;
    std::vector<int> pred;

    // 1. Calculate initial SSSP
    std::cout << "--- Initial SSSP Calculation ---" << std::endl;
    dijkstra(start_node, graph, dist, pred);
    print_paths(start_node, dist, pred);

    // --- METIS Partitioning ---
    std::cout << "\n--- METIS Graph Partitioning ---" << std::endl;
    idx_t metis_nvtxs;
    std::vector<idx_t> metis_xadj;
    std::vector<idx_t> metis_adjncy;
    std::vector<idx_t> metis_adjwgt; // Using unit weights

    convert_to_metis_csr(graph, metis_nvtxs, metis_xadj, metis_adjncy, metis_adjwgt);

    if (metis_nvtxs > 0) {
        idx_t ncon = 1; // Number of balancing constraints (typically 1 for vertex weights)
        idx_t nparts = 2; // Desired number of partitions (e.g., for 2 MPI processes later)
                         // Ensure nparts <= metis_nvtxs

        if (nparts > metis_nvtxs) {
            std::cout << "Warning: Number of partitions (" << nparts
                      << ") is greater than number of vertices (" << metis_nvtxs
                      << "). Setting nparts to " << metis_nvtxs << std::endl;
            nparts = metis_nvtxs;
        }
        if (nparts <=0 && metis_nvtxs > 0) {
             std::cout << "Warning: Invalid number of partitions. Setting to 1." << std::endl;
             nparts = 1;
        }


        std::vector<idx_t> part(metis_nvtxs); // Stores partition ID for each vertex
        idx_t objval; // Stores the edge-cut of the partitioning

        // METIS options
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // Minimize edge-cut
        options[METIS_OPTION_NUMBERING] = 0; // C-style 0-based numbering

        // If adjwgt is empty, pass NULL (METIS treats as unweighted)
        idx_t* adjwgt_ptr = metis_adjwgt.empty() ? NULL : metis_adjwgt.data();

        int ret = METIS_PartGraphKway(&metis_nvtxs, &ncon, metis_xadj.data(), metis_adjncy.data(),
                                      NULL, NULL, adjwgt_ptr, // vwgt, vsize, adjwgt
                                      &nparts, NULL, NULL,    // tpwgts (target partition weights), ubvec (load imbalance tolerance)
                                      options, &objval, part.data());

        if (ret == METIS_OK) {
            std::cout << "METIS partitioning successful." << std::endl;
            std::cout << "  Number of vertices: " << metis_nvtxs << std::endl;
            std::cout << "  Number of partitions: " << nparts << std::endl;
            std::cout << "  Edge-cut: " << objval << std::endl;
            for (idx_t i = 0; i < metis_nvtxs; ++i) {
                std::cout << "  Vertex " << i << " -> Partition " << part[i] << std::endl;
            }
        } else {
            std::cerr << "METIS_PartGraphKway failed with error code: " << ret << std::endl;
            if (ret == METIS_ERROR_INPUT) std::cerr << "  (Likely input error. Check CSR arrays, nvtxs, nparts)" << std::endl;
            if (ret == METIS_ERROR_MEMORY) std::cerr << "  (Memory allocation error within METIS)" << std::endl;
        }
    } else {
        std::cout << "Graph has no vertices, skipping METIS partitioning." << std::endl;
    }


    // --- Example SSSP Update (Optional, for testing update logic after partitioning info) ---
    std::cout << "\n--- Simulating Edge Weight Decrease (0 -> 1) from 10 to 1 ---" << std::endl;
    update_sssp(0, 1, 1.0, start_node, graph, dist, pred); // graph modified here
    print_paths(start_node, dist, pred);
    // Note: If graph is modified, METIS partitioning should ideally be re-done or
    // dynamic partitioning strategies employed for a fully dynamic system.
    // For now, we partition once at the beginning.

    return 0;
}