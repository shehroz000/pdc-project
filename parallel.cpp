#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <string>
#include <algorithm> // For std::remove_if, std::find_if
#include <cmath>     // For std::abs
#include <map>       // For message buffers
#include <vector>    // Ensure vector is included
#include <stddef.h>  // Required for offsetof

// --- METIS Header ---
#include <metis.h>
// --- MPI Header ---
#include <mpi.h>


const double INF = std::numeric_limits<double>::infinity();
const int NO_PREDECESSOR = -1;

// Represents an edge: target vertex and weight (for SSSP)
struct Edge {
    int to;
    double weight;
    // Add a default constructor for vector resize
    Edge() : to(0), weight(0.0) {}
    Edge(int t, double w) : to(t), weight(w) {}
};
using Graph = std::vector<std::vector<Edge>>;
using PrioQueueEntry = std::pair<double, int>;


// --- SSSP Functions (dijkstra, print_paths) ---
void dijkstra(int start_node, const Graph& graph, std::vector<double>& dist, std::vector<int>& pred) {
    int num_vertices = graph.size();
    if (start_node < 0 || start_node >= num_vertices) { 
        if (num_vertices > 0) { 
             if (MPI::COMM_WORLD.Get_rank() == 0) { // Only rank 0 prints errors for clarity
                std::cerr << "Error: Start node " << start_node << " is out of bounds for graph with " << num_vertices << " vertices." << std::endl;
             }
        } else {
            if (MPI::COMM_WORLD.Get_rank() == 0) {
                std::cerr << "Error: Graph is empty." << std::endl;
            }
        }
        dist.assign(num_vertices, INF);
        pred.assign(num_vertices, NO_PREDECESSOR);
        return;
    }
    dist.assign(num_vertices, INF);
    pred.assign(num_vertices, NO_PREDECESSOR);

    if (num_vertices == 0) return; 

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
        if (u < 0 || u >= graph.size()) continue; // Bounds check for u

        for (const Edge& edge : graph[u]) {
            int v = edge.to;
            double weight = edge.weight;
            if (v < 0 || v >= num_vertices) continue; 

            if (dist[u] != INF && dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                pred[v] = u;
                pq.push({dist[v], v});
            }
        }
    }
}

void print_paths(int rank, int start_node, int num_vertices, const std::vector<double>& dist, const std::vector<int>& pred) {
    if (rank != 0) return; 

    std::cout << "Rank " << rank << ": Distances and paths from node " << start_node << ":" << std::endl;
    if (num_vertices == 0) {
        std::cout << "  Graph is empty." << std::endl;
        return;
    }
    for (int i = 0; i < num_vertices; ++i) {
        std::cout << "  Vertex " << i << ": Dist = ";
        if (i < dist.size() && dist[i] == INF) {
            std::cout << "INF";
        } else if (i < dist.size()) {
            std::cout << dist[i];
        } else {
            std::cout << "N/A (out of bounds)";
        }
        std::cout << ", Pred = ";
        if (i < pred.size()) {
            std::cout << pred[i];
        } else {
            std::cout << "N/A (out of bounds)";
        }
        std::cout << std::endl;
    }
}


// --- METIS Integration Functions ---
void convert_to_metis_csr(const Graph& directed_graph,
                          idx_t& nvtxs,
                          std::vector<idx_t>& xadj,
                          std::vector<idx_t>& adjncy,
                          std::vector<idx_t>& adjwgt) {
    nvtxs = directed_graph.size();
    xadj.clear();
    adjncy.clear();
    adjwgt.clear();

    if (nvtxs == 0) return;
    xadj.resize(nvtxs + 1);

    std::vector<std::vector<std::pair<int, int>>> undirected_adj(nvtxs);
    for (idx_t u = 0; u < nvtxs; ++u) {
        if (u >= directed_graph.size()) continue; // Bounds check
        for (const auto& edge : directed_graph[u]) {
            idx_t v = edge.to;
            if (v < 0 || v >= nvtxs) continue; 

            bool found_uv = false;
            if (u < undirected_adj.size()) { // Bounds check
                for(const auto& p : undirected_adj[u]) if(p.first == v) found_uv = true;
                if(!found_uv) undirected_adj[u].push_back({(int)v, 1});
            }


            bool found_vu = false;
            if (v < undirected_adj.size()) { // Bounds check
                for(const auto& p : undirected_adj[v]) if(p.first == u) found_vu = true;
                if(!found_vu) undirected_adj[v].push_back({(int)u, 1});
            }
        }
    }
    
    for (idx_t i = 0; i < nvtxs; ++i) {
        if (i < undirected_adj.size()) // Bounds check
            std::sort(undirected_adj[i].begin(), undirected_adj[i].end());
    }

    idx_t current_edge_idx = 0;
    xadj[0] = 0;
    for (idx_t i = 0; i < nvtxs; ++i) {
        if (i < undirected_adj.size()) { // Bounds check
            for (const auto& edge_pair : undirected_adj[i]) {
                adjncy.push_back(edge_pair.first); 
                adjwgt.push_back(edge_pair.second); 
                current_edge_idx++;
            }
        }
        xadj[i + 1] = current_edge_idx;
    }
}


// --- MPI Specific Data Structures & Functions ---
struct UpdateMessage {
    int target_vertex;
    double new_distance;
    int predecessor_vertex;
};


void update_sssp_mpi(
    int u_changed, int v_changed, double new_weight, 
    int start_node,
    Graph& graph,                     
    std::vector<double>& dist,        
    std::vector<int>& pred,           
    const std::vector<idx_t>& part,   
    int world_rank, int world_size,
    MPI_Comm comm)
{
    int num_vertices = graph.size();
    if (num_vertices == 0) return;

    bool change_applied_locally = false;
    if (u_changed >= 0 && u_changed < num_vertices) {
        if (new_weight == INF) { // Deletion
            auto& u_edges = graph[u_changed];
            auto old_size = u_edges.size();
            u_edges.erase(
                std::remove_if(u_edges.begin(), u_edges.end(),
                               [v_changed](const Edge& e){ return e.to == v_changed; }),
                u_edges.end()
            );
            if (u_edges.size() < old_size) change_applied_locally = true;
        } else { // Update or Insertion
            bool found = false;
            for(auto& edge : graph[u_changed]) {
                if (edge.to == v_changed) {
                    edge.weight = new_weight;
                    found = true;
                    change_applied_locally = true;
                    break;
                }
            }
            if (!found) { 
                graph[u_changed].push_back({v_changed, new_weight});
                change_applied_locally = true;
            }
        }
    }


    std::priority_queue<PrioQueueEntry, std::vector<PrioQueueEntry>, std::greater<PrioQueueEntry>> local_pq;

    // Initial push to PQ: This logic needs to be robust and based on the paper's method
    // for identifying affected nodes. Here's a basic start:
    if (u_changed >= 0 && u_changed < num_vertices &&
        v_changed >= 0 && v_changed < num_vertices) {
        if (dist[u_changed] != INF && dist[u_changed] + new_weight < dist[v_changed]) {
            if (part[v_changed] == world_rank) {
                dist[v_changed] = dist[u_changed] + new_weight;
                pred[v_changed] = u_changed;
                local_pq.push({dist[v_changed], v_changed});
            }
            // If v_changed is remote, its owner will get this via message propagation
        }
        // More comprehensive: if the original edge (u,v) was part of SSSP and is now worse/deleted,
        // the owner of v might need to re-evaluate v.
        // This initial seeding is critical and paper-dependent.
    }


    bool active_globally = true;
    int iter_count = 0;
    const int MAX_ITERATIONS = num_vertices * world_size * 2; // Adjusted safety break

    std::vector<std::vector<UpdateMessage>> outgoing_messages_map(world_size); // Using map for clarity

    // Define MPI Datatype for UpdateMessage (do this once)
    MPI_Datatype mpi_update_message_type;
    int blocklengths[] = {1, 1, 1}; 
    MPI_Aint displacements[] = {offsetof(UpdateMessage, target_vertex), offsetof(UpdateMessage, new_distance), offsetof(UpdateMessage, predecessor_vertex)};
    MPI_Datatype types[] = {MPI_INT, MPI_DOUBLE, MPI_INT};
    MPI_Type_create_struct(3, blocklengths, displacements, types, &mpi_update_message_type);
    MPI_Type_commit(&mpi_update_message_type);

    bool forced_termination_signal = false;

    while(active_globally) {
        iter_count++;
        bool made_local_update_in_iteration = false;

        while(!local_pq.empty()) {
            int u = local_pq.top().second;
            double d_u = local_pq.top().first;
            local_pq.pop();

            if (u < 0 || u >= num_vertices) continue; 
            if (d_u > dist[u]) continue; 

            if (u >= graph.size()) continue; // Check u against actual graph size
            for (const auto& edge : graph[u]) {
                int v = edge.to;
                double weight_uv = edge.weight;
                if (v < 0 || v >= num_vertices) continue; 

                if (dist[u] != INF && dist[u] + weight_uv < dist[v]) {
                    dist[v] = dist[u] + weight_uv; 
                    pred[v] = u;
                    made_local_update_in_iteration = true;

                    if (v >= part.size()) { /* error or resize part */ continue;}
                    if (part[v] == world_rank) { 
                        local_pq.push({dist[v], v});
                    } else if (part[v] < world_size && part[v] >= 0) { // Check valid rank
                        outgoing_messages_map[part[v]].push_back({v, dist[v], u});
                    }
                }
            }
        }

        std::vector<int> send_counts(world_size, 0);
        std::vector<int> recv_counts(world_size, 0);
        std::vector<int> sdispls(world_size, 0);
        std::vector<int> rdispls(world_size, 0);
        std::vector<UpdateMessage> send_buffer;
        
        int current_sdispl = 0;
        for(int i=0; i<world_size; ++i) {
            send_counts[i] = outgoing_messages_map[i].size();
            sdispls[i] = current_sdispl;
            send_buffer.insert(send_buffer.end(), outgoing_messages_map[i].begin(), outgoing_messages_map[i].end());
            current_sdispl += outgoing_messages_map[i].size();
            outgoing_messages_map[i].clear(); 
        }

        MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm);

        std::vector<UpdateMessage> recv_buffer;
        int current_rdispl = 0;
        int total_to_receive = 0;
        for(int i=0; i<world_size; ++i) {
            rdispls[i] = current_rdispl;
            current_rdispl += recv_counts[i];
            total_to_receive += recv_counts[i];
        }
        recv_buffer.resize(total_to_receive);
        
        MPI_Alltoallv(send_buffer.data(), send_counts.data(), sdispls.data(), mpi_update_message_type,
                      recv_buffer.data(), recv_counts.data(), rdispls.data(), mpi_update_message_type, comm);

        for (const auto& msg : recv_buffer) {
            if (msg.target_vertex < 0 || msg.target_vertex >= num_vertices) continue;

            if (msg.new_distance < dist[msg.target_vertex]) {
                dist[msg.target_vertex] = msg.new_distance;
                pred[msg.target_vertex] = msg.predecessor_vertex;
                made_local_update_in_iteration = true;
                if (msg.target_vertex >= part.size()) { /* error */ continue; }
                if (part[msg.target_vertex] == world_rank) { // Should always be true
                    local_pq.push({dist[msg.target_vertex], msg.target_vertex});
                }
            }
        }
        
        bool current_process_active = made_local_update_in_iteration || !local_pq.empty();
        bool any_process_active = false;
        MPI_Allreduce(&current_process_active, &any_process_active, 1, MPI_CXX_BOOL, MPI_LOR, comm);
        active_globally = any_process_active;

        // Check for forced termination signal from rank 0
        if (world_rank == 0 && iter_count > MAX_ITERATIONS && active_globally) {
            if (world_rank == 0) std::cout << "Rank 0: Warning - exceeded max iterations (" << iter_count << "). Forcing termination." << std::endl;
            forced_termination_signal = true;
        }
        MPI_Bcast(&forced_termination_signal, 1, MPI_CXX_BOOL, 0, comm);
        if (forced_termination_signal) {
            active_globally = false;
        }
    }
    MPI_Type_free(&mpi_update_message_type);
    if (world_rank == 0) std::cout << "Rank 0: MPI Update finished after " << iter_count << " iterations." << std::endl;
}


// --- Main Function ---
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int num_vertices = 0;
    Graph graph; 
    std::vector<double> dist;
    std::vector<int> pred;
    std::vector<idx_t> part_assignment; // Changed name for clarity
    int start_node = 0;


    if (world_rank == 0) {
        std::cout << "MPI SSSP Update running with " << world_size << " processes." << std::endl;
        num_vertices = 5; 
        graph.resize(num_vertices);
        graph[0].push_back({1, 10.0}); graph[0].push_back({3, 5.0});
        graph[1].push_back({2, 1.0}); graph[1].push_back({3, 2.0});
        graph[2].push_back({4, 4.0});
        graph[3].push_back({1, 3.0}); graph[3].push_back({2, 9.0}); graph[3].push_back({4, 2.0});

        idx_t metis_nvtxs_check; // Use a different name to avoid conflict if idx_t is not int
        std::vector<idx_t> metis_xadj, metis_adjncy, metis_adjwgt;
        convert_to_metis_csr(graph, metis_nvtxs_check, metis_xadj, metis_adjncy, metis_adjwgt);
        
        part_assignment.resize(num_vertices);
        if (metis_nvtxs_check > 0) {
            idx_t ncon = 1;
            idx_t nparts_metis = world_size; // Use a different name for metis nparts
            if (nparts_metis > metis_nvtxs_check) nparts_metis = metis_nvtxs_check;
            if (nparts_metis <= 0 && metis_nvtxs_check > 0 ) nparts_metis = 1;
            else if (nparts_metis <=0 && metis_nvtxs_check == 0) nparts_metis =0;


            idx_t objval;
            idx_t options[METIS_NOPTIONS];
            METIS_SetDefaultOptions(options);
            options[METIS_OPTION_NUMBERING] = 0;
            idx_t* adjwgt_ptr = metis_adjwgt.empty() ? NULL : metis_adjwgt.data();

            if (nparts_metis > 0) { // Only call METIS if there's something to partition
                int ret = METIS_PartGraphKway(&metis_nvtxs_check, &ncon, metis_xadj.data(), metis_adjncy.data(),
                                            NULL, NULL, adjwgt_ptr, &nparts_metis, NULL, NULL,
                                            options, &objval, part_assignment.data());
                if (ret != METIS_OK) {
                    std::cerr << "Rank 0: METIS failed! Error: " << ret << std::endl;
                    for(size_t i=0; i < part_assignment.size(); ++i) part_assignment[i] = 0; // Default
                } else {
                    std::cout << "Rank 0: METIS Partitioning Complete. Edge-cut: " << objval << std::endl;
                    for(int i=0; i<num_vertices; ++i) std::cout << "  Vertex " << i << " -> Partition " << part_assignment[i] << std::endl;
                }
            } else {
                 for(size_t i=0; i < part_assignment.size(); ++i) part_assignment[i] = i % world_size; // Simple distribution if METIS not called
            }

        } else { // metis_nvtxs_check is 0 or graph is empty
            for(size_t i=0; i < part_assignment.size(); ++i) part_assignment[i] = i % world_size; 
        }
        
        dist.resize(num_vertices);
        pred.resize(num_vertices);
        dijkstra(start_node, graph, dist, pred);
        std::cout << "Rank 0: Initial SSSP calculated." << std::endl;
        print_paths(world_rank, start_node, num_vertices, dist, pred); 
    }

    MPI_Bcast(&num_vertices, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (world_rank != 0) { 
        graph.resize(num_vertices); 
        for(int i=0; i < num_vertices; ++i) graph[i].clear(); // Clear before receiving
        part_assignment.resize(num_vertices);
        dist.resize(num_vertices);
        pred.resize(num_vertices);
    }

    for (int i = 0; i < num_vertices; ++i) {
        int num_edges_for_i;
        if (world_rank == 0) {
            num_edges_for_i = (i < graph.size()) ? graph[i].size() : 0;
        }
        MPI_Bcast(&num_edges_for_i, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (world_rank != 0) {
            if (i < graph.size()) graph[i].resize(num_edges_for_i);
        }
        for (int j = 0; j < num_edges_for_i; ++j) {
             if (i < graph.size() && j < graph[i].size()) { // Bounds check for safety
                MPI_Bcast(&(graph[i][j].to), 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Bcast(&(graph[i][j].weight), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
        }
    }

    // Assuming idx_t is compatible with int. If idx_t is long, use MPI_LONG.
    // Check your metis.h for how idx_t is defined if MPI_INT causes issues.
    MPI_Bcast(part_assignment.data(), num_vertices, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast(dist.data(), num_vertices, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(pred.data(), num_vertices, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD); 

    if (world_rank == 0) {
        std::cout << "\n--- Simulating Parallel Edge Weight Decrease (0 -> 1) from 10 to 1 ---" << std::endl;
    }
    int u_update = 0;
    int v_update = 1;
    double new_w = 1.0;
    
    update_sssp_mpi(u_update, v_update, new_w, start_node, graph, dist, pred, part_assignment, world_rank, world_size, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD); 

    if (world_rank == 0) {
        std::cout << "\n--- Updated SSSP Results (Rank 0) ---" << std::endl;
        print_paths(world_rank, start_node, num_vertices, dist, pred);
    }
    
    if (world_rank == 0) {
        std::cout << "\n--- Simulating Parallel Edge Weight Increase (3 -> 4) from 2 to 10 ---" << std::endl;
    }
    u_update = 3; v_update = 4; new_w = 10.0; 
    
    update_sssp_mpi(u_update, v_update, new_w, start_node, graph, dist, pred, part_assignment, world_rank, world_size, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank == 0) {
        std::cout << "\n--- Updated SSSP Results after Increase (Rank 0) ---" << std::endl;
        print_paths(world_rank, start_node, num_vertices, dist, pred);
    }

    MPI_Finalize();
    return 0;
}