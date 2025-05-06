#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>

#include <metis.h>
#include <mpi.h>

#define INF INT_MAX
#define NO_PREDECESSOR -1

typedef struct {
    int to;
    double weight;
} Edge;

typedef struct {
    int target_vertex;
    double new_distance;
    int predecessor_vertex;
} UpdateMessage;

void dijkstra(int start_node, const Edge **graph, int num_vertices, double *dist, int *pred) {
    if (start_node < 0 || start_node >= num_vertices) {
        if (num_vertices > 0) {
            int world_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
            if (world_rank == 0) {
                fprintf(stderr, "Error: Start node %d is out of bounds for graph with %d vertices.\n", start_node, num_vertices);
            }
        } else {
            int world_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
            if (world_rank == 0) {
                fprintf(stderr, "Error: Graph is empty.\n");
            }
        }
        for (int i = 0; i < num_vertices; ++i) {
            dist[i] = INF;
            pred[i] = NO_PREDECESSOR;
        }
        return;
    }

    for (int i = 0; i < num_vertices; ++i) {
        dist[i] = INF;
        pred[i] = NO_PREDECESSOR;
    }

    if (num_vertices == 0) return;

    dist[start_node] = 0.0;
    int *pq = (int *)malloc(num_vertices * sizeof(int));
    double *pq_dist = (double *)malloc(num_vertices * sizeof(double));
    int pq_size = 0;

    pq[pq_size] = start_node;
    pq_dist[pq_size++] = 0.0;

    while (pq_size > 0) {
        double d = pq_dist[0];
        int u = pq[0];

        // Remove from priority queue (simple min extraction)
        for (int i = 0; i < pq_size - 1; ++i) {
            pq[i] = pq[i + 1];
            pq_dist[i] = pq_dist[i + 1];
        }
        pq_size--;

        if (d > dist[u]) {
            continue;
        }

        if (u < 0 || u >= num_vertices) continue;

        for (int i = 0; graph[u] != NULL && graph[u][i].to != -1; ++i) {
            Edge edge = graph[u][i];
            int v = edge.to;
            double weight = edge.weight;
            if (v < 0 || v >= num_vertices) continue;

            if (dist[u] != INF && dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                pred[v] = u;

                // Insert into priority queue (simple append)
                pq[pq_size] = v;
                pq_dist[pq_size++] = dist[v];
            }
        }

        // Simple priority queue re-ordering (bubble sort style)
        for (int i = 0; i < pq_size - 1; ++i) {
            for (int j = 0; j < pq_size - i - 1; ++j) {
                if (pq_dist[j] > pq_dist[j + 1]) {
                    // Swap
                    double temp_dist = pq_dist[j];
                    pq_dist[j] = pq_dist[j + 1];
                    pq_dist[j + 1] = temp_dist;

                    int temp_node = pq[j];
                    pq[j] = pq[j + 1];
                    pq[j + 1] = temp_node;
                }
            }
        }
    }

    free(pq);
    free(pq_dist);
}

void print_paths(int rank, int start_node, int num_vertices, const double *dist, const int *pred) {
    if (rank != 0) return;

    printf("Rank %d: Distances and paths from node %d:\n", rank, start_node);
    if (num_vertices == 0) {
        printf("  Graph is empty.\n");
        return;
    }
    for (int i = 0; i < num_vertices; ++i) {
        printf("  Vertex %d: Dist = ", i);
        if (i < num_vertices && dist[i] == INF) {
            printf("INF");
        } else if (i < num_vertices) {
            printf("%f", dist[i]);
        } else {
            printf("N/A (out of bounds)");
        }
        printf(", Pred = ");
        if (i < num_vertices) {
            printf("%d", pred[i]);
        } else {
            printf("N/A (out of bounds)");
        }
        printf("\n");
    }
}

void convert_to_metis_csr(const Edge ***directed_graph,
                          idx_t *nvtxs,
                          idx_t **xadj,
                          idx_t **adjncy,
                          idx_t **adjwgt) {
    *nvtxs = 0;
    if (directed_graph == NULL) return;
    while (directed_graph[*nvtxs] != NULL) (*nvtxs)++;

    *xadj = (idx_t *)malloc(((*nvtxs) + 1) * sizeof(idx_t));
    *adjncy = NULL;
    *adjwgt = NULL;
    idx_t total_edges = 0;

    // First pass to count edges (since C doesn't have dynamic vector resizing)
    for (idx_t i = 0; i < *nvtxs; ++i) {
        if (directed_graph[i] == NULL) continue;
        for (int j = 0; directed_graph[i][j].to != -1; ++j) {
            total_edges++;
        }
    }

    *adjncy = (idx_t *)malloc(2 * total_edges * sizeof(idx_t)); // *2 for undirected
    *adjwgt = (idx_t *)malloc(2 * total_edges * sizeof(idx_t));
    idx_t current_edge_idx = 0;
    (*xadj)[0] = 0;

    for (idx_t u = 0; u < *nvtxs; ++u) {
        if (directed_graph[u] == NULL) {
            (*xadj)[u + 1] = current_edge_idx;
            continue;
        }

        for (int j = 0; directed_graph[u][j].to != -1; ++j) {
            idx_t v = directed_graph[u][j].to;
            int weight = 1; // In the C++ code, all weights were 1 for METIS
            (*adjncy)[current_edge_idx] = v;
            (*adjwgt)[current_edge_idx] = weight;
            current_edge_idx++;

            // Add reverse edge for undirected graph
            (*adjncy)[current_edge_idx] = u;
            (*adjwgt)[current_edge_idx] = weight;
            current_edge_idx++;
        }
        (*xadj)[u + 1] = current_edge_idx;
    }
}

void free_csr(idx_t *xadj, idx_t *adjncy, idx_t *adjwgt) {
    if (xadj) free(xadj);
    if (adjncy) free(adjncy);
    if (adjwgt) free(adjwgt);
}

void update_sssp_mpi(
    int u_changed, int v_changed, double new_weight,
    int start_node,
    Edge ***graph,
    double *dist,
    int *pred,
    const idx_t *part,
    int world_rank, int world_size,
    MPI_Comm comm) {

    int num_vertices = 0;
    if (graph != NULL) {
        while (graph[num_vertices] != NULL) num_vertices++;
    }
    if (num_vertices == 0) return;

    bool change_applied_locally = false;
    if (u_changed >= 0 && u_changed < num_vertices) {
        if (new_weight == INF) {
            //Remove edge
            int old_size = 0;
            if (graph[u_changed] != NULL) {
                while (graph[u_changed][old_size].to != -1) old_size++;
            }
            Edge *new_edges = (Edge *)malloc((old_size) * sizeof(Edge)); // One less edge
            int new_index = 0;
            if (graph[u_changed] != NULL) {
                for (int i = 0; i < old_size; ++i) {
                    if (graph[u_changed][i].to != v_changed) {
                        new_edges[new_index++] = graph[u_changed][i];
                    } else {
                        change_applied_locally = true;
                    }
                }
                free(graph[u_changed]);
            }
            new_edges[new_index].to = -1; // Sentinel
            graph[u_changed] = new_edges;

        } else {
            // Update or add edge
            bool found = false;
            if (graph[u_changed] != NULL) {
                for (int i = 0; graph[u_changed][i].to != -1; ++i) {
                    if (graph[u_changed][i].to == v_changed) {
                        graph[u_changed][i].weight = new_weight;
                        found = true;
                        change_applied_locally = true;
                        break;
                    }
                }
            }
            if (!found) {
                int old_size = 0;
                if (graph[u_changed] != NULL) {
                    while (graph[u_changed][old_size].to != -1) old_size++;
                }
                Edge *new_edges = (Edge *)malloc((old_size + 2) * sizeof(Edge)); // One more edge + sentinel
                if (graph[u_changed] != NULL) {
                    for (int i = 0; i < old_size; ++i) {
                        new_edges[i] = graph[u_changed][i];
                    }
                    free(graph[u_changed]);
                }
                new_edges[old_size].to = v_changed;
                new_edges[old_size].weight = new_weight;
                new_edges[old_size + 1].to = -1; // Sentinel
                graph[u_changed] = new_edges;
                change_applied_locally = true;
            }
        }
    }

    int *local_pq = (int *)malloc(num_vertices * sizeof(int));
    double *local_pq_dist = (double *)malloc(num_vertices * sizeof(double));
    int local_pq_size = 0;

    if (u_changed >= 0 && u_changed < num_vertices &&
        v_changed >= 0 && v_changed < num_vertices) {
        if (dist[u_changed] != INF && dist[u_changed] + new_weight < dist[v_changed]) {
            if (part[v_changed] == world_rank) {
                dist[v_changed] = dist[u_changed] + new_weight;
                pred[v_changed] = u_changed;
                local_pq[local_pq_size] = v_changed;
                local_pq_dist[local_pq_size++] = dist[v_changed];
            }
        }
    }

    bool active_globally = true;
    int iter_count = 0;
    const int MAX_ITERATIONS = num_vertices * world_size * 2;

    UpdateMessage *outgoing_messages_map[world_size];
    for (int i = 0; i < world_size; ++i) {
        outgoing_messages_map[i] = (UpdateMessage *)malloc(num_vertices * sizeof(UpdateMessage)); // Max size
    }
    int outgoing_messages_count[world_size];
    for (int i = 0; i < world_size; ++i) outgoing_messages_count[i] = 0;


    MPI_Datatype mpi_update_message_type;
    int blocklengths[] = {1, 1, 1};
    MPI_Aint displacements[] = {offsetof(UpdateMessage, target_vertex), offsetof(UpdateMessage, new_distance), offsetof(UpdateMessage, predecessor_vertex)};
    MPI_Datatype types[] = {MPI_INT, MPI_DOUBLE, MPI_INT};
    MPI_Type_create_struct(3, blocklengths, displacements, types, &mpi_update_message_type);
    MPI_Type_commit(&mpi_update_message_type);

    bool forced_termination_signal = false;

    while (active_globally) {
        iter_count++;
        bool made_local_update_in_iteration = false;

        while (local_pq_size > 0) {
            int u = local_pq[0];
            double d_u = local_pq_dist[0];

            // Remove from local priority queue
            for (int i = 0; i < local_pq_size - 1; ++i) {
                local_pq[i] = local_pq[i + 1];
                local_pq_dist[i] = local_pq_dist[i + 1];
            }
            local_pq_size--;

            if (u < 0 || u >= num_vertices) continue;
            if (d_u > dist[u]) continue;

            if (graph[u] == NULL) continue;
            for (int i = 0; graph[u][i].to != -1; ++i) {
                int v = graph[u][i].to;
                double weight_uv = graph[u][i].weight;
                if (v < 0 || v >= num_vertices) continue;

                if (dist[u] != INF && dist[u] + weight_uv < dist[v]) {
                    dist[v] = dist[u] + weight_uv;
                    pred[v] = u;
                    made_local_update_in_iteration = true;

                    if (part[v] == world_rank) {
                        local_pq[local_pq_size] = v;
                        local_pq_dist[local_pq_size++] = dist[v];
                    } else if (part[v] < world_size && part[v] >= 0) {
                        outgoing_messages_map[part[v]][outgoing_messages_count[part[v]]].target_vertex = v;
                        outgoing_messages_map[part[v]][outgoing_messages_count[part[v]]].new_distance = dist[v];
                        outgoing_messages_map[part[v]][outgoing_messages_count[part[v]]].predecessor_vertex = u;
                        outgoing_messages_count[part[v]]++;
                    }
                }
            }
        }

        // Local priority queue re-ordering (bubble sort style)
        for (int i = 0; i < local_pq_size - 1; ++i) {
            for (int j = 0; j < local_pq_size - i - 1; ++j) {
                if (local_pq_dist[j] > local_pq_dist[j + 1]) {
                    // Swap
                    double temp_dist = local_pq_dist[j];
                    local_pq_dist[j] = local_pq_dist[j + 1];
                    local_pq_dist[j + 1] = temp_dist;

                    int temp_node = local_pq[j];
                    local_pq[j] = local_pq[j + 1];
                    local_pq[j + 1] = temp_node;
                }
            }
        }


        int send_counts[world_size];
        int recv_counts[world_size];
        int sdispls[world_size];
        int rdispls[world_size];

        UpdateMessage *send_buffer = (UpdateMessage *)malloc(num_vertices * world_size * sizeof(UpdateMessage)); // Max size
        UpdateMessage *recv_buffer = NULL;

        int current_sdispl = 0;
        for (int i = 0; i < world_size; ++i) {
            send_counts[i] = outgoing_messages_count[i];
            sdispls[i] = current_sdispl;
            for (int j = 0; j < outgoing_messages_count[i]; ++j) {
                send_buffer[current_sdispl + j] = outgoing_messages_map[i][j];
            }
            current_sdispl += outgoing_messages_count[i];
            outgoing_messages_count[i] = 0;
        }

        MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);

        int total_to_receive = 0;
        current_rdispl = 0;
        for (int i = 0; i < world_size; ++i) {
            rdispls[i] = current_rdispl;
            total_to_receive += recv_counts[i];
            current_rdispl += recv_counts[i];
        }

        if (total_to_receive > 0) {
            recv_buffer = (UpdateMessage *)malloc(total_to_receive * sizeof(UpdateMessage));
            MPI_Alltoallv(send_buffer, send_counts, sdispls, mpi_update_message_type,
                          recv_buffer, recv_counts, rdispls, mpi_update_message_type, comm);

            for (int i = 0; i < total_to_receive; ++i) {
                UpdateMessage msg = recv_buffer[i];
                if (msg.target_vertex < 0 || msg.target_vertex >= num_vertices) continue;

                if (msg.new_distance < dist[msg.target_vertex]) {
                    dist[msg.target_vertex] = msg.new_distance;
                    pred[msg.target_vertex] = msg.predecessor_vertex;
                    made_local_update_in_iteration = true;
                    if (part[msg.target_vertex] == world_rank) {
                        local_pq[local_pq_size] = msg.target_vertex;
                        local_pq_dist[local_pq_size++] = dist[msg.target_vertex];
                    }
                }
            }
            free(recv_buffer);
        }
        free(send_buffer);


        // Local priority queue re-ordering (bubble sort style)
        for (int i = 0; i < local_pq_size - 1; ++i) {
            for (int j = 0; j < local_pq_size - i - 1; ++j) {
                if (local_pq_dist[j] > local_pq_dist[j + 1]) {
                    // Swap
                    double temp_dist = local_pq_dist[j];
                    local_pq_dist[j] = local_pq_dist[j + 1];
                    local_pq_dist[j + 1] = temp_dist;

                    int temp_node = local_pq[j];
                    local_pq[j] = local_pq[j + 1];
                    local_pq[j + 1] = temp_node;
                }
            }
        }


        bool current_process_active = made_local_update_in_iteration || local_pq_size > 0;
        bool any_process_active = false;
        MPI_Allreduce(&current_process_active, &any_process_active, 1, MPI_C_BOOL, MPI_LOR, comm);
        active_globally = any_process_active;

        if (world_rank == 0 && iter_count > MAX_ITERATIONS && active_globally) {
            printf("Rank 0: Warning - exceeded max iterations (%d). Forcing termination.\n", iter_count);
            forced_termination_signal = true;
        }
        MPI_Bcast(&forced_termination_signal, 1, MPI_C_BOOL, 0, comm);
        if (forced_termination_signal) {
            active_globally = false;
        }
    }

    MPI_Type_free(&mpi_update_message_type);
    free(local_pq);
    free(local_pq_dist);
    for (int i = 0; i < world_size; ++i) {
        free(outgoing_messages_map[i]);
    }
    if (world_rank == 0) printf("Rank 0: MPI Update finished after %d iterations.\n", iter_count);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int num_vertices = 0;
    Edge ***graph = NULL;
    double *dist = NULL;
    int *pred = NULL;
    idx_t *part_assignment = NULL;
    int start_node = 0;

    if (world_rank == 0) {
        printf("MPI SSSP Update running with %d processes.\n", world_size);

        // Example graph initialization (replace with your graph loading)
        num_vertices = 5;
        graph = (Edge ***)malloc((num_vertices + 1) * sizeof(Edge **));
        for (int i = 0; i < num_vertices; ++i) {
            graph[i] = (Edge **)malloc(5 * sizeof(Edge *)); // Assume max 4 edges per vertex + sentinel
            for (int j = 0; j < 5; ++j) {
                graph[i][j] = (Edge *)malloc(sizeof(Edge));
                graph[i][j]->to = -1; // Sentinel value
            }
        }

        // Initialize the graph
        graph[0][0]->to = 1; graph[0][0]->weight = 10.0;
        graph[0][1]->to = 3; graph[0][1]->weight = 5.0;
        graph[0][2]->to = -1;

        graph[1][0]->to = 2; graph[1][0]->weight = 1.0;
        graph[1][1]->to = 3; graph[1][1]->weight = 2.0;
        graph[1][2]->to = -1;

        graph[2][0]->to = 4; graph[2][0]->weight = 4.0;
        graph[2][1]->to = -1;

        graph[3][0]->to = 1; graph[3][0]->weight = 3.0;
        graph[3][1]->to = 2; graph[3][1]->weight = 9.0;
        graph[3][2]->to = 4; graph[3][2]->weight = 2.0;
        graph[3][3]->to = -1;

        graph[4][0]->to = 0; graph[4][0]->weight = 7.0;
        graph[4][1]->to = 2; graph[4][1]->weight = 6.0;
        graph[4][2]->to = -1;


        dist = (double *)malloc(num_vertices * sizeof(double));
        pred = (int *)malloc(num_vertices * sizeof(int));
        part_assignment = (idx_t *)malloc(num_vertices * sizeof(idx_t));

        // Example partitioning (replace with METIS)
        for (int i = 0; i < num_vertices; ++i) {
            part_assignment[i] = i % world_size;
        }
    }

    MPI_Bcast(&num_vertices, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (world_rank != 0) {
        graph = (Edge ***)malloc((num_vertices + 1) * sizeof(Edge **));
        for (int i = 0; i < num_vertices; ++i) {
            graph[i] = (Edge **)malloc(5 * sizeof(Edge *));
            for (int j = 0; j < 5; ++j) {
                graph[i][j] = (Edge *)malloc(sizeof(Edge));
                graph[i][j]->to = -1;
            }
        }
        dist = (double *)malloc(num_vertices * sizeof(double));
        pred = (int *)malloc(num_vertices * sizeof(int));
        part_assignment = (idx_t *)malloc(num_vertices * sizeof(idx_t));
    }

    // Broadcast graph, dist, pred, part_assignment
    for (int i = 0; i < num_vertices; ++i) {
        for (int j = 0; j < 5; ++j) {
            MPI_Bcast(graph[i][j], sizeof(Edge), MPI_BYTE, 0, MPI_COMM_WORLD);
        }
    }
    MPI_Bcast(dist, num_vertices, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(pred, num_vertices, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(part_assignment, num_vertices, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);


    dijkstra(start_node, (const Edge **)graph, num_vertices, dist, pred);
    print_paths(world_rank, start_node, num_vertices, dist, pred);

    // Example update (replace with your update logic)
    if (world_rank == 0) {
        printf("Performing update: Removing edge 0->1, Adding edge 2->1 with weight 2.0\n");
        update_sssp_mpi(0, 1, INF, start_node, graph, dist, pred, part_assignment, world_rank, world_size, MPI_COMM_WORLD);
        update_sssp_mpi(2, 1, 2.0, start_node, graph, dist, pred, part_assignment, world_rank, world_size, MPI_COMM_WORLD);
    } else {
        update_sssp_mpi(-1, -1, 0.0, start_node, graph, dist, pred, part_assignment, world_rank, world_size, MPI_COMM_WORLD);
    }

    dijkstra(start_node, (const Edge **)graph, num_vertices, dist, pred);
    print_paths(world_rank, start_node, num_vertices, dist, pred);


    // Cleanup
    if (graph != NULL) {
        for (int i = 0; i < num_vertices; ++i) {
            if (graph[i] != NULL) {
                for (int j = 0; j < 5; ++j) {
                    free(graph[i][j]);
                }
                free(graph[i]);
            }
        }
        free(graph);
    }
    if (dist) free(dist);
    if (pred) free(pred);
    if (part_assignment) free(part_assignment);


    MPI_Finalize();
    return 0;
}