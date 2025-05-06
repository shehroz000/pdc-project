#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <metis.h>
#include <cstdlib> // For EXIT_FAILURE
#include <string>   // For std::string

// Function to read the graph from the custom file format and get the number of nodes
int read_graph_num_nodes(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    int num_nodes = 0;
    if (std::getline(file, line)) {
        std::stringstream ss_first(line);
        ss_first >> num_nodes;
    }
    return num_nodes;
}

// Function to read the graph from the custom file format and convert to CSR
void read_graph_and_convert_to_csr(const std::string& filename, std::vector<int>& xadj, std::vector<int>& adjncy, std::vector<int>& adjwgt) {
    std::ifstream file(filename);
    std::string line;
    int num_nodes = 0;

    // Read the first line: num_nodes num_edges source_node (source_node is not used here)
    if (std::getline(file, line)) {
        std::stringstream ss_first(line);
        ss_first >> num_nodes;
        int num_edges;
        ss_first >> num_edges;
        int temp_source_node;
        ss_first >> temp_source_node;
    }

    // Initialize CSR structures
    xadj.resize(num_nodes + 1, 0);
    std::vector<std::vector<std::pair<int, int>>> adj_list(num_nodes); // Temporary adjacency list

    // Read the rest of the lines to build the adjacency list
    for (int i = 0; i < num_nodes; ++i) {
        std::getline(file, line);
        std::stringstream ss(line);
        int node_id = i;

        int neighbor;
        int weight;
        while (ss >> neighbor >> weight) {
            adj_list[node_id].push_back({neighbor - 1, weight}); // Store neighbor and weight (0-based indexing)
            xadj[node_id + 1]++; // Increment the number of neighbors for this node
        }
    }

    // Calculate prefix sum for xadj (CSR row pointers)
    for (int i = 1; i <= num_nodes; ++i) {
        xadj[i] += xadj[i - 1];
    }

    // Fill adjncy and adjwgt from the adjacency list
    int total_edges = xadj[num_nodes];
    adjncy.resize(total_edges);
    adjwgt.resize(total_edges);
    int current_edge_index = 0;
    for (int i = 0; i < num_nodes; ++i) {
        for (const auto& neighbor : adj_list[i]) {
            adjncy[current_edge_index] = neighbor.first;
            adjwgt[current_edge_index] = neighbor.second;
            current_edge_index++;
        }
    }
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: ./metis_partition <graph_file>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string graph_file = argv[1];
    int num_nodes_global = read_graph_num_nodes(graph_file);
    std::vector<int> xadj, adjncy, adjwgt;

    // Read the graph from the file and convert it to CSR
    read_graph_and_convert_to_csr(graph_file, xadj, adjncy, adjwgt);

    // Partition the graph using METIS
    std::vector<int> part(num_nodes_global);
    idx_t ncon = 1;
    int num_parts = 2; // You can change this to the desired number of partitions
    int edgecut;

    idx_t metis_num_nodes = num_nodes_global;
    idx_t metis_ncon = ncon;
    int metis_result = METIS_PartGraphKway(&metis_num_nodes, &metis_ncon, &xadj[0], &adjncy[0], nullptr, nullptr, &adjwgt[0], &num_parts, nullptr, nullptr, nullptr, &edgecut, &part[0]);

    if (metis_result != METIS_OK) {
        std::cerr << "METIS partitioning failed with error: " << metis_result << std::endl;
        return EXIT_FAILURE;
    }

    // Create separate files for each partition
    std::vector<std::ofstream> partition_files(num_parts);
    for (int i = 0; i < num_parts; ++i) {
        partition_files[i].open("partition_" + std::to_string(i) + ".txt");
        if (!partition_files[i].is_open()) {
            std::cerr << "Error opening output file for partition " << i << std::endl;
            return EXIT_FAILURE;
        }
    }

    // Write the nodes belonging to each partition into their respective files
    for (int i = 0; i < num_nodes_global; ++i) {
        partition_files[part[i]] << i + 1 << std::endl; // Write the original node number (1-based)
    }

    // Close all the output files
    for (int i = 0; i < num_parts; ++i) {
        partition_files[i].close();
    }

    std::cout << "Graph partitioned into " << num_parts << " parts. Edge cut: " << edgecut << std::endl;
    std::cout << "Partition information saved to partition_0.txt, partition_1.txt, ..." << std::endl;

    return 0;
}