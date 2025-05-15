#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <sstream>
#include <limits>
#include <string>
#include <metis.h>
#include<mpi.h>
#include <unordered_map>
#include <omp.h>
#include <algorithm>

using namespace std;

const int INF = numeric_limits<int>::max();


using AffectedNodesMap = std::unordered_map<int, int>;

// Function to update the partition file
void update_partition_file(int rank, const std::pair<int, int>& affected_node) {
    //cout<<"hello from func";
    int affected_id = affected_node.first;
    int new_dist = affected_node.second;
    //cout<<"affected vertex: "<<affected_id<<"  and its new dist: "<<new_dist<<"\n\n";
    std::string partition_filename = "partition_" + std::to_string(rank) + ".txt";
    //cout<<rank<<endl;
    std::ifstream in_partition(partition_filename);
    if (!in_partition) {
        std::cerr << "Failed to open file: " << partition_filename << std::endl;
        return;
    }

    std::vector<std::string> updated_lines;
    std::string line;

    while (std::getline(in_partition, line)) {
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;

        while (iss >> token) {
            tokens.push_back(token);
        }

        if (tokens.empty()) continue;

        int source = std::stoi(tokens[0]);
        int source_dist = std::stoi(tokens[1]);


        for (size_t i = 2; i + 1 <= tokens.size(); i += 2) {
            int neighbor = std::stoi(tokens[i]);
            int weight = std::stoi(tokens[i + 1]);

            //if(neighbor==10)cout<<"10 found10 found10 found10 found while affected_id is: "<<affected_id<<"\n\n"; 
            //cout<<"neighbour: "<<neighbor;
            if (neighbor == affected_id) {
                //cout<<"neighbour Matched"<<neighbor<<"\n\n";
                int path_dist = source_dist + weight;
                //cout<<"old path_dist= "<<path_dist<<"\n\n";
                //cout<<"new: "<<new_dist<<"\n\n";

                if (path_dist > new_dist) {
                    
                    //cout<<"Removing pair and using: "<<new_dist<<"\n\n";
                    // Remove the pair <neighbor> <weight>
                    tokens.erase(tokens.begin() + i, tokens.begin() + i + 2);
                    i -= 2; // Adjust index
                }
            }
        }

        // Reconstruct line
        std::ostringstream oss;
        for (size_t i = 0; i < tokens.size(); ++i) {
            if (i > 0) oss << " ";
            oss << tokens[i];
        }

        updated_lines.push_back(oss.str());
    }

    in_partition.close();

    // CRITICAL: Overwrite the file with the updated content
   
    std::ofstream out_partition(partition_filename);
    for (const auto& updated_line : updated_lines) {
        out_partition << updated_line << "\n";
    }
    out_partition.close();


}






// Structure for message to send between ranks
struct Message {
    int node_id;
    int tentative_distance;
};


struct Edge {
    int to;
    int weight;
};

void insertEdgeInOutput(const string& filename, int src, int dest, int weight) {
    // Open the file to read and modify the output
    ifstream infile(filename);
    string line;
    vector<string> lines;
    
    // Read all lines into memory
    while (getline(infile, line)) {
        lines.push_back(line);
    }
    
    infile.close();
    
    // Parse the lines into a structure
    vector<vector<int>> graph;
    vector<int> nodes;
    
    for (const string& l : lines) {
        stringstream ss(l);
        int node;
        vector<int> row;
        
        while (ss >> node) {
            row.push_back(node);
        }
        graph.push_back(row);
        nodes.push_back(row[0]); // First number is the node identifier
    }
    
    // Check if the src node exists in the graph
    bool srcExists = false;
    bool destExists = false;
    int srcIdx = -1, destIdx = -1;
    
    // Find the indices of src and dest nodes if they exist
    for (int i = 0; i < graph.size(); ++i) {
        if (graph[i][0] == src) {
            srcExists = true;
            srcIdx = i;
        }
        if (graph[i][0] == dest) {
            destExists = true;
            destIdx = i;
        }
    }
    
    // If the src node exists, add the new edge to its row
    if (srcExists) {
        graph[srcIdx].push_back(dest);
        graph[srcIdx].push_back(weight);
    } else {
        // If src doesn't exist, add a new row for src with only the src node
        graph.push_back({src, -1, dest, weight}); // -1 signifies that no connections exist yet
    }
    
    // If the destination node does not exist, add it to the graph
    if (!destExists) {
        graph.push_back({dest, -1, src, weight});
    }

    // Write the updated graph to output.txt
    ofstream outfile("inserted.txt");
    for (const auto& row : graph) {
        for (size_t i = 0; i < row.size(); ++i) {
            outfile << row[i];
            if (i < row.size() - 1) {
                outfile << " ";
            }
        }
        outfile << endl;
    }
    outfile.close();
}



void runDijkstra(const string& filename) {
    ifstream infile(filename);
    ofstream outfile("output.txt");

    if (!infile || !outfile) {
        cerr << "Error opening file for Dijkstra.\n";
        return;
    }

    int V, E, start;
    infile >> V >> E >> start;

    vector<vector<Edge>> adj(V + 1);
    string line;
    getline(infile, line); // consume rest of the first line

    for (int i = 1; i <= V; ++i) {
        getline(infile, line);
        stringstream ss(line);
        int neighbor, weight;
        while (ss >> neighbor >> weight) {
            adj[i].push_back({neighbor, weight});
        }
    }

    vector<int> dist(V + 1, INF);
    dist[start] = 0;
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;
    pq.push({0, start});

    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dist[u]) continue;
        for (auto &edge : adj[u]) {
            int v = edge.to, w = edge.weight;
            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                pq.push({dist[v], v});
            }
        }
    }

    for (int i = 1; i <= V; ++i) {
        outfile << i << " ";
        if (dist[i] == INF) outfile << -1;
        else outfile << dist[i];
        for (const auto &edge : adj[i]) {
            outfile << " " << edge.to << " " << edge.weight;
        }
        outfile << "\n";
    }

    infile.close();
    outfile.close();
    cout << "Dijkstra output written to output.txt\n";
}

void runMetis(const string& inputFile, int nparts = 4) {
    ifstream infile(inputFile);
    if (!infile.is_open()) {
        cerr << "Error: Cannot open input file for METIS.\n";
        return;
    }

    idx_t nvtxs, nedges, startVertex;
    infile >> nvtxs >> nedges >> startVertex;
    string dummy;
    getline(infile, dummy); // skip rest of first line

    vector<idx_t> xadj = {0};
    vector<idx_t> adjncy;
    vector<idx_t> adjwgt;

    string line;
    idx_t edgeCount = 0;
    while (getline(infile, line)) {
        istringstream iss(line);
        idx_t neighbor, weight;
        while (iss >> neighbor >> weight) {
            adjncy.push_back(neighbor - 1);
            adjwgt.push_back(weight);
            edgeCount++;
        }
        xadj.push_back(edgeCount);
    }
    infile.close();

    idx_t ncon = 1;
    idx_t objval;
    vector<idx_t> part(nvtxs);

    idx_t* vwgt = NULL;
    idx_t* vsize = NULL;
    real_t* tpwgts = NULL;
    real_t* ubvec = NULL;

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;

    int status = METIS_PartGraphKway(
        &nvtxs, &ncon, xadj.data(), adjncy.data(),
        vwgt, vsize, adjwgt.data(),
        &nparts, tpwgts, ubvec,
        options, &objval, part.data()
    );

    if (status != METIS_OK) {
        cerr << "METIS_PartGraphKway failed.\n";
        return;
    }

    string partFileName = inputFile + ".part." + to_string(nparts);
    ofstream outfile(partFileName);
    if (!outfile.is_open()) {
        cerr << "Error writing to " << partFileName << "\n";
        return;
    }

    for (idx_t i = 0; i < nvtxs; ++i) {
        outfile << part[i] << "\n";
    }
    outfile.close();
    cout << "METIS partition written to " << partFileName << "\n";
}

// ---------------------------
// Integration Step
// ---------------------------
void writePartitionFiles(const string& metisPartitionFile, const string& dijkstraOutputFile, int numPartitions) {
    static int callCount = 0;
    callCount++;
    //cout << "[DEBUG] writePartitionFiles() called " << callCount << " time(s)\n";

    ifstream metisIn(metisPartitionFile);
    ifstream dijkstraIn(dijkstraOutputFile);

    if (!metisIn.is_open() || !dijkstraIn.is_open()) {
        cerr << "Error: Cannot open metis or dijkstra output file.\n";
        return;
    }

    // Open output files for each partition
    vector<ofstream> outFiles(numPartitions);
    for (int i = 0; i < numPartitions; ++i) {
        outFiles[i].open("partition_" + to_string(i) + ".txt", std::ios::out | std::ios::trunc);

        if (!outFiles[i].is_open()) {
            cerr << "Error creating partition_" << i << ".txt\n";
            return;
        }
    }

    string line;
    int lineNumber = 0;
    vector<int> writeCounts(numPartitions, 0); // Track lines written to each file

    while (getline(dijkstraIn, line)) {
        string partStr;
        if (!getline(metisIn, partStr)) {
            cerr << "Mismatch between Metis and Dijkstra lines!\n";
            return;
        }

        int partition = stoi(partStr);
        if (partition < 0 || partition >= numPartitions) {
            cerr << "Invalid partition number: " << partition << "\n";
            return;
        }

        outFiles[partition] << line << "\n";
        writeCounts[partition]++;
        lineNumber++;
    }

    // Close all files
    for (int i = 0; i < numPartitions; ++i) {
        outFiles[i].close();
        //cout << "partition_" << i << ".txt was written with " << writeCounts[i] << " lines.\n";
    }

    //cout << "Total lines processed: " << lineNumber << "\n";
    //cout << "Partitioned output written to partition_*.txt files.\n";
}


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    double global_start = MPI_Wtime();  // Start timing


    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int omp_threads = omp_get_max_threads();  // Total threads used by OpenMP

    string filename = "t1.txt";
    int nparts = 2;
    int numVertices=7;
    int half=numVertices/2;
    int temp=0;
std::vector<std::pair<int, int>> copied_nodes;

if (rank == 0) {
    for (int i = 0; i < nparts; ++i) {
        string fname = "partition_" + to_string(i) + ".txt";
        ofstream clearFile(fname, ios::out | ios::trunc); // force clearing the file
        clearFile.close();
        //cout << "[DEBUG] Cleared file: " << fname << "\n";
    }
}


    if(rank == 0){
        runDijkstra(filename);
        runMetis(filename, nparts);

        // After both are done, now write partitioned outputs
        string metisOut = filename + ".part." + to_string(nparts);
         int src = 2;
        int dest = 6;
        int weight = 1;

        // Call the function to insert the edge into output.txt
        insertEdgeInOutput("output.txt", src, dest, weight);
        //cout<<"parting at"<<rank<<endl<<endl;
        writePartitionFiles(metisOut, "inserted.txt", nparts);
        
    }

MPI_Barrier(MPI_COMM_WORLD);
////////////////////////////////////////////////////////////////////////                        finding affected nodes                               /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    

    while (true) {
        bool local_change = false;
        filename = "partition_" + to_string(rank) + ".txt";
        ifstream infile(filename);
        ofstream affected_out("affected_" + to_string(rank) + ".txt", ios::trunc);


        unordered_map<int, int> nodeDistances;
        vector<tuple<int, int, vector<pair<int, int>>>> lines;

        string line;
        while (getline(infile, line)) {
            stringstream ss(line);
            int node_id, dist;
            ss >> node_id >> dist;
            nodeDistances[node_id] = dist;

            vector<pair<int, int>> neighbors;
            int neighbor_id, weight;
            while (ss >> neighbor_id >> weight) {
                neighbors.emplace_back(neighbor_id, weight);
            }
            lines.emplace_back(node_id, dist, neighbors);
        }
        infile.close();

        for (auto& [node_id, dist, neighbors] : lines) {
            for (auto& [neighbor_id, weight] : neighbors) {
                int tentative = dist + weight;

                if (nodeDistances.find(neighbor_id) != nodeDistances.end()) {
                    if (tentative < nodeDistances[neighbor_id]) {
                        //cout<<"the following vertex was affected: "<<neighbor_id<<"\n\n\n";
                        affected_out << "0 " << neighbor_id << " " << tentative << endl;
                        local_change = true;
                    }
                } else {
                    //cout<<"sending message\n\n";
                    Message msg{neighbor_id, tentative};
                    int dest = (rank == 0) ? 1 : 0;
                    MPI_Send(&msg, sizeof(Message), MPI_BYTE, dest, 0, MPI_COMM_WORLD);
                }
            }
        }

MPI_Barrier(MPI_COMM_WORLD);

        affected_out.close();

        // Probe for messages
        MPI_Status status;
        int flag;
        ofstream affected_out_msg("affected_" + to_string(rank) + ".txt", ios::app);
        while (true) {
            MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
            if (!flag) break;

            Message msg;
            MPI_Recv(&msg, sizeof(Message), MPI_BYTE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

            int neighbor_id = msg.node_id;
            int tentative = msg.tentative_distance;
            //cout<<"received neighbor id and tentative distance as follows: "<<msg.node_id<<" "<<msg.tentative_distance<<"\n\n\n\n";

            if (nodeDistances.find(neighbor_id) != nodeDistances.end()) {
                if (tentative < nodeDistances[neighbor_id]) {
                    //cout<<"the following vertex was affected: "<<neighbor_id<<"\n\n\n";
                    affected_out_msg << "0 " << neighbor_id << " " << tentative << endl;
                    local_change = true;
                }
            }
        }
        affected_out_msg.close();

        MPI_Barrier(MPI_COMM_WORLD);

        // Read and update affected nodes
        ifstream affected_file("affected_" + to_string(rank) + ".txt");
        vector<pair<int, int>> affected_nodes;
        vector<string> updated_affected_lines;

        if (affected_file.is_open()) {
            string line;
            while (getline(affected_file, line)) {
                istringstream iss(line);
                int flag, node_id, new_distance;
                if (!(iss >> flag >> node_id >> new_distance)) continue;
                if (flag == 0) {
                    affected_nodes.emplace_back(node_id, new_distance);
                    updated_affected_lines.push_back("1 " + to_string(node_id) + " " + to_string(new_distance));
                } else {
                    updated_affected_lines.push_back(line);
                }
            }
            affected_file.close();
        }

        string partition_filename = "partition_" + to_string(rank) + ".txt";
        ifstream partition_in(partition_filename);
        vector<string> updated_lines;

        if (partition_in.is_open()) {
            while (getline(partition_in, line)) {
                istringstream iss(line);
                int node;
                iss >> node;

                bool is_affected = false;
                int new_distance = -1;

                for (auto& p : affected_nodes) {
                    if (p.first == node) {
                        is_affected = true;
                        new_distance = p.second;
                        break;
                    }
                }

                if (is_affected) {
                    ostringstream oss;
                    oss << node << " " << new_distance;

                    int old_distance;
                    iss >> old_distance;
                    string rest_of_line;
                    getline(iss, rest_of_line);
                    oss << rest_of_line;

                    updated_lines.push_back(oss.str());
                } else {
                    updated_lines.push_back(line);
                }
            }
            partition_in.close();
        }

        ofstream partition_out(partition_filename, ios::trunc);
        for (const auto& line : updated_lines) {
            partition_out << line << endl;
        }
        partition_out.close();

        ofstream affected_out_updated("affected_" + to_string(rank) + ".txt", ios::trunc);
        for (const auto& line : updated_affected_lines) {
            affected_out_updated << line << endl;
        }
        affected_out_updated.close();
        




MPI_Barrier(MPI_COMM_WORLD);
      /* Print only if non-empty
      if (!affected_nodes.empty()) {
          cout << "Rank " << rank << " has affected_nodes:\n";
          for (size_t i = 0; i < affected_nodes.size(); ++i) {
              cout << "  [" << i << "] = (" << affected_nodes[i].first << ", " << affected_nodes[i].second << ")\n";
          }
}*/
      

      // Flatten regardless of size
      std::vector<int> flat_data;
      for (auto &p : affected_nodes) {
          flat_data.push_back(p.first);
          flat_data.push_back(p.second);
      }
      int local_size = flat_data.size();

      int world_size;
      MPI_Comm_size(MPI_COMM_WORLD, &world_size);

      // Gather sizes
      std::vector<int> recv_counts(world_size);
      MPI_Allgather(&local_size, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

      // Compute displacements
      std::vector<int> displs(world_size, 0);
      int total_size = recv_counts[0];
      for (int i = 1; i < world_size; ++i) {
          displs[i] = displs[i - 1] + recv_counts[i - 1];
          total_size += recv_counts[i];
      }

      // Allgatherv
      std::vector<int> all_flat_data(total_size);
      MPI_Allgatherv(flat_data.data(), local_size, MPI_INT,
                     all_flat_data.data(), recv_counts.data(), displs.data(), MPI_INT,
                     MPI_COMM_WORLD);

      // Reconstruct copied_nodes
      copied_nodes.clear();
      for (size_t i = 0; i + 1 < all_flat_data.size(); i += 2) {
          copied_nodes.push_back({all_flat_data[i], all_flat_data[i + 1]});
      }


/*cout << "Rank " << rank << " has copied_nodes: in whileeeeeeeeee:\n";
for (size_t i = 0; i < copied_nodes.size(); ++i) {
    cout << "  [" << i << "] = (" << copied_nodes[i].first << ", " << copied_nodes[i].second << ")\n";
}*/

if(temp==0){temp=copied_nodes.size();}

if(copied_nodes.size()>=half){
        if(temp>copied_nodes.size()){
        cout<<temp<<" i.e more(or equal) than 50% of the total nodes got affected, hence running djikstra after insertion\n\n";temp=-1;}else{temp=copied_nodes.size();}
        runDijkstra("inserted.txt");
break;

}

else{
if(copied_nodes.size()!=0){
if(temp>copied_nodes.size()){
cout<<temp<<" i.e less than 50% of the total nodes got affected, hence not running djikstra(performing updates for affected nodes)\n\n";}else{temp=copied_nodes.size();}}
#pragma omp parallel for
for (int i = 0; i < copied_nodes.size(); ++i) {
    
    update_partition_file(rank, copied_nodes[i]);  // Now takes one node
}
}

      
////////////          

        
        int global_change;
        int local_change_int = local_change ? 1 : 0;
        MPI_Allreduce(&local_change_int, &global_change, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);

        if (global_change == 0){
        

      //cout << "Rank " << rank << " has affected_nodes where needed:\n";
          for (size_t i = 0; i < affected_nodes.size(); ++i) {
              //cout << "  [" << i << "] = (" << affected_nodes[i].first << ", " << affected_nodes[i].second << ")\n";
          }
     



          
          
        break;}
    }

//cout << "Rank " << rank << " has copied_nodes:\n";
for (size_t i = 0; i < copied_nodes.size(); ++i) {
    //cout << "  [" << i << "] = (" << copied_nodes[i].first << ", " << copied_nodes[i].second << ")\n";
}


for (int i = 0; i < copied_nodes.size(); ++i) {
    
    update_partition_file(rank, copied_nodes[i]);  // Now takes one node
}



  






MPI_Barrier(MPI_COMM_WORLD);

double global_end = MPI_Wtime();
double exec_time = global_end - global_start;

if (rank == 0) {
    ofstream perf_log("perf_log.csv", ios::app);
    if (!perf_log) {
        cerr << "Warning: Could not open perf_log.csv\n";
    } else {
        perf_log << size << "," << omp_threads << "," << exec_time << "\n";
        perf_log.close();
        cout << "[Rank 0] Performance logged to perf_log.csv\n";
    }
}



  
  MPI_Finalize();

    return 0;
}