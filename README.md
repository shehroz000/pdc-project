# pdc-project
Project Overview: Parallel Dynamic Single-Source Shortest Path (SSSP)
The core goal of this project is to implement and evaluate a parallel algorithm for updating Single-Source Shortest Paths (SSSP) in large-scale dynamic networks. This is based on the research paper you provided:

Research Paper: "A Parallel Algorithm Template for Updating Single-Source Shortest Paths in Large-Scale Dynamic Networks" by Arindam Khanda, Sriram Srinivasan, Sanjukta Bhowmick, Boyana Norris, and Sajal K. Das (IEEE Transactions on Parallel and Distributed Systems, Vol. 33, No. 4, April 2022).   
The paper proposes a framework to efficiently update SSSP when network structures change (edges added/deleted, weights modified), aiming to be faster than recomputing from scratch. The key idea involves identifying the affected portion of the network and updating a relevant tree data structure.

Key Technologies & Requirements for the Project:

Parallel Programming: MPI (Message Passing Interface) for distributed memory parallelism.   
Node-Level Parallelism (Future Step): OpenMP or OpenCL for shared-memory or GPU parallelism within MPI processes.
Graph Partitioning: METIS library to divide the graph for MPI processing.
Testing & Evaluation: Use multiple datasets (public datasets were specified), evaluate scalability (strong/weak), compare performance (Sequential vs. MPI vs. MPI+Hybrid), and visualize results.
Implemented Solutions (Step-by-Step)
Here's a summary of the steps we've worked through and the solutions implemented:

Step 1: Environment Setup and Prerequisites

This foundational step involved ensuring the necessary software components were ready and understanding the core algorithm.

MPI (Message Passing Interface):
Verification: Confirmed how to check for an MPI installation (e.g., mpicc --version) and test it with a simple "Hello World" MPI program.
Purpose: Essential for distributed memory parallel processing across multiple nodes or cores.
Compiler with OpenMP Support:
Verification: Checked for a C/C++ compiler (like GCC, Clang) supporting OpenMP (e.g., using the -fopenmp flag) and tested with a basic OpenMP "Hello Threads" program.
Purpose: For future implementation of shared-memory parallelism within an MPI process.
OpenCL (Optional):
Consideration: Discussed checking for OpenCL SDKs and drivers, noting its complexity compared to OpenMP. This was marked as an optional path for node-level parallelism.
METIS Library:
Setup: Guided on downloading METIS from its official source and the typical compilation process (make config, make). Emphasized noting the library and header file paths for linking.
Purpose: To partition the graph into a specified number of parts, aiming to minimize edge cuts between partitions for efficient parallel processing.   
Understanding the Algorithm:
Action: Critically emphasized the need to thoroughly read and understand the specific SSSP update algorithm detailed in the provided research paper. This understanding is paramount for a correct implementation.
Step 2: Implement the Sequential Algorithm (Baseline)

This step focused on creating a correct sequential version of the SSSP update algorithm, which serves as a baseline for correctness and performance comparison.

Graph Representation: Used a C++ std::vector<std::vector<Edge>> (adjacency list) where Edge is a struct containing to (target vertex) and weight.
Initial SSSP: Implemented a standard Dijkstra's algorithm (dijkstra function) to compute the initial shortest paths from a source node on the entire graph.
SSSP Update Function (update_sssp):
Provided a C++ function structure for update_sssp.
Crucial Caveat: Since the exact details of the paper's update algorithm were initially unknown (before the PDF was processed), the provided C++ code for this function was a generic placeholder/template. It included basic logic for:
Modifying the graph structure (updating edge weight, handling rudimentary insertion/deletion).
A Dijkstra-like propagation mechanism using a priority queue, initiated if an edge weight decrease directly offered a shorter path.
Explicitly stated that the user MUST replace this placeholder logic with the specific, detailed algorithm from their research paper, especially for handling weight increases and deletions, and for accurately identifying the "affected portion" of the graph. The paper's method is expected to be more sophisticated.
Output: The functions would output updated distance (dist) and predecessor (pred) arrays.
Step 3: Integrate METIS for Graph Partitioning

This step involved incorporating the METIS library to partition the graph, preparing it for distributed processing.

CSR Format Conversion (convert_to_metis_csr):
Implemented a C++ function to convert the project's Graph representation (adjacency list of structs) into the Compressed Sparse Row (CSR) format (xadj, adjncy, adjwgt arrays) required by METIS.
Undirected Assumption for METIS: For partitioning purposes, the graph was treated as undirected (i.e., if an edge u -> v exists, METIS effectively sees an edge between u and v).
Unit Weights for METIS: The adjwgt array for METIS was populated with unit weights (1) for simplicity, as partitioning often relies more on connectivity than precise SSSP edge weights.
Calling METIS:
In the main function, after graph creation, the CSR conversion function was called.
METIS_PartGraphKway was invoked with the CSR data to partition the graph into a specified number of parts (nparts, hardcoded for sequential testing, intended to be world_size in MPI).
METIS options were set to minimize edge-cut and use 0-based numbering.
Output: The partition assignment for each vertex (stored in a part array, where part[i] is the partition ID for vertex i) and the resulting edge-cut were printed.
Purpose: This step provides the crucial part array, which dictates how vertices (and their associated data/computation) will be distributed among MPI processes.
Step 4: Implement the Parallel Algorithm using MPI

This was the most complex step, focusing on distributing the SSSP update process using MPI.

MPI Initialization: Standard MPI_Init, MPI_Comm_rank, MPI_Comm_size.   
Graph Loading & Initial Distribution (Rank 0 as Coordinator):
Rank 0 loads/defines the graph, performs METIS partitioning (using world_size for nparts).
Broadcasts global information: num_vertices, the part array (partition assignments), the entire graph structure (adjacency lists), and the initial global dist and pred arrays (computed by Rank 0 using sequential Dijkstra).
Simplification Note: Broadcasting the entire graph structure was chosen for implementation simplicity in this phase. A more optimized approach for very large graphs would involve Rank 0 sending only the relevant subgraph and ghost node information to each process.
Parallel Update Function (update_sssp_mpi):
Input: The global graph change (e.g., u_changed, v_changed, new_weight), graph data, current dist/pred arrays, part array, and MPI communicators.
Local Data: Each process works with its copy of the dist and pred arrays. It "owns" vertices assigned to its rank by METIS.
Local Priority Queue: Each process maintains its own std::priority_queue<PrioQueueEntry> for vertices it owns that require processing.
Graph Modification: The triggering edge change is applied to each process's local copy of the graph structure.
Initial PQ Seeding: Logic to seed the local PQs based on the update. This part was highlighted as needing to closely follow the paper's method for identifying the "affected portion."
Iterative Update Loop:
Local Processing: Processes pop vertices from their local PQ. For each relaxed edge (u,v):
If v is local (owned by the current process): dist[v] and pred[v] are updated, and v is added back to the local PQ if improved.
If v is remote (owned by another process P_owner): An UpdateMessage (containing target_vertex, new_distance, predecessor_vertex) is buffered for P_owner.
Communication Phase:
A custom MPI_Datatype (mpi_update_message_type) was created for the UpdateMessage struct to allow sending mixed data types.
MPI_Alltoallv is used for efficient, collective exchange of these update messages. Each process sends its buffered messages and receives messages destined for its owned vertices.
Processing Incoming Messages: Received updates are processed. If a distance to an owned vertex is improved, its dist/pred are updated, and it's added to the local PQ.
Termination Detection:
A boolean flag made_local_update_in_iteration (true if a distance was changed or the local PQ is not empty) is used.
MPI_Allreduce with MPI_LOR aggregates these flags across all processes. If no process is active, the global loop terminates.
A safety break (max iterations) was included.
Error Correction: Corrected an error where MPI_IDX_T (not a standard MPI type) was used for broadcasting the METIS part array, changing it to MPI_INT (assuming idx_t from METIS is compatible with int).
Output: Rank 0 prints the final updated SSSP distances and predecessors.
Summary of Where We Are
We have successfully:

Set up the conceptual environment.
Implemented a sequential SSSP update algorithm (with the understanding that its core logic for updates needs to be precisely aligned with the research paper).
Integrated METIS to partition the graph.
Implemented a parallel SSSP update algorithm using MPI, which includes data distribution (simplified), message passing for updates across partitions, and global termination detection.
The MPI implementation provides a framework for distributing the computation. The effectiveness of the update part heavily relies on how accurately the paper's specific logic for "identifying the affected portion" and handling increases/decreases is translated into the initial PQ seeding and the relaxation rules within the update_sssp_mpi function.   
