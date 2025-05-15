#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <metis.h>

int main() {
    std::string inputFile = "t.txt";
    int nparts = 4;

    std::ifstream infile(inputFile);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open input file.\n";
        return 1;
    }

    idx_t nvtxs, nedges, startVertex;
    infile >> nvtxs >> nedges >> startVertex;
    std::string dummy;
    std::getline(infile, dummy); // skip rest of first line

    std::vector<idx_t> xadj = {0};
    std::vector<idx_t> adjncy;
    std::vector<idx_t> adjwgt;

    std::string line;
    idx_t edgeCount = 0;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        idx_t neighbor, weight;
        while (iss >> neighbor >> weight) {
            adjncy.push_back(neighbor - 1); // Convert to 0-based indexing
            adjwgt.push_back(weight);
            edgeCount++;
        }
        xadj.push_back(edgeCount);
    }
    infile.close();

    // METIS parameters
    idx_t ncon = 1;
    idx_t objval;
    std::vector<idx_t> part(nvtxs);

    // Optional parameters
    idx_t* vwgt = NULL;
    idx_t* vsize = NULL;
    real_t* tpwgts = NULL;
    real_t* ubvec = NULL;

    // Options
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0; // 0-based

    // Call METIS
    int status = METIS_PartGraphKway(
        &nvtxs, &ncon, xadj.data(), adjncy.data(),
        vwgt, vsize, adjwgt.data(),
        &nparts, tpwgts, ubvec,
        options, &objval, part.data()
    );

    if (status != METIS_OK) {
        std::cerr << "METIS_PartGraphKway failed.\n";
        return 1;
    }

    // Write result to t.txt.part.4
    std::string partFileName = inputFile + ".part." + std::to_string(nparts);
    std::ofstream outfile(partFileName);
    if (!outfile.is_open()) {
        std::cerr << "Error writing to " << partFileName << "\n";
        return 1;
    }

    for (idx_t i = 0; i < nvtxs; ++i) {
        outfile << part[i] << "\n";
    }
    outfile.close();

    std::cout << "Partitioning completed. Result written to " << partFileName << "\n";
    return 0;
}

