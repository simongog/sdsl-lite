//
// Created by d056848 on 6/7/16.
//

#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <complex>
#include <sdsl/k2_tree_comp.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sys/times.h>
#include <sdsl/k2_tree_utility.hpp>
#include "sdsl/k2_tree_utility.hpp"

using std::ifstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using ::std::ofstream;
using namespace sdsl;
using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

void print(std::vector<uint32_t> &to_print, uint32_t source_node) {
    std::cout << "Reachable nodes from " << source_node << std::endl;
    for (auto &link: to_print) {
        std::cout << link << "\t";
    }
    std::cout << std::endl;
}

int main(int argc, char *argv[]) {




/* end Time meassuring */

    /*Graph Format:
    <NumNodes>
    Source  Target
    Source  Target


    First lines of eu2005
     862664
     0       240
     0       420
     0       620
     0       630
    */
    if (argc < 3) {
        fprintf(stderr, "USAGE: %s <GRAPH> <Output-File> [Hash-Size]\n <GRAPH> has to be in the .ladrabin format", argv[0]);
        return (-1);
    }

    std::string file_name(argv[1]);
    std::string output_file_name(argv[2]);

    const uint8_t k = 4;
    typedef k2_tree_comp<k, bit_vector, bit_vector> tested_type;
    //typedef k2_tree_hybrid<4,6,2,8, bit_vector, bit_vector> k2_rrr;
    //typedef k2_tree_partitioned<8, k2_rrr> tested_type;
    // Initialize treap with a vector of (x,y,weight) elements
    //construct_im(k2treap, coordinates, numberOfNodes - 1);

    tested_type k2tree;

    k2tree.load_from_ladrabin(file_name);
    store_to_file(output_file_name, k2tree);
}

