//
// Created by d056848 on 6/7/16.
//

#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <complex>
#include <sdsl/k2_tree.hpp>
#include <sdsl/k2_tree_algorithm.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sys/times.h>

using std::ifstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using ::std::ofstream;
using namespace sdsl;
using namespace std;

/* Time meassuring */
double ticks;
struct tms t1, t2;

void start_clock() {
    times(&t1);
}

double stop_clock() {
    times(&t2);
    return (t2.tms_utime - t1.tms_utime) / ticks;
}

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
    const uint8_t k = 4;
    //typedef k2_tree<k, bit_vector, bit_vector, true> k2_rrr;
    //typedef k2_tree_partitioned<8, k2_rrr, true> k2_part;
    typedef k2_tree_hybrid<4,5,2,8, bit_vector, bit_vector, true> k2_rrr;
    // Initialize treap with a vector of (x,y,weight) elements
    //construct_im(k2treap, coordinates, numberOfNodes - 1);

    uint64_t hash_size = 0;
    if (argc > 3){
        hash_size = stoull(argv[3]);
    }
    k2_rrr k2;
    k2.load_from_ladrabin(file_name, hash_size, true);

    std::string output_file_name(argv[2]);
    store_to_file(k2, output_file_name);

    write_structure<HTML_FORMAT>(k2, output_file_name + "k_" + std::to_string(k) + ".html");
}

