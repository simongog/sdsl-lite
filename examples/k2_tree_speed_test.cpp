#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <complex>
#include "sdsl/k2_tree.hpp"
#include "sdsl/k2_tree_algorithm.hpp"
#include <sdsl/bit_vectors.hpp>
#include <chrono>
#include <sys/times.h>
#include "sdsl/k2_tree_utility.hpp"

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

using namespace sdsl;

int main(int argc, char *argv[]) {

    if (argc < 3) {
        fprintf(stderr,
                "USAGE: %s <path to serialized k tree (set templates correctly!)> query_file (binary format) [use access shortcut]\n",
                argv[0]);
        return (-1);
    }

    //char *filename = (char *)malloc(sizeof(char)*20);
    const uint8_t k = 4;
    //typedef k2_tree_hybrid<k,k,k,k, bit_vector, bit_vector,true> k2_rrr;
//    typedef k2_tree<k, bit_vector, bit_vector> k2_rrr;
    //const uint8_t k = 4;

    //typedef k2_tree_partitioned<4, k2_rrr, true> k2_part;

//    typedef k2_tree_hybrid<4,5,2,8, bit_vector, bit_vector, false> k2_rrr;
//    typedef k2_tree_partitioned<4, k2_rrr, true> k2_part;
    //typedef k2_tree<k, bit_vector, bit_vector> tested_type;
    typedef k2_tree_hybrid<4,6,2,8, bit_vector, bit_vector> k2_rrr;
    typedef k2_tree_partitioned<8, k2_rrr> tested_type;
    // Initialize treap with a vector of (x,y,weight) elements
    //construct_im(k2treap, coordinates, numberOfNodes - 1);

    bool use_shortcut = argc > 3;

    tested_type k2tree;
    std::string fileName = argv[1];
    load_from_file(k2tree, fileName);
    std::string queryFile = argv[2];

    access_times times = perform_speed_test(queryFile, k2tree, use_shortcut);

    if (use_shortcut){
        //Construction Time	Compressed Size (Byte)	Bpe	Direct Short (ns)	Direct (ns)	Inverse Short (ns)	Inverse (ns)	Check S (ns)	Check (ns)
        std::cout << "Hereyougo:" << times.direct_short_time <<","<< times.direct_time <<","<< times.inverse_short_time <<","<< times.inverse_time <<","<< times.check_short_time <<","<< times.check_time << std::endl;
    } else {
        //Construction Time	Compressed Size (Byte)	Bpe	Direct (ns)	Inverse (ns)	Check (ns)
        std::cout << "Hereyougo:" << times.direct_time << "," << times.inverse_time << "," << times.check_time << std::endl;
    }
    return 0;
}
