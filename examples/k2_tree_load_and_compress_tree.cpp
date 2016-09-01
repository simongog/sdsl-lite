#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <complex>
#include <sdsl/k2_tree.hpp>
#include <sdsl/k2_tree_algorithm.hpp>
#include <sdsl/bit_vectors.hpp>
#include <chrono>
#include <sys/times.h>

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
    typedef k2_tree<k, bit_vector, bit_vector> k2_rrr;
    //const uint8_t k = 4;
    //typedef k2_tree_hybrid<4,5,2,8, bit_vector, bit_vector,false> k2_rrr;
    typedef k2_tree_partitioned<8, k2_rrr> tested_type;

    //typedef k2_tree_partitioned<4, k2_rrr, true> k2_part;

//    typedef k2_tree_hybrid<4,5,2,8, bit_vector, bit_vector, false> k2_rrr;
//    typedef k2_tree_partitioned<4, k2_rrr, true> k2_part;

    tested_type k2tree;
    std::string fileName = argv[1];
    load_from_file(k2tree, fileName);

    k2tree.compress_leaves();

    std::string output_file_name(argv[2]);
    store_to_file(k2, output_file_name);

    write_structure<HTML_FORMAT>(k2, output_file_name + + "("+ k2.get_type_string()+ ")" + ".html");



    return 0;
}
