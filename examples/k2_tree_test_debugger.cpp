//
// Created by d056848 on 17.08.16.
//

#include <string>
#include <ios>
#include <iostream>
#include <vector>
#include <sdsl/vectors.hpp>
#include <sdsl/k2_tree.hpp>
#include <sdsl/k2_tree_partitioned.hpp>

using namespace std;
using namespace sdsl;

string test_file;
string temp_file;
bool in_memory;

void performTest(){
    int_vector_buffer<> buf_x(test_file + ".x", std::ios::in);
    int_vector_buffer<> buf_y(test_file + ".y", std::ios::in);
    typedef k2_tree<2, bit_vector, bit_vector, false> k2;
    k2_tree_partitioned<2, k2, true> failingK(buf_x, buf_y, false);
    k2 compare(buf_x, buf_y, false);

    std::vector<uint> result;
    compare.direct_links2((uint) 20184, result);

    std::cout << "results compare" << std::endl;
    for (auto res : result) {
        std::cout << res << "\t";
    }
    std::cout << std::endl;
    failingK.direct_links2((uint) 20184, result);

    std::cout << "results failing" << std::endl;
    for (auto res : result) {
        std::cout << res << "\t";
    }
    std::cout << std::endl;
}




int main(int argc, char **argv) {
    if (argc < 3) {
        // LCOV_EXCL_START
        cout << "Usage: " << argv[0] << " file temp_file [in-memory]" << endl;
        cout << " (1) Generates a k2-treap out of file.x, file.y, and file.w." << endl;
        cout << "     Result is stored in temp_file." << endl;
        cout << "     If `in-memory` is specified, the in-memory construction is tested." << endl;
        cout << " (2) Performs tests." << endl;
        cout << " (3) Deletes temp_file." << endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    test_file = argv[1];
    temp_file = argv[2];
    in_memory = argc > 3;
    if (in_memory) {
        auto load_and_store_in_mem = [&](string suf) {
            int_vector<> data;
            string file = temp_file + suf;
            load_vector_from_file(data, file);
            string ram_file = ram_file_name(file);
            store_to_file(data, ram_file);
        };
        load_and_store_in_mem("x");
        load_and_store_in_mem("y");
        temp_file = ram_file_name(temp_file);
    }
    performTest();

    exit(0);
}