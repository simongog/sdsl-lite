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

// Writes primitive-typed variable t to stream out
template<class T>
size_t _write_member(const T &t, std::ostream &out) {
    out.write((char *) &t, sizeof(t));
    size_t written_bytes = sizeof(t);
    return written_bytes;
}

// Writes array-typed variable t of length length to stream out
template<class T>
size_t _write_member(const T *t, size_t length, std::ostream &out) {
    out.write((char *) t, length * sizeof(T));
    size_t written_bytes = length * sizeof(T);
    return written_bytes;
}

// Reads primitive-typed variable t from stream om
template<class T>
void _read_member(T &t, std::istream &in) {
    in.read((char *) &t, sizeof(t));
}

// Reads array-typed variable t of length length from stream in
template<class T>
void _read_member(T **t, size_t length, std::istream &in) {
    *t = (T *) malloc(sizeof(T) * length);
    in.read((char *) *t, length * sizeof(T));
}

bool hasEnding (std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
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
        fprintf(stderr, "USAGE: %s <GRAPH> S (tm)\n <GRAPH> has to be in the .ladrabin format", argv[0]);
        return (-1);
    }

    std::string fileName(argv[1]);
    if(!hasEnding(fileName, ".ladrabin")){
        fileName.append(".ladrabin");
        std::cout << "Appending .graph-txt to filename as file has to be in .ladrabin format" << std::endl;
    }
    std::fstream fileStream(fileName, std::ios_base::in);
    if (fileStream.is_open()) {

        uint number_of_nodes;
        ulong number_of_edges;
        _read_member(number_of_nodes, fileStream);
        _read_member(number_of_edges, fileStream);

        uint s = stoul(argv[2]);
        uint nodes_read = 0;
        uint source_id;
        int target_id;

/*
        uint res = 0;
        while (res <= 64 and pow(k0, res) <= number_of_nodes) { ++res; }
        if (res == 65) {
            throw std::logic_error("Maximal element of input is too big.");
        }
*/
        /*uint64_t matrix_size = pow(k0, res);
        std::cout << matrix_size << std::endl;*/
        uint64_t part_matrix_size = 1<<s;
        std::cout << part_matrix_size << std::endl;

        //std::vector<std::pair<uint, uint>> coords(number_of_edges);

        uint number_of_parts = (number_of_nodes/part_matrix_size);
        if (number_of_nodes%part_matrix_size != 0){
            number_of_parts++;
        }
        std::vector<uint> counter(number_of_parts*number_of_parts);
        for (uint64_t i = 0; i < number_of_nodes + number_of_edges; i++) {
            _read_member(target_id, fileStream);
            if (target_id < 0) {
                nodes_read++;
            } else {
                source_id = nodes_read - 1;
                uint p1 = source_id / part_matrix_size;
                uint p2 = target_id / part_matrix_size;
                uint corresponding_matrix = p1 * number_of_parts + p2;
                counter[corresponding_matrix]++;
            }
        }

        uint number_of_ones = 0;
        uint number_of_zeros = 0;
        for (uint i = 0; i < number_of_parts; ++i) {
            for (uint j = 0; j < number_of_parts; ++j) {
                if (counter[i*number_of_parts+j] > 0){
                    cout << "\033[1;31m \u25A0 \033[0m";
                    number_of_ones++;
                } else {
                    cout << "\033[1;32m \u25A0 \033[0m";
                    number_of_zeros++;
                }
            }
            std::cout << std::endl;
        }

        std::cout << "Number of ones: " << number_of_ones << std::endl;
        std::cout << "Number of zeros: " << number_of_zeros << std::endl;

        fileStream.close();
    } else {
        throw "Could not load file";
    }
}