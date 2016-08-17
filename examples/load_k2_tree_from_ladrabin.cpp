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
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <sys/times.h>

using std::ifstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using ::std::ofstream;
using namespace sdsl;
using namespace std;
using namespace boost;

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
        fprintf(stderr, "USAGE: %s <GRAPH> <Output-File>\n <GRAPH> has to be in the .graph-txt format", argv[0]);
        return (-1);
    }

    std::string fileName(argv[1]);
    if (!boost::algorithm::ends_with(fileName, ".ladrabin")) {
        fileName.append(".ladrabin");
        std::cout << "Appending .graph-txt to filename as file has to be in .ladrabin format" << std::endl;
    }
    std::fstream fileStream(fileName, std::ios_base::in);
    if (fileStream.is_open()) {

        double t2 = 0;
        ticks = (double) sysconf(_SC_CLK_TCK);
        start_clock();

        uint number_of_nodes;
        ulong number_of_edges;
        _read_member(number_of_nodes, fileStream);
        _read_member(number_of_edges, fileStream);

        uint nodes_read = 0;
        uint edges_read = 0;
        uint source_id;
        int target_id;

        std::vector<std::pair<uint, uint>> coords(number_of_edges);
        for (uint64_t i = 0; i < number_of_nodes + number_of_edges; i++) {
            _read_member(target_id, fileStream);
            if (target_id < 0) {
                nodes_read++;
            } else {
                source_id = nodes_read - 1;
                coords.push_back(std::make_pair(source_id, target_id));
                edges_read++;
            }
        }

        fileStream.close();

        std::cerr << "Finished Reading File " << std::endl;
        std::cerr << "Amount of edges: " << coords.size() << std::endl;
        /*
        uint8_t k = (uint8_t) atoi(argv[3]);
        uint8_t shortcut_size = (uint8_t) atoi(argv[4]);
        std::stringstream ss(argv[5]);
        bool compress_leaves;

        if(!(ss >> std::boolalpha >> compress_leaves)) {
            throw std::runtime_error("Could not convert argv[5] to bool, please use true/false");

        }

        std::shared_ptr<k2_tree> k2 = get_k2_tree(k, shortcut_size, compress_leaves, coords, number_of_nodes-1);*/
        const uint8_t k = 4;
        typedef k2_tree_hybrid<4,5,2,8, bit_vector, bit_vector, false> k2_rrr;
        typedef k2_tree_partitioned<4, k2_rrr, true> k2_part;
        //typedef k2_tree<k, bit_vector, bit_vector, true, 4> k2_rrr;
        // Initialize treap with a vector of (x,y,weight) elements
        //construct_im(k2treap, coordinates, numberOfNodes - 1);
        k2_part k2tree("", false, coords, number_of_nodes - 1);
        coords.clear();

        std::string output_file_name(argv[2]);
        store_to_file(k2tree, output_file_name);

        t2 += stop_clock();
        t2 *= 1000; // to milliseconds
        fprintf(stderr, "Initialization time (ms): %f\n", t2);

        write_structure<HTML_FORMAT>(k2tree, output_file_name + "k_" + std::to_string(k) + ".html");
    } else {
        throw "Could not load file";
    }
}

