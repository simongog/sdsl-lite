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

using std::ifstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using ::std::ofstream;
using namespace sdsl;
using namespace std;

template <typename t_x>
void load_and_convert_ladrabin(ifstream& input_file, ofstream& output_file, t_x){
    typedef typename std::make_signed<t_x>::type t_x_signed;
    t_x number_of_nodes;
    ulong number_of_edges;
    read_member(number_of_nodes, input_file);
    read_member(number_of_edges, input_file);

    int node_identifier_size = sizeof(t_x)*8; //Use 32 Bit Node Identifiers
    write_member(node_identifier_size, output_file);
    write_member(number_of_nodes, output_file);
    write_member(number_of_edges, output_file);

    t_x nodes_read = 0;
    t_x source_id;
    t_x_signed target_id;

    vector<pair<t_x,t_x>> buffer;
    buffer.reserve(10000000);

    for (uint64_t i = 0; i < number_of_nodes + number_of_edges; i++) {
        read_member(target_id, input_file);
        if (target_id < 0) {
            nodes_read++;
            if (nodes_read % 10000000 == 0){
                cout << "Read " << nodes_read << endl;
                write_member(&buffer[0], buffer.size(), output_file);
                buffer.clear();
            }
        } else {
            source_id = nodes_read - 1;
            buffer.push_back(make_pair(source_id, target_id));
        }
    }
    write_member(&buffer[0], buffer.size(), output_file);
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
        fprintf(stderr, "USAGE: %s <Graph> <Output-File>\n <GRAPH> has to be in the .ladrabin format", argv[0]);
        return (-1);
    }

    std::string fileName(argv[1]);
    std::string outputFileName(argv[2]);
    bool use_64_bits = false;
    if (has_ending(fileName, ".ladrabin64")){
        use_64_bits = true;
    }

    if(!has_ending(fileName, ".ladrabin") && !has_ending(fileName, ".ladrabin64")){
        fileName.append(".ladrabin");
        std::cout << "Appending .ladrabin to filename and asuming 32 bit node identifiers as file has to be in .ladrabin format" << std::endl;
    }

    if(!has_ending(outputFileName, ".edge")){
        outputFileName.append(".edge");
        std::cout << "Appending .edge to output filename" << std::endl;
    }
    std::ifstream input_file(fileName, std::ios_base::in);
    std::ofstream output_file(outputFileName, std::ios_base::binary);
    if (!input_file.is_open()) {
        std::cout << "Could not open input file" << std::endl;
        return -1;
    }

    if (!output_file.is_open()) {
        std::cout << "Could not open output file" << std::endl;
        return -2;
    }

    if (use_64_bits){
        uint64_t type = 0;
        load_and_convert_ladrabin(input_file, output_file, type);
    } else {
        uint32_t type = 0;
        load_and_convert_ladrabin(input_file, output_file, type);
    }

    input_file.close();
    output_file.close();
}