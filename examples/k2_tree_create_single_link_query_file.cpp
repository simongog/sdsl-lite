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
std::vector<std::vector<t_x>> load_and_convert_ladrabin(ifstream& input_file, t_x){
    typedef typename std::make_signed<t_x>::type t_x_signed;
    t_x number_of_nodes;
    ulong number_of_edges;
    read_member(number_of_nodes, input_file);
    read_member(number_of_edges, input_file);

    t_x nodes_read = 0;
    t_x source_id;
    t_x_signed target_id;

    std::vector<std::vector<t_x>> adjList(number_of_nodes);

    for (uint64_t i = 0; i < number_of_nodes + number_of_edges; i++) {
        read_member(target_id, input_file);
        if (target_id < 0) {
            nodes_read++;
        } else {
            source_id = nodes_read - 1;
            adjList[source_id].push_back(target_id);
        }
    }

    return adjList;
}

template <typename t_x>
void write_query_file(vector<vector<t_x>>& adjList, ofstream& ofstream, int count){
    double present_link_factor = 0.5;

    long random_links = count * present_link_factor;
    long present_links = count - random_links;

    srand(0);

    write_member(count, ofstream);
    for (int i = 0; i < random_links; ++i) {
        int source = rand() % (adjList.size() - 1);
        int target = rand() % (adjList.size() - 1);
        write_member(source, ofstream);
        write_member(target, ofstream);

    }

    for (int i = 0; i < present_links; ++i) {
        int node = rand() % (adjList.size() - 1);
        while  (adjList[node].size() <= 1)
            node = rand() % (adjList.size() - 1);

        int targetIndex = rand() % (adjList[node].size() - 1);

        int target = adjList[node][targetIndex];
        write_member(node, ofstream);
        write_member(target, ofstream);
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
    if (argc < 4) {
        fprintf(stderr, "USAGE: %s <Graph> <Output-File> <Query_Count>\n <GRAPH> has to be in the .ladrabin format", argv[0]);
        return (-1);
    }

    std::string fileName(argv[1]);
    std::string outputFileName(argv[2]);
    int queryCount = atoi(argv[3]);
    bool use_64_bits = false;
    if (has_ending(fileName, ".ladrabin64")){
        use_64_bits = true;
    }

    if(!has_ending(fileName, ".ladrabin") && !has_ending(fileName, ".ladrabin64")){
        fileName.append(".ladrabin");
        std::cout << "Appending .ladrabin to filename and asuming 32 bit node identifiers as file has to be in .ladrabin format" << std::endl;
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
        auto adjList = load_and_convert_ladrabin(input_file, type);
        write_query_file(adjList, output_file, queryCount);
    } else {
        uint32_t type = 0;
        auto adjList = load_and_convert_ladrabin(input_file, type);
        write_query_file(adjList, output_file, queryCount);
    }

    input_file.close();
    output_file.close();
}