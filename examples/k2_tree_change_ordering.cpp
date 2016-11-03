#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <complex>
#include <sdsl/k2_tree.hpp>
#include <sdsl/bit_vectors.hpp>
#include <chrono>
#include <sys/times.h>
#include <sdsl/k2_tree_hybrid.hpp>
#include <sdsl/k2_tree_partitioned.hpp>

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

using namespace sdsl;

int main(int argc, char *argv[]) {

    unsigned int number_of_nodes;
    unsigned long number_of_edges;
    
    if(argc<4){
        fprintf(stderr,"USAGE: %s <GRAPH> <map file> <out file>\n",argv[0]);
        return(-1);
    }

    std::fstream graph(argv[1], std::ios_base::in);

    if(!graph.is_open()){
        std::cout << "could not open graph file in ladrabin format" << std::endl;
        return 1;
    }

    read_member(number_of_nodes, graph);
    read_member(number_of_edges, graph);

    std::vector<uint> mapping(number_of_nodes);

    std::string read_buffer;
    std::fstream mappingFile(argv[2], std::ios_base::in);
    uint64_t counter = 0;

    if(!mappingFile.is_open()){
        std::cout << "could not open graph file in ladrabin format" << std::endl;
        return 1;
    }

    while (std::getline(mappingFile, read_buffer)) {
        std::istringstream ss(read_buffer);
        ss >> mapping[counter];
        counter++;
	if (counter % 10000000 == 0){
		std::cout << "Read " << counter << " from map" << std::endl;
	}
    }
    mappingFile.close();

    std::vector<std::vector<uint>> adjList(number_of_nodes);
    std::vector<uint> buffer;

    //load ladrabin format and apply mapping
    int nodes_read = -1;
    int target_id;

    for (uint64_t i = 0; i < (number_of_nodes + number_of_edges); i++) {
        read_member(target_id, graph);
        if (target_id < 0) {
            if(nodes_read != -1)//skip first
                adjList[mapping[nodes_read]].swap(buffer);
            buffer.clear();
            nodes_read++;
	    if (nodes_read % 10000000 == 0){
		std::cout << "Read " << nodes_read << " from adjlist" << std::endl;
	    }
        } else {
            buffer.push_back(target_id);
        }
    }

    //cover leftovers
    adjList[mapping[nodes_read]].swap(buffer);

    //serialize adjList as ladrabin
    std::fstream remapped_graph(argv[3], std::ios_base::out);

    write_member(number_of_nodes, remapped_graph);
    write_member(number_of_edges, remapped_graph);

    if(adjList.size() > INT_MAX){
        std::cout << "Currently only signed integer node ids are supported in ladrabin format" << std::endl;
    }

    for(int i = 0; i < adjList.size(); ++i){
        write_member(-(i+1), remapped_graph);
        write_member(&adjList[i][0],adjList[i].size(), remapped_graph);
    }

    remapped_graph.close();

    return 0;
}
