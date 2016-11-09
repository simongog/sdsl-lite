#include <iostream>
#include <vector>
#include <tuple>
#include <complex>
#include <fstream>
#include <sdsl/io.hpp>
#include <sdsl/k2_tree_utility.hpp>

using std::ifstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using ::std::ofstream;

using namespace std;
using namespace sdsl;

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
    if (argc < 5) {
        fprintf(stderr, "USAGE: %s <GRAPH> <Output-File> <Number of nodes> <Number of edges>\n <GRAPH> has to be a csv edge list (separated by tab)", argv[0]);
        return (-1);
    }

    std::string file_name(argv[1]);
    std::fstream fileStream(file_name, std::ios_base::in);

    //simple string dictionary
    if (fileStream.is_open()) {
        cout << "Reading file " << argv[1] << endl;
        //std::vector<std::tuple<unsigned int, unsigned int>> idBuffer;

        std::string read_buffer;
        uint64_t number_of_nodes = strtoull(argv[3], NULL, 10);
        uint64_t number_of_edges = strtoull(argv[4], NULL, 10);

        std::string output_file_name(argv[2]);

        if (!has_ending(output_file_name, ".ladrabin64"))
            output_file_name = output_file_name + ".ladrabin64";

        std::ofstream out;
        out.open(output_file_name);
        if (!out)
            throw std::runtime_error("Can not open file \"" + output_file_name + "\" for writing.");

        write_member(number_of_nodes, out);
        write_member(number_of_edges, out);

        uint64_t edge_counter = 0;
        std::vector<int64_t> coords;
        uint64_t source_id, target_id;
        uint64_t last_written_source_id = 0;
        uint64_t previous_source_id = 0;
        int64_t negative_node_index = -1;

        if (std::getline(fileStream, read_buffer)){
            istringstream ss(read_buffer);
            ss >> source_id >> target_id;

            previous_source_id = source_id;
            last_written_source_id = -1;
            coords.push_back(target_id);
            edge_counter++;
        }

        while (std::getline(fileStream, read_buffer)) {
            //tokenizer<escaped_list_separator<char> > tok(readBuffer);
            istringstream ss(read_buffer);
            ss >> source_id >> target_id;

            if (previous_source_id != source_id){
		if (source_id % 10000000 == 0) std::cout << "Wrote " << source_id << " nodes" << std::endl;
                for (uint i = 0; i < previous_source_id - last_written_source_id; ++i) {
                    write_member(negative_node_index, out);
                    negative_node_index--;
                }
                last_written_source_id = previous_source_id;
                previous_source_id = source_id;
                write_member(&coords[0], coords.size(), out);

                coords.clear();
            }

            coords.push_back(target_id);
            edge_counter++;
        }

        for (uint i = 0; i < previous_source_id - last_written_source_id; ++i) {
            write_member(negative_node_index, out);
            negative_node_index--;
        }
        
        write_member(&coords[0], coords.size(), out);
    	coords.clear();

        for (uint i = 0; i < number_of_nodes - previous_source_id -1; ++i) {
            write_member(negative_node_index, out);
            negative_node_index--;
        }

        fileStream.close();
        out.close();

        if (edge_counter != number_of_edges){
            std::cerr << "Specified edge number was incorrect, correct amount: " << edge_counter << std::endl;
        }
    } else {
        throw "Could not load file";
    }

    return 0;
}

