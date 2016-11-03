#include <iostream>
#include <vector>
#include <tuple>
#include <complex>
#include <fstream>
#include <sdsl/io.hpp>

using std::ifstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using ::std::ofstream;

using namespace std;
using namespace sdsl;

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
    if (argc < 4) {
        fprintf(stderr, "USAGE: %s <GRAPH> <Output-File> <Number of edges>\n <GRAPH> has to be in the .graph-txt format", argv[0]);
        return (-1);
    }

    std::string file_name(argv[1]);
    if(!hasEnding(file_name, ".graph-txt")){
        file_name.append(".graph-txt");
        std::cout << "Appending .graph-txt to filename as file has to be in .graph-txt format" << std::endl;
    }

    std::fstream fileStream(file_name, std::ios_base::in);

    //simple string dictionary
    if (fileStream.is_open()) {
        cout << "Reading file " << argv[1] << endl;
        //std::vector<std::tuple<unsigned int, unsigned int>> idBuffer;

        std::string read_buffer;
        std::getline(fileStream, read_buffer);
        uint number_of_nodes = stoul(read_buffer);


        std::string output_file_name(argv[2]);
        output_file_name = output_file_name + ".ladrabin";

        std::ofstream out;
        out.open(output_file_name);
        if (!out)
            throw std::runtime_error("Can not open file \"" + output_file_name + "\" for writing.");

        write_member(number_of_nodes, out);

        ulong number_of_edges = stoull(argv[3]);
        ulong edge_counter = 0;
        write_member(number_of_edges, out);

        std::vector<int> coords;
        uint source_id, target_id;
        uint last_written_source_id = 0;
        uint previous_source_id = 0;
        int negative_node_index = -1;

        if (std::getline(fileStream, read_buffer)){
            istringstream ss(read_buffer);
            ss >> source_id >> target_id;
            target_id /= 10; //0 delimieted format

            previous_source_id = source_id;
            last_written_source_id = -1;
            coords.push_back(target_id);
            edge_counter++;
        }

        while (std::getline(fileStream, read_buffer)) {
            //tokenizer<escaped_list_separator<char> > tok(readBuffer);
            istringstream ss(read_buffer);
            ss >> source_id >> target_id;
            target_id /= 10; //0 delimieted format



            if (previous_source_id != source_id){
                for (int i = 0; i < previous_source_id - last_written_source_id; ++i) {
                    write_member(negative_node_index, out);
                    negative_node_index--;
                }
                last_written_source_id = previous_source_id;
                previous_source_id = source_id;
                for (auto target : coords) {
                    write_member(target, out);
                }
                coords.clear();
            }

            coords.push_back(target_id);
            edge_counter++;
        }

        for (int i = 0; i < previous_source_id - last_written_source_id; ++i) {
            write_member(negative_node_index, out);
            negative_node_index--;
        }
        write_member(&coords[0], coords.size(), out);
        coords.clear();

        for (int i = 0; i < number_of_nodes - previous_source_id -1; ++i) {
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

