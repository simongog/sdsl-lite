//
// Created by d056848 on 6/7/16.
//

#include <iostream>
#include <vector>
#include <tuple>
#include <complex>
#include <sdsl/k2_tree_algorithm.hpp>
#include <sys/times.h>
#include <sdsl/k2_tree_hybrid.hpp>

using std::ifstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using ::std::ofstream;
using namespace sdsl;
using namespace std;

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
        fprintf(stderr, "USAGE: %s <GRAPH> <Output-File>\n <GRAPH> has to be in the .graph-txt format", argv[0]);
        return (-1);
    }

    std::string fileName(argv[1]);
    if(!hasEnding(fileName, ".graph-txt")){
        fileName.append(".graph-txt");
        std::cout << "Appending .graph-txt to filename as file has to be in .graph-txt format" << std::endl;
    }

    std::fstream fileStream(fileName, std::ios_base::in);

    //simple string dictionary
    if (fileStream.is_open()) {
        cout << "Reading file " << argv[1] << endl;
        //std::vector<std::tuple<unsigned int, unsigned int>> idBuffer;

        std::string readBuffer;
        //read lines

        std::getline(fileStream, readBuffer);
        uint numberOfNodes = stoul(readBuffer);

        std::string output_file_name(argv[2]);

        int_vector_buffer<32> x_buffer(output_file_name+".x", std::ios::out);
        int_vector_buffer<32> y_buffer(output_file_name+".y", std::ios::out);

        uint source_id, target_id;
        while (std::getline(fileStream, readBuffer)) {
            //tokenizer<escaped_list_separator<char> > tok(readBuffer);
            stringstream ss(readBuffer);
            ss >> source_id >> target_id;

            target_id /= 10; //0 delimieted format

            x_buffer.push_back(source_id);
            y_buffer.push_back(target_id);
        }

    } else {
        throw "Could not load file";
    }

    return 0;
}

