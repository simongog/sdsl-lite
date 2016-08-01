//
// Created by d056848 on 6/7/16.
//

#include <iostream>
#include <vector>
#include <tuple>
#include <complex>
#include <sdsl/k2_tree_algorithm.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
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
using namespace boost;

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
    if(!boost::algorithm::ends_with(fileName, ".graph-txt")){
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

        while (std::getline(fileStream, readBuffer)) {
            //tokenizer<escaped_list_separator<char> > tok(readBuffer);
            vector <string> sourceNodeAndTargets;
            boost::split(sourceNodeAndTargets, readBuffer, boost::is_any_of("\t"));

            uint sourceId = std::stoul(sourceNodeAndTargets[0].c_str());
            uint targetId = std::stoul(sourceNodeAndTargets[1].c_str()) / 10;//Strange 0 delimited format
            //cout << "adding " << sourceId << "," << targetId << endl;
            x_buffer.push_back(sourceId);
            y_buffer.push_back(sourceId);
        }

    } else {
        throw "Could not load file";
    }

    return 0;
}

