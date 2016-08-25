//
// Created by d056848 on 6/7/16.
//

#include <iostream>
#include <vector>
#include <tuple>
#include <complex>
#include <sdsl/k2_tree_hybrid.hpp>
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
        vector<pair<uint32_t, uint32_t>> coordinates;

        uint source_id, target_id;
        while (std::getline(fileStream, readBuffer)) {
            //tokenizer<escaped_list_separator<char> > tok(readBuffer);
            stringstream ss(readBuffer);
            ss >> source_id >> target_id;

            target_id /= 10; //0 delimieted format

            //cout << "adding " << sourceId << "," << targetId << endl;
            coordinates.push_back(make_pair(source_id, target_id));
        }

        fileStream.close();

        std::cerr << "Finished Reading File " << std::endl;
        std::cerr << "Amount of edges: " << coordinates.size() << std::endl;
        //typedef k2_tree<2,rrr_vector<63>> k2_part;
        const uint8_t k0 = 2;
        const uint8_t k = 2;
        typedef k2_tree_hybrid<4, 5, 2, 4, bit_vector, rrr_vector<63>,true> k2_rrr;
        typedef k2_tree_partitioned<k0,k2_rrr> k2_part;

        double t2 = 0;
        ticks = (double) sysconf(_SC_CLK_TCK);
        start_clock();
        // Initialize treap with a vector of (x,y,weight) elements
        k2_part k2tree_part;

    	construct_im_bottom_up(k2tree_part, coordinates, numberOfNodes - 1);

        t2 += stop_clock();
        t2 *= 1000; // to milliseconds

        coordinates.clear();

        fprintf(stderr, "Initialization time (ms): %f\n", t2);


        std::cerr << "Processing direct neighbour queries\n" << std::endl;

        uint i;

        srand(0);
        std::vector<uint> queries;
        uint queryCount = 10000;
        for (i = 0; i < queryCount; i++) {
            queries.push_back(rand() % 862664);
        }


        double t = 0;
        start_clock();
        uint count = 0;
        std::vector<uint32_t> result;
        for (auto query : queries){
            k2tree_part.direct_links2(query, result);
            count+= result.size();
        }

        t += stop_clock();
        t *= 1000; // to milliseconds

        std::cerr << "Total time(ms): "<< t <<"\n";
        std::cerr << "Recovered Nodes: "<< count << "\n";
        std::cerr << "Queries: " << queryCount << "\n";
        std::cerr << "Total time(ms): " << t << "\n";
        std::cerr << "Time per query: " << t/queryCount << "\n";
        std::cerr << "Time per link: " << t/count << "\n";

        std::cerr << "Processing reverse neighbour queries" << std::endl;

        double t3 = 0;
        start_clock();
        uint count2 = 0;
        for (auto query : queries){
            k2tree_part.inverse_links2(query, result);
            count2 += result.size();
        }

        //uint count = 342045;
        t3 += stop_clock();
        t3 *= 1000; // to milliseconds

        std::cerr << "Total time(ms): "<< t3 <<"\n";
        std::cerr << "Recovered Nodes: "<< count2 << "\n";
        std::cerr << "Queries: " << queryCount << "\n";
        std::cerr << "Total time(ms): " << t3 << "\n";
        std::cerr << "Time per query: " << t3/queryCount << "\n";
        std::cerr << "Time per link: " << t3/count2 << "\n";

        std::string output_file_name(argv[2]);
        store_to_file(k2tree_part,output_file_name);
        write_structure<HTML_FORMAT>(k2tree_part,output_file_name+"k0_"+std::to_string(k0)+"k_"+std::to_string(k)+".html");
    } else {
        throw "Could not load file";
    }
}

