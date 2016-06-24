//
// Created by d056848 on 6/7/16.
//

#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <complex>
#include <sdsl/k2_treap.hpp>
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
    fileName.append(".graph-txt");
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

        while (std::getline(fileStream, readBuffer)) {
            //tokenizer<escaped_list_separator<char> > tok(readBuffer);
            vector <string> sourceNodeAndTargets;
            boost::split(sourceNodeAndTargets, readBuffer, boost::is_any_of("\t"));

            uint sourceId = std::stoul(sourceNodeAndTargets[0].c_str());
            uint targetId = std::stoul(sourceNodeAndTargets[1].c_str()) / 10;//Strange 0 delimited format
            //cout << "adding " << sourceId << "," << targetId << endl;
            coordinates.push_back(make_pair(sourceId, targetId));
        }

        fileStream.close();

        std::cerr << "Finished Reading File " << std::endl;
        std::cerr << "Amount of edges: " << coordinates.size() << std::endl;
        typedef k2_treap<4, rrr_vector<63>> k2_rrr;
        k2_rrr k2treap;


        double t2 = 0;
        ticks = (double) sysconf(_SC_CLK_TCK);
        start_clock();
        // Initialize treap with a vector of (x,y,weight) elements
        construct_im(k2treap, coordinates);

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
            k2treap.direct_links(query, result);
            count+= result.size();
        }

        t += stop_clock();
        t *= 1000; // to milliseconds

        fprintf(stderr, "Total time(ms): %f\n", t);
        fprintf(stderr,"Recovered Nodes: %lu\n",count);
        fprintf(stderr,"Queries: %d\n",queryCount);
        fprintf(stderr,"Total time(ms): %f",t);
        fprintf(stderr,"Time per query: %f\n",t/queryCount);
        fprintf(stderr,"Time per link: %f\n",t/count);



        std::cerr << "Processing reverse neighbour queries\n" << std::endl;

        double t3 = 0;
        start_clock();
        uint count2 = 0;
        for (auto query : queries){
            k2treap.inverse_links(query, result);
            count2 += result.size();
        }

        //uint count = 342045;
        t3 += stop_clock();
        t3 *= 1000; // to milliseconds

        fprintf(stderr, "Total time(ms): %f\n", t3);
        fprintf(stderr,"Recovered Nodes: %lu\n",count2);
        fprintf(stderr,"Queries: %d\n",queryCount);
        fprintf(stderr,"Total time(ms): %f",t3);
        fprintf(stderr,"Time per query: %f\n",t3/queryCount);
        fprintf(stderr,"Time per link: %f\n",t3/count2);


        store_to_file(k2treap,argv[2]);
        write_structure<HTML_FORMAT>(k2treap,"structure.json");
    } else {
        throw "Could not load file";
    }
}

