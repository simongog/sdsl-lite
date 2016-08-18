#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <complex>
#include <sdsl/k2_tree.hpp>
#include <sdsl/k2_tree_algorithm.hpp>
#include <sdsl/bit_vectors.hpp>
#include <boost/algorithm/string.hpp>
#include <chrono>
#include <sys/times.h>

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

using namespace sdsl;

int main(int argc, char* argv[]){

    if(argc<3){
        fprintf(stderr,"USAGE: %s <path to serialized k tree (set templates correctly!)> query_file (binary format) [use access shortcut]\n",argv[0]);
        return(-1);
    }

    //char *filename = (char *)malloc(sizeof(char)*20);

    const uint8_t k = 4;
    //typedef k2_tree_hybrid<k,k,k,k, bit_vector, bit_vector,true> k2_rrr;
         typedef k2_tree<k, bit_vector, bit_vector, true, 4> k2_rrr;
        //const uint8_t k = 4;
        //typedef k2_tree_hybrid<4,5,2,8, bit_vector, bit_vector> k2_rrr;
        //typedef k2_tree_partitioned<4, k2_rrr, true> k2_part;

//    typedef k2_tree_hybrid<4,5,2,8, bit_vector, bit_vector, false> k2_rrr;
//    typedef k2_tree_partitioned<4, k2_rrr, true> k2_part;

    bool use_shortut = argc > 3;

    k2_rrr k2tree;
    std::string fileName = argv[1];
    load_from_file(k2tree, fileName);

    char * list_file = argv[2];
    FILE * list_fp = fopen(list_file,"r");
    uint queries;
    fread(&queries, sizeof(uint), 1, list_fp);
    uint *qry = (uint *) malloc(sizeof(uint)*queries);
    fread(qry,sizeof(uint),queries,list_fp);
    fclose(list_fp);
    {
        uint i;
        ulong recovered = 0;
        std::vector<uint32_t> result;
        if (use_shortut) {
            std::cout << "Direct with access shortcut" << std::endl;
            auto start = timer::now();
            for (i = 0; i < queries; i++) {
                k2tree.direct_links_shortcut(qry[i], result);
                recovered += result.size();
            }
            auto stop = timer::now();

            std::cout <<  "Recovered Nodes:" << recovered << "\n";
            std::cout <<  "Queries:" << queries<< "\n";
            std::cout << "Total time(ms): "<<  duration_cast<nanoseconds>(stop-start).count() << "\n";
            std::cout << "Time per query: "<< duration_cast<nanoseconds>(stop-start).count()/queries << "\n";
            std::cout << "Time per link: "<< duration_cast<nanoseconds>(stop-start).count()/recovered << "\n";
        }
    }

    {
        uint i;
        ulong recovered = 0;
        std::vector<uint32_t> result;
        if (use_shortut) {
            std::cout << "Direct without access shortcut" << std::endl;
            auto start = timer::now();
            for (i = 0; i < queries; i++) {
                k2tree.direct_links2(qry[i], result);
                recovered += result.size();
            }
            auto stop = timer::now();

            std::cout <<  "Recovered Nodes:" << recovered << "\n";
            std::cout <<  "Queries:" << queries<< "\n";
            std::cout << "Total time(ms): "<<  duration_cast<nanoseconds>(stop-start).count() << "\n";
            std::cout << "Time per query: "<< duration_cast<nanoseconds>(stop-start).count()/queries << "\n";
            std::cout << "Time per link: "<< duration_cast<nanoseconds>(stop-start).count()/recovered << "\n";
        }
    }

    {
        uint i;
        ulong recovered = 0;
        std::vector<uint32_t> result;
        if (use_shortut) {
            std::cout << "Inverse with access shortcut" << std::endl;
            auto start = timer::now();
            for (i = 0; i < queries; i++) {
                k2tree.inverse_links_shortcut(qry[i], result);
                recovered += result.size();
            }
            auto stop = timer::now();

            std::cout <<  "Recovered Nodes:" << recovered << "\n";
            std::cout <<  "Queries:" << queries<< "\n";
            std::cout << "Total time(ms): "<<  duration_cast<nanoseconds>(stop-start).count() << "\n";
            std::cout << "Time per query: "<< duration_cast<nanoseconds>(stop-start).count()/queries << "\n";
            std::cout << "Time per link: "<< duration_cast<nanoseconds>(stop-start).count()/recovered << "\n";
        }
    }

    {
        uint i;
        ulong recovered = 0;
        std::vector<uint32_t> result;
        if (use_shortut) {
            std::cout << "Inverse without access shortcut" << std::endl;
            auto start = timer::now();
            for (i = 0; i < queries; i++) {
                k2tree.inverse_links2(qry[i], result);
                recovered += result.size();
            }
            auto stop = timer::now();

            std::cout <<  "Recovered Nodes:" << recovered << "\n";
            std::cout <<  "Queries:" << queries<< "\n";
            std::cout << "Total time(ms): "<<  duration_cast<nanoseconds>(stop-start).count() << "\n";
            std::cout << "Time per query: "<< duration_cast<nanoseconds>(stop-start).count()/queries << "\n";
            std::cout << "Time per link: "<< duration_cast<nanoseconds>(stop-start).count()/recovered << "\n";
        }
    }


    return 0;
}
