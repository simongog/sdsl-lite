#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <complex>
#include <sdsl/k2_tree.hpp>
#include <sdsl/k2_tree_algorithm.hpp>
#include <sdsl/bit_vectors.hpp>
#include <chrono>
#include <sys/times.h>

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

using namespace sdsl;

int main(int argc, char *argv[]) {

    if (argc < 3) {
        fprintf(stderr,
                "USAGE: %s <path to serialized k tree (set templates correctly!)> query_file (binary format) [use access shortcut]\n",
                argv[0]);
        return (-1);
    }

    //char *filename = (char *)malloc(sizeof(char)*20);

    const uint8_t k = 4;
    //typedef k2_tree_hybrid<k,k,k,k, bit_vector, bit_vector,true> k2_rrr;
    //typedef k2_tree<k, bit_vector, bit_vector, true, 4> k2_rrr;
    //const uint8_t k = 4;
    typedef k2_tree_hybrid<4,5,2,8, bit_vector, bit_vector> k2_rrr;
    //typedef k2_tree_partitioned<4, k2_rrr, true> k2_part;

//    typedef k2_tree_hybrid<4,5,2,8, bit_vector, bit_vector, false> k2_rrr;
//    typedef k2_tree_partitioned<4, k2_rrr, true> k2_part;


    auto direct_short_time = 0;
    auto direct_time = 0;
    auto inverse_short_time = 0;
    auto inverse_time = 0;
    auto check_short_time= 0;
    auto check_time= 0;

    bool use_shortut = argc > 3;

    k2_rrr k2tree;
    std::string fileName = argv[1];
    load_from_file(k2tree, fileName);

    char *list_file = argv[2];
    FILE *list_fp = fopen(list_file, "r");
    uint queries;
    fread(&queries, sizeof(uint), 1, list_fp);
    uint *qry = (uint *) malloc(sizeof(uint) * queries);
    fread(qry, sizeof(uint), queries, list_fp);
    fclose(list_fp);
    /*{
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

            std::cout << "Recovered Nodes:" << recovered << "\n";
            std::cout << "Queries:" << queries << "\n";
            std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
            std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / queries << "\n";
            std::cout << "Time per link: " << duration_cast<nanoseconds>(stop - start).count() / recovered << "\n";

            direct_short_time = duration_cast<nanoseconds>(stop - start).count();
        }
    }*/

    {
        uint i;
        ulong recovered = 0;
        std::vector<uint32_t> result;
        if (use_shortut) {
            std::cout << "Direct with access shortcut2" << std::endl;
            auto start = timer::now();
            for (i = 0; i < queries; i++) {
                k2tree.direct_links_shortcut_2(qry[i], result);
                recovered += result.size();
            }
            auto stop = timer::now();

            std::cout << "Recovered Nodes:" << recovered << "\n";
            std::cout << "Queries:" << queries << "\n";
            std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
            std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / queries << "\n";
            std::cout << "Time per link: " << duration_cast<nanoseconds>(stop - start).count() / recovered << "\n";

            direct_short_time = duration_cast<nanoseconds>(stop - start).count() /recovered;
        }
    }

    {
        uint i;
        ulong recovered = 0;
        std::vector<uint32_t> result;
        std::cout << "Direct without access shortcut" << std::endl;
        auto start = timer::now();
        for (i = 0; i < queries; i++) {
            k2tree.direct_links2(qry[i], result);
            recovered += result.size();
        }
        auto stop = timer::now();

        std::cout << "Recovered Nodes:" << recovered << "\n";
        std::cout << "Queries:" << queries << "\n";
        std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
        std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / queries << "\n";
        std::cout << "Time per link: " << duration_cast<nanoseconds>(stop - start).count() / recovered << "\n";

        direct_time = duration_cast<nanoseconds>(stop - start).count() /recovered;
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

            std::cout << "Recovered Nodes:" << recovered << "\n";
            std::cout << "Queries:" << queries << "\n";
            std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
            std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / queries << "\n";
            std::cout << "Time per link: " << duration_cast<nanoseconds>(stop - start).count() / recovered << "\n";

            inverse_short_time = duration_cast<nanoseconds>(stop - start).count() / recovered;
        }
    }

    {
        uint i;
        ulong recovered = 0;
        std::vector<uint32_t> result;
        std::cout << "Inverse without access shortcut" << std::endl;
        auto start = timer::now();
        for (i = 0; i < queries; i++) {
            k2tree.inverse_links2(qry[i], result);
            recovered += result.size();
        }
        auto stop = timer::now();

        std::cout << "Recovered Nodes:" << recovered << "\n";
        std::cout << "Queries:" << queries << "\n";
        std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
        std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / queries << "\n";
        std::cout << "Time per link: " << duration_cast<nanoseconds>(stop - start).count() / recovered << "\n";

        inverse_time = duration_cast<nanoseconds>(stop - start).count() / recovered;
    }

    srand(0);
    uint link_query_count = 1000000;
    std::vector<std::pair<uint, uint>> check_link_queries(link_query_count);
    for (uint i = 0; i < link_query_count; i++) {
        check_link_queries.push_back(
                std::make_pair(rand() % (k2tree.m_max_element - 1), rand() % (k2tree.m_max_element - 1)));
    }

    {
        if (use_shortut) {
            std::cout << "Performing single link with shortcut" << std::endl;
            auto start = timer::now();
            for (auto pair: check_link_queries) {
                k2tree.check_link_shortcut(pair);
            }
            auto stop = timer::now();

            std::cout << "Queries:" << link_query_count << "\n";
            std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
            std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / link_query_count
                      << "\n";

            check_short_time = duration_cast<nanoseconds>(stop - start).count() / link_query_count;
        }
    }

    {
        std::cout << "Performing single link without shortcut" << std::endl;
        auto start = timer::now();
        for (auto pair: check_link_queries) {
            k2tree.check_link(pair);
        }
        auto stop = timer::now();

        std::cout << "Queries:" << link_query_count << "\n";
        std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
        std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / link_query_count << "\n";

        check_time = duration_cast<nanoseconds>(stop - start).count()/ link_query_count;
    }


    if (use_shortut){
        //Construction Time	Compressed Size (Byte)	Bpe	Direct Short (ns)	Direct (ns)	Inverse Short (ns)	Inverse (ns)	Check S (ns)	Check (ns)
        std::cout << "Hereyougo:" << direct_short_time <<","<< direct_time <<","<< inverse_short_time <<","<< inverse_time <<","<< check_short_time <<","<< check_time << std::endl;
    } else {
        //Construction Time	Compressed Size (Byte)	Bpe	Direct (ns)	Inverse (ns)	Check (ns)
        std::cout << "Hereyougo:" << direct_time <<"," << inverse_time <<"," << check_time << std::endl;
    }




    return 0;
}
