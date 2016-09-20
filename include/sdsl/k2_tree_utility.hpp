//
// Created by d056848 on 9/12/16.
//

#ifndef SDSL_K2_TREE_UTILITY_H
#define SDSL_K2_TREE_UTILITY_H


#include <string>
#include <sdsl/mem_monitor.hpp>
#include <sdsl/config.hpp>
#include <sdsl/io.hpp>

namespace  sdsl {

    using namespace std::chrono;
    using timer = std::chrono::high_resolution_clock;

    struct access_times {
        long direct_short_time = 0;
        long direct_time = 0;
        long inverse_short_time = 0;
        long inverse_time = 0;
        long check_short_time= 0;
        long check_time= 0;
    };

    template <typename k2tree>
    void store_to_file(std::string file_name, k2tree& tree){
        store_to_file(tree, file_name);
        write_structure<HTML_FORMAT>(tree, file_name + "(" + tree.get_type_string() + ")" + ".html");
    }

    template <typename k2tree>
    access_times perform_speed_test(std::string query_file, k2tree& tree, bool use_shortcut = false){
        if(!has_ending(query_file, ".queries")){
            query_file.append(".queries");
            std::cout << "Appending .queries to filename as file has to be in .queries format (binary format uint n number of queries followed by n uints representing the node id to query" << std::endl;
        }
        FILE *list_fp = fopen(query_file.c_str(), "r");
        uint queries;
        fread(&queries, sizeof(uint), 1, list_fp);
        uint *qry = (uint *) malloc(sizeof(uint) * queries);
        fread(qry, sizeof(uint), queries, list_fp);
        fclose(list_fp);

        struct access_times times;

        {
            uint i;
            ulong recovered = 0;
            std::vector<uint32_t> result;
            if (use_shortcut) {
                std::cout << "Direct with access shortcut2" << std::endl;
                auto start = timer::now();
                for (i = 0; i < queries; i++) {
                    tree.direct_links_shortcut_2(qry[i], result);
                    recovered += result.size();
                }
                auto stop = timer::now();

                std::cout << "Recovered Nodes:" << recovered << "\n";
                std::cout << "Queries:" << queries << "\n";
                std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
                std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / queries << "\n";
                std::cout << "Time per link: " << duration_cast<nanoseconds>(stop - start).count() / recovered << "\n";

                times.direct_short_time = duration_cast<nanoseconds>(stop - start).count() /recovered;
            }
        }

        {
            uint i;
            ulong recovered = 0;
            std::vector<uint32_t> result;
            std::cout << "Direct without access shortcut" << std::endl;
            auto start = timer::now();
            for (i = 0; i < queries; i++) {
                tree.direct_links2(qry[i], result);
                recovered += result.size();
            }
            auto stop = timer::now();

            std::cout << "Recovered Nodes:" << recovered << "\n";
            std::cout << "Queries:" << queries << "\n";
            std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
            std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / queries << "\n";
            std::cout << "Time per link: " << duration_cast<nanoseconds>(stop - start).count() / recovered << "\n";

            times.direct_time = duration_cast<nanoseconds>(stop - start).count() /recovered;
        }

        {
            uint i;
            ulong recovered = 0;
            std::vector<uint32_t> result;
            if (use_shortcut) {
                std::cout << "Inverse with access shortcut" << std::endl;
                auto start = timer::now();
                for (i = 0; i < queries; i++) {
                    tree.inverse_links_shortcut(qry[i], result);
                    recovered += result.size();
                }
                auto stop = timer::now();

                std::cout << "Recovered Nodes:" << recovered << "\n";
                std::cout << "Queries:" << queries << "\n";
                std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
                std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / queries << "\n";
                std::cout << "Time per link: " << duration_cast<nanoseconds>(stop - start).count() / recovered << "\n";

                times.inverse_short_time = duration_cast<nanoseconds>(stop - start).count() / recovered;
            }
        }

        {
            uint i;
            ulong recovered = 0;
            std::vector<uint32_t> result;
            std::cout << "Inverse without access shortcut" << std::endl;
            auto start = timer::now();
            for (i = 0; i < queries; i++) {
                tree.inverse_links2(qry[i], result);
                recovered += result.size();
            }
            auto stop = timer::now();

            std::cout << "Recovered Nodes:" << recovered << "\n";
            std::cout << "Queries:" << queries << "\n";
            std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
            std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / queries << "\n";
            std::cout << "Time per link: " << duration_cast<nanoseconds>(stop - start).count() / recovered << "\n";

            times.inverse_time = duration_cast<nanoseconds>(stop - start).count() / recovered;
        }

        srand(0);
        uint link_query_count = 1000000;
        std::vector<std::pair<uint, uint>> check_link_queries(link_query_count);
        for (uint i = 0; i < link_query_count; i++) {
            check_link_queries.push_back(
                    std::make_pair(rand() % (tree.m_max_element - 1), rand() % (tree.m_max_element - 1)));
        }

        {
            if (use_shortcut) {
                std::cout << "Performing single link with shortcut" << std::endl;
                auto start = timer::now();
                for (auto pair: check_link_queries) {
                    tree.check_link_shortcut(pair);
                }
                auto stop = timer::now();

                std::cout << "Queries:" << link_query_count << "\n";
                std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
                std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / link_query_count
                          << "\n";

                times.check_short_time = duration_cast<nanoseconds>(stop - start).count() / link_query_count;
            }
        }

        {
            std::cout << "Performing single link without shortcut" << std::endl;
            auto start = timer::now();
            for (auto pair: check_link_queries) {
                tree.check_link(pair);
            }
            auto stop = timer::now();

            std::cout << "Queries:" << link_query_count << "\n";
            std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
            std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / link_query_count << "\n";

            times.check_time = duration_cast<nanoseconds>(stop - start).count()/ link_query_count;
        }

        return times;
    }
}
#endif //SDSL_K2_TREE_UTILITY_H
