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
        long direct_short_recovered = 0;
        long direct_time = 0;
        long direct_recovered = 0;
        long inverse_short_time = 0;
        long inverse_short_recovered = 0;
        long inverse_time = 0;
        long inverse_recovered = 0;
        long check_short_time= 0;
        long check_short_present=0;
        long check_time= 0;
        long check_present=0;
    };

    template <typename k2tree>
    void store_to_file(std::string file_name, k2tree& tree){
        store_to_file(tree, file_name);
        write_structure<R_FORMAT>(tree, file_name + "(" + tree.get_type_string() + ")" + ".html");
    }

    bool has_ending(std::string const &fullString, std::string const &ending) {
        if (fullString.length() >= ending.length()) {
            return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
        } else {
            return false;
        }
    }

    void load_single_link_query_file(std::string basic_string, vector<std::pair<uint, uint>>& vector);

    template <typename k2tree>
    access_times perform_speed_test(std::string query_file, k2tree& tree, std::string single_link_query_file = "", bool use_shortcut = false){
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
                //std::cout << "Direct with access shortcut2" << std::endl;
                auto start = timer::now();
                for (i = 0; i < queries; i++) {
                    tree.direct_links_shortcut_2(qry[i], result);
                    recovered += result.size();
                }
                auto stop = timer::now();

                times.direct_short_recovered = recovered;
                times.direct_short_time = duration_cast<nanoseconds>(stop - start).count() /recovered;
            }
        }

        {
            uint i;
            ulong recovered = 0;
            std::vector<uint32_t> result;
            //std::cout << "Direct without access shortcut" << std::endl;
            auto start = timer::now();
            for (i = 0; i < queries; i++) {
                tree.direct_links2(qry[i], result);
                recovered += result.size();
            }
            auto stop = timer::now();

            times.direct_recovered = recovered;
            times.direct_time = duration_cast<nanoseconds>(stop - start).count() /recovered;
        }

        {
            uint i;
            ulong recovered = 0;
            std::vector<uint32_t> result;
            if (use_shortcut) {
                //std::cout << "Inverse with access shortcut" << std::endl;
                auto start = timer::now();
                for (i = 0; i < queries; i++) {
                    tree.inverse_links_shortcut(qry[i], result);
                    recovered += result.size();
                }
                auto stop = timer::now();

                times.inverse_short_recovered = recovered;
                times.inverse_short_time = duration_cast<nanoseconds>(stop - start).count() / recovered;
            }
        }

        {
            uint i;
            ulong recovered = 0;
            std::vector<uint32_t> result;
            //std::cout << "Inverse without access shortcut" << std::endl;
            auto start = timer::now();
            for (i = 0; i < queries; i++) {
                tree.inverse_links2(qry[i], result);
                recovered += result.size();
            }
            auto stop = timer::now();

            times.inverse_recovered = recovered;
            times.inverse_time = duration_cast<nanoseconds>(stop - start).count() / recovered;
        }

        srand(0);
        std::vector<std::pair<uint, uint>> check_link_queries;
        uint link_query_count = 1000000;
        if (!single_link_query_file.empty()){
            load_single_link_query_file(single_link_query_file, check_link_queries);
        }

        if (check_link_queries.empty()){
            std::cout << "Using random link queries" << std::endl;
            check_link_queries.resize(link_query_count);
            for (uint i = 0; i < link_query_count; i++) {
                check_link_queries.push_back(
                        std::make_pair(rand() % (tree.m_max_element - 1), rand() % (tree.m_max_element - 1)));
            }
        }
        link_query_count = check_link_queries.size();

        {
            if (use_shortcut) {
                //std::cout << "Performing single link with shortcut" << std::endl;
                uint present_links = 0;
                auto start = timer::now();
                for (auto pair: check_link_queries) {
                    present_links += tree.check_link_shortcut(pair);
                }
                auto stop = timer::now();

                times.check_short_time = duration_cast<nanoseconds>(stop - start).count() / link_query_count;
                times.check_short_present = present_links;
            }
        }

        {
            uint present_links = 0;
            //std::cout << "Performing single link without shortcut" << std::endl;
            auto start = timer::now();
            for (auto pair: check_link_queries) {
                present_links += tree.check_link(pair);
            }
            auto stop = timer::now();

            times.check_time = duration_cast<nanoseconds>(stop - start).count()/ link_query_count;
            times.check_present = present_links;
        }

        return times;
    }

    void load_single_link_query_file(std::string query_file, vector<std::pair<uint, uint>>& queries) {
	std::cout << "Attempting to load " << query_file << std::endl;
        std::ifstream input_file(query_file, std::ios_base::in);
        if (!input_file.is_open()) {
            return;
        }

       int source, target, count = 0;
       read_member(count, input_file);

        for (int i = 0; i < count; ++i) {
            read_member(source, input_file);
            read_member(target, input_file);
            queries.push_back(std::make_pair(source, (uint)target));
        }
    }
}
#endif //SDSL_K2_TREE_UTILITY_H
