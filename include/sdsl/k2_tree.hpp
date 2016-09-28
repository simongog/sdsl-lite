/* sdsl - succinct data structures library
    Copyright (C) 2014 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file k2_tree.hpp
    \brief k2_tree.hpp contains a compact k^2-tree.
    \author Jan Broß, based on the k2 treap code of Simon Gog
*/
#ifndef INCLUDED_SDSL_K2_TREE
#define INCLUDED_SDSL_K2_TREE

#include "k2_tree_base.hpp"
#include <gtest/gtest_prod.h>

//! Namespace for the succinct data structure library.
namespace sdsl {

/*! A hybrid k2 tree implementation
 *  \par References
 *       [1] Nieves R. Brisaboa, Susana Ladra, and Gonzalo Navarro:
 *           Compact representation of Web graphs with extended functionality.
 *           Inf. Syst. 39 (January 2014), 152-174. DOI=http://dx.doi.org/10.1016/j.is.2013.08.003
 */
    template<uint8_t t_k,
            typename t_lev=bit_vector,
            typename t_leaf=bit_vector,
            typename t_rank=typename t_lev::rank_1_type>

    class k2_tree : public k2_tree_base<t_k, t_k, t_lev, t_leaf, t_rank> {
        static_assert(t_k <= 8, "t_k can at most be 8 because of the current dac compression implementation.");//int_vectors support only 64Bit
        static_assert(t_k == 2 || t_k == 4 || t_k == 8 || t_k == 16, "t_k_l_1 has to one of 2,4,8,16");

        FRIEND_TEST(K2TreeInternalTest, testZOrderSort);
        FRIEND_TEST(K2TreeInternalTest, testZOrderSort2);
        FRIEND_TEST(K2TreeInternalTest, testZOrder100);
        FRIEND_TEST(K2TreeInternalTest, test_calculate_subtree_number_and_new_relative_coordinates);
        FRIEND_TEST(K2TreeInternalTest, test_access_shortcut);

    public:
        typedef int_vector<>::size_type size_type;

        typedef stxxl::VECTOR_GENERATOR<std::pair<uint32_t, uint32_t>>::result stxxl_32bit_pair_vector;
        typedef stxxl::VECTOR_GENERATOR<std::pair<uint64_t, uint64_t>>::result stxxl_64bit_pair_vector;

        using node_type = k2_treap_ns::node_type;
        using point_type = k2_treap_ns::point_type;
        using t_p = k2_treap_ns::t_p;

        using k2_tree_base<t_k, t_k, t_lev, t_leaf, t_rank>::operator=;
        using k2_tree_base<t_k, t_k, t_lev, t_leaf, t_rank>::operator==;

        enum {
            k = t_k
        };

        k2_tree() = default;

        /*
        k2_tree(const k2_tree &tr)
                : k2_tree_base<t_k, t_k, t_lev, t_leaf, t_rank>(tr) {
            *this = tr;
        }

        k2_tree(k2_tree &&tr)
        : k2_tree_base<t_k, t_k, t_lev, t_leaf, t_rank>(tr) {
                *this = std::move(tr);
        }*/

        template<typename t_vector>
        k2_tree(t_vector &v, construction_algorithm construction_algo, uint64_t max_hint = 0, std::string temp_file_prefix="") {

            using namespace k2_treap_ns;
            if (v.size() > 0) {
                if (max_hint == 0) {
                    max_hint = get_maximum(v);
                }
                this->m_tree_height = get_tree_height(max_hint);
                construct(v, construction_algo, temp_file_prefix);
                this->post_init();
            }
        }

        k2_tree(int_vector_buffer<> &buf_x,
                int_vector_buffer<> &buf_y, construction_algorithm construction_algo, uint64_t max_hint = 0) {
            using namespace k2_treap_ns;

            if (buf_x.size() == 0) {
                return;
            }

            typedef int_vector_buffer<> *t_buf_p;
            if (buf_x.size() == 0) {
                return;
            }

            std::vector<t_buf_p> bufs = {&buf_x, &buf_y};

            uint64_t max;
            if (max_hint != 0) {
                max = max_hint;//temporarily set
            } else {
                auto max_element = [](int_vector_buffer<> &buf) {
                    uint64_t max_val = 0;
                    for (auto val : buf) {
                        max_val = std::max((uint64_t) val, max_val);
                    }
                    return max_val;
                };

                auto max_buf_element = [&]() {
                    uint64_t max_v = 0;
                    for (auto buf : bufs) {
                        uint64_t _max_v = max_element(*buf);
                        max_v = std::max(max_v, _max_v);
                    }
                    return max_v;
                };

                max = max_buf_element();
            };

            uint8_t res = 0;
            while (res <= 64 and precomp<t_k>::exp(res) <= max) { ++res; }
            if (res == 65) {
                throw std::logic_error("Maximal element of input is too big.");
            }

            this->m_tree_height = res;
            this->m_max_element = precomp<t_k>::exp(res);

            if (precomp<t_k>::exp(res) <= std::numeric_limits<uint32_t>::max()) {
                auto v = read<uint32_t, uint32_t>(bufs);
                construct(v, construction_algo, buf_x.filename());
            } else {
                auto v = read<uint64_t, uint64_t>(bufs);
                construct(v, construction_algo, buf_x.filename());
            }

            this->post_init();
        }
        /* Accesor methods for links, duplicated because no virtual template funtions possible */
        bool check_link(std::pair<uint, uint> link) const override {
            return this->check_link_internal(link, &precomp<t_k>::divexp, &precomp<t_k>::modexp);
        }
        bool check_link_shortcut(std::pair<uint, uint> link) const override {
            return this->check_link_shortcut_internal(link, &precomp<t_k>::divexp, &precomp<t_k>::modexp);
        }
        void inverse_links2(uint target_id, std::vector<uint> &result) const override {
            this->inverse_links2_internal(target_id, result, &precomp<t_k>::divexp, &precomp<t_k>::modexp, &precomp<t_k>::multexp);
        }
        void inverse_links_shortcut(uint target_id, std::vector<uint> &result) const override {
            this->inverse_links_shortcut_internal(target_id, result, &precomp<t_k>::divexp, &precomp<t_k>::modexp, &precomp<t_k>::multexp);
        }
        void direct_links_shortcut_2(uint source_id, std::vector<uint> &result) const override {
            this->direct_links_shortcut_2_internal(source_id, result, &precomp<t_k>::divexp, &precomp<t_k>::modexp, &precomp<t_k>::multexp);
        }
        void direct_links_shortcut(uint source_id, std::vector<uint> &result) const override {
            this->direct_links_shortcut_internal(source_id, result, &precomp<t_k>::divexp, &precomp<t_k>::modexp, &precomp<t_k>::multexp);
        }
        void direct_links2(uint source_id, std::vector<uint> &result) const override {
            this->direct_links2_internal(source_id, result, &precomp<t_k>::divexp, &precomp<t_k>::modexp, &precomp<t_k>::multexp);
        }

        bool check_link(std::pair<uint64_t, uint64_t> link) const override {
            return this->check_link_internal(link, &precomp<t_k>::divexp, &precomp<t_k>::modexp);
        }
        bool check_link_shortcut(std::pair<uint64_t, uint64_t> link) const override {
            return this->check_link_shortcut_internal(link, &precomp<t_k>::divexp, &precomp<t_k>::modexp);
        }
        void inverse_links2(uint64_t target_id, std::vector<uint64_t> &result) const override {
            this->inverse_links2_internal(target_id, result, &precomp<t_k>::divexp, &precomp<t_k>::modexp, &precomp<t_k>::multexp);
        }
        void inverse_links_shortcut(uint64_t target_id, std::vector<uint64_t> &result) const override {
            this->inverse_links_shortcut_internal(target_id, result, &precomp<t_k>::divexp, &precomp<t_k>::modexp, &precomp<t_k>::multexp);
        }
        void direct_links_shortcut_2(uint64_t source_id, std::vector<uint64_t> &result) const override {
            this->direct_links_shortcut_2_internal(source_id, result, &precomp<t_k>::divexp, &precomp<t_k>::modexp, &precomp<t_k>::multexp);
        }
        void direct_links_shortcut(uint64_t source_id, std::vector<uint64_t> &result) const override {
            this->direct_links_shortcut_internal(source_id, result, &precomp<t_k>::divexp, &precomp<t_k>::modexp, &precomp<t_k>::multexp);
        }
        void direct_links2(uint64_t source_id, std::vector<uint64_t> &result) const override {
            this->direct_links2_internal(source_id, result, &precomp<t_k>::divexp, &precomp<t_k>::modexp, &precomp<t_k>::multexp);
        }


        inline uint8_t get_k(uint8_t) const {
            return t_k;
        }

        uint word_size() const override  {
            return div_ceil((uint) t_k * t_k, kUcharBits);
        }

        size_t words_count() const override {
            if (this->m_tree_height == 0) {
                return 0;
            }

            return this->m_leaves.size() / t_k / t_k;
        }

        void load(std::istream &in) override {
            k2_tree_base<t_k, t_k, t_lev, t_leaf, t_rank>::load(in);
            if (this->m_tree_height > 0) {
                if (this->m_access_shortcut_size > 0) {
                    this->perform_access_shortcut_precomputations();
                }
            }
        }

        std::string get_type_string_without_compression() const {
            return "k2_tree<"+std::to_string(t_k)+">";
        }

        std::string get_type_string() const {
            return "k2_tree<"+std::to_string(t_k)+","+get_compression_name(this->m_used_compression)+">";
        }

        //hack a the moment, because construct cannot be virtual
        void load_from_ladrabin(std::string fileName, construction_algorithm construction_algo = COUNTING_SORT, uint8_t access_shortcut_size = 0, std::string temp_file_prefix = "") {
            using namespace k2_treap_ns;
            if (!has_ending(fileName, ".ladrabin")) {
                fileName.append(".ladrabin");
                std::cout << "Appending .ladrabin to filename as file has to be in .ladrabin format" << std::endl;
            }

            std::fstream fileStream(fileName, std::ios_base::in);

            if (fileStream.is_open()) {
                uint number_of_nodes;
                ulong number_of_edges;

                read_member(number_of_nodes, fileStream);
                read_member(number_of_edges, fileStream);

                this->m_max_element = number_of_nodes - 1;
                this->m_size = number_of_edges;
                this->m_tree_height = get_tree_height(this->m_max_element);

                uint nodes_read = 0;
                uint source_id;
                int target_id;

                std::vector<std::pair<uint, uint>> coords;
                //stxxl_32bit_pair_vector coords;
                coords.reserve(number_of_edges);
                for (uint64_t i = 0; i < number_of_nodes + number_of_edges; i++) {
                    read_member(target_id, fileStream);
                    if (target_id < 0) {
                        nodes_read++;
                    } else {
                        source_id = nodes_read - 1;
                        coords.push_back(std::make_pair(source_id, target_id));
                    }
                }
                fileStream.close();

                if (coords.size() > 0) {
                    construct(coords, construction_algo, temp_file_prefix);
                    coords.clear();

                    this->m_access_shortcut_size = access_shortcut_size;
                    this->post_init();
                }
            } else {
                throw std::runtime_error("Could not load ladrabin file");
            }
        }

    private:
        //FIXME: declared here and in k2_tree as a workaround, because virtual template methods are not possible
        template<typename t_vector>
        void construct(t_vector &v, const construction_algorithm construction_algo,
                       const std::string &temp_file_prefix = 0) {
            std::cout << "Constructing using " << get_construction_name(construction_algo) << std::endl;
            switch (construction_algo){
                case COUNTING_SORT:
                    this->construct_counting_sort(v);
                    break;
                case PARTITIONBASED:
                    this->construct(v, temp_file_prefix);
                    break;
                case ZORDERSORT:
                    construct_by_z_order_sort_internal(v, temp_file_prefix);
                    break;
            }
        }

        template<typename t_vector>
        void construct_counting_sort(t_vector &links) {
            this->construct_counting_sort_internal(links, &precomp<t_k>::divexp, &precomp<t_k>::exp);
        }

        template<typename t_vector>
        void construct(t_vector &links, std::string temp_file_prefix = "") {
            this->construct_internal(links, &precomp<t_k>::divexp, temp_file_prefix);
        }

        /**
         * Constructs the tree corresponding to the points in the links vector inpace by performing a z order sort and subsequently constructing the tree top down
         * @param links
         * @param temp_file_prefix
         */
        template<typename t_vector>
        void
        construct_by_z_order_sort_internal(t_vector &links, std::string temp_file_prefix = "") {
            using namespace k2_treap_ns;
            typedef decltype(links[0].first) t_x;
            typedef decltype(links[0].second) t_y;

            using t_e = std::pair<t_x, t_y>;

            this->m_size = links.size();

            std::string id_part = util::to_string(util::pid())
                                  + "_" + util::to_string(util::id());

            std::cout << "Sorting By Z Order" << std::endl;
            __gnu_parallel::sort(links.begin(), links.end(), [&](const t_e &lhs, const t_e &rhs) {
                return interleave<k>::bits(lhs) < interleave<k>::bits(rhs);
            });

            std::cout << "Sorting Finished, Constructing Bitvectors" << std::endl;

            /*for (int m = 0; m < links.size(); ++m) {
                std::cout << links[m].first << "," << links[m].second << std::endl;
            }*/
            
            
            std::vector<t_x> previous_subtree_number(this->m_tree_height, -1);

            {
                t_x subtree_distance;
                bool fill_to_k2_entries = false; //begin extra case!
                std::vector<uint> gap_to_k2(this->m_tree_height, k * k);
                bool firstLink = true;
                t_x current_subtree_number = 0;

                std::vector<int_vector_buffer<1>> level_buffers = this->create_level_buffers(temp_file_prefix, id_part);

                std::pair<t_x, t_y> previous_link;
                for (auto &current_link: links) {
                    auto tmp = std::make_pair(current_link.first, current_link.second);
                    for (uint current_level = 0; current_level < this->m_tree_height; ++current_level) {
                        current_subtree_number = calculate_subtree_number_and_new_relative_coordinates(current_link,
                                                                                                       current_level);
                        subtree_distance = current_subtree_number - previous_subtree_number[current_level];

                        if (subtree_distance > 0) {
                            //invalidate previous subtree numbers as new relative frame
                            for (uint i = current_level + 1; i < this->m_tree_height; ++i) {
                                previous_subtree_number[i] = -1;
                            }

                            if (fill_to_k2_entries && current_level != 0) {
                                for (uint j = 0; j < gap_to_k2[current_level]; ++j) {
                                    level_buffers[current_level].push_back(0);
                                }
                                gap_to_k2[current_level] = k * k;
                            }

                            for (uint j = 0; j < subtree_distance - 1; ++j) {
                                level_buffers[current_level].push_back(0);
                                gap_to_k2[current_level]--;
                            }

                            level_buffers[current_level].push_back(1);
                            gap_to_k2[current_level]--;

                            if (!firstLink)
                                fill_to_k2_entries = true;
                        } else if (subtree_distance == 0) {
                            fill_to_k2_entries = false;
                        } else {
                            std::string error_message(
                                    "negative subtree_distance after z_order sort is not possible, somethings wrong current_level=" +
                                    std::to_string(current_level) + " subtree_distance=" +
                                    std::to_string(subtree_distance) +
                                    " current_subtree_number=" + std::to_string(current_subtree_number) +
                                    " previous_subtree_number[current_level]=" +
                                    std::to_string(previous_subtree_number[current_level]) + "current_link=" +
                                    std::to_string(current_link.first) + "," + std::to_string(current_link.second) +
                                    "previous_link=" + std::to_string(previous_link.first) + "," +
                                    std::to_string(previous_link.second));
                            throw std::logic_error(error_message);
                        }
                        //std::cout << "Setting previous_subtree_number[" << current_level << "] = "<< current_subtree_number << std::endl;
                        previous_subtree_number[current_level] = current_subtree_number;
                    }
                    //FIXME: special case treatment for last level (doesn't need to be sorted --> set corresponding bit, but don't append)
                    firstLink = false;
                    previous_link = tmp;
                }

                //fill rest with 0s
                for (uint l = 0; l < gap_to_k2.size(); ++l) {
                    for (uint i = 0; i < gap_to_k2[l]; ++i) {
                        level_buffers[l].push_back(0);
                    }
                }
            }

            this->load_vectors_from_file(temp_file_prefix, id_part);
        }

        /**
         * Calculates the corresponding subtree of link on a given level as well as
         * the new relative coordinates (relative to the upper left corner of the corresponding submatrix)
         * of link on the next level
         *
         * @param link
         *      in:  current coordinates
         *      out: relative coordinates on level "level"
         * @param level
         * @return
         */
        template<typename t_x, typename t_y>
        t_x inline calculate_subtree_number_and_new_relative_coordinates(std::pair<t_x, t_y> &link, int level) {
            using namespace k2_treap_ns;
            t_x exponent = this->m_tree_height - level - 1;
            t_x result = k * precomp<t_k>::divexp(link.first, exponent) + precomp<t_k>::divexp(link.second, exponent);
            link.first = precomp<t_k>::modexp(link.first, exponent);
            link.second = precomp<t_k>::modexp(link.second, exponent);

            return result;
        }

        uint8_t get_tree_height(const uint64_t max) {
            uint8_t res = 0;
            while (precomp<t_k>::exp(res) <= max) { ++res; }

            this->m_max_element = precomp<t_k>::exp(res);
            return res;
        }
    };
}
#endif