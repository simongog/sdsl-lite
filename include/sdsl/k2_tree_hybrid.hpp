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
/*! \file hybrid_k2_tree.hpp
    \brief hybrid_k2_tree.hpp contains a compact hybrid k^2-tree.
    \author Jan Broß, based on the k2 treap code of Simon Gog, leaf compression is based on the libk2tree implemenetation which uses DACs impelemented by Brisaboa, Ladra et.al.
*/
#ifndef INCLUDED_SDSL_HYBRID_K2_TREE
#define INCLUDED_SDSL_HYBRID_K2_TREE

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
    template<uint8_t t_k_l_1,
            uint8_t t_k_l_1_size,
            uint8_t t_k_l_2,
            uint8_t t_k_leaves,
            typename t_lev=bit_vector,
            typename t_leaf=bit_vector,
            typename t_rank=typename t_lev::rank_1_type>

    class k2_tree_hybrid : public k2_tree_base<t_k_l_1, t_k_leaves, t_lev, t_leaf, t_rank> {
        static_assert(t_k_l_1 == 2    || t_k_l_1 == 4    || t_k_l_1 == 8    || t_k_l_1 == 16,    "t_k_l_1 has to one of 2,4,8,16");
        static_assert(t_k_l_2 == 2    || t_k_l_2 == 4    || t_k_l_2 == 8    || t_k_l_2 == 16,    "t_k_l_2 has to one of 2,4,8,16");
        static_assert(t_k_leaves == 2 || t_k_leaves == 4 || t_k_leaves == 8 || t_k_leaves == 16, "t_k_leaves has to one of 2,4,8,16");
        static_assert(t_k_leaves >= t_k_l_1,
                      "t_k_leaves has to be larger than t_k_l_1,  otherwise this could lead to different word sizes and thus to a problem for the k2part approach"); //if smaller than t_k_l_1 it could be that t_k_leaves is not used
        static_assert(t_k_leaves <= 8, "t_k can at most be 8 because of the current dac compression implementation.");

    private:

        std::vector<uint8_t> m_k_for_level;
        std::vector<uint8_t> m_shift_table; //for fast div and mod operations by submatrix size

    public:
        typedef stxxl::VECTOR_GENERATOR<std::pair<uint32_t, uint32_t>>::result stxxl_32bit_pair_vector;
        typedef stxxl::VECTOR_GENERATOR<std::pair<uint64_t, uint64_t>>::result stxxl_64bit_pair_vector;
        typedef int_vector<>::size_type size_type;

        k2_tree_hybrid() = default;

        k2_tree_hybrid(const k2_tree_hybrid &tr)
                : k2_tree_base<t_k_l_1, t_k_leaves, t_lev, t_leaf, t_rank>(tr) {
            *this = tr;
        }

        k2_tree_hybrid(k2_tree_hybrid &&tr)
                : k2_tree_base<t_k_l_1, t_k_leaves, t_lev, t_leaf, t_rank>(tr) {
            *this = std::move(tr);
        }

        template<typename t_vector>
        k2_tree_hybrid(t_vector &v, construction_algorithm construction_algo, uint64_t max_hint = 0, std::string temp_file_prefix = "") {

            using namespace k2_treap_ns;
            if (v.size() > 0) {
                if (max_hint == 0) {
                    max_hint = get_maximum(v);
                }
                this->m_tree_height = get_tree_height(max_hint);
                initialize_shift_table();

                contruct(v, construction_algo, temp_file_prefix);

                this->post_init();
            }
        }

        k2_tree_hybrid(int_vector_buffer<> &buf_x,
                       int_vector_buffer<> &buf_y, construction_algorithm construction_algo, uint64_t max_hint = 0)  {
            using namespace k2_treap_ns;
            typedef int_vector_buffer<> *t_buf_p;

            if (buf_x.size() == 0){
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

            this->m_tree_height = get_tree_height(max);
            initialize_shift_table();

            if (this->m_max_element <= std::numeric_limits<uint32_t>::max()) {
                auto v = read<uint32_t, uint32_t>(bufs);
                contruct(v, construction_algo, buf_x.filename());
            } else {
                auto v = read<uint64_t, uint64_t>(bufs);
                contruct(v, construction_algo, buf_x.filename());
            }

            this->post_init();
        }

        /* Accesor methods for links, duplicated because no virtual template funtions possible */
        bool check_link(std::pair<uint, uint> link) const override {
            return this->check_link_internal(link,
                                      [&](uint64_t x, uint8_t l){return divexp(x,l);},
                                      [&](uint64_t x, uint8_t l){return modexp(x, l);});
        }
        bool check_link_shortcut(std::pair<uint, uint> link) const override {
            return this->check_link_shortcut_internal(link,
                                               [&](uint x, uint8_t l){return divexp(x,l);},
                                               [&](uint x, uint8_t l){return modexp(x, l);});
        }
        void inverse_links2(uint target_id, std::vector<uint> &result) const override {
            this->inverse_links2_internal(target_id, result,
                                          [&](uint x, uint8_t l){return divexp(x,l);},
                                          [&](uint x, uint8_t l){return modexp(x, l);},
                                          [&](uint x, uint8_t l){return multexp(x, l);});
        }
        void inverse_links_shortcut(uint target_id, std::vector<uint> &result) const override {
            this->inverse_links_shortcut_internal(target_id, result,
                                                  [&](uint x, uint8_t l){return divexp(x,l);},
                                                  [&](uint x, uint8_t l){return modexp(x, l);},
                                                  [&](uint x, uint8_t l){return multexp(x, l);});
        }
        void direct_links_shortcut_2(uint source_id, std::vector<uint> &result) const override {
            this->direct_links_shortcut_2_internal(source_id, result,
                                                   [&](uint x, uint8_t l){return divexp(x,l);},
                                                   [&](uint x, uint8_t l){return modexp(x, l);},
                                                   [&](uint x, uint8_t l){return multexp(x, l);});
        }
        void direct_links_shortcut(uint source_id, std::vector<uint> &result) const override {
            this->direct_links_shortcut_internal(source_id, result,
                                                 [&](uint x, uint8_t l){return divexp(x,l);},
                                                 [&](uint x, uint8_t l){return modexp(x, l);},
                                                 [&](uint x, uint8_t l){return multexp(x, l);});
        }
        void direct_links2(uint source_id, std::vector<uint> &result) const override {
            this->direct_links2_internal(source_id, result,
                                         [&](uint x, uint8_t l){return divexp(x,l);},
                                         [&](uint x, uint8_t l){return modexp(x, l);},
                                         [&](uint x, uint8_t l){return multexp(x, l);});
        }

        bool check_link(std::pair<uint64_t, uint64_t> link) const override {
            return this->check_link_internal(link,
                                      [&](uint64_t x, uint8_t l){return divexp(x,l);},
                                      [&](uint64_t x, uint8_t l){return modexp(x, l);});
        }
        bool check_link_shortcut(std::pair<uint64_t, uint64_t> link) const override {
            return this->check_link_shortcut_internal(link,
                                               [&](uint64_t x, uint8_t l){return divexp(x,l);},
                                               [&](uint64_t x, uint8_t l){return modexp(x, l);});
        }
        void inverse_links2(uint64_t target_id, std::vector<uint64_t> &result) const override {
            this->inverse_links2_internal(target_id, result,
                                          [&](uint64_t x, uint8_t l){return divexp(x,l);},
                                          [&](uint64_t x, uint8_t l){return modexp(x, l);},
                                          [&](uint64_t x, uint8_t l){return multexp(x, l);});
        }
        void inverse_links_shortcut(uint64_t target_id, std::vector<uint64_t> &result) const override {
            this->inverse_links_shortcut_internal(target_id, result,
                                                  [&](uint64_t x, uint8_t l){return divexp(x,l);},
                                                  [&](uint64_t x, uint8_t l){return modexp(x, l);},
                                                  [&](uint64_t x, uint8_t l){return multexp(x, l);});
        }
        void direct_links_shortcut_2(uint64_t source_id, std::vector<uint64_t> &result) const override {
            this->direct_links_shortcut_2_internal(source_id, result,
                                                   [&](uint64_t x, uint8_t l){return divexp(x,l);},
                                                   [&](uint64_t x, uint8_t l){return modexp(x, l);},
                                                   [&](uint64_t x, uint8_t l){return multexp(x, l);});
        }
        void direct_links_shortcut(uint64_t source_id, std::vector<uint64_t> &result) const override {
            this->direct_links_shortcut_internal(source_id, result,
                                                 [&](uint64_t x, uint8_t l){return divexp(x,l);},
                                                 [&](uint64_t x, uint8_t l){return modexp(x, l);},
                                                 [&](uint64_t x, uint8_t l){return multexp(x, l);});
        }
        void direct_links2(uint64_t source_id, std::vector<uint64_t> &result) const override {
            this->direct_links2_internal(source_id, result,
                                         [&](uint64_t x, uint8_t l){return divexp(x,l);},
                                         [&](uint64_t x, uint8_t l){return modexp(x, l);},
                                         [&](uint64_t x, uint8_t l){return multexp(x, l);});
        }

        size_type serialize(std::ostream &out, structure_tree_node *v, std::string name) const override {
            return k2_tree_base<t_k_l_1, t_k_leaves, t_lev, t_leaf, t_rank>::serialize(out, v, name);
        }

        uint word_size() const override {
            return div_ceil((uint) t_k_leaves * t_k_leaves, kUcharBits);
        }


        size_t words_count() const override {
            return this->m_leaves.size() / t_k_leaves / t_k_leaves;
        }

        void load(std::istream &in) override {
            k2_tree_base<t_k_l_1, t_k_leaves, t_lev, t_leaf, t_rank>::load(in);
            if (this->m_tree_height > 0) {
                for (int i = 1; i <= std::min(t_k_l_1_size, (uint8_t) (this->m_tree_height - 1)); ++i) {
                    m_k_for_level.push_back(t_k_l_1);
                }

                for (int j = t_k_l_1_size + 1; j <= (this->m_tree_height - 1); ++j) {
                    m_k_for_level.push_back(t_k_l_2);
                }

                m_k_for_level.push_back(t_k_leaves);

                initialize_shift_table();

                if (this->m_access_shortcut_size > 0) {
                    this->perform_access_shortcut_precomputations();
                }
            }
        }

        //hack a the moment, because construct cannot be virtual
        void load_from_ladrabin(std::string fileName, construction_algorithm construction_algo = COUNTING_SORT, uint8_t access_shortcut_size = 0, std::string temp_file_prefix = "") {
            using namespace k2_treap_ns;
            if (!has_ending(fileName, ".ladrabin")) {
                fileName.append(".ladrabin");
                std::cout << "Appending .graph-txt to filename as file has to be in .ladrabin format" << std::endl;
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
                initialize_shift_table();

                uint nodes_read = 0;
                uint source_id;
                int target_id;

                std::vector<std::pair<uint, uint>> coords;
                coords.reserve(number_of_edges);
                /*typedef stxxl::VECTOR_GENERATOR<pair<uint32_t, uint32_t>>::result stxxl_pair_vector;
                stxxl_pair_vector coords(number_of_nodes);*/
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
                    contruct(coords, construction_algo, temp_file_prefix);

                    std::cout << "Finished Construction" << std::endl;

                    coords.clear();

                    this->m_access_shortcut_size = access_shortcut_size;
                    this->post_init();
                }
            } else {
                throw std::runtime_error("Could not load ladrabin file");
            }
        }

        k2_tree_hybrid &operator=(k2_tree_hybrid &&tr) {
            k2_tree_base<t_k_l_1, t_k_leaves, t_lev, t_leaf, t_rank>::operator=(tr);
            if (this != &tr) {
                m_k_for_level = tr.m_k_for_level;
                m_shift_table = tr.m_shift_table;
            }
            return *this;
        }

        k2_tree_hybrid &operator=(const k2_tree_hybrid &tr) {
            k2_tree_base<t_k_l_1, t_k_leaves, t_lev, t_leaf, t_rank>::operator=(tr);
            if (this != &tr) {
                m_k_for_level = tr.m_k_for_level;
            }
            return *this;
        }

        bool operator==(const k2_tree_hybrid &tr) const {
            if (!k2_tree_base<t_k_l_1, t_k_leaves, t_lev, t_leaf, t_rank>::operator==(tr)) {
                return false;
            }

            if (m_k_for_level.size() !=
                tr.m_k_for_level.size()) {//must be the same for the same template parameters and data
                return false;
            }

            if (m_shift_table.size() !=
                tr.m_shift_table.size()) {//must be the same for the same template parameters and data
                return false;
            }

            return true;
        }

        void swap(k2_tree_hybrid &tr) {
            k2_tree_base<t_k_l_1, t_k_leaves, t_lev, t_leaf, t_rank>::swap(tr);
            std::swap(m_k_for_level, tr.m_k_for_level);
            std::swap(m_shift_table, tr.m_shift_table);
        }

        virtual void construct_access_shortcut(uint8_t access_shortcut_size) override {
            if (this->m_access_shortcut_size <= t_k_l_1_size){
                k2_tree_base<t_k_l_1, t_k_leaves, t_lev, t_leaf, t_rank>::construct_access_shortcut(access_shortcut_size);
            } else {
                throw std::runtime_error("access shortcut is only supported over levels with the same k values, therefore t_k_l_1_size must be >= access_shortcut_size");
            }
        }

        std::string get_type_string_without_compression() const {
            return "k2_tree_hybrid<"+std::to_string(t_k_l_1)+","+std::to_string(t_k_l_1_size)+","+std::to_string(t_k_l_2)+","+std::to_string(t_k_leaves)+">";
        }

        std::string get_type_string() const {
            return "k2_tree_hybrid<"+std::to_string(t_k_l_1)+","+std::to_string(t_k_l_1_size)+","+std::to_string(t_k_l_2)+","+std::to_string(t_k_leaves)+","+get_compression_name(this->m_used_compression)+">";
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

            if (this->m_size == 0) {
                return;
            }

            std::string id_part = util::to_string(util::pid())
                                  + "_" + util::to_string(util::id());

            /*amount of levels with a k value of t_k_l_1 might differ from t_k_l_1_size as
             * it is enforced that the leaf level has k=t_k_leaves therefore if
             * m_tree_height <= t_k_l_1_size, the actual amount of levels with k=t_k_l_1
             * is smaller than t_k_l_1_size
             * for details look at get_tree_height, which in case of a hybrid tree calculates the tree height
             * considering above constraints
             * */
            uint8_t levels_with_k1 = 0;
            uint8_t levels_with_k2 = 0;
            uint8_t levels_with_k_leaves = 1;

            int ctr = 0;
            while (ctr < (this->m_tree_height-1)){
                auto k = get_k(ctr);
                if (k == t_k_l_1){
                    levels_with_k1++;
                } else if (k == t_k_l_2) {
                    levels_with_k2++;
                }
                ctr++;
            }

            const t_x bitsToInterleaveForK1 = bits::hi(t_k_l_1) * levels_with_k1;
            const t_x bitsToInterleaveForK2 = bits::hi(t_k_l_2) * levels_with_k2;
            const t_x bitsToInterleaveForKLeaves = bits::hi(t_k_leaves) * levels_with_k_leaves;

            t_x bitsOfMaximalValue = bitsToInterleaveForK1+bitsToInterleaveForK2+bitsToInterleaveForKLeaves;

            //set to one between 2*(bitsToInterleaveForK2+bitsToInterleaveForKLeaves) and 2*bitsOfMaximalValue
            t_x k1_bitmask = createBitmask(2*(bitsToInterleaveForK2+bitsToInterleaveForKLeaves), 2*bitsOfMaximalValue);
            t_x k2_bitmask = createBitmask(2*(bitsToInterleaveForKLeaves), 2*(bitsToInterleaveForK2+bitsToInterleaveForKLeaves));
            t_x k_leaves_bitmask = createBitmask(t_x(0), 2*(bitsToInterleaveForKLeaves));
            
            std::cout << "Sorting By Z Order" << std::endl;
            __gnu_parallel::sort(links.begin(), links.end(), [&](const t_e &lhs, const t_e &rhs) {
                auto lhs_interleaved = ((interleave<t_k_l_1>::bits(lhs) & k1_bitmask) | (interleave<t_k_l_2>::bits(lhs) & k2_bitmask) | (interleave<t_k_leaves>::bits(lhs) & k_leaves_bitmask));
                auto rhs_interleaved = ((interleave<t_k_l_1>::bits(rhs) & k1_bitmask) | (interleave<t_k_l_2>::bits(rhs) & k2_bitmask) | (interleave<t_k_leaves>::bits(rhs) & k_leaves_bitmask));
                return lhs_interleaved < rhs_interleaved;
                //return (t_k_leaves * divexp(lhs.first, 0) + divexp(lhs.second, 0) < (t_k_leaves * divexp(rhs.first, 0) + divexp(rhs.second, 0)));
            });

            std::cout << "Sorting Finished, Constructing Bitvectors" << std::endl;

            /*for (int m = 0; m < links.size(); ++m) {
                std::cout << links[m].first << "," << links[m].second << std::endl;
            }*/


            std::vector<int64_t> previous_subtree_number(this->m_tree_height, -1);

            {
                int64_t subtree_distance;
                bool fill_to_k2_entries = false; //begin extra case!
                std::vector<uint> gap_to_k2(this->m_tree_height);
                for (uint i = 0; i < gap_to_k2.size(); ++i) {
                    gap_to_k2[i] = get_k(i) * get_k(i);
                }
                bool firstLink = true;
                uint current_subtree_number = 0;

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
                                gap_to_k2[current_level] = get_k(current_level) * get_k(current_level);
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

    private:
        //FIXME: declared here and in k2_tree as a workaround, because virtual template methods are not possible
        template<typename t_vector>
        void contruct(t_vector &v, const construction_algorithm construction_algo,
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
            auto k = get_k(level);
            t_x result = k * divexp(link.first, exponent) + divexp(link.second, exponent);
            link.first = modexp(link.first, exponent);
            link.second = modexp(link.second, exponent);

            return result;
        }

        template <typename t_x>
        t_x createBitmask(t_x start, t_x end){
            t_x result = 0;
            for (auto i = start; i < end ; ++i) {
                result |= (1 << i);
            }
            return result;
        }

        template<typename t_vector>
        void construct_counting_sort(t_vector &links) {
            this->construct_counting_sort_internal(links, [&](uint64_t x, uint8_t l){return divexp(x,l);}, [&](uint8_t l){return exp(l);});
        }

        template<typename t_vector>
        void construct(t_vector &links, std::string temp_file_prefix = "") {
            this->construct_internal(links, [&](uint64_t x, uint8_t l){return divexp(x,l);}, temp_file_prefix);
        }

        inline uint64_t exp(uint8_t l) const {
            assert(l >= 0 && l <= m_tree_height);
            return 1Ull << m_shift_table[l];
        }

        template<typename t_x>
        inline t_x divexp(t_x x, uint8_t l) const {
            assert(l >= 0 && l <= m_tree_height);
            return x >> m_shift_table[l];
        }

        template<typename t_x>
        inline t_x modexp(t_x x, uint8_t l) const {
            assert(l >= 0 && l <= m_tree_height);
            return x & bits::lo_set[m_shift_table[l]];
        }

        template<typename t_x>
        inline t_x multexp(t_x x, uint8_t l) const {
            assert(l >= 0 && l <= m_tree_height);
            return x << m_shift_table[l];
        }


        void initialize_shift_table() {
            m_shift_table.resize(this->m_tree_height + 1);

            m_shift_table[0] = 0;
            for (uint i = 0; i < this->m_tree_height; ++i) {
                m_shift_table[i + 1] = m_shift_table[i] + bits::hi(m_k_for_level[this->m_tree_height-i-1]);
            }
        }

        inline uint8_t get_k(uint8_t level) const {
            return m_k_for_level[level];
        }

        uint8_t get_tree_height(const uint64_t max) {
            uint8_t res = 1;
            if (t_k_l_1 <= max) {
                this->m_max_element = t_k_l_1;
                m_k_for_level.push_back(t_k_l_1);
            } else {
                //in case only one level use k_leaves (has to be higher than t_k_l_1 (compile time-checked)
                this->m_max_element = t_k_leaves;
                m_k_for_level.push_back(t_k_leaves);
            }

            while (this->m_max_element <= max) {
                if ((this->m_max_element * t_k_leaves) > max) {
                    this->m_max_element = this->m_max_element * t_k_leaves;
                    m_k_for_level.push_back(t_k_leaves);
                } else if (res < t_k_l_1_size) {
                    this->m_max_element = this->m_max_element * t_k_l_1;
                    m_k_for_level.push_back(t_k_l_1);
                } else {
                    this->m_max_element = this->m_max_element * t_k_l_2;
                    m_k_for_level.push_back(t_k_l_2);
                }
                ++res;
            }

            m_k_for_level.resize(res);
            if (res <= t_k_l_1_size) {

                std::cerr
                        << "The tree height is smaller than t_k_l_1_size using t_k_l_1 up to m_tree_height-1 and then t_k_leaves"
                        << std::endl;
            } else if (res == t_k_l_1_size + 1) {
                std::cerr << "The tree equals t_k_l_1_size+1, only using t_k_l_1 and t_k_leaves" << std::endl;
            }

            return res;
        }
    };
}
#endif

