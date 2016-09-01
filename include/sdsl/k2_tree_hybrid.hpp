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

#include "vectors.hpp"
#include "bits.hpp"
#include "k2_tree_base.hpp"
#include "k2_tree_helper.hpp"
#include "k2_tree_vocabulary.hpp"
#include <tuple>
#include <algorithm>
#include <stxxl/vector>
#include <climits>
#include <vector>
#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/rank_support_v.hpp>
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
            bool t_comp = false,
            uint8_t t_access_shortcut_size = 0,
            typename t_rank=typename t_lev::rank_1_type>

    class k2_tree_hybrid : public k2_tree_base<t_k_l_1, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank> {
        static_assert(t_k_l_1 > 1, "t_k has to be larger than 1.");
        static_assert(t_k_l_1 <= 16, "t_k has to be smaller than 17.");
        static_assert(t_k_l_2 > 1, "t_k has to be larger than 1.");
        static_assert(t_k_l_2 <= 16, "t_k has to be smaller than 17.");
        static_assert(t_k_leaves > 1, "t_k has to be larger than 1.");
        static_assert(t_k_leaves >= t_k_l_1,
                      "t_k_leaves has to be larger than t_k_l_1,  otherwise this could lead to different word sizes and thus to a problem for the k2part approach"); //if smaller than t_k_l_1 it could be that t_k_leaves is not used
        static_assert(t_k_leaves <= 16, "t_k has to be smaller than 17.");
        static_assert(
                t_access_shortcut_size == 0 || (t_access_shortcut_size > 0 && t_access_shortcut_size <= t_k_l_1_size),
                "when using the access shortcut, the the levels up to the access shortcut need to have the same k value");


    public:

        std::vector<uint8_t> m_k_for_level;

        typedef stxxl::VECTOR_GENERATOR<std::pair<uint32_t, uint32_t>>::result stxxl_32bit_pair_vector;
        typedef stxxl::VECTOR_GENERATOR<std::pair<uint64_t, uint64_t>>::result stxxl_64bit_pair_vector;
        typedef int_vector<>::size_type size_type;

        k2_tree_hybrid() = default;

        k2_tree_hybrid(const k2_tree_hybrid &tr)
                : k2_tree_base<t_k_l_1, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank>(tr) {
            *this = tr;
        }

        k2_tree_hybrid(k2_tree_hybrid &&tr)
                : k2_tree_base<t_k_l_1, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank>(tr) {
            *this = std::move(tr);
        }

        template<typename t_vector>
        k2_tree_hybrid(std::string temp_file_prefix, bool use_counting_sort, t_vector &v, uint64_t max_hint = 0,
                       uint64_t hash_size = 0) {

            using namespace k2_treap_ns;
            if (v.size() > 0) {
                if (max_hint == 0) {
                    max_hint = get_maximum(v);
                }
                this->m_tree_height = get_tree_height(max_hint);

                if (use_counting_sort) {
                    k2_tree_base<t_k_l_1, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank>::template construct_counting_sort(
                            v, temp_file_prefix);
                    //construct_bottom_up(v, temp_file_prefix);
                } else {
                    construct(v, temp_file_prefix);
                }

                this->m_access_shortcut_size = t_access_shortcut_size;
                if (t_access_shortcut_size > 0) {
                    this->construct_access_shortcut();
                }

                if (t_comp) {
                    this->compress_leaves(hash_size);
                }
            }
        }

        k2_tree_hybrid(int_vector_buffer<> &buf_x,
                       int_vector_buffer<> &buf_y, bool use_counting_sort = false, uint64_t max_hint = 0,
                       uint64_t hash_size = 0) {
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
            if (this->m_max_element <= std::numeric_limits<uint32_t>::max()) {
                auto v = k2_tree_base<t_k_l_1, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank>::template read<uint32_t, uint32_t>(
                        bufs);
                if (use_counting_sort) {
                    k2_tree_base<t_k_l_1, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank>::template construct_counting_sort(
                            v, buf_x.filename());
                } else {
                    construct(v, buf_x.filename());
                }

            } else {
                auto v = k2_tree_base<t_k_l_1, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank>::template read<uint64_t, uint64_t>(
                        bufs);
                if (use_counting_sort) {
                    k2_tree_base<t_k_l_1, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank>::template construct_counting_sort(
                            v, buf_x.filename());
                } else {
                    construct(v, buf_x.filename());
                }
            }

            if (this->m_tree_height > 0) {
                this->m_access_shortcut_size = t_access_shortcut_size;
                if (t_access_shortcut_size > 0) {
                    this->construct_access_shortcut();
                }

                if (t_comp) {
                    this->compress_leaves(hash_size);
                }
            }
        }

        virtual size_type serialize(std::ostream &out, structure_tree_node *v, std::string name) const override {
            return k2_tree_base<t_k_l_1, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank>::serialize(out, v,
                                                                                                           name);
        }

        inline uint8_t get_k(uint8_t level) const {
            return m_k_for_level[level];
        }

        uint word_size() const {
            return div_ceil((uint) t_k_leaves * t_k_leaves, kUcharBits);
        }


        size_t words_count() const {
            return this->m_leaves.size() / t_k_leaves / t_k_leaves;
        }

        void load(std::istream &in) override {
            k2_tree_base<t_k_l_1, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank>::load(in);
            if (this->m_tree_height > 0) {
                for (int i = 1; i <= std::min(t_k_l_1_size, (uint8_t) (this->m_tree_height - 1)); ++i) {
                    m_k_for_level.push_back(t_k_l_1);
                }

                for (int j = t_k_l_1_size + 1; j <= (this->m_tree_height - 1); ++j) {
                    m_k_for_level.push_back(t_k_l_2);
                }

                m_k_for_level.push_back(t_k_leaves);

                if (t_access_shortcut_size > 0) {
                    this->perform_access_shortcut_precomputations();
                }
            }
        }

        //hack a the moment, because construct cannot be virtual
        void load_from_ladrabin(std::string fileName, uint64_t hash_size = 0, bool use_counting_sort = false,
                                std::string temp_file_prefix = "") {
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
                    if (use_counting_sort) {
                        k2_tree_base<t_k_l_1, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank>::template construct_counting_sort(
                                coords, temp_file_prefix);
                        //construct_bottom_up(v, temp_file_prefix);
                    } else {
                        construct(coords, temp_file_prefix);
                    }

                    std::cout << "Finished Construction" << std::endl;

                    coords.clear();

                    if (t_comp) {
                        this->compress_leaves(hash_size);
                    }

                    this->m_access_shortcut_size = t_access_shortcut_size;
                    if (t_access_shortcut_size > 0) {
                        this->construct_access_shortcut();
                    }
                }
            } else {
                throw std::runtime_error("Could not load ladrabin file");
            }
        }

        k2_tree_hybrid &operator=(k2_tree_hybrid &&tr) {
            k2_tree_base<t_k_l_1, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank>::operator=(tr);
            if (this != &tr) {
                m_k_for_level = tr.m_k_for_level;
            }
            return *this;
        }

        k2_tree_hybrid &operator=(const k2_tree_hybrid &tr) {
            k2_tree_base<t_k_l_1, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank>::operator=(tr);
            if (this != &tr) {
                m_k_for_level = tr.m_k_for_level;
            }
            return *this;
        }

        bool operator==(const k2_tree_hybrid &tr) const {
            if (!k2_tree_base<t_k_l_1, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank>::operator==(tr)) {
                return false;
            }

            if (m_k_for_level.size() !=
                tr.m_k_for_level.size()) {//must be the same for the same template parameters and data
                return false;
            }

            return true;
        }

        void swap(k2_tree_hybrid &tr) {
            k2_tree_base<t_k_l_1, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank>::swap(tr);
            std::swap(m_k_for_level, tr.m_k_for_level);
        }

    private:

        /**
         * Constructs the tree corresponding to the points in the links vector by partitioning the input multiple times
         * @param links
         * @param temp_file_prefix
         */
        template<typename t_vector>
        void construct(t_vector &links, std::string temp_file_prefix = "") {
            using namespace k2_treap_ns;
            typedef decltype(links[0].first) t_x;
            typedef decltype(links[0].second) t_y;
            using t_e = std::pair<t_x, t_y>;

            this->m_size = links.size();

            if (this->m_tree_height == 0) {//might occur in k2part
                return;
            }

            std::string id_part = util::to_string(util::pid())
                                  + "_" + util::to_string(util::id());

            {
                std::vector<int_vector_buffer<1>> level_buffers = this->create_level_buffers(temp_file_prefix, id_part);

                auto end = std::end(links);
                //uint64_t last_level_bits = 0;
                //uint64_t level_bits = 0;

                //recursively partition that stuff
                uint64_t submatrix_size = this->m_max_element;
                for (int l = this->m_tree_height; l + 1 > 0; --l) {

                    //std::cout << "Processing Level " << l << std::endl;
                    //level_bits = 0;

                    uint8_t k = 0;
                    uint64_t current_submatrix_size = 0;
                    if (l > 0) {
                        k = get_k(this->m_tree_height - l);
                        current_submatrix_size = submatrix_size / k;
                    }

                    auto sp = std::begin(links);
                    for (auto ep = sp; ep != end;) {

                        //Iterator which only returns the nodes within a certain subtree
                        ep = std::find_if(sp, end, [&submatrix_size, &sp, &l](const t_e &e) {
                            auto x1 = std::get<0>(*sp);
                            auto y1 = std::get<1>(*sp);
                            auto x2 = std::get<0>(e);
                            auto y2 = std::get<1>(e);

                            bool in_sub_tree = (x1 / submatrix_size) != (x2 / submatrix_size)
                                               or (y1 / submatrix_size) != (y2 / submatrix_size);

                            return in_sub_tree;
                        });

                        if (l > 0) {
                            auto _sp = sp;

                            for (uint8_t i = 0; i < k; ++i) {
                                auto _ep = ep;
                                if (i + 1 < k) {  //partition t_k -1 times vertically (1 in the case of k=2)
                                    _ep = std::partition(_sp, _ep,
                                                         [& current_submatrix_size, &k, &i, &l](const t_e &e) {
                                                             return (std::get<0>(e) / current_submatrix_size) % k <= i;
                                                         });
                                }
                                auto __sp = _sp;

                                for (uint8_t j = 0;
                                     j < k; ++j) { //partition the t_k vertical partitions t_k -1 times horizontally
                                    auto __ep = _ep;
                                    if (j + 1 < k) {
                                        __ep = std::partition(__sp, _ep,
                                                              [& current_submatrix_size, &k, &j, &l](const t_e &e) {
                                                                  return (std::get<1>(e) / current_submatrix_size) %
                                                                         k <= j;
                                                              });
                                    }
                                    bool not_empty = __ep > __sp;
                                    level_buffers[this->m_tree_height - l].push_back(not_empty);
                                    //level_bits++;
                                    __sp = __ep;
                                }
                                _sp = _ep;
                            }
                        }
                        if (ep != end) {
                            //++ep;
                            sp = ep;
                        }
                    }
                    //last_level_bits = level_bits;
                    submatrix_size = current_submatrix_size;
                }
            }

            this->load_vectors_from_file(temp_file_prefix, id_part);

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

