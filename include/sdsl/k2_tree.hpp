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

#include "vectors.hpp"
#include "bits.hpp"
#include "k2_tree_helper.hpp"
#include "k2_tree_base.hpp"
#include <stxxl/vector>
#include <tuple>
#include <algorithm>
#include <climits>
#include <vector>
#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/rank_support_v.hpp>
#include <gtest/gtest_prod.h>
#include <parallel/algorithm>

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
            bool t_comp=false,
            typename t_rank=typename t_lev::rank_1_type>

    class k2_tree : public k2_tree_base<t_lev, t_leaf, t_rank, t_comp> {
        static_assert(t_k > 1, "t_k has to be larger than 1.");
        static_assert(t_k <= 16, "t_k has to be smaller than 17.");

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

        using k2_tree_base<t_lev, t_leaf, t_rank, t_comp>::operator=;
        using k2_tree_base<t_lev, t_leaf, t_rank, t_comp>::operator==;

        enum {
            k = t_k
        };

        k2_tree() = default;

        /*k2_tree(const k2_tree &tr) {
            *this = tr;
        }

        k2_tree(k2_tree &&tr) {
            *this = std::move(tr);
        }*/

        template<typename t_vector>
        k2_tree(std::string temp_file_prefix, bool use_counting_sort, uint access_shortcut_size, t_vector &v, uint64_t max_hint=0){
            this->m_access_shortcut_size = access_shortcut_size;
            this->m_tree_height = get_tree_height(v, max_hint);

            if (v.size() > 0) {
                if (use_counting_sort){
                    k2_tree_base<t_lev, t_leaf, t_rank, t_comp>::template construct_counting_sort(v, temp_file_prefix);
                    //construct_bottom_up(v, temp_file_prefix);
                } else  {
                    construct(v, temp_file_prefix);
                }
            }

            if (this->m_access_shortcut_size > 0){
                construct_access_shortcut();
            }

            if (t_comp){
                this->compress_leaves();
            }
        }

        k2_tree(int_vector_buffer<> &buf_x,
                 int_vector_buffer<> &buf_y, bool use_counting_sort = false, uint access_shortcut_size = 0) {
            using namespace k2_treap_ns;
            typedef int_vector_buffer<> *t_buf_p;
            std::vector<t_buf_p> bufs = {&buf_x, &buf_y};

            this->m_access_shortcut_size = access_shortcut_size;

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

            uint64_t x = max_buf_element();
            uint8_t res = 0;
            while (res <= 64 and precomp<t_k>::exp(res) <= x) { ++res; }
            if (res == 65) {
                throw std::logic_error("Maximal element of input is too big.");
            }

            this->m_tree_height = res;
            this->m_max_element = precomp<t_k>::exp(res);

            if (precomp<t_k>::exp(res) <= std::numeric_limits<uint32_t>::max()) {
                auto v = k2_tree_base<t_lev, t_leaf, t_rank, t_comp>::template read< uint32_t, uint32_t>(bufs);
                if (use_counting_sort){
                    k2_tree_base<t_lev, t_leaf, t_rank, t_comp>::template construct_counting_sort<std::vector<std::pair<uint32_t, uint32_t>>>(v, buf_x.filename());
                } else  {
                    construct(v, buf_x.filename());
                }

            } else {
                auto v = k2_tree_base<t_lev, t_leaf, t_rank, t_comp>::template read< uint64_t, uint64_t>(bufs);
                if (use_counting_sort){
                    k2_tree_base<t_lev, t_leaf, t_rank, t_comp>::template construct_counting_sort<std::vector<std::pair<uint64_t, uint64_t>>>(v, buf_x.filename());
                } else  {
                    construct(v, buf_x.filename());
                }
            }

            if (this->m_access_shortcut_size > 0){
                construct_access_shortcut();
            }

            if (t_comp){
                this->compress_leaves();
            }
        }

        inline uint8_t get_k(uint8_t) const{
            return t_k;
        }

        /**
         * Checks whether link from p = link.first to q = link.second is present i.e. matrix entry a_pq = 1
         */
        template<typename t_x, typename t_y>
        bool check_link(std::pair<t_x,t_y> link) const{

            //Patological case happening e.g. when using k2part
            if (this->m_leaves.size() == 0){
                return false;
            }
            if (t_comp){
                return check_link_internal(0, link.first, link.second, 0, [this](int64_t pos, uint8_t leafK){
                    return this->is_leaf_bit_set_comp(pos, leafK);
                });
            } else {
                return check_link_internal(0, link.first, link.second, 0, [this](int64_t pos, uint8_t leafK){
                    return this->is_leaf_bit_set(pos, leafK);
                });
            }
        }


    private:
        /**
         * Checks wether the edge p-->q exists recursively
         * @param level
         *  current level, initialy 0
         * @param p
         *  source_node
         * @param q
         *  target_node
         * @param index
         *  contains the index of the first child of the previous node, initially set to 0
         * @return
         */
        template<typename t_x, typename t_y, typename Function>
        bool check_link_internal(int level, t_x p, t_y q, int64_t index, Function check_leaf) const {
            using namespace k2_treap_ns;

            uint64_t current_submatrix_size = precomp<t_k>::exp(this->m_tree_height-level-1);
            int64_t y = index + t_k * (p/current_submatrix_size) + (q/current_submatrix_size);

            if (this->is_leaf_level(level)){
                return check_leaf(y,k);
            } else if (this->m_levels[level][y]) {
                return check_link_internal(level+1, p % current_submatrix_size, q % current_submatrix_size, this->get_child_index(0, y, level), check_leaf);
            } else {
                return false;
            }
        }

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

            if (this->m_tree_height == 0){//might occur in k2part
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
                for (int l = this->m_tree_height; l + 1 > 0; --l) {

                    //std::cout << "Processing Level " << l << std::endl;
                    //level_bits = 0;

                    auto sp = std::begin(links);
                    for (auto ep = sp; ep != end;) {

                        //Iterator which only returns the nodes within a certain subtree
                        ep = std::find_if(sp, end, [=, &sp, &l](const t_e &e) {
                            auto x1 = std::get<0>(*sp);
                            auto y1 = std::get<1>(*sp);
                            auto x2 = std::get<0>(e);
                            auto y2 = std::get<1>(e);
                            bool in_sub_tree = precomp<t_k>::divexp(x1, l) != precomp<t_k>::divexp(x2, l)
                                               or precomp<t_k>::divexp(y1, l) != precomp<t_k>::divexp(y2, l);

                            return in_sub_tree;
                        });


                        if (l > 0) {
                            auto _sp = sp;

                            for (uint8_t i = 0; i < k; ++i) {
                                auto _ep = ep;
                                if (i + 1 < k) {  //partition t_k -1 times vertically (1 in the case of k=2)
                                    _ep = std::partition(_sp, _ep, [=, &i, &l](const t_e &e) {
                                        return precomp<t_k>::divexp(std::get<0>(e), l - 1) % k <= i;
                                    });
                                }
                                auto __sp = _sp;

                                for (uint8_t j = 0; j < k; ++j) { //partition the t_k vertical partitions t_k -1 times horizontally
                                        auto __ep = _ep;
                                        if (j + 1 < k) {
                                            __ep = std::partition(__sp, _ep, [=, &j, &l](const t_e &e) {
                                                return precomp<t_k>::divexp(std::get<1>(e), l - 1) % k <= j;
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
                }
            }

            this->load_vectors_from_file(temp_file_prefix, id_part);

        }

        /**
         * Constructs the tree corresponding to the points in the links vector inpace by performing a z order sort and subsequently constructing the tree top down
         * @param links
         * @param temp_file_prefix
         */
        template<typename t_x, typename t_y>
        void construct_by_z_order_sort_internal(std::vector<std::pair<t_x, t_y>> &links, std::string temp_file_prefix = "") {
            using namespace k2_treap_ns;
            using t_e = std::pair<t_x, t_y>;

            this->m_size = links.size();

            std::string id_part = util::to_string(util::pid())
                                  + "_" + util::to_string(util::id());

            std::cout << "Sorting By Z Order" << std::endl;
            __gnu_parallel::sort(links.begin(), links.end(), [&](const t_e& lhs, const t_e& rhs){
                return sort_by_z_order(lhs, rhs);
            });

            std::cout << "Sorting Finished, Constructing Bitvectors" << std::endl;
            std::vector<int> previous_subtree_number(this->m_tree_height,-1);
            uint64_t total_size = 0;

            {
                int subtree_distance;
                bool fill_to_k2_entries = false; //begin extra case!
                std::vector<uint> gap_to_k2(this->m_tree_height,k*k);
                bool firstLink = true;
                uint current_subtree_number = 0;

                std::vector<int_vector_buffer<1>> level_buffers = this->create_level_buffers(temp_file_prefix, id_part);

                std::pair<t_x,t_y> previous_link;
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

                            for (int j = 0; j < subtree_distance - 1; ++j) {
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
                            std::string error_message("negative subtree_distance after z_order sort is not possible, somethings wrong current_level="+
                                                      std::to_string(current_level)+" subtree_distance="+std::to_string(subtree_distance)+
                                                      " current_subtree_number="+std::to_string(current_subtree_number)+" previous_subtree_number[current_level]="+
                                                      std::to_string(previous_subtree_number[current_level])+"current_link="+std::to_string(current_link.first)+","+std::to_string(current_link.second)+
                                                      "previous_link="+std::to_string(previous_link.first)+","+std::to_string(previous_link.second));
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
         * Comparision function for z_order_sort.
         * @param lhs
         * @param rhs
         * @return
         * True if lhs is smaller in z order, false otherwise
         */
        template<typename t_x, typename t_y>
        bool sort_by_z_order(const std::pair<t_x, t_y> lhs, const std::pair<t_x, t_y> rhs) {
            using namespace k2_treap_ns;
            if (lhs.first <= rhs.first && lhs.second <= rhs.second) {
                return true;
            } else if (lhs.first >= rhs.first && lhs.second >= rhs.second) {
                return false;
            } else if (lhs.first < rhs.first && lhs.second > rhs.second) {

                t_x lhsFirst = lhs.first;
                t_x lhsSecond = lhs.second;
                t_x rhsFirst = rhs.first;
                t_x rhsSecond = rhs.second;


                t_x lhsFirstDiv;
                t_x lhsSecondDiv;

                t_x rhsFirstDiv;
                t_x rhsSecondDiv;

                for (int i = this->m_tree_height; i > 0; --i) {
                    lhsFirstDiv = precomp<t_k>::divexp(lhsFirst, i);
                    lhsSecondDiv = precomp<t_k>::divexp(lhsSecond, i);
                    rhsFirstDiv = precomp<t_k>::divexp(rhsFirst, i);
                    rhsSecondDiv = precomp<t_k>::divexp(rhsSecond, i);

                    if (lhsFirstDiv < rhsFirstDiv) {
                        return true;
                    } else if (lhsFirstDiv == rhsFirstDiv && lhsSecondDiv > rhsSecondDiv) {
                        return false;
                    }
                }

                return true;

            } else { //lhs.first > rhs.first && lhs.second < rhs.second

                t_x lhsFirst = lhs.first;
                t_x lhsSecond = lhs.second;
                t_x rhsFirst = rhs.first;
                t_x rhsSecond = rhs.second;

                t_x lhsFirstDiv;
                t_x lhsSecondDiv;

                t_x rhsFirstDiv;
                t_x rhsSecondDiv;

                for (int i = this->m_tree_height; i > 0; --i) {
                    lhsFirstDiv = precomp<t_k>::divexp(lhsFirst, i);
                    lhsSecondDiv = precomp<t_k>::divexp(lhsSecond, i);
                    rhsFirstDiv = precomp<t_k>::divexp(rhsFirst, i);
                    rhsSecondDiv = precomp<t_k>::divexp(rhsSecond, i);

                    if (lhsFirstDiv > rhsFirstDiv) {
                        return false;
                    } else if (lhsFirstDiv == rhsFirstDiv && lhsSecondDiv < rhsSecondDiv) {
                        return true;
                    }
                }
                return false;
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
        uint inline calculate_subtree_number_and_new_relative_coordinates(std::pair<t_x, t_y>& link, int level) {
            using namespace k2_treap_ns;
            t_x exponent = this->m_tree_height-level-1;
            t_x result = k*precomp<t_k>::divexp(link.first,exponent)+precomp<t_k>::divexp(link.second,exponent);
            link.first = precomp<t_k>::modexp(link.first, exponent);
            link.second = precomp<t_k>::modexp(link.second, exponent);

            return result;
        }

        /**
        * Used for accelerating the check whether a certain link exists by skipping m_access_shortcut_size levels
        *
        * @param p Identifier of first object.
        * @param q Identifier of second object.
        *
        * @return Returns the subtree on level m_access_shortcut_size if present and nullptr otherwise
        */
        template<typename t_x, typename t_y>
        node_type *check_link_shortcut(t_x p, t_y q) const {
            using namespace k2_treap_ns;

            //FIXME: height if k_L tree!, it depends as we're not only targeting the last level anymore
            //FIXME: check which points are in the same tree and only fetch once
            //how to get corresponding subtree on level x of a point efficiently? (for k=2^x, interleave x-bitwise the top h bits
            //implement subtree calculation in general and for 2^x special-cases manually, think about precomp in the case of k=3

            //z = interleaved first h-1 set bits of p,q
            uint8_t real_size = 0;
            while (precomp<2>::exp(real_size) <= this->m_max_element) { ++real_size; }

            uint z = access_shortcut_helper<k>::corresponding_subtree(p, q, real_size, this->m_access_shortcut_size);
            //y = zth 1 via rank on B_
            uint y = this->m_access_shortcut_select_1_support(z+1);
            //check if exists and if B_[y-1] == 0 otherwise no link
            if (y < 0 || this->m_levels[this->m_access_shortcut_size][y-1] == true){
                return nullptr;
            }
            //rank 01 pattern on B[0,p] to find out how many non-empty trees are there until p
            //directly get corresponding data from leaf array
            uint number_of_present_trees_searched_value_is_in = this->m_access_shortcut_rank_01_support(y);

            //hack to get corresponding coordinates, might not be necessary later on
            uint field_size = precomp<t_k>::exp(this->m_tree_height - this->m_access_shortcut_size -1);
            uint upper_left_corner_x = 0;
            while (upper_left_corner_x <= p){
                upper_left_corner_x += field_size;
            }
            upper_left_corner_x -= field_size;

            uint upper_left_corner_y = 0;
            while (upper_left_corner_y <= q){
                upper_left_corner_y += field_size;
            }
            upper_left_corner_y -= field_size;

            uint index = number_of_present_trees_searched_value_is_in*k*k;
            node_type* result =  new node_type(this->m_access_shortcut_size, t_p(upper_left_corner_x, upper_left_corner_y), index);
            return result;
        }

        template<typename t_tv>
        uint8_t get_tree_height(const t_tv &v, uint64_t max_hint) {
            using namespace k2_treap_ns;
            if (v.size() == 0) {
                return 0;
            }

            uint64_t max;
            if (max_hint != 0){
                max = max_hint;
            } else  {
                using t_e = typename t_tv::value_type;
                auto tupmax = [](t_e a) {
                    return std::max(a.first, a.second);
                };
                auto max_it = std::max_element(std::begin(v), std::end(v), [&](t_e a, t_e b) {
                    return tupmax(a) < tupmax(b);
                });
                max = tupmax(*max_it);
            }

            uint8_t res = 0;
            while (precomp<t_k>::exp(res) <= max) { ++res; }

            this->m_max_element = precomp<t_k>::exp(res);
            return res;
        }

        /**
         * Hier noch eine Idee um den k^2-tree zu beschleunigen: Um nicht erst durch h Levels zu navigieren kann man sich erst ein bit_vector B bauen, der aus 4^h Einsen und höchstens 4^h Nullen besteht. Für jeden der 4^h Teilbäume schreibt man eine Eins; für nichtleere Teilbäume zusätzlich eine Null vor der entsprechenden Eins.
         *  Also für das Beispiel in Abb.1.3 (in Jans Bericht) mit h=2:
         *
         *       0 12 345678 901 2345
         *  B = 010110111111011101111
         *  P = 3 4  5      8   1
         *                      2
         *  Zum Teilbaum der Koordinate (x,y) kommt man indem man
         *   (1) die oberen h Bits von x und y interleaved; nennen wir das z
         *   (2) die Position p der (z+1)te Eins selektieren
         *   (3) Falls p=0 oder B[p-1]=1 ist, so ist der Teilbaum leer. Andernfalls
         *        ist der Teilbaum nicht leer und durch das Ergebnis r einer
         *       rank Operation auf das Bitpattern ,01' in B[0,p]
         *       adressieren wir einen Array P, der die Präfixsummer der
         *       Subbaumgrößen enthält. Der Eintrag P[r] kann als Pointer auf die
         *       k^2 Repräsentation des nichtleeren Subtrees dienen.
         *
         *  Das ist praktikabel für h=8. Der Bitvektor würde höchstens 16kBytes
         *  benötigen und P im worst case 2MB. Da nur ein rank und ein select
         *  gemacht werden sollte das deutlich schneller sein als die 8 ranks
         *  in der vorherigen Implementierung.

         */
        void construct_access_shortcut() {
            using namespace k2_treap_ns;

            //maximal size of shortcut is tree height
            if (this->m_access_shortcut_size > this->m_tree_height -2) {
                std::cout << "Reducing shortcut size to tree height -2" << std::endl;
                this->m_access_shortcut_size = this->m_tree_height - 2;
            }

            //Use 1 to code empty trees in level height-1 and 01 to code non-empty trees, height has to be calculated with kL_ as height of hybrid tree is different
            //coresponds to the amount of non-empty trees in level h-1
            //amount of Zeros = actually existent number of trees --> (level_begin_idx[level+2] - level_begin_idx[level+1])/k² or rank l(evel_begin_idx[level], level_begin_idx[level+1])

            uint64_t amountOfZeros = this->m_levels[this->m_access_shortcut_size].size() / (k*k); //spares rank of comp. level
            //corresponds to the theoretical amount of trees in level m_access_shortcut_size (round up (in case not divisible by k^2)
            uint64_t amountOfOnes = precomp<t_k*k>::exp((uint8_t) (this->m_access_shortcut_size-1));
            bit_vector access_shortcut(amountOfOnes+amountOfZeros, 1);

            //BitArray<uint> B(amountOfOnes + amountOfZeros);
            uint counter = 0;
            construct_access_shortcut_by_dfs(access_shortcut, this->root(), 0, counter);
            this->m_access_shortcut.swap(access_shortcut);

            sdsl::util::init_support(this->m_access_shortcut_rank_01_support, &this->m_access_shortcut);
            sdsl::util::init_support(this->m_access_shortcut_select_1_support, &this->m_access_shortcut);
        }

        /**
         * Constructs the bitvector m_access_shortcut used for speeding up tree traversal/link checks
         */
        void construct_access_shortcut_by_dfs(bit_vector& access_shortcut, node_type root, uint current_level, uint& counter) {
            using namespace k2_treap_ns;
            uint64_t rank = this->m_levels_rank[current_level](root.idx);
            auto x = std::real(root.p);
            auto y = std::imag(root.p);

            for (size_t i = 0; i < k; ++i) {
                for (size_t j = 0; j < k; ++j) {
                    // get_int better for compressed bitvectors
                    // or introduce cache for bitvectors
                    if (current_level == this->m_access_shortcut_size){
                        if (this->m_levels[current_level][root.idx + k * i + j]) { //if subtree present
                            access_shortcut[counter] = 0;//save 01 at counter position (m_access_shortcut gets initialised with 1s)
                            counter++;
                        }
                        counter++;


                    } else { //continue dfs tree traversal
                        if (this->m_levels[current_level][root.idx + k * i + j]) { //if subtree present
                            ++rank;
                            auto _x = x + i * precomp<t_k>::exp(root.t - 1);
                            auto _y = y + j * precomp<t_k>::exp(root.t - 1);

                            node_type subtree_root(root.t - 1, t_p(_x, _y), rank * k * k);
                            construct_access_shortcut_by_dfs(access_shortcut, subtree_root, current_level+1, counter);
                        } else {
                            counter += precomp<t_k>::exp(2 * (this->m_access_shortcut_size - current_level));
                        }
                    }
                }
            }
        }
    };
}
#endif

