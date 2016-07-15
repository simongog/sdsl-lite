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
/*! \file k2_treap.hpp
    \brief k2_treap.hpp contains a compact k^2-treap.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_K2_TREAP
#define INCLUDED_SDSL_K2_TREAP

#include "vectors.hpp"
#include "bits.hpp"
#include "k2_treap_helper.hpp"
#include "k2_treap_algorithm.hpp"
#include <tuple>
#include <algorithm>
#include <climits>
#include <vector>
#include <iostream>

//! Namespace for the succinct data structure library.
namespace sdsl {

//! A k^2-treap.
/*! A k^2-treap is an indexing structure for a set of weighted points. The set
 *  consists of triples (x,y,w), where the first two components x and y are
 *  the coordinates of the point and w is the point's weight.
 *
 *  The k^2 treap supports 4-sided range count queries and 4-sided prioritized
 *  range queries in 2d. Using the latter functionality it is also possible to
 *  support 6-sided range queries in 3d. An example can be found in
 *  examples/k2_treap_in_mem.cpp .
 *
 *  The k^2-treap constructed in-place. The construct method expects either
 *  a vector of std::array<X,3> elements (each array represent a tuple x,y,w)
 *  or a file prefix FILE. In the latter case three serialized int_vector<>
 *  have to be present at FILE.x, FILE.y, and FILE.w. One int_vector<> per
 *  component.
 *
 *  \par References
 *       [1] N. Brisaboa, G. de Bernardo, R. Konow, and G. Navarro:
 *           ,,$K^2$-Treaps: Range Top-$k$ Queries in Compact Space,
 *           Proceedings of SPIRE 2014.
 */
    template<uint8_t t_k,
            typename t_bv=bit_vector,
            typename t_rank=typename t_bv::rank_1_type>
    class k2_treap {
        static_assert(t_k > 1, "t_k has to be larger than 1.");
        static_assert(t_k <= 16, "t_k has to be smaller than 17.");

    public:
        typedef int_vector<>::size_type size_type;
        using node_type = k2_treap_ns::node_type;
        using point_type = k2_treap_ns::point_type;
        using t_p = k2_treap_ns::t_p;

        enum {
            k = t_k
        };

    private:
        uint8_t m_tree_height = 0;
        t_bv m_bp;
        t_rank m_bp_rank;
        int_vector<64> m_level_begin_idx;
        size_type m_size = 0;

        template<typename t_tv>
        uint8_t get_tree_height(const t_tv &v) {
            using namespace k2_treap_ns;
            if (v.size() == 0) {
                return 0;
            }
            using t_e = typename t_tv::value_type;
            auto tupmax = [](t_e a) {
                return std::max(a.first, a.second);
            };
            auto max_it = std::max_element(std::begin(v), std::end(v), [&](t_e a, t_e b) {
                return tupmax(a) < tupmax(b);
            });
            uint64_t x = tupmax(*max_it);
            uint8_t res = 0;
            while (precomp<t_k>::exp(res) <= x) { ++res; }
            return res;
        }


        //maybe change to table lookup version and extract B,S
        uint inline interleaveLowerBits2(uint x, uint y) {
            static const unsigned int B[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF};
            static const unsigned int S[] = {1, 2, 4, 8};

            // Interleave lower 16 bits of x and y, so the bits of x
            // are in the even positions and bits from y in the odd;
            unsigned int z; // z gets the resulting 32-bit Morton Number.
            // x and y must initially be less than 65536.


            x = (x | (x << S[3])) & B[3];
            x = (x | (x << S[2])) & B[2];
            x = (x | (x << S[1])) & B[1];
            x = (x | (x << S[0])) & B[0];

            y = (y | (y << S[3])) & B[3];
            y = (y | (y << S[2])) & B[2];
            y = (y | (y << S[1])) & B[1];
            y = (y | (y << S[0])) & B[0];

            z = x | (y << 1);

            return z;
        }

    public:
        uint8_t &t = m_tree_height;

        k2_treap() = default;

        k2_treap(const k2_treap &tr) {
            *this = tr;
        }

        k2_treap(k2_treap &&tr) {
            *this = std::move(tr);
        }

        //! Move assignment operator
        k2_treap &operator=(k2_treap &&tr) {
            if (this != &tr) {
                m_tree_height = tr.m_tree_height;
                m_size = tr.m_size;
                m_bp = std::move(tr.m_bp);
                m_bp_rank = std::move(tr.m_bp_rank);
                m_bp_rank.set_vector(&m_bp);
                m_level_begin_idx = std::move(tr.m_level_begin_idx);
            }
            return *this;
        }

        //! Assignment operator
        k2_treap &operator=(k2_treap &tr) {
            if (this != &tr) {
                m_tree_height = tr.m_tree_height;
                m_size = tr.m_size;
                m_bp = tr.m_bp;
                m_bp_rank = tr.m_bp_rank;
                m_bp_rank.set_vector(&m_bp);
                m_level_begin_idx = tr.m_level_begin_idx;
            }
            return *this;
        }

        //! Assignment operator
        bool operator==(const k2_treap &tr) const {
            if (m_tree_height != tr.m_tree_height)
                return false;
            if (m_size != tr.m_size)
                return false;
            if (m_bp.size() != tr.m_bp.size())
                return false;

            for (uint i = 0; i < m_bp.size(); ++i) {
                if (m_bp[i] != tr.m_bp[i]){
                    std::cout << "m_bp vectors differ at " << i << std::endl;
                    return false;
                }

            }

            if (m_level_begin_idx != tr.m_level_begin_idx)
                return false;

            return true;
        }

        //! Number of points in the 2^k treap
        size_type
        size() const {
            return m_size;
        }

        //! Swap operator
        void swap(k2_treap &tr) {
            if (this != &tr) {
                std::swap(m_tree_height, tr.m_tree_height);
                std::swap(m_size, tr.m_size);
                m_bp.swap(tr.m_bp);
                util::swap_support(m_bp_rank, tr.m_bp_rank, &m_bp, &(tr.m_bp));
                m_level_begin_idx.swap(tr.m_level_begin_idx);
            }
        }

        k2_treap(int_vector_buffer<> &buf_x,
                 int_vector_buffer<> &buf_y, bool bottom_up) {
            using namespace k2_treap_ns;
            typedef int_vector_buffer<> *t_buf_p;
            std::vector<t_buf_p> bufs = {&buf_x, &buf_y};

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

            if (precomp<t_k>::exp(res) <= std::numeric_limits<uint32_t>::max()) {
                auto v = read < uint32_t, uint32_t>(bufs);
                if (bottom_up){
                    construct_counting_sort(v, buf_x.filename());
                } else  {
                    construct(v, buf_x.filename());
                }

            } else {
                auto v = read < uint64_t, uint64_t>(bufs);
                if (bottom_up){
                    construct_counting_sort(v, buf_x.filename());
                } else  {
                    construct(v, buf_x.filename());
                }
            }
        }

        template<typename t_x=uint64_t, typename t_y=uint64_t>
        std::vector<std::pair<t_x, t_y>>
        read(std::vector<int_vector_buffer<> * > &bufs) {
            typedef std::vector<std::pair<t_x, t_y>> t_tuple_vec;
            t_tuple_vec v = t_tuple_vec(bufs[0]->size());
            for (uint64_t j = 0; j < v.size(); ++j) {
                std::get<0>(v[j]) = (*(bufs[0]))[j];
            }
            for (uint64_t j = 0; j < v.size(); ++j) {
                std::get<1>(v[j]) = (*(bufs[1]))[j];
            }

            return v;
        }


        template<typename t_x, typename t_y>
        k2_treap(std::vector<std::pair<t_x, t_y>> &v, std::string temp_file_prefix = "", bool bottom_up = false) {
            if (v.size() > 0) {
                if (bottom_up){
                    construct_counting_sort(v, temp_file_prefix);
                } else  {
                    construct(v, temp_file_prefix);
                }
            }
        }

        template<typename t_x, typename t_y>
        void construct(std::vector<std::pair<t_x, t_y>> &v, std::string temp_file_prefix = "") {
            using namespace k2_treap_ns;
            using t_e = std::pair<t_x, t_y>;

            m_size = v.size();
            m_tree_height = get_tree_height(v);
            uint64_t M = precomp<t_k>::exp(t);
            t_e MM = t_e(M, M);

            std::string id_part = util::to_string(util::pid())
                                  + "_" + util::to_string(util::id());

            m_level_begin_idx = int_vector<64>(1 + t, 0);

            std::string bp_file = temp_file_prefix + "_bp_" + id_part
                                  + ".sdsl";

            {
                int_vector_buffer<1> bp_buf(bp_file, std::ios::out);

                auto end = std::end(v);
                uint64_t last_level_bits = 0;
                uint64_t level_bits = 0;

                //recursively partition that stuff
                for (uint64_t l = t; l + 1 > 0; --l) {

                    if (l > 0) {
                        m_level_begin_idx[l - 1] = m_level_begin_idx[l] + last_level_bits;
                        std::cout << "Setring m_level_begin_idx["<<l-1<<"] =" << m_level_begin_idx[l] + last_level_bits << std::endl;
                    }

                    //std::cout << "Processing Level " << l << std::endl;

                    level_bits = 0;

                    auto sp = std::begin(v);
                    for (auto ep = sp; ep != end;) {

                        //Iterator which only returns the nodes within a certain subtree
                        ep = std::find_if(sp, end, [&sp, &l](const t_e &e) {
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

                            for (uint8_t i = 0; i < t_k; ++i) {
                                auto _ep = ep;
                                if (i + 1 < t_k) {  //partition t_k -1 times vertically (1 in the case of k=2)
                                    _ep = std::partition(_sp, _ep, [&i, &l](const t_e &e) {
                                        return precomp<t_k>::divexp(std::get<0>(e), l - 1) % t_k <= i;
                                    });
                                }
                                auto __sp = _sp;
                                for (uint8_t j = 0;
                                     j < t_k; ++j) { //partition the t_k vertical partitions t_k -1 times horizontally
                                    auto __ep = _ep;
                                    if (j + 1 < t_k) {
                                        __ep = std::partition(__sp, _ep, [&j, &l](const t_e &e) {
                                            return precomp<t_k>::divexp(std::get<1>(e), l - 1) % t_k <= j;
                                        });
                                    }
                                    bool not_empty = __ep > __sp;
                                    bp_buf.push_back(not_empty);
                                    level_bits++;
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
                    last_level_bits = level_bits;
                }
            }

            bit_vector bp;
            load_from_file(bp, bp_file);
            {
                bit_vector _bp;
                _bp.swap(bp);
                m_bp = t_bv(_bp);
                /*std::cout << "m_bp: \t";
                for (auto i = 0; i < m_bp.size(); ++i) {
                    std::cout << m_bp[i];
                }
                std::cout << std::endl;*/
            }

            util::init_support(m_bp_rank, &m_bp);
            sdsl::remove(bp_file);
        }

        template<typename t_x, typename t_y>
        void construct_counting_sort(std::vector<std::pair<t_x, t_y>> &links, std::string temp_file_prefix = "") {
            using namespace k2_treap_ns;
            using t_e = std::pair<t_x, t_y>;

            m_size = links.size();
            m_tree_height = get_tree_height(links);
            uint64_t M = precomp<t_k>::exp(t);
            t_e MM = t_e(M, M);

            std::string id_part = util::to_string(util::pid())
                                  + "_" + util::to_string(util::id());

            m_level_begin_idx = int_vector<64>(1 + t, 0);

            std::string bp_file = temp_file_prefix + "_bp_" + id_part
                                  + ".sdsl";


            {
                int_vector_buffer<1> bp_buf(bp_file, std::ios::out);

                //                  upper left          lower right                 interval in links   level
                typedef std::tuple<std::pair<t_x,t_y>,std::pair<uint64_t,uint64_t>, std::pair<t_x,t_y>, uint> t_queue;
                std::queue <t_queue> queue;

                //partition recursively until reaching the leaves
                uint64_t matrix_size = precomp<t_k>::exp(m_tree_height); //could be bigger than 32 bit although biggest value in links is 32 bit
                queue.push(std::make_tuple(std::make_pair<t_x,t_y>(0,0), std::make_pair(matrix_size-1, matrix_size-1), std::make_pair<t_x,t_x>(0, links.size()),0));

                uint64_t number_of_bits = 0; //for speed comparison purposes of different k

                m_level_begin_idx[m_tree_height+1] = 0;
                //std::cout << "Setring m_level_begin_idx["<<m_tree_height-1<<"] =" << 0 << std::endl;
                while (!queue.empty()) {
                    auto upper_left = std::get<0>(queue.front());
                    auto lower_right = std::get<1>(queue.front());
                    auto links_interval = std::get<2>(queue.front());
                    auto current_level = std::get<3>(queue.front());

                    queue.pop();

                    auto submatrix_size =
                            lower_right.first - upper_left.first + 1;//precomp<k>::exp(m_tree_height-level);
                    std::vector<t_x> intervals(k * k + 1);

                    //do counting sort
                    auto x1 = upper_left.first;
                    auto y1 = upper_left.second;
                    auto subDivK = (submatrix_size / k);
                    for (uint64_t j = links_interval.first; j < links_interval.second; ++j) {
                        auto x = links[j].first;
                        auto y = links[j].second;
                        auto p1 = (x - x1) / subDivK;
                        auto p2 = (y - y1) / subDivK;
                        uint corresponding_matrix = p1 * k + p2;
                            intervals[corresponding_matrix + 1]++;//offset corresponding matrix by one to allow for more efficient in interval comparision
                    }

                    intervals[0] = 0;
                    //append bits to level_vectors[level] based on result



                    for (uint i = 1; i < intervals.size(); ++i) {
                        if (intervals[i] > 0) {
                            bp_buf.push_back(1);
                            //std::cout << "1";
                            number_of_bits++;
                        } else {
                            bp_buf.push_back(0);
                            //std::cout << "0";
                            number_of_bits++;
                        }
                        intervals[i] += intervals[i - 1]; //build prefix sum
                    }


                    if (submatrix_size > k) {

                        //std::cout << std::endl;
                        m_level_begin_idx[m_tree_height-current_level-2] = number_of_bits;
                        //std::cout << "Setring m_level_begin_idx["<<m_tree_height-current_level-2<<"] =" << number_of_bits << std::endl;

                        std::vector<t_x> offset(k * k);
                        offset[0] = 0;
                        for (size_t l = 1; l < offset.size(); ++l) {
                            offset[l] = intervals[l] + 1;
                        }

                        auto begin = links.begin() + links_interval.first;
                        auto it = begin;

                        uint64_t index = 0;
                        while (it != links.begin() + links_interval.second){
                            auto x = (*it).first;
                            auto y = (*it).second;
                            auto p1 = (x - x1) / subDivK;
                            auto p2 = (y - y1) / subDivK;
                            uint corresponding_matrix = p1 * k + p2;

                            if (index >= intervals[corresponding_matrix] &&
                                index < intervals[corresponding_matrix + 1]) {
                                //element is at correct position
                                offset[corresponding_matrix]++;

                                it++;
                                index++;
                            } else {
                                //swap to correct position either swapping a match or a non match back
                                //there are at most m swaps if all m elements are at the wrong position
                                //every value is at most read twice
                                std::iter_swap(it, begin+offset[corresponding_matrix]-1);
                                offset[corresponding_matrix]++;
                            }
                        }


                        /*
                        //swap entries to the correct place based on prefix sum
                        uint64_t index = 0;
                        for (auto it = links.begin() + links_interval.first;
                             it != links.begin() + links_interval.second; ++it) {
                            auto x = (*it).first;
                            auto y = (*it).second;
                            uint corresponding_matrix =
                                    floor((x - x1) / (submatrix_size / k)) * k + floor((y - y1) / (submatrix_size / k));

                            if (index >= intervals[corresponding_matrix] &&
                                index < intervals[corresponding_matrix + 1]) {
                                //element is at correct position
                                offset[corresponding_matrix]++;
                            } else {
                                //search for a matching swap

                                uint submatrix_corresponding_to_index = 0;
                                while (intervals[submatrix_corresponding_to_index] <= index) {
                                    ++submatrix_corresponding_to_index;
                                }
                                submatrix_corresponding_to_index -= 1;

                                //think about more efficient implementation
                                if (nearest_point_belonging_to_matrix[corresponding_matrix] > index) {
                                    auto it2 = links.begin() + links_interval.first +
                                               nearest_point_belonging_to_matrix[corresponding_matrix];
                                    std::iter_swap(it, it2);
                                    offset[corresponding_matrix]++;
                                    offset[submatrix_corresponding_to_index]++;
                                } else {
                                    uint64_t index2 = 0;
                                    for (auto it2 = links.begin() + links_interval.first + index + 1;
                                         it2 != links.begin() + links_interval.second; ++it2) {
                                        auto x = (*it2).first;
                                        auto y = (*it2).second;
                                        uint corresponding_matrix2 = floor((x - x1) / (submatrix_size / k)) * k +
                                                                     floor((y - y1) / (submatrix_size / k));

                                        if (corresponding_matrix2 == submatrix_corresponding_to_index) {
                                            std::iter_swap(it, it2);
                                            offset[submatrix_corresponding_to_index]++;
                                            break;
                                        } else {
                                            //buffer not set                                            or      invalid
                                            if (nearest_point_belonging_to_matrix[corresponding_matrix] == 0 ||
                                                nearest_point_belonging_to_matrix[corresponding_matrix] < index) {
                                                nearest_point_belonging_to_matrix[corresponding_matrix] =
                                                        index + 1 + index2;
                                            }
                                        }
                                        index2++;
                                    }
                                }
                            }
                            index++;
                        }*/
                        auto new_submatrix_size = submatrix_size / k;
                        for (uint x = 0; x < k; ++x) {
                            for (uint y = 0; y < k; ++y) {
                                auto new_interval = std::make_pair(intervals[x * k + y]+links_interval.first, intervals[x * k + y + 1]+links_interval.first);
                                if (new_interval.first != new_interval.second) {
                                    auto new_upper_left = std::make_pair<t_x, t_y>(x * new_submatrix_size  + upper_left.first,
                                                                                   y * new_submatrix_size  + upper_left.second);
                                    auto new_lower_right = std::make_pair<t_x, t_y>((x + 1) * new_submatrix_size - 1 + upper_left.first,
                                                                                    (y + 1) * new_submatrix_size - 1 + upper_left.second);
                                    queue.push(std::make_tuple(new_upper_left, new_lower_right, new_interval,current_level+1));
                                }
                            }
                        }
                    }
                }
            }

            bit_vector bp;
            load_from_file(bp, bp_file);
            {
                bit_vector _bp;
                _bp.swap(bp);
                m_bp = t_bv(_bp);
                /*std::cout << "m_bp: \t";
                for (auto i = 0; i < m_bp.size(); ++i) {
                    std::cout << m_bp[i];
                }
                std::cout << std::endl;*/
            }

            util::init_support(m_bp_rank, &m_bp);
            sdsl::remove(bp_file);
        }

        template<typename t_x, typename t_y>
        void recursively_partition_with_counting_sort(std::vector<std::pair<t_x, t_y>> &links, std::pair<t_x, t_y> upper_left,
                                                      std::pair<t_x, t_y> lower_right, std::pair<t_x, t_x> corresponding_interval,
                                                      std::vector<int_vector_buffer<>>& level_vectors, uint level) {
            using namespace k2_treap_ns;
            if (level == m_tree_height || corresponding_interval.first == corresponding_interval.second){ //think about m_tree_height
                return;
            }
        }


        template<typename t_x, typename t_y>
        uint inline calculate_subtree_number_and_new_relative_coordinates(std::pair<t_x, t_y>& link, int level) {
            using namespace k2_treap_ns;
            t_x exponent = m_tree_height-level-1;
            t_x result = k*precomp<t_k>::divexp(link.first,exponent)+precomp<t_k>::divexp(link.second,exponent);
            link.first = precomp<t_k>::modexp(link.first, exponent);
            link.second = precomp<t_k>::modexp(link.second, exponent);

            return result;
        }


        /*template <size_t k>
        struct sort_by_z_order {
            //static const int MortonTable256[256];


            uint inline interleaveLowerBits(ushort x, ushort y) {
                unsigned int z;   // z gets the resulting 32-bit Morton Number.

                z = MortonTable256[y >> 8] << 17 |
                    MortonTable256[x >> 8] << 16 |
                    MortonTable256[y & 0xFF] << 1 |
                    MortonTable256[x & 0xFF];

                return z;
            }*/


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

                for (int i = m_tree_height; i > 0; --i) {
                    lhsFirstDiv = precomp<k>::divexp(lhsFirst, i);
                    lhsSecondDiv = precomp<k>::divexp(lhsSecond, i);
                    rhsFirstDiv = precomp<k>::divexp(rhsFirst, i);
                    rhsSecondDiv = precomp<k>::divexp(rhsSecond, i);

                    if (lhsFirstDiv < rhsFirstDiv) {
                        return true;
                    } else if (lhsFirstDiv == rhsFirstDiv && lhsSecondDiv > rhsSecondDiv) {
                        return false;
                    }
                }

                return true;

                //throw std::logic_error("Shouldn't be possible");//maybe true

            } else { //lhs.first > rhs.first && lhs.second < rhs.second

                t_x lhsFirst = lhs.first;
                t_x lhsSecond = lhs.second;
                t_x rhsFirst = rhs.first;
                t_x rhsSecond = rhs.second;

                t_x lhsFirstDiv;
                t_x lhsSecondDiv;

                t_x rhsFirstDiv;
                t_x rhsSecondDiv;

                for (int i = m_tree_height; i > 0; --i) {
                    lhsFirstDiv = precomp<k>::divexp(lhsFirst, i);
                    lhsSecondDiv = precomp<k>::divexp(lhsSecond, i);
                    rhsFirstDiv = precomp<k>::divexp(rhsFirst, i);
                    rhsSecondDiv = precomp<k>::divexp(rhsSecond, i);

                    if (lhsFirstDiv > rhsFirstDiv) {
                        return false;
                    } else if (lhsFirstDiv == rhsFirstDiv && lhsSecondDiv < rhsSecondDiv) {
                        return true;
                    }
                }
                return false;
                //throw std::logic_error("Shouldn't be possible");//maybe false
            }
        }


            /*
            inline bool operator()(const std::pair<uint32_t, uint32_t> lhs, const std::pair<uint32_t, uint32_t> rhs) {
                uint lhsInterleaved = interleaveLowerBits(lhs.second >> 16, lhs.first >> 16);
                uint rhsInterleaved = interleaveLowerBits(rhs.second >> 16, rhs.first >> 16);

                if (lhsInterleaved < rhsInterleaved) {
                    return true;
                } else if (lhsInterleaved > rhsInterleaved) {
                    return false;
                } else {
                    //interleaveLowerBits (shift away upper 16 Bits
                    lhsInterleaved = interleaveLowerBits((lhs.second << 16) >> 16, (lhs.first << 16) >> 16);
                    rhsInterleaved = interleaveLowerBits((rhs.second << 16) >> 16, (rhs.first << 16) >> 16);
                    if (lhsInterleaved < rhsInterleaved) {
                        return true;
                    } else if (lhsInterleaved > rhsInterleaved) {
                        return false;
                    } else {
                        return true;
                    }
                }
            }

            inline bool operator()(const std::pair<uint64_t, uint64_t>& lhs, const std::pair<uint64_t, uint64_t>& rhs) {
                const uint64_t lhsFirst = lhs.first;
                const uint64_t lhsSecond = lhs.second;
                const uint64_t rhsFirst = rhs.first;
                const uint64_t rhsSecond = rhs.second;
                if (lhsFirst <= std::numeric_limits<uint32_t>::max() && lhsSecond <= std::numeric_limits<uint32_t>::max() && rhsFirst <= std::numeric_limits<uint32_t>::max() && rhsSecond <= std::numeric_limits<uint32_t>::max()){
                    return this->operator()((std::pair<uint32_t, uint32_t>) lhs, (std::pair<uint32_t, uint32_t>) rhs);
                } else {
                    //can't call 32 bit as it returns true in the equal case for the upper 32 bit
                    uint32_t lhsFirst1 = (uint32_t)(lhsFirst >>32);
                    uint32_t rhsFirst1 = (uint32_t)(rhsFirst >>32);
                    uint32_t lhsSecond1 = (uint32_t)(lhsSecond >>32);
                    uint32_t rhsSecond1 = (uint32_t)(rhsSecond >>32);

                    if (lhsFirst1 == rhsFirst1 && lhsSecond1 == rhsSecond1){
                        //if first 32 bit of each of the two numbers of each pair are the same
                        //only look at lower 32 bit by implicit conversion of 64 bit integer to 32 bit
                        return this->operator()((std::pair<uint32_t, uint32_t>) lhs, (std::pair<uint32_t, uint32_t>) rhs);
                    } else {
                        return this->operator()(std::make_pair(lhsFirst1, lhsSecond1), std::make_pair(rhsFirst1, rhsSecond1));
                    }
                }
            }

        };*/

        /*
        //TODO: add more preshifted Morton Tables
        const int sort_by_z_order::MortonTable256[256] =
                {
                        0x0000, 0x0001, 0x0004, 0x0005, 0x0010, 0x0011, 0x0014, 0x0015,
                        0x0040, 0x0041, 0x0044, 0x0045, 0x0050, 0x0051, 0x0054, 0x0055,
                        0x0100, 0x0101, 0x0104, 0x0105, 0x0110, 0x0111, 0x0114, 0x0115,
                        0x0140, 0x0141, 0x0144, 0x0145, 0x0150, 0x0151, 0x0154, 0x0155,
                        0x0400, 0x0401, 0x0404, 0x0405, 0x0410, 0x0411, 0x0414, 0x0415,
                        0x0440, 0x0441, 0x0444, 0x0445, 0x0450, 0x0451, 0x0454, 0x0455,
                        0x0500, 0x0501, 0x0504, 0x0505, 0x0510, 0x0511, 0x0514, 0x0515,
                        0x0540, 0x0541, 0x0544, 0x0545, 0x0550, 0x0551, 0x0554, 0x0555,
                        0x1000, 0x1001, 0x1004, 0x1005, 0x1010, 0x1011, 0x1014, 0x1015,
                        0x1040, 0x1041, 0x1044, 0x1045, 0x1050, 0x1051, 0x1054, 0x1055,
                        0x1100, 0x1101, 0x1104, 0x1105, 0x1110, 0x1111, 0x1114, 0x1115,
                        0x1140, 0x1141, 0x1144, 0x1145, 0x1150, 0x1151, 0x1154, 0x1155,
                        0x1400, 0x1401, 0x1404, 0x1405, 0x1410, 0x1411, 0x1414, 0x1415,
                        0x1440, 0x1441, 0x1444, 0x1445, 0x1450, 0x1451, 0x1454, 0x1455,
                        0x1500, 0x1501, 0x1504, 0x1505, 0x1510, 0x1511, 0x1514, 0x1515,
                        0x1540, 0x1541, 0x1544, 0x1545, 0x1550, 0x1551, 0x1554, 0x1555,
                        0x4000, 0x4001, 0x4004, 0x4005, 0x4010, 0x4011, 0x4014, 0x4015,
                        0x4040, 0x4041, 0x4044, 0x4045, 0x4050, 0x4051, 0x4054, 0x4055,
                        0x4100, 0x4101, 0x4104, 0x4105, 0x4110, 0x4111, 0x4114, 0x4115,
                        0x4140, 0x4141, 0x4144, 0x4145, 0x4150, 0x4151, 0x4154, 0x4155,
                        0x4400, 0x4401, 0x4404, 0x4405, 0x4410, 0x4411, 0x4414, 0x4415,
                        0x4440, 0x4441, 0x4444, 0x4445, 0x4450, 0x4451, 0x4454, 0x4455,
                        0x4500, 0x4501, 0x4504, 0x4505, 0x4510, 0x4511, 0x4514, 0x4515,
                        0x4540, 0x4541, 0x4544, 0x4545, 0x4550, 0x4551, 0x4554, 0x4555,
                        0x5000, 0x5001, 0x5004, 0x5005, 0x5010, 0x5011, 0x5014, 0x5015,
                        0x5040, 0x5041, 0x5044, 0x5045, 0x5050, 0x5051, 0x5054, 0x5055,
                        0x5100, 0x5101, 0x5104, 0x5105, 0x5110, 0x5111, 0x5114, 0x5115,
                        0x5140, 0x5141, 0x5144, 0x5145, 0x5150, 0x5151, 0x5154, 0x5155,
                        0x5400, 0x5401, 0x5404, 0x5405, 0x5410, 0x5411, 0x5414, 0x5415,
                        0x5440, 0x5441, 0x5444, 0x5445, 0x5450, 0x5451, 0x5454, 0x5455,
                        0x5500, 0x5501, 0x5504, 0x5505, 0x5510, 0x5511, 0x5514, 0x5515,
                        0x5540, 0x5541, 0x5544, 0x5545, 0x5550, 0x5551, 0x5554, 0x5555
                };
           */




        template<typename t_x, typename t_y>
        void construct_bottom_up(std::vector<std::pair<t_x, t_y>> &links, std::string temp_file_prefix = "") {
            using namespace k2_treap_ns;
            using t_e = std::pair<t_x, t_y>;

            m_size = links.size();
            m_tree_height = get_tree_height(links);
            uint64_t M = precomp<t_k>::exp(t);
            t_e MM = t_e(M, M);

            std::string id_part = util::to_string(util::pid())
                                  + "_" + util::to_string(util::id());

            m_level_begin_idx = int_vector<64>(1 + t, 0);

            std::cout << "Sorting By Z Order" << std::endl;
            std::sort(links.begin(), links.end(), [&](const std::pair<t_x, t_y>& lhs, const std::pair<t_x, t_y>& rhs){
                return sort_by_z_order(lhs, rhs);
            });

            std::cout << "Sorting Finished, Constructing Bitvectors" << std::endl;
            std::vector<int> previous_subtree_number(m_tree_height,-1);
            uint64_t total_size = 0;

            {
                int subtree_distance;
                bool fill_to_k2_entries = false; //begin extra case!
                std::vector<uint> gap_to_k2(m_tree_height,k*k);
                bool firstLink = true;
                uint current_subtree_number = 0;

                std::vector<int_vector_buffer<1>> level_vectors;
                for (int i = 0; i < m_tree_height; i++) {
                    std::string bp_file = temp_file_prefix + "_bp_" + id_part + "_" + std::to_string(i) + ".sdsl";
                    level_vectors.push_back(int_vector_buffer<1>(bp_file, std::ios::out));
                }

                std::pair<t_x,t_y> previous_link;
                for (auto &current_link: links) {
                    auto tmp = std::make_pair(current_link.first, current_link.second);
                    for (uint current_level = 0; current_level < m_tree_height; ++current_level) {
                        current_subtree_number = calculate_subtree_number_and_new_relative_coordinates(current_link,
                                                                                                       current_level);
                        subtree_distance = current_subtree_number - previous_subtree_number[current_level];

                        if (subtree_distance > 0) {
                            //invalidate previous subtree numbers as new relative frame
                            for (uint i = current_level + 1; i < m_tree_height; ++i) {
                                previous_subtree_number[i] = -1;
                            }

                            if (fill_to_k2_entries && current_level != 0) {
                                for (uint j = 0; j < gap_to_k2[current_level]; ++j) {
                                    level_vectors[current_level].push_back(0);
                                }
                                gap_to_k2[current_level] = k * k;
                            }

                            for (int j = 0; j < subtree_distance - 1; ++j) {
                                level_vectors[current_level].push_back(0);
                                gap_to_k2[current_level]--;
                            }

                            level_vectors[current_level].push_back(1);
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
                        level_vectors[l].push_back(0);
                    }
                }

                for (uint m = 0; m < level_vectors.size(); ++m) {
                    total_size+= level_vectors[m].size();
                }
            }

            std::cout << "Construction finished, concatenating bit vectors" << std::endl;

            bit_vector concat = bit_vector(total_size,0);

            uint64_t level_begin_offset = 0;
            for (uint i = 0; i < m_tree_height; i++) {
                std::string bp_file = temp_file_prefix + "_bp_" + id_part + "_" + std::to_string(i) + ".sdsl";
                //std::cout << "Reading from " << temp_file_prefix << "_bp_" << id_part << "_" << std::to_string(i) << ".sdsl" << std::endl;
                {
                    bit_vector bp;
                    load_from_file(bp, bp_file);

                    m_level_begin_idx[m_tree_height-i-1] = level_begin_offset;
                    //copy values using old total length as offset
                    //FIXME: Probably pre calculate vector size and/or introduce append for int_vectors in sdsl
                    for (uint j = 0; j < bp.size(); ++j) {
                        //std::cout << bp[j] << "\t";
                        concat[j+level_begin_offset] = bp[j];
                    }
                    //std::cout << std::endl;
                    level_begin_offset += bp.size();
                    //std::cout << "Bp size" << bp.size() << std::endl;
                }
            }

            std::cout << "Concatenation Finished" << std::endl;

            {
                bit_vector _bp;
                _bp.swap(concat);
                m_bp = t_bv(_bp);
            }

            std::cout << "m_bp set" << std::endl;

            util::init_support(m_bp_rank, &m_bp);
            for (int i = 0; i < m_tree_height; i++) {
                std::string bp_file = temp_file_prefix + "_bp_" + id_part + "_" + std::to_string(i) + ".sdsl";
                sdsl::remove(bp_file);
            }

            std::cout << "initialized rank support" << std::endl;
        }


        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream &out, structure_tree_node *v = nullptr,
                            std::string name = "") const {
            structure_tree_node *child = structure_tree::add_child(
                    v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_tree_height, out, child, "t");
            written_bytes += write_member(m_size, out, child, "s");
            written_bytes += m_bp.serialize(out, child, "bp");
            written_bytes += m_bp_rank.serialize(out, child, "bp_rank");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream &in) {
            read_member(m_tree_height, in);
            read_member(m_size, in);
            m_bp.load(in);
            m_bp_rank.load(in);
            m_bp_rank.set_vector(&m_bp);
        }

        bool
        is_leaf(const node_type &v) const {
            //FIXME: check if this works
            return v.idx >= m_level_begin_idx[0];
        }

        node_type root() const {
            return node_type(t, t_p(0, 0), 0);
        }

        template<typename t_x>
        void direct_links(t_x source_id, std::vector<t_x> &result) const {
            result.clear();
            traverse_tree<t_x, DirectImpl>(this->root(), source_id, result);
        }

        template<typename t_x>
        void inverse_links(t_x source_id, std::vector<t_x> &result) const {
            result.clear();
            traverse_tree<t_x, InverseImpl>(this->root(), source_id, result);
        }

        template<typename t_x, class Impl>
        void traverse_tree(node_type root, t_x source_id, std::vector<t_x> &result) const {
            using namespace k2_treap_ns;
            if (!is_leaf(root)) {
                uint64_t rank = m_bp_rank(root.idx);
                auto x = std::real(root.p);
                auto y = std::imag(root.p);

                for (size_t i = 0; i < t_k; ++i) {
                    for (size_t j = 0; j < t_k; ++j) {
                        // get_int better for compressed bitvectors
                        // or introduce cache for bitvectors
                        if (m_bp[root.idx + t_k * i + j]) { //if subtree present
                            ++rank;
                            auto _x = x + i * precomp<t_k>::exp(root.t - 1);
                            auto _y = y + j * precomp<t_k>::exp(root.t - 1);

                            node_type subtree_root(root.t - 1, t_p(_x, _y), rank * t_k * t_k);
                            if (Impl::is_relevant_subtree(source_id, subtree_root)) {
                                traverse_tree<t_x, Impl>(subtree_root, source_id, result);
                            }
                        }
                    }
                }
            } else {
                //add corresponding values to result
                auto x = std::real(root.p);
                auto y = std::imag(root.p);

                for (size_t i = 0; i < t_k; ++i) {
                    auto _x = (x + i * precomp<t_k>::exp(root.t - 1));
                    for (size_t j = 0; j < t_k; ++j) {
                        if (m_bp[root.idx + t_k * i + j]) {
                            Impl::add_to_result_if_relevant(source_id, y, _x, j, root, result);
                        }
                    }
                }
            }
        }

        struct DirectImpl {
            template<typename t_x>
            inline static bool is_relevant_subtree(t_x row_id, node_type subtree_root) {
                using namespace k2_treap_ns;
                uint64_t d = precomp<t_k>::exp(subtree_root.t) - 1;

                return row_id >= real(subtree_root.p) and row_id <= real(subtree_root.p) + d;
            }

            template<typename t_x>
            inline static void add_to_result_if_relevant(t_x source_id, point_type::value_type y, point_type::value_type _x,
                                                         size_t j, node_type root, std::vector<t_x> &result) {
                using namespace k2_treap_ns;
                if (source_id == _x) { //if bit set and leaf part of correct row, add to result
                    auto _y = y + j * precomp<t_k>::exp(root.t - 1);
                    result.push_back(_y);
                }
            }
        };

        struct InverseImpl {
            template<typename t_x>
            inline static bool is_relevant_subtree(t_x column_id, node_type subtree_root) {
                using namespace k2_treap_ns;
                uint64_t d = precomp<t_k>::exp(subtree_root.t) - 1;

                return column_id >= imag(subtree_root.p) and column_id <= imag(subtree_root.p) + d;
            }

            template<typename t_x>
            inline static void add_to_result_if_relevant(t_x source_id, point_type::value_type y, point_type::value_type _x,
                                                         size_t j, node_type root, std::vector<t_x> &result) {
                using namespace k2_treap_ns;
                auto _y = y + j * precomp<t_k>::exp(root.t - 1);
                if (source_id == _y) { //if bit set and leaf part of correct row, add to result
                    result.push_back(_x);
                }
            }
        };

        //use only for testing purposes
        void set_height(uint height){
            m_tree_height = height;
        }
    };
}
#endif