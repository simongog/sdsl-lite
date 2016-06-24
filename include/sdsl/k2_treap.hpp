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
                m_tree_height = tr.m_size;
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
                m_tree_height = tr.m_size;
                m_bp = tr.m_bp;
                m_bp_rank = tr.m_bp_rank;
                m_bp_rank.set_vector(&m_bp);
                m_level_begin_idx = tr.m_level_begin_idx;
            }
            return *this;
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
                 int_vector_buffer<> &buf_y) {
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
                construct(v, buf_x.filename());
            } else {
                auto v = read < uint64_t, uint64_t>(bufs);
                construct(v, buf_x.filename());
            }
        }

        template<typename t_x=uint64_t, typename t_y=uint64_t>
        std::vector<std::pair<t_x, t_y>>
        read(std::vector<int_vector_buffer<> *> &bufs) {
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
        k2_treap(std::vector<std::pair<t_x, t_y>> &v, std::string temp_file_prefix = "") {
            if (v.size() > 0) {
                construct(v, temp_file_prefix);
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

                auto begin = std::begin(v);
                auto end = std::end(v);
                uint64_t last_level_bits = 0;
                uint64_t level_bits = 0;

                //recursively partition that stuff
                for (uint64_t l = t; l + 1 > 0; --l) {

                    if (l > 0) {
                        //std::cout << "Setting: " << "m_level_begin_idx[" << (l - 1) << "] = "
                        //         << m_level_begin_idx[l] + last_level_bits << std::endl;
                        m_level_begin_idx[l - 1] = m_level_begin_idx[l] + last_level_bits;
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

                            /*if (asd)
                            {
                                std::cout << "ep at " << e.first << "," << e.second << std::endl;
                            }*/
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
                                    /*std::cout << "After partition 0: \t";
                                    for (auto &pair: v) {
                                        std::cout << pair.first << "," << pair.second << "\t";
                                    }
                                    std::cout << std::endl;

                                    if (_ep < end) {
                                        std::cout << "Split at " << _ep.operator*().first << ","
                                                  << _ep.operator*().second << std::endl;
                                    } else {
                                        std::cout << "Split at end" << std::endl;
                                    }
                                    */

                                }
                                auto __sp = _sp;
                                for (uint8_t j = 0;
                                     j < t_k; ++j) { //partition the t_k vertical partitions t_k -1 times horizontally
                                    auto __ep = _ep;
                                    if (j + 1 < t_k) {
                                        __ep = std::partition(__sp, _ep, [&j, &l](const t_e &e) {
                                            return precomp<t_k>::divexp(std::get<1>(e), l - 1) % t_k <= j;
                                        });
                                        /*std::cout << "After partition 1: \t";
                                        for (auto &pair: v) {
                                            std::cout << pair.first << "," << pair.second << "\t";
                                        }
                                        std::cout << std::endl;
                                        if (_ep < end) {
                                            std::cout << "Split at " << __ep.operator*().first << ","
                                                      << __ep.operator*().second << std::endl;
                                        } else {
                                            std::cout << "Split at end" << std::endl;
                                        }*/
                                    }
                                    bool not_empty = __ep > __sp;
                                    //std::cout << "Pushing " << not_empty << " to bp_buf" << std::endl;
                                    bp_buf.push_back(not_empty);
                                    level_bits++;
                                    __sp = __ep;
                                }
                                _sp = _ep;
                            }
                        }
                        if (ep!= end){
                            //++ep;
                            sp = ep;
                        }
                    }
                    last_level_bits = level_bits;
                    //std::cout << "Last Level Nodes: " << last_level_nodes << std::endl;
                }
            }

            bit_vector bp;
            load_from_file(bp, bp_file);
            {
                bit_vector _bp;
                _bp.swap(bp);
                m_bp = t_bv(_bp);
            }

            /*
            std::cout << "m_bp";
            for (size_t m = 0; m < m_bp.size(); ++m) {
                std::cout << m_bp[m];
            }
            std::cout << std::endl;
            */

            util::init_support(m_bp_rank, &m_bp);
            sdsl::remove(bp_file);
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
            return node_type(t, t_p(0,0), 0);
        }

        template<typename t_x>
        void direct_links(t_x source_id, std::vector<t_x>& result) const {
            result.clear();
            traverse_tree<t_x,DirectImpl>(this->root(), source_id, result);
        }

        template<typename t_x>
        void inverse_links(t_x source_id, std::vector<t_x>& result) const {
            result.clear();
            traverse_tree<t_x,InverseImpl>(this->root(), source_id, result);
        }

        template<typename t_x, class Impl>
        void traverse_tree(node_type root, t_x source_id, std::vector<t_x>& result) const{
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
                                traverse_tree<t_x, Impl>(subtree_root,source_id,result);
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
                        if (m_bp[root.idx + t_k * i + j]){
                            Impl::add_to_result_if_relevant(source_id, x, y, _x, i, j, root, result);
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
            inline static void add_to_result_if_relevant(t_x source_id, point_type::value_type x, point_type::value_type y, point_type::value_type _x, size_t i, size_t j, node_type root, std::vector<t_x>& result) {
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
            inline static void add_to_result_if_relevant(t_x source_id, point_type::value_type x, point_type::value_type y, point_type::value_type _x, size_t i, size_t j, node_type root, std::vector<t_x>& result) {
                using namespace k2_treap_ns;
                auto _y = y + j * precomp<t_k>::exp(root.t - 1);
                if (source_id == _y) { //if bit set and leaf part of correct row, add to result
                    result.push_back(_x);
                }
            }
        };

    };
}
#endif