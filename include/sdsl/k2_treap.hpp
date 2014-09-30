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

#include "sdsl/vectors.hpp"
#include "sdsl/bits.hpp"
#include "sdsl/k2_treap_helper.hpp"
#include "sdsl/k2_treap_algorithm.hpp"
#include <tuple>
#include <algorithm>
#include <climits>
#include <vector>

//! Namespace for the succinct data structure library.
namespace sdsl
{

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
template<uint8_t  t_k,
         typename t_bv=bit_vector,
         typename t_rank=typename t_bv::rank_1_type,
         typename t_max_vec=dac_vector<>>
class k2_treap
{
        static_assert(t_k>1, "t_k has to be larger than 1.");
        static_assert(t_k<=16, "t_k has to be smaller than 17.");

    public:
        typedef int_vector<>::size_type size_type;
        using node_type = k2_treap_ns::node_type;
        using point_type = k2_treap_ns::point_type;
        using t_p = k2_treap_ns::t_p;

        enum { k = t_k };

    private:
        uint8_t                   m_t = 0;
        t_bv                      m_bp;
        t_rank                    m_bp_rank;
        t_max_vec                 m_maxval;
        std::vector<int_vector<>> m_coord;
        int_vector<64>            m_level_idx;

        template<typename t_tv>
        uint8_t get_t(const t_tv& v)
        {
            using namespace k2_treap_ns;
            if (v.size() == 0) {
                return 0;
            }
            using t_e = typename t_tv::value_type;
            auto tupmax = [](t_e a) {
                return std::max(std::get<0>(a),std::get<1>(a));
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
        uint8_t& t = m_t;

        k2_treap() = default;

        k2_treap(const k2_treap& tr)
        {
            *this = tr;
        }

        k2_treap(k2_treap&& tr)
        {
            *this = std::move(tr);
        }

        //! Move assignment operator
        k2_treap& operator=(k2_treap&& tr)
        {
            if (this != &tr) {
                m_t = tr.m_t;
                m_bp = std::move(tr.m_bp);
                m_bp_rank = std::move(tr.m_bp_rank);
                m_bp_rank.set_vector(&m_bp);
                m_maxval = std::move(tr.m_maxval);
                m_coord = std::move(tr.m_coord);
                m_level_idx = std::move(tr.m_level_idx);
            }
            return *this;
        }

        //! Assignment operator
        k2_treap& operator=(k2_treap& tr)
        {
            if (this != &tr) {
                m_t = tr.m_t;
                m_bp = tr.m_bp;
                m_bp_rank = tr.m_bp_rank;
                m_bp_rank.set_vector(&m_bp);
                m_maxval = tr.m_maxval;
                m_coord = tr.m_coord;
                m_level_idx = tr.m_level_idx;
            }
            return *this;
        }

        //! Number of points in the 2^k treap
        size_type
        size() const
        {
            return m_maxval.size();
        }

        //! Swap operator
        void swap(k2_treap& tr)
        {
            if (this != &tr) {
                std::swap(m_t, tr.m_t);
                m_bp.swap(tr.m_bp);
                util::swap_support(m_bp_rank, tr.m_bp_rank, &m_bp, &(tr.m_bp));
                m_maxval.swap(tr.m_maxval);
                m_coord.swap(tr.m_coord);
                m_level_idx.swap(tr.m_level_idx);
            }
        }

        k2_treap(int_vector_buffer<>& buf_x,
                 int_vector_buffer<>& buf_y,
                 int_vector_buffer<>& buf_w)
        {
            using namespace k2_treap_ns;
            typedef int_vector_buffer<>* t_buf_p;
            std::vector<t_buf_p> bufs = {&buf_x, &buf_y, &buf_w};

            auto max_element = [](int_vector_buffer<>& buf) {
                uint64_t max_val = 0;
                for (auto val : buf) {
                    max_val = std::max((uint64_t)val, max_val);
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
                auto v = read<uint32_t,uint32_t,uint32_t>(bufs);
                construct(v, buf_x.filename());
            } else {
                auto v = read<uint64_t,uint64_t,uint64_t>(bufs);
                construct(v, buf_x.filename());
            }
        }

        template<typename t_x=uint64_t, typename t_y=uint64_t, typename t_w=uint64_t>
        std::vector<std::tuple<t_x, t_y, t_w>>
                                            read(std::vector<int_vector_buffer<>*>& bufs)
        {
            typedef std::vector<std::tuple<t_x, t_y, t_w>> t_tuple_vec;
            t_tuple_vec v = t_tuple_vec(bufs[0]->size());
            for (uint64_t j=0; j<v.size(); ++j) {
                std::get<0>(v[j]) = (*(bufs[0]))[j];
            }
            for (uint64_t j=0; j<v.size(); ++j) {
                std::get<1>(v[j]) = (*(bufs[1]))[j];
            }
            for (uint64_t j=0; j<v.size(); ++j) {
                std::get<2>(v[j]) = (*(bufs[2]))[j];
            }
            return v;
        }


        template<typename t_x, typename t_y, typename t_w>
        k2_treap(std::vector<std::tuple<t_x, t_y, t_w>>& v, std::string temp_file_prefix="")
        {
            if (v.size() > 0) {
                construct(v, temp_file_prefix);
            }
        }

        template<typename t_x, typename t_y, typename t_w>
        void construct(std::vector<std::tuple<t_x, t_y, t_w>>& v, std::string temp_file_prefix="")
        {
            using namespace k2_treap_ns;
            using t_e = std::tuple<t_x, t_y, t_w>;
            m_t = get_t(v);
            uint64_t M = precomp<t_k>::exp(t);
            t_e MM = t_e(M,M,M);

            std::string id_part = util::to_string(util::pid())
                                  + "_" + util::to_string(util::id());

            m_coord.resize(t);
            m_level_idx = int_vector<64>(1+t, 0);

            std::string val_file =  temp_file_prefix + "_k2_treap_"
                                    + id_part + ".sdsl";
            std::string bp_file  = temp_file_prefix + "_bp_" + id_part
                                   + ".sdsl";

            {
                int_vector_buffer<> val_buf(val_file, std::ios::out);
                int_vector_buffer<1> bp_buf(bp_file, std::ios::out);

                auto end = std::end(v);
                uint64_t last_level_nodes = 1;
                uint64_t level_nodes;
                for (uint64_t l=t, cc=0; l+1 > 0; --l) {
                    if (l > 0) {
                        m_level_idx[l-1] = m_level_idx[l] + last_level_nodes;
                        m_coord[l-1] = int_vector<>(2*last_level_nodes,0, bits::hi(precomp<t_k>::exp(l))+1);
                    }
                    level_nodes = 0;
                    cc = 0;
                    auto sp = std::begin(v);
                    for (auto ep = sp; ep != end;) {
                        ep = std::find_if(sp, end, [&sp,&l](const t_e& e) {
                            auto x1 = std::get<0>(*sp);
                            auto y1 = std::get<1>(*sp);
                            auto x2 = std::get<0>(e);
                            auto y2 = std::get<1>(e);
                            return    precomp<t_k>::divexp(x1,l) != precomp<t_k>::divexp(x2,l)
                                      or precomp<t_k>::divexp(y1,l) != precomp<t_k>::divexp(y2,l);
                        });
                        auto max_it = std::max_element(sp, ep, [](t_e a, t_e b) {
                            if (std::get<2>(a) != std::get<2>(b))
                                return std::get<2>(a) < std::get<2>(b);
                            else if (std::get<0>(a) != std::get<0>(b))
                                return std::get<0>(a) > std::get<0>(b);
                            return std::get<1>(a) > std::get<1>(b);
                        });
                        if (l > 0) {
                            m_coord[l-1][2*cc]   = precomp<t_k>::modexp(std::get<0>(*max_it), l);
                            m_coord[l-1][2*cc+1] = precomp<t_k>::modexp(std::get<1>(*max_it), l);
                            ++cc;
                        }

                        val_buf.push_back(std::get<2>(*max_it));
                        *max_it = MM;
                        --ep;
                        std::swap(*max_it, *ep);
                        if (l > 0) {
                            auto _sp = sp;

                            for (uint8_t i=0; i < t_k; ++i) {
                                auto _ep = ep;
                                if (i+1 < t_k) {
                                    _ep = std::partition(_sp, ep, [&i,&l](const t_e& e) {
                                        return precomp<t_k>::divexp(std::get<0>(e),l-1)%t_k <= i;
                                    });
                                }
                                auto __sp = _sp;
                                for (uint8_t j=0; j < t_k; ++j) {
                                    auto __ep = _ep;
                                    if (j+1 < t_k) {
                                        __ep = std::partition(__sp, _ep, [&j,&l](const t_e& e) {
                                            return precomp<t_k>::divexp(std::get<1>(e),l-1)%t_k <= j;
                                        });
                                    }
                                    bool not_empty = __ep > __sp;
                                    bp_buf.push_back(not_empty);
                                    level_nodes += not_empty;
                                    __sp = __ep;
                                }
                                _sp = _ep;
                            }
                        }
                        ++ep;
                        sp = ep;
                    }
                    end = std::remove_if(begin(v), end, [&](t_e e) {
                        return e == MM;
                    });
                    last_level_nodes = level_nodes;
                }
            }
            bit_vector bp;
            load_from_file(bp, bp_file);
            {
                int_vector_buffer<> val_rw(val_file, std::ios::in | std::ios::out);
                int_vector_buffer<> val_r(val_file, std::ios::in);
                uint64_t bp_idx = bp.size();
                uint64_t r_idx = m_level_idx[0];
                uint64_t rw_idx = val_rw.size();
                while (bp_idx > 0) {
                    --r_idx;
                    for (size_t i=0; i < t_k*t_k; ++i) {
                        if (bp[--bp_idx]) {
                            --rw_idx;
                            val_rw[rw_idx] = val_r[r_idx] - val_rw[rw_idx];
                        }
                    }
                }
            }
            {
                int_vector_buffer<> val_r(val_file);
                m_maxval = t_max_vec(val_r);
            }
            {
                bit_vector _bp;
                _bp.swap(bp);
                m_bp = t_bv(_bp);
            }
            util::init_support(m_bp_rank, &m_bp);
            sdsl::remove(bp_file);
            sdsl::remove(val_file);
        }


        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr,
                            std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(
                                             v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_t, out, child, "t");
            written_bytes += m_bp.serialize(out, child, "bp");
            written_bytes += m_bp_rank.serialize(out, child, "bp_rank");
            written_bytes += serialize_vector(m_coord, out, child, "coord");
            written_bytes += m_maxval.serialize(out, child, "maxval");
            written_bytes += m_level_idx.serialize(out, child, "level_idx");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in)
        {
            read_member(m_t, in);
            m_bp.load(in);
            m_bp_rank.load(in);
            m_bp_rank.set_vector(&m_bp);
            m_coord.resize(t);
            load_vector(m_coord, in);
            m_maxval.load(in);
            m_level_idx.load(in);
        }

        node_type
        root() const
        {
            return node_type(t, t_p(0,0), 0, m_maxval[0],
                             t_p(m_coord[t-1][0], m_coord[t-1][1]));
        }

        bool
        is_leaf(const node_type& v) const
        {
            return v.idx >= m_bp.size();
        }

        std::vector<node_type>
        children(const node_type& v) const
        {
            using namespace k2_treap_ns;
            std::vector<node_type> res;
            if (!is_leaf(v)) {
                uint64_t rank = m_bp_rank(v.idx);
                auto x = std::real(v.p);
                auto y = std::imag(v.p);

                for (size_t i=0; i<t_k; ++i) {
                    for (size_t j=0; j<t_k; ++j) {
                        // get_int better for compressed bitvectors
                        // or introduce cache for bitvectors
                        if (m_bp[v.idx+t_k*i+j]) {
                            ++rank;

                            auto _x = x + i*precomp<t_k>::exp(v.t-1);
                            auto _y = y + j*precomp<t_k>::exp(v.t-1);

                            auto _max_v = v.max_v - m_maxval[rank];
                            auto _max_p = t_p(_x, _y);
                            if (v.t > 1) {
                                auto y = rank-m_level_idx[v.t-1];
                                _max_p = t_p(_x+m_coord[v.t-2][2*y],
                                             _y+m_coord[v.t-2][2*y+1]);
                            }
                            res.emplace_back(v.t-1, t_p(_x,_y), rank*t_k*t_k,
                                             _max_v, _max_p);
                        }
                    }
                }
            }
            return res;
        }

};

}
#endif
