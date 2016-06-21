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
/*! \file k2_treap_algorithm.hpp
    \brief k2_treap_algorithm.hpp contains k^2-treap algorithms.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_K2_TREAP_ALGORITHM
#define INCLUDED_SDSL_K2_TREAP_ALGORITHM

#include "vectors.hpp"
#include "bits.hpp"
#include "k2_treap_helper.hpp"
#include <tuple>
#include <algorithm>
#include <iterator>
#include <climits>
#include <vector>
#include <complex>
#include <queue>
#include <array>

//! Namespace for the succinct data structure library.
namespace sdsl {

    namespace k2_treap_ns {

//! Check if point x is contained in the rectangle (p1,p2)
/*! \param p Point.
 *  \param Lower left corner of the rectangle.
 *  \param Upper right corner of the rectangle.
 */
        bool
        contained(const point_type p, const point_type &p1, const point_type &p2) {
            return real(p) >= real(p1) and real(p) <= real(p2) and
                   imag(p) >= imag(p1) and imag(p) <= imag(p2);
        }

//! Check if the rectangle of node v is contained in the rectangle (p1,p2)
        template<uint8_t t_k>
        bool
        contained(const point_type &p1, const point_type &p2, const node_type &v) {
//    uint64_t d = (1ULL << v.t)-1;
//    uint64_t d = (1ULL << v.t)-1;
            uint64_t d = precomp<t_k>::exp(v.t) - 1;
            return real(p1) <= real(v.p) and real(p2) >= real(v.p) + d and
                   imag(p1) <= imag(v.p) and imag(p2) >= imag(v.p) + d;
        }

//! Check if rectangle (p1,p2) and the area of node v overlap
        template<uint8_t t_k>
        bool
        overlap(const point_type &p1, const point_type &p2, const node_type &v) {
//    uint64_t d = (1ULL << v.t)-1;
            uint64_t d = precomp<t_k>::exp(v.t) - 1;
            return real(p1) <= real(v.p) + d and real(p2) >= real(v.p) and
                   imag(p1) <= imag(v.p) + d and imag(p2) >= imag(v.p);
        }

        template<typename t_k2_treap>
        class range_iterator {
        public:
            typedef void(*t_mfptr)();

            typedef point_type t_point_val;

        private:
            typedef k2_treap_ns::node_type node_type;
            typedef std::pair<node_type, bool> t_nt_b;

            const t_k2_treap *m_treap = nullptr;
            std::priority_queue <t_nt_b> m_pq;
            t_point_val m_point_val;
            point_type m_p1;
            point_type m_p2;
            bool m_valid = false;

            void pq_emplace(node_type v, bool b) {
                m_pq.emplace(v, b);
            }

        public:
            range_iterator() = default;

            range_iterator(const range_iterator &) = default;

            range_iterator(range_iterator &&) = default;

            range_iterator &operator=(const range_iterator &) = default;

            range_iterator &operator=(range_iterator &&) = default;

            range_iterator(const t_k2_treap &treap, point_type p1, point_type p2) :
                    m_treap(&treap), m_p1(p1), m_p2(p2), m_valid(treap.size() > 0) {
                if (m_treap->size() > 0) {
                    pq_emplace(m_treap->root(), false);
                    ++(*this);
                }
            }

            //! Prefix increment of the iterator
            range_iterator &operator++() {
                m_valid = false;
                while (!m_pq.empty()) {
                    auto v = std::get<0>(m_pq.top());
                    auto is_contained = std::get<1>(m_pq.top());
                    m_pq.pop();
                    if (is_contained) {
                        auto nodes = m_treap->children(v);
                        for (auto node : nodes){
                            pq_emplace(node, true);
                        }

                        m_valid = true;
                        m_point_val = t_point_val(v.max_p);
                        break;
                    } else {
                        if (contained<t_k2_treap::k>(m_p1, m_p2, v)) {
                            m_pq.emplace(v, true);
                        } else if (overlap<t_k2_treap::k>(m_p1, m_p2, v)) {
                            auto nodes = m_treap->children(v);
                            for (auto node : nodes)
                                pq_emplace(node, false);
                            if (contained(v.max_p, m_p1, m_p2)) {
                                m_valid = true;
                                m_point_val = t_point_val(v.max_p);
                                break;
                            }
                        }
                    }
                }
                return *this;
            }

            //! Postfix increment of the iterator
            range_iterator operator++(int) {
                range_iterator it = *this;
                ++(*this);
                return it;
            }

            t_point_val operator*() const {
                return m_point_val;
            }

            //! Cast to a member function pointer
            // Test if there are more elements
            operator t_mfptr() const {
                return (t_mfptr) (m_valid);
            }
        };

    } // end namespace k2_treap_ns

//! Get iterator for all points in rectangle (p1,p2) with weights in range
/*! \param treap k2-treap
 *  \param p1    Lower left corner of the rectangle
 *  \param p2    Upper right corner of the rectangle
 *  \param range Range {w1,w2}.
 *  \return Iterator to list of all points in the range.
 *  \pre real(p1) <= real(p2) and imag(p1)<=imag(p2)
 *       real(range) <= imag(range)
 */
    template<typename t_k2_treap>
    k2_treap_ns::range_iterator<t_k2_treap>
    range_3d(const t_k2_treap &t,
             k2_treap_ns::point_type p1,
             k2_treap_ns::point_type p2) {
        return k2_treap_ns::range_iterator<t_k2_treap>(t, p1, p2);
    }


// forward declaration
    template<typename t_k2_treap>
    uint64_t __count(const t_k2_treap &, typename t_k2_treap::node_type);

// forward declaration
    template<typename t_k2_treap>
    uint64_t _count(const t_k2_treap &, k2_treap_ns::point_type,
                    k2_treap_ns::point_type, typename t_k2_treap::node_type);

//! Count how many points are in the rectangle (p1,p2)
/*! \param treap k2-treap
 *  \param p1    Lower left corner of the rectangle.
 *  \param p2    Upper right corner of the rectangle.
 *  \return The number of points in rectangle (p1,p2).
 *  \pre real(p1) <= real(p2) and imag(p1)<=imag(p2)
 */
    template<typename t_k2_treap>
    uint64_t
    count(const t_k2_treap &treap,
          k2_treap_ns::point_type p1,
          k2_treap_ns::point_type p2) {
        if (treap.size() > 0) {
            return _count(treap, p1, p2, treap.root());
        }
        return 0;
    }


    template<typename t_k2_treap>
    uint64_t
    _count(const t_k2_treap &treap,
           k2_treap_ns::point_type p1,
           k2_treap_ns::point_type p2,
           typename t_k2_treap::node_type v) {
        using namespace k2_treap_ns;
        if (contained<t_k2_treap::k>(p1, p2, v)) {
            return __count(treap, v);
        } else if (overlap<t_k2_treap::k>(p1, p2, v)) {
            uint64_t res = contained(v.max_p, p1, p2);
            auto nodes = treap.children(v);
            for (auto node : nodes) {
                res += _count(treap, p1, p2, node);
            }
            return res;
        }
        return 0;
    }


    template<typename t_k2_treap>
    uint64_t
    __count(const t_k2_treap &treap,
            typename t_k2_treap::node_type v) {
        uint64_t res = 1; // count the point at the node
        auto nodes = treap.children(v);
        for (auto node : nodes)
            res += __count(treap, node);
        return res;
    }


// forward declaration
    template<uint8_t t_k,
            typename t_bv,
            typename t_rank>
    class k2_treap;


//! Specialized version of method ,,construct'' for k2_treaps.
    template<uint8_t t_k,
            typename t_bv,
            typename t_rank>
    void
    construct(k2_treap<t_k, t_bv, t_rank> &idx, std::string file) {
        int_vector_buffer<> buf_x(file + ".x", std::ios::in);
        int_vector_buffer<> buf_y(file + ".y", std::ios::in);
        k2_treap<t_k, t_bv, t_rank> tmp(buf_x, buf_y);
        tmp.swap(idx);
    }

//! Specialized version of method ,,construct_im'' for k2_treaps.
    template<uint8_t t_k,
            typename t_bv,
            typename t_rank>
    void
    construct_im(k2_treap<t_k, t_bv, t_rank> &idx, std::vector <std::array<uint64_t, 2>> data) {
        std::string tmp_prefix = ram_file_name("k2_treap_");
        std::vector <std::pair<uint64_t, uint64_t>> d;
        for (auto x : data) {
            d.push_back(std::make_pair(x[0], x[1]));
        }
        k2_treap<t_k, t_bv, t_rank> tmp(d, tmp_prefix);
        tmp.swap(idx);
    }

//! Specialized version of method ,,construct_im'' for k2_treaps.
    template<uint8_t t_k,
            typename t_bv,
            typename t_rank>
    void
    construct_im(k2_treap<t_k, t_bv, t_rank> &idx, std::vector <std::pair<uint64_t, uint64_t>> data) {
        std::string tmp_prefix = ram_file_name("k2_treap_");
        tmp.swap(idx);
        k2_treap<t_k, t_bv, t_rank> tmp(data, tmp_prefix);
    }

}
#endif
