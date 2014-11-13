/* sdsl - succinct data structures library
    Copyright (C) 2013 Simon Gog

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
/*! \file wt_topk.hpp
    \brief wt_topk.hpp contains a class for the wavelet tree plus rmq data structure.
            This data structure is based on the solution for solving top-k queries on grids by
            G.Navarro, Y. Nekrich and L. Russo, Space-Efficient Data-Analysis Queries on Grids.
            Theoretical Computer Science 482:60-72, 2013

    \author Simon Gog, Roberto Konow
*/

#ifndef INCLUDED_SDSL_WT_TOPK
#define INCLUDED_SDSL_WT_TOPK

#include "sdsl/vectors.hpp"
#include "sdsl/wavelet_trees.hpp"
#include "sdsl/rmq_support.hpp"
#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>
#include <complex>

using namespace std;
using namespace sdsl;

namespace sdsl
{



template<typename t_wt=wt_int<>,
         typename t_rmq=rmq_succinct_sct<false>,
         typename t_weight_vec=dac_vector<>
         >
class wt_topk
{
    public:
        typedef std::complex<uint64_t> point_type;
        typedef int_vector<>::size_type size_type;
        typedef tuple<size_type, size_type, size_type> range_type;
        typedef pair<size_type,size_type> result_pair;
        typedef priority_queue<result_pair, vector<result_pair>, std::greater<result_pair>> result_type;

        class top_k_iterator
        {
            public:
                typedef void(*t_mfptr)();
                typedef std::pair<point_type, uint64_t> t_point_val;
                // state consists of weight, range, index of minimum in range
                typedef std::tuple<uint64_t, range_type, uint64_t> t_state;

            private:
                const wt_topk*                m_topk;
                const priority_queue<t_state> m_pq;
                t_point_val                   m_point_val;
                bool                          m_valid = false;
            public:
                top_k_iterator() = default;
                top_k_iterator(const top_k_iterator&) = default;
                top_k_iterator(top_k_iterator&&) = default;
                top_k_iterator& operator=(const top_k_iterator&) = default;
                top_k_iterator& operator=(top_k_iterator&&) = default;
                top_k_iterator(const wt_topk& topk, point_type p1, point_type p2) :
                    m_topk(&topk), m_valid(topk.size()>0)
                {
                    // get all range from topk->wt
                    // for each range compute the position of the maximum weight;
                    // get the weight weight and push the tuple
                    // (maximum weight, range, pos of maximum weight) to the
                    // priority queue

                    // call ++iterator
                    if (topk.size() > 0) {
                        ++(*this);
                    }
                }

                //! Prefix increment of the iterator
                top_k_iterator& operator++()
                {
                    m_valid = false;
                    if (!m_pq.empty()) {
                        // pop top element
                        // calc x pos and store result in m_point_val
                        // calculate the two subranges
                        // if not empty: get position of maximum weight
                        // and the weight and push new tuple to the queue
                    }
                    return *this;
                }

                //! Postfix increment of the iterator
                top_k_iterator operator++(int)
                {
                    top_k_iterator it = *this;
                    ++(*this);
                    return it;
                }

                t_point_val operator*() const
                {
                    return m_point_val;
                }

                operator t_mfptr() const
                {
                    return (t_mfptr)(m_valid);
                }
        };

    private:
        t_wt                m_wt;       // wavelet tree over y-values
        t_weight_vec        m_weights;  // array for weights
        t_rmq               m_rmq;      // range maximum query structure on top of the weights

    public:
        wt_topk() = default;

        wt_topk(const wt_topk& tr) = default;
        wt_topk(wt_topk&& tr) = default;
        wt_topk& operator=(const wt_topk& tr) = default;
        wt_topk& operator=(wt_topk&& tr) = default;

        //! Number of points in the grid
        size_type size() const
        {
            return m_wt.size();
        }

        void swap(wt_topk& tr)
        {
            if (this != &tr) {
                m_wt.swap(tr.m_wt);
                m_rmq.swap(tr.m_rmq);
                m_weights.swap(tr.m_weights);
            }
        }

        wt_topk(int_vector_buffer<>&,
                int_vector_buffer<>& buf_y,
                int_vector_buffer<>& buf_w)
        {
            std::string temp_file = buf_w.filename() +
                                    + "_wt_topk_" + to_string(util::pid())
                                    + "_" + to_string(util::id());
            // (1) Calculate permuted weight vector
            {
                int_vector<> perm = sorted_perm(buf_y);
                int_vector<> weights;
                load_from_file(weights, buf_w.filename());

                int_vector_buffer<> permuted_weights(temp_file, std::ios::out);
                for (size_type i=0; i<weights.size(); ++i) {
                    permuted_weights[i] = weights[perm[i]];
                }
            }
            // (2) Construct range maximum query structure over permuted weights
            {
                int_vector<> permuted_weights;
                load_from_file(permuted_weights, temp_file);
                m_rmq = t_rmq(&permuted_weights);
            }
            // (3) Construct the WT over y values
            m_wt = t_wt(buf_y, buf_y.size()); // construct the wavelet trees over the y-values
            // (4) Load the permuted weights and delete temporary file
            {
                int_vector_buffer<> permuted_weights(temp_file);
                m_weights = t_weight_vec(permuted_weights);
            }
            sdsl::remove(temp_file);
        }

        top_k_iterator topk(point_type p1, point_type p2) const
        {
            return top_k_iterator(*this, p1, p2);
        }


        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr,
                            std::string name="") const
        {
            structure_tree_node* child = structure_tree::add_child(
                                             v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_wt.serialize(out, child, "wt");
            written_bytes += m_rmq.serialize(out, child, "rmq");
            written_bytes += m_weights.serialize(out, child, "weights");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in)
        {
            m_wt.load(in);
            m_rmq.load(in);
            m_weights.load(in);
        }

        // main topk function, stores the result into the m_result class variable.
        void topk(const unsigned int k, const range_type& x_range, const range_type& y_range)
        {
            /*
                        m_result = result_type(); // clear the m_result;
                         // range_type = <left, right, symbol>
                        queue<range_type> ranges;
                        auto mts_it = map_to_sorted_sequence(m_wt, x_range, y_range);
                        size_t l, r = 0;

                        while (mts_it) { // maybe this can be improved by not storing the initial ranges into the queue?
                            l = get<0>(get<1>(*mts_it));
                            r = get<1>(get<1>(*mts_it));
                            ranges.emplace(l, r, get<0>(*mts_it));
                            ++mts_it;
                        }

                        size_t pos, w_value, x_value = 0;
                        size_t left_l, left_r, right_l, right_r = 0;
                        while (!ranges.empty()) {
                            range_type r = ranges.front();
                            ranges.pop();

                            pos = m_rmq(get<0>(r), get<1>(r));
                            w_value = m_weights[pos];
                            x_value = m_wt.select(pos-get<0>(r), get<2>(r));

                            if (m_result.size() < k) {
                                m_result.emplace(w_value, x_value);
                            } else {
                                if (get<0>(m_result.top()) < w_value) {
                                    m_result.pop();
                                    m_result.emplace(w_value, x_value);
                                }

                            }

                            left_l = get<0>(r);
                            left_r = pos - 1;
                            right_l = pos + 1;
                            right_r = get<1>(r);

                            if (left_l <= left_r and
                                    left_r != (size_t) -1) {
                                ranges.emplace(left_l, left_r, get<2>(r));
                            }

                            if (right_l <= right_r and
                                    right_r < m_wt.size()) {
                                ranges.emplace(right_l, right_r, get<2>(r));
                            }
                        }
            */
        }
        /*
            void print_results() {
                while (!m_result.empty()) {
                    cout << "( " << m_result.top().first << " , " << m_result.top().second << " ) " << endl;
                    m_result.pop();
                }
            }

            void get_result_reference(result_type &result) const {
                result = m_result;
            }
        */
        void print_info() const
        {
            std::cout<<"m_wt     ="<<m_wt<<std::endl;
            std::cout<<"m_weights="<<m_weights<<std::endl;
        }

    private:

        //! Returns a permutation P which is stable sorted according to data_buf
        /*!
         * I.e. data_buf[P[i]] <= data_buf[P[i+1]] and if
         * data_buf[P[i]] == data_buf[P[i+1]] then P[i] < P[i+1].
         */
        int_vector<> sorted_perm(int_vector_buffer<>& data_buf)
        {
            using value_type = int_vector<>::value_type;
            int_vector<> v;
            load_from_file(v, data_buf.filename());
            int_vector<> perm(v.size(), 0, bits::hi(v.size())+1);
            util::set_to_id(perm);
            sort(perm.begin(), perm.end(), [&v](const value_type& i,
            const value_type& j) {
                auto x = v[i], y = v[j];
                if (x == y)
                    return i < j;
                return x < y;
            });
            return std::move(perm);
        }

};
}

#endif
