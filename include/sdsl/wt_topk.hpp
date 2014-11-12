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

using namespace std;
using namespace sdsl;

namespace sdsl {

template<typename t_wt=wt_int,
        typename t_rmq=rmq_succinct_sct<>,
        typename t_weight_vec=dac_vector<>
>
class wt_topk
{
    public:
        typedef int_vector<>::size_type size_type;
        typedef pair<size_t, size_t> range_type;
        typedef priority_queue<range_type, vector<range_type>, std::greater<range_type>> result_type;

    private:
//        typedef wt_int<> t_wt;
//        typedef rmq_succinct_sct<> t_rmq;
//        typedef dac_vector<> t_weight_vec;
        t_weight_vec        m_weights;
        t_wt                m_wt;
        result_type         m_result; // to be deleted
        t_rmq               m_rmq;


    public:
        wt_topk() = default;

        wt_topk(wt_topk& tr)
        {
            *this = tr;
        }

        wt_topk(wt_topk&& tr)
        {
            *this = std::move(tr);
        }

        //! Move assignment operator
        wt_topk& operator=(wt_topk&& tr)
        {
            if (this != &tr) {
                m_weights = std::move(tr.m_weights);
                m_wt = std::move(tr.m_wt);
                m_rmq = std::move(tr.m_rmq);
            }
            return *this;
        }

        wt_topk& operator=(wt_topk& tr)
        {
            if (this != &tr) {
                m_weights = tr.m_weights;
                m_wt = tr.m_wt;
                m_rmq = tr.m_rmq;
            }
            return *this;
        }


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

        wt_topk(int_vector_buffer<>& buf_x,
                int_vector_buffer<>& buf_y,
                int_vector_buffer<>& buf_w)
        {
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
            if (x <= std::numeric_limits<uint32_t>::max()) {
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
        wt_topk(std::vector<std::tuple<t_x, t_y, t_w>>& v, std::string temp_file_prefix="")
        {
            if (v.size() > 0) {
                construct(v, temp_file_prefix);
            }
        }

        template<typename t_x, typename t_y, typename t_w>
        void construct(std::vector<std::tuple<t_x, t_y, t_w>>& v, std::string temp_file_prefix="") {
            using t_e = std::tuple<t_x, t_y, t_w>;

            std::string id_part = util::to_string(util::pid())
                    + "_" + util::to_string(util::id());

            std::string val_file = temp_file_prefix + "_wt_topk_"
                    + id_part + ".sdsl";
            {
                int_vector_buffer<> val_buf(val_file, std::ios::out);
                // TODO: efficient construction here (not sure how to do this)
                int_vector<> depths(v.size());
                int_vector<> weights(v.size());
                size_type i = 0;
                for (const auto& sp : v) {
                    depths[i] = std::get<1>(sp);
                    weights[i] = std::get<2>(sp);
                    i++;
                }
                construct_im(m_wt, depths);
                sort(depths, weights, depths.size());
                m_weights = dac_vector<>(weights);
                m_rmq = rmq_succinct_sct<>(&weights);
            }
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
            read_member(m_size, in);
            m_wt.load(in);
            m_rmq.load(in);
            m_weights.load(in);
        }

        // main topk function, stores the result into the m_result class variable.

        void topk(const unsigned int k, const range_type &x_range, const range_type &y_range) {
            m_result = result_type(); // clear the m_result;
            queue<range_type> ranges;
            auto mts_it = map_to_sorted_sequence(m_wt, x_range, y_range);
            size_t l, r = 0;

            while (mts_it) { // maybe this can be improved by not storing the initial ranges into the queue?
                l = get<0>(get<1>(*mts_it));
                r = get<1>(get<1>(*mts_it));
                ranges.emplace(l, r);
                ++mts_it;
            }

            size_t pos, w_value, d_value = 0;
            size_t left_l, left_r, right_l, right_r = 0;
            while (!ranges.empty()) {

                range_type r = ranges.front();
                ranges.pop();

                pos = m_rmq(get<0>(r), get<1>(r));
                w_value = m_weights[pos];
                d_value = m_documents[pos];

                if (m_result.size() < k) {
                    m_result.emplace(w_value, d_value);
                } else {
                    if (get<0>(m_result.top()) < w_value) {
                        m_result.pop();
                        m_result.emplace(w_value, d_value);
                    }

                }

                left_l = get<0>(r);
                left_r = pos - 1;
                right_l = pos + 1;
                right_r = get<1>(r);

                if (left_l <= left_r and
                        left_r != (size_t) -1) {
                    ranges.emplace(left_l, left_r);
                }

                if (right_l <= right_r and
                        right_r < m_size) {
                    ranges.emplace(right_l, right_r);
                }
            }
        }

//    // performs a self_check by doing a exhaustive and brute-force traversal and compare
//    // the results against the ones obtained by the wt operations.
//    static bool self_check(const unsigned int k, int_vector<> &depths, int_vector<> &weights, int_vector<> &documents,
//            const range_type &x_range, const range_type &y_range) {
//
//        result_type result;
//        for (size_t x = get<0>(x_range); x <= get<1>(x_range); x++) {
//            if (depths[x] >= get<0>(y_range) and
//                    depths[x] <= get<1>(y_range)) {
//
//                if (result.size() < k) {
//                    result.emplace(weights[x], documents[x]);
//                }
//                else {
//                    if (get<0>(result.top()) < weights[x]) {
//                        result.pop();
//                        result.emplace(weights[x], documents[x]);
//                    }
//                }
//            }
//        }
//
//        wt_topk topkwt_test(depths, weights, documents);
//        topkwt_test.topk(k, x_range, y_range);
//        result_type other_result;
//        topkwt_test.get_result_reference(other_result);
//        assert(other_result.size() == result.size());
//        size_t i = 0;
//        while (!result.empty()) {
//            if (result.top().first != other_result.top().first) {
//                cerr << "ERROR AT " << i << endl;
//                cerr << "weight = \t " << result.top().first << " \t  other_result = " << other_result.top().first << endl;
//                return false;
//            }
//            i++;
//            result.pop();
//            other_result.pop();
//        }
//        return true;
//    }

    void print_results() {
        while (!m_result.empty()) {
            cout << "( " << m_result.top().first << " , " << m_result.top().second << " ) " << endl;
            m_result.pop();
        }
    }

    void get_result_reference(result_type &result) const {
        result = m_result;
    }


protected:
    // auxiliary function for the sort
    void permute_values(int_vector<> &array, int_vector<> &tmp, vector<pair<size_t, size_t> > &v, size_t n) {
        for (size_t i = 0; i < n; i++) {
            tmp[i] = array[v[i].second];
        }
        for (size_t i = 0; i < n; i++) {
            array[i] = tmp[i];
        }
    }
    // SG: replace by one stable sort using the right comparison function
        // RK: Done
    // SG: done keep all things in main memory during construction
        // RK: Do you have something for external sorting in sdsl already?

    // performs a sort of array and permutes array2 and array3 to be according to the sorted array.
    void sort(int_vector<> &array, int_vector<> &array2, size_t n) {
        vector<pair<size_t, size_t> > vec(n);
        for (size_t i = 0; i < n; i++) {
            vec[i] = make_pair(array[i], i);
        }
        std::stable_sort(vec.begin(), vec.end(),
                [](pair<size_t,size_t> x, pair<size_t,size_t> y)
                {
                    return (get<0>(x) < get<0>(y));
                });
        int_vector<> tmp(n);
        permute_values(array, tmp, vec, n);
        permute_values(array2, tmp, vec, n);
    }
};
}

#endif