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
/*! \file wt_topk_algorithm.hpp
    \brief wt_topk_algorithm.hpp contains a class for the wavelet tree plus rmq data structure.
            This data structure is based on the solution for solving top-k queries on grids by
            G.Navarro, Y. Nekrich and L. Russo, Space-Efficient Data-Analysis Queries on Grids.
            Theoretical Computer Science 482:60-72, 2013

    \author Simon Gog, Roberto Konow
*/

#ifndef INCLUDED_SDSL_WT_TOPK_ALGORITHM
#define INCLUDED_SDSL_WT_TOPK_ALGORITHM

namespace sdsl
{
// forward declaration
template<typename t_wt,
        typename t_rmq,
        typename t_weight_vec
>
class wt_topk;

//! Specialized version of method ,,construct'' for wt_topk.
template<typename t_wt,
        typename t_rmq,
        typename t_weight_vec
>
void
construct(wt_topk<t_wt, t_rmq, t_weight_vec>& idx, std::string file)
{
    int_vector_buffer<> buf_x(file+".x", std::ios::in);
    int_vector_buffer<> buf_y(file+".y", std::ios::in);
    int_vector_buffer<> buf_w(file+".w", std::ios::in);
    wt_topk<t_wt, t_rmq, t_weight_vec> tmp(buf_x, buf_y, buf_w);
    tmp.swap(idx);
}

//! Specialized version of method ,,construct_im'' for k2_treaps.
template<typename t_wt,
        typename t_rmq,
        typename t_weight_vec
>
void
construct_im(wt_topk<t_wt, t_rmq, t_weight_vec>& idx, std::vector<std::array<uint64_t, 3>> data)
{
    std::string tmp_prefix = ram_file_name("wt_topk_");
    std::vector<std::tuple<uint64_t,uint64_t,uint64_t>> d;
    for (auto x : data) {
        d.push_back(std::make_tuple(x[0],x[1],x[2]));
    }
    wt_topk<t_wt, t_rmq, t_weight_vec> tmp(d, tmp_prefix);
    tmp.swap(idx);
}
}

#endif
