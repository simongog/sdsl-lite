/* sdsl - succinct data structures library
    Copyright (C) 2016 Francisco Montoto

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
    \author Francisco Montoto
*/
#ifndef INCLUDED_SDSL_K2_TREE
#define INCLUDED_SDSL_K2_TREE

#include <sdsl/bit_vectors.hpp>
#include "sdsl/k2_tree_helper.hpp"


//! Namespace for the succint data structure library
namespace sdsl
{
//! A k^2-tree
/*! A k^2-tree is a compact tree structure to represent a web graph. The
 *  structure takes advantage of large empty areas of the adjacency matrix of
 *  the graph.
 *
 *  \par References
 *      [1] Brisaboa, N. R., Ladra, S., & Navarro, G. (2009, August):
 *          k2-trees for compact web graph representation. In International
 *          Symposium on String Processing and Information Retrieval
 *          (pp. 18-30). Springer Berlin Heidelberg.
 */

template<uint8_t k,
         typename t_bv=bit_vector,
         typename t_rank=typename t_bv::rank_1_type>
class k2_tree
{
    private:
        //! Bit array to store all the bits of the tree, except those in the
        //! last level.
        t_bv    k_t;
        //! Bit array to store the last level of the tree.
        t_bv    k_l;

        t_rank  k_t_rank;

    public:
        typedef k2_tree_ns::idx_type idx_type;
        typedef k2_tree_ns::size_type size_type;

        k2_tree() = default;

        // TODO There is a build process proportional to the number of 1s in
        // the paper. Implement it!
        // TODO just for testing purposes. Do not keep this constructor.
        //! Constructor
        k2_tree(std::vector<std::vector <int>>& matrix)
        {
            std::vector<bit_vector> t;
            build_from_matrix(matrix);

            k_t_rank = t_rank(&k_t);
        }

        void build_from_matrix(const std::vector<std::vector <int>>& matrix)
        {
            // height = log_k n
            int height = std::ceil(
                    std::log(matrix.size())/static_cast<double>(std::log(k)));
            // printf("k:: %d\n", k);
            // printf("n: %d\n", matrix.size());
            // printf("Height: %d\n", height);
            // Makes the size a power of k.
            int simulated_size = std::pow(k, height);
            // printf("simulated_size: %d\n", simulated_size);
            std::vector<std::deque<t_bv>> acc(height + 1);
            // Check vectors have the same size
            k2_tree_ns::_build_from_matrix<t_bv>(matrix, k, simulated_size,
                                                 height, 1, 0, 0, acc);

            size_type t_size = 0;
            size_type l_size = 0;
            for(int i = 1; i < height; i++)
                for(auto it = acc[i].begin(); it != acc[i].end(); it++)
                    t_size += (*it).size();

            for(auto it = acc[height].begin(); it != acc[height].end(); it++)
                l_size += (*it).size();

            k_t = t_bv(t_size, 0);
            k_l = t_bv(l_size, 0);

            int n = 0;
            for(int j = 1; j < height; j++)
                for(auto it = acc[j].begin(); it != acc[j].end(); it++)
                    // TODO erase
                    for(unsigned i = 0; i < (*it).size(); i++) {
                        k_t.set_int(n, (*it).get_int(i, 1), 1);
                        n++;
                    }
            n = 0;
            for(auto it = acc[height].begin(); it != acc[height].end(); it++)
                // TODO erase
                for(unsigned i = 0; i < (*it).size(); i++) {
                    k_l.set_int(n * 1, (*it).get_int(i, 1), 1);
                    n++;
                }
        }


        k2_tree(const k2_tree& tr)
        {
            *this = tr;
        }

        k2_tree(k2_tree&& tr)
        {
            *this = std::move(tr);
        }

        //! Move assignment operator
        k2_tree& operator=(k2_tree&& tr)
        {
            if(this != &tr) {
                k_t = std::move(tr.k_t);
                k_l = std::move(tr.k_l);
                k_t_rank = std::move(tr.k_t_rank);
                k_t_rank.set_vector(&k_t);
            }
            return *this;
        }

        //! Assignment operator
        k2_tree& operator=(k2_tree& tr)
        {
            if(this != &tr) {
                k_t = tr.k_t;
                k_l = tr.k_l;
                k_t_rank = tr.k_t_rank;
                k_t_rank.set_vector(&k_t);
            }
            return *this;
        }

        //! Swap operator
        void swap(k2_tree& tr)
        {
            if(this != &tr) {
                std::swap(k_t, tr.k_t);
                util::swap_support(k_t_rank, tr.k_t_rank, &k_t, &(tr.k_t));
            }
        }

        t_bv get_t()
        {
            return k_t;
        }

        t_bv get_l()
        {
            return k_l;
        }

};
}

#endif
