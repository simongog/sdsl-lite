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
    public:
        typedef k2_tree_ns::idx_type idx_type;
        typedef k2_tree_ns::size_type size_type;

    private:
        //! Bit array to store all the bits of the tree, except those in the
        //! last level.
        t_bv        k_t;
        //! Bit array to store the last level of the tree.
        t_bv        k_l;

        t_rank      k_t_rank;

        uint8_t     k_k;
        uint16_t    k_height;

    protected:

        void build_from_matrix(const std::vector<std::vector <int>>& matrix)
        {
            // Makes the size a power of k.
            int simulated_size = std::pow(k, k_height);
            std::vector<std::deque<t_bv>> acc(k_height + 1);

            k2_tree_ns::_build_from_matrix<t_bv>(matrix, k, simulated_size,
                                                 k_height, 1, 0, 0, acc);

            size_type t_size = 0;
            size_type l_size = 0;
            for(int i = 1; i < k_height; i++)
                for(auto it = acc[i].begin(); it != acc[i].end(); it++)
                    t_size += (*it).size();

            for(auto it = acc[k_height].begin(); it != acc[k_height].end(); it++)
                l_size += (*it).size();

            k_t = t_bv(t_size, 0);
            k_l = t_bv(l_size, 0);

            int n = 0;
            for(int j = 1; j < k_height; j++)
                for(auto it = acc[j].begin(); it != acc[j].end(); it++)
                    // TODO erase
                    for(unsigned i = 0; i < (*it).size(); i++) {
                        k_t.set_int(n, (*it).get_int(i, 1), 1);
                        n++;
                    }
            n = 0;
            for(auto it = acc[k_height].begin(); it != acc[k_height].end(); it++)
                // TODO erase
                for(unsigned i = 0; i < (*it).size(); i++) {
                    k_l.set_int(n * 1, (*it).get_int(i, 1), 1);
                    n++;
                }
        }

        void _neigh(size_type n, idx_type row, idx_type col, size_type level,
                    std::vector<idx_type>& acc) const
        {
            // TODO we are only using n/k_k fix this, there is a better way.
            if(level >= k_t.size()) { // Last level
                if(k_l[level - k_t.size()] == 1)
                    acc.push_back(col);
                return;
            }

            if(k_t[level] == 1) {
                idx_type y = k_t_rank(level + 1) * std::pow(k_k, 2) +
                        k_k * std::floor(row/static_cast<double>(n/k_k));
                for(unsigned j = 0; j < k_k; j++)
                    _neigh(n/k_k, row % (n/k_k), col + (n/k_k) * j, y + j, acc);
            }
        }

        void _reverse_neigh(size_type n, idx_type row, idx_type col,
                            size_type level, std::vector<idx_type>& acc) const
        {
            if(level >= k_t.size()) { // Last level
                if(k_l[level - k_t.size()] == 1) {
                    acc.push_back(row);
                }
                return;
            }

            if(k_t[level] == 1) {
                idx_type y = k_t_rank(level + 1) * std::pow(k_k, 2) +
                        std::floor(col/static_cast<double>(n/k_k));
                for(unsigned j = 0; j < k_k; j++)
                    _reverse_neigh(n/k_k, row + (n/k_k) * j, col % (n/k_k),
                                   y + j * k_k, acc);
            }
        }

    public:

        k2_tree() = default;

        // TODO There is a build process proportional to the number of 1s in
        // the paper. Implement it!
        //! Constructor
        /*! This constructos takes the graph adjacency matrix.
         *  The time complexity for this constructor is linear in the matrix
         *  size
         *  \param matrix Adjacency matrix of the graph. It must be a binary
         *      square matrix.
         */
        k2_tree(std::vector<std::vector <int>>& matrix)
        {
            if(matrix.size() < 1) {
                throw std::logic_error("Matrix has no elements");
            }
            // TODO Assert matrix is an square matrix?
            std::vector<bit_vector> t;
            k_k = k;
            if(matrix.size() < k_k)
                k_height = 1;
            else // height = log_k n
                k_height = std::ceil(std::log(matrix.size())/std::log(k));

            build_from_matrix(matrix);

            k_t_rank = t_rank(&k_t);
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

        //! Indicates wheter node j is adjacent to node i or not.
        /*!
         *  \param i Node i.
         *  \param j Node j.
         *  \returns true if there is an edge going from node i to node j,
         *           false otherwise.
         */
        bool adj(idx_type i, idx_type j) const
        {
            size_type n = std::pow(k_k, k_height - 1);
            size_type k_2 = std::pow(k_k, 2);
            idx_type col, row;

            // This is duplicated to avoid an extra if at the loop. As idx_type
            // is unsigned and rank has an offset of one, is not possible to run
            // k_t_rank with zero as parameter at the first iteration.
            row = std::floor(i/static_cast<double>(n));
            col = std::floor(j/static_cast<double>(n));
            i = i % n;
            j = j % n;
            idx_type level = k_k * row + col;
            n = n/k_k;
            idx_type y;

            while(level < k_t.size()) {
                if(k_t[level] == 0)
                    return false;
                row = std::floor(i/static_cast<double>(n));
                col = std::floor(j/static_cast<double>(n));
                i = i % n;
                j = j % n;
                level = k_t_rank(level + 1) * k_2 + k_k * row + col;
                n = n/k_k;
            }

            return k_l[level - k_t.size()] == 1;
        }

        //! Returns a list of neighbors of node i.
        /*!
         *  \param i Node to get neighbors from.
         *  \returns A list of neighbors of node i.
         */
        std::vector<idx_type>neigh(idx_type i) const
        {
            std::vector<idx_type> acc{};
            size_type n = std::pow(k_k, k_height);
            idx_type y = k_k * std::floor(i/static_cast<double>(n/k_k));
            for(unsigned j = 0; j < k_k; j++)
                _neigh(n/k_k, i % (n/k_k), (n/k_k) * j, y + j, acc);
            return acc;
        }

        //! Returns a list of reverse neighbors of node i.
        /*!
         *  \param i Node to get reverse neighbors from.
         *  \returns A list of reverse neighbors of node i.
         */
        std::vector<idx_type> reverse_neigh(idx_type i) const
        {
            std::vector<idx_type> acc{};
            size_type n = std::pow(k_k, k_height);
            idx_type y = k_k * std::floor(i/static_cast<double>(n/k_k));
            for(unsigned j = 0; j < k_k; j++)
                _reverse_neigh(n/k_k, (n/k_k) * j, i % (n/k), y + j * k_k, acc);

            return acc;
        }

};
}

#endif
