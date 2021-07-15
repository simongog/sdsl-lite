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

#include <queue>
#include <stdexcept>
#include <tuple>
#include <cmath>
#include <memory>
#include "sdsl/bit_vectors.hpp"
#include "sdsl/k2_tree_helper.hpp"
#include "sdsl/k2_tree_iterator.hpp"
#include "sdsl/int_vector_buffer.hpp"


#include <ctime>
// #include <boost/circular_buffer.hpp>
// #include <boost/fusion/container/deque.hpp>
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

template <uint16_t k,
          typename t_bv = bit_vector,
          typename t_rank = typename t_bv::rank_1_type>
class k2_tree
{
public:
    typedef k2_tree_ns::idx_type idx_type;
    typedef k2_tree_ns::size_type size_type;
    using edg_iterator = edge_iterator<k2_tree<k, t_bv, t_rank>>;
    using nod_iterator = node_iterator<k2_tree<k, t_bv, t_rank>>;
    using neigh_iterator = neighbour_iterator<k2_tree<k, t_bv, t_rank>>;

    t_bv k_t;
    uint k_t_size = 0;
    uint k_l_size = 0;
    //! Bit array to store the last level of the tree.
    t_bv k_l;

    t_rank k_t_rank;

    uint16_t k_k;
    uint16_t k_height = 0;

protected:
    //! Bit array to store all the bits of the tree, except those in the
    //! last level.
    uint64_t n_marked_edges = 0;
    uint64_t n_edges = 0;
    size_t n_vertices = 0;

    std::vector<uint64_t> pointerL;
    std::vector<uint64_t> div_level_table;
    uint16_t max_level;
    uint last_level_rank = 0;

    edg_iterator it_e_begin, it_e_end;
    neigh_iterator it_neigh_begin, it_neigh_end = neigh_iterator().end();

    void build_from_matrix(const std::vector<std::vector<int>> &matrix)
    {
        // Makes the size a power of k.
        int simulated_size = std::pow(k, k_height);
        std::vector<std::deque<bit_vector>> acc(k_height + 1);

        k2_tree_ns::_build_from_matrix<bit_vector>(matrix, k,
                                                   simulated_size, k_height,
                                                   1, 0, 0, acc);

        size_type t_size = 0;
        size_type l_size = 0;
        for (int i = 1; i < k_height; i++)
            for (auto it = acc[i].begin(); it != acc[i].end(); it++)
                t_size += (*it).size();

        for (auto it = acc[k_height].begin(); it != acc[k_height].end(); it++)
            l_size += (*it).size();

        bit_vector k_t_(t_size, 0);
        bit_vector k_l_(l_size, 0);

        int n = 0;
        for (int j = 1; j < k_height; j++)
            for (auto it = acc[j].begin(); it != acc[j].end(); it++)
                for (unsigned i = 0; i < (*it).size(); i++)
                {
                    // TODO there should be a better way to do this
                    k_t_.set_int(n, (*it).get_int(i, 1), 1);
                    n++;
                }
        n = 0;
        for (auto it = acc[k_height].begin(); it != acc[k_height].end(); it++)
            for (unsigned i = 0; i < (*it).size(); i++)
            {
                // TODO there should be a better way to do this
                k_l_.set_int(n * 1, (*it).get_int(i, 1), 1);
                n++;
                if((*it).get_int(i, 1) == 1)
                    ++n_edges;
            }

        k2_tree_ns::build_template_vector<t_bv>(k_t_, k_l_, k_t, k_l);
        k_t_rank = t_rank(&k_t);
        
        k_t_size = k_t.size();
        k_l_size = k_l.size();
    }

    /*! Recursive function to retrieve list of neighbors.
         *
         *  \param n Size of the submatrix in the next recursive step.
         *  \param row Row of interest in the current submatrix, this is the
         *      row corresponding the node we are looking neighbors for.
         *  \param col Column offset of the current submatrix in the global
         *      matrix.
         *  \param level Position in k_t:k_l (k_l appended to k_t) of the node
         *      or leaf being processed at this step.
         *  \param acc Accumulator to store the neighbors found.
         */
    void _neigh(size_type n, idx_type row, idx_type col, size_type level,
                std::vector<idx_type> &acc) const
    {
        if (level >= k_t.size())
        { // Last level
            if (k_l[level - k_t.size()] == 1)
                acc.push_back(col);
            return;
        }

        if (k_t[level] == 1)
        {
            idx_type y = k_t_rank(level + 1) * std::pow(k_k, 2) +
                         k_k * std::floor(row / static_cast<double>(n));
            for (unsigned j = 0; j < k_k; j++)
                _neigh(n / k_k, row % n, col + n * j, y + j, acc);
        }
    }

    void _neigh_rec(size_type n, idx_type row, idx_type col, size_type level,
                std::function<void(uint64_t)> func) const {
        if (level >= k_t.size())
        { // Last level
            if (k_l[level - k_t.size()] == 1)
                func(col);
            return;
        }

        if (k_t[level] == 1)
        {
            idx_type y = k_t_rank(level + 1) * std::pow(k_k, 2) +
                         k_k * std::floor(row / static_cast<double>(n));
            for (unsigned j = 0; j < k_k; j++)
                _neigh_rec(n / k_k, row % n, col + n * j, y + j, func);
        }
    }

    void _edge_it_rec(uint64_t dp, uint64_t dq, int64_t x, int16_t l, std::function<void(uint64_t, uint64_t)> func) {

        if(l == max_level) {
            if(k_l[x] == 1) {
                func(dp, dq);
            }
        }

        if(((l == max_level-1) && (x != -1) && k_t[x] == 1)) {
            uint64_t y = pointerL[l+1];
            pointerL[l+1] += k_k*k_k;

            for(uint i = 0; i < k_k; i++) {
                for(uint j = 0; j < k_k; j++) {
                    _edge_it_rec(dp+i, dq+j, y+k_k*i+j, l+1, func);
                }
            }
        }

        if((x == -1) || ((l < max_level-1) && (k_t[x] == 1) )) {
            uint64_t y = pointerL[l+1];
            pointerL[l+1] += k_k*k_k;

            uint64_t div_level = div_level_table[l+1];
            for(uint i = 0; i < k_k; i++) {
                for(uint j = 0; j < k_k; j++) {
                    _edge_it_rec(dp+div_level*i, dq+div_level*j, y+k_k*i+j, l+1, func);
                }
            }
        }
    }

    /*! Recursive function to retrieve list of reverse neighbors.
         *
         *  \param n Size of the submatrix in the next recursive step.
         *  \param row Row offset of the current submatrix in the global matrix.
         *  \param col Column of interest in the current submatrix, this is the
         *      column corresponding the node we are looking reverse neighbors
         *      for.
         *  \param level Position in k_t:k_l (k_l appended to k_t) of the node
         *      or leaf being processed at this step.
         *  \param acc Accumulator to store the neighbors found.
         */
    void _reverse_neigh(size_type n, idx_type row, idx_type col,
                        size_type level, std::vector<idx_type> &acc) const
    {
        if (level >= k_t.size())
        { // Last level
            if (k_l[level - k_t.size()] == 1)
            {
                acc.push_back(row);
            }
            return;
        }

        if (k_t[level] == 1)
        {
            idx_type y = k_t_rank(level + 1) * std::pow(k_k, 2) +
                         std::floor(col / static_cast<double>(n));
            for (unsigned j = 0; j < k_k; j++)
                _reverse_neigh(n / k_k, row + n * j, col % n,
                               y + j * k_k, acc);
        }
    }

    //! Build a tree from an edges collection
    /*! This method takes a vector of edges describing the graph
         *  and the graph size. And takes linear time over the amount of
         *  edges to build the k_2 representation.
         *  \param edges A vector with all the edges of the graph, it can
         *               not be empty.
         *  \param size Size of the graph, all the nodes in edges must be
         *              within 0 and size ([0, size[).
         */
    void build_from_edges(std::vector<std::tuple<idx_type, idx_type>> &edges,
                          const size_type size)
    {

        typedef std::tuple<idx_type, idx_type, size_type, idx_type,
                           idx_type>
            t_part_tuple;

        k_k = k;
        k_height = std::ceil(std::log(size) / std::log(k_k));
        k_height = k_height > 1 ? k_height : 1; // If size == 0
        size_type k_2 = std::pow(k_k, 2);
        bit_vector k_t_ = bit_vector(k_2 * k_height * edges.size(), 0);
        bit_vector k_l_;

        std::queue<t_part_tuple> q;
        idx_type t = 0, last_level = 0;
        idx_type i, j, r_0, c_0, it, c, r;
        size_type l = std::pow(k_k, k_height - 1);
        std::vector<idx_type> pos_by_chunk(k_2 + 1, 0);

        q.push(t_part_tuple(0, edges.size(), l, 0, 0));

        while (!q.empty())
        {
            std::vector<idx_type> amount_by_chunk(k_2, 0);
            std::tie(i, j, l, r_0, c_0) = q.front();
            q.pop();
            // Get size for each chunk
            for (it = i; it < j; it++)
                amount_by_chunk[k2_tree_ns::get_chunk_idx(
                    std::get<0>(edges[it]), std::get<1>(edges[it]),
                    c_0, r_0, l, k_k)] += 1;
            if (l == 1)
            {
                if (last_level == 0)
                {
                    last_level = t;
                    k_l_ = bit_vector(k_t_.size() - last_level, 0);
                    k_t_.resize(last_level);
                    last_level = 1; // if t was 0
                    t = 0;          // Restart counter as we're storing at k_l_ now.
                }
                for (it = 0; it < k_2; it++, t++)
                    if (amount_by_chunk[it] != 0) {
                        k_l_[t] = 1;
                        ++n_edges;
                    }
                // At l == 1 we do not put new elements at the queue.
                continue;
            }

            // Set starting position in the vector for each chunk
            pos_by_chunk[0] = i;
            for (it = 1; it < k_2; it++)
                pos_by_chunk[it] =
                    pos_by_chunk[it - 1] + amount_by_chunk[it - 1];
            // To handle the last case when it = k_2 - 1
            pos_by_chunk[k_2] = j;
            // Push to the queue every non zero elements chunk
            for (it = 0; it < k_2; it++, t++)
                // If not empty chunk, set bit to 1
                if (amount_by_chunk[it] != 0)
                {
                    r = it / k_k;
                    c = it % k_k;
                    k_t_[t] = 1;
                    q.push(t_part_tuple(pos_by_chunk[it],
                                        pos_by_chunk[it + 1],
                                        l / k_k,
                                        r_0 + r * l,
                                        c_0 + c * l));
                }
            idx_type chunk;

            // Sort edges' vector
            for (unsigned ch = 0; ch < k_2; ch++)
            {
                idx_type be = ch == 0 ? i : pos_by_chunk[ch - 1];
                for (it = pos_by_chunk[ch]; it < be + amount_by_chunk[ch];)
                {
                    chunk = k2_tree_ns::get_chunk_idx(
                        std::get<0>(edges[it]), std::get<1>(edges[it]),
                        c_0, r_0, l, k_k);

                    if (pos_by_chunk[chunk] != it)
                        std::iter_swap(edges.begin() + it,
                                       edges.begin() + pos_by_chunk[chunk]);
                    else
                        it++;
                    pos_by_chunk[chunk]++;
                }
            }
        }
        k_l_.resize(t);
        k2_tree_ns::build_template_vector<t_bv>(k_t_, k_l_, k_t, k_l);

        k_t_rank = t_rank(&k_t);
        k_t_size = k_t.size();
        k_l_size = k_l.size();
    }

public:
    k2_tree()
    {
        k_k = k;
        k_t = bit_vector(0, 0);
        k_l = bit_vector(0, 0);
        k_t_rank = t_rank(&k_t);
    }

    k2_tree(uint n_vertices) : n_vertices(n_vertices)
    {
        k_k = k;
        k_t = bit_vector(0, 0);
        k_l = bit_vector(0, 0);
        k_t_rank = t_rank(&k_t);
    }

    //! Constructor
    /*! This constructos takes the graph adjacency matrix.
         *  The time complexity for this constructor is linear in the matrix
         *  size
         *  \param matrix Adjacency matrix of the graph. It must be a binary
         *      square matrix.
         */
    k2_tree(const std::vector<std::vector<int>> &matrix)
    {
        if (matrix.size() < 1)
        {
            throw std::logic_error("Matrix has no elements");
        }
        std::vector<bit_vector> t;
        k_k = k;
        if (matrix.size() < k_k)
            k_height = 1;
        else // height = log_k n
            k_height = std::ceil(std::log(matrix.size()) / std::log(k_k));

        build_from_matrix(matrix);
        this->n_vertices = matrix.size(); 

        max_level = floor(log(n_vertices)/log(k_k));
        if(max_level != 0 && floor(log(n_vertices)/log(k_k)) == (log(n_vertices)/log(k_k)))
            max_level = max_level-1;

        div_level_table = std::vector<uint64_t>(max_level+1);
        for(int64_t i = 0; i <= max_level; i++)
            div_level_table[i] = exp_pow(k_k, max_level-i);
    }

    //! Constructor
    /*! This constructos takes a vector of edges describing the graph
         *  and the graph size. And takes linear time over the amount of
         *  edges to build the k_2 representation.
         *  \param edges A vector with all the edges of the graph, it can
         *               not be empty.
         *  \param size Size of the graph, all the nodes in edges must be
         *              within 0 and size ([0, size[).
         */
    k2_tree(std::vector<std::tuple<idx_type, idx_type>> &edges,
            const size_type size)
    { 
        assert(size > 0);
        assert(edges.size() > 0);

        build_from_edges(edges, size);
        this->n_vertices = size;

        max_level = floor(log(n_vertices)/log(k_k));
        if(max_level != 0 && floor(log(n_vertices)/log(k_k)) == (log(n_vertices)/log(k_k)))
            max_level = max_level-1;

        div_level_table = std::vector<uint64_t>(max_level+1);
        for(int64_t i = 0; i <= max_level; i++)
            div_level_table[i] = exp_pow(k_k, max_level-i);
    }

    //! Constructor
    /*! This constructos expects a filename prefix. Two serialized
         *  int_vectors have to be present at filename.x and filename.y.
         *  Each pair x,y describes an edge of the graph, from the node x
         *  to the node y.
         *  \param filename String with the prefix of the files filename.x,
         *					filename.y each of them containing a serialized
         *					int_vector<>.
         *  \param size Size of the graph, all the nodes in the edges defined
         *				by the files must be within 0 and size ([0, size[). If
         *				size==0, the size will be taken as the max node
         *				in the edges.
         */
    k2_tree(std::string filename, size_type size = 0)
    {
        int_vector_buffer<> buf_x(filename + ".x", std::ios::in);
        int_vector_buffer<> buf_y(filename + ".y", std::ios::in);

        assert(buf_x.size() == buf_y.size());
        assert(buf_x.size() > 0);

        std::vector<std::tuple<idx_type, idx_type>> edges;
        edges.reserve(buf_x.size());

        if (size == 0)
        {
            size_type max = 0;
            for (auto v : buf_x)
                max = std::max(static_cast<size_type>(v), max);
            for (auto v : buf_y)
                max = std::max(static_cast<size_type>(v), max);
            size = max + 1;
        }

        for (uint64_t i = 0; i < buf_x.size(); i++)
            edges.push_back(
                std::tuple<idx_type, idx_type>{buf_x[i], buf_y[i]});

        build_from_edges(edges, size);
        n_vertices = size;

        max_level = floor(log(n_vertices)/log(k_k));
        if(max_level != 0 && floor(log(n_vertices)/log(k_k)) == (log(n_vertices)/log(k_k)))
            max_level = max_level-1;

        div_level_table = std::vector<uint64_t>(max_level+1);
        for(int64_t i = 0; i <= max_level; i++)
            div_level_table[i] = exp_pow(k_k, max_level-i);
    }

    k2_tree(k2_tree &tr)
    {
        *this = tr;
    }

    k2_tree(k2_tree &&tr)
    {
        *this = std::move(tr);
    }

    size_t t_size() const {
        return k_t_size;
    }

    size_t l_size() const {
        return k_l_size;
    }

    size_t get_int_t(size_t i) {
        assert(i < k_t_size);
        return k_t[i];
    }

    size_t get_int_l(size_t i) {
        assert(i < k_l_size);
        return k_l[i];
    }

    bool t_empty() {
        return k_t.empty();
    }

    bool l_empty() {
        return k_l.empty();
    }

    uint64_t get_marked_edges() const
    {
        return n_marked_edges;
    }

    uint64_t get_number_edges() const
    {
        return n_edges - n_marked_edges;
    }

    size_type get_number_nodes() const
    {
        return n_vertices;
    }

    typedef struct union_node {
        uint16_t level;
        uint8_t rA, rB;
    } union_node;

    //! Union Operation
    /*! Performs the union operation between two tree. This operations requires both 
         * trees to have the same number of nodes.
         *  \param k2_B a k2_tree with the same number of nodes
         *  \par References
         *      [2] Brisaboa, Nieves R., et al. "Efficient Set Operations over 
         *      k2-Trees." 2015 Data Compression Conference. IEEE, 2015.
         */
    void unionOp(shared_ptr<k2_tree> k2_B)
    {
        if(k2_B == nullptr)
            return;
        
        if(k2_B->get_number_edges() == 0)
            return;

        if (get_number_edges() == 0) {
            *this = *k2_B;
            return;
        }
        assert(k_k == k2_B->k_k);

        if (n_vertices != k2_B->n_vertices)
            throw std::logic_error("Trees must have the same number of nodes.");
        if (k_height != k2_B->k_height)
            throw std::logic_error("Trees must have the same height.");

        const uint16_t max_height = k_height;

        const uint64_t t_size_A = k_t.size();
        const uint64_t t_size_B = k2_B->k_t.size();

        const uint64_t l_size_A = k_l.size();
        const uint64_t l_size_B = k2_B->k_l.size();

        // C Initialization
        uint64_t calc = 1;
        uint64_t max_bits = 0;
        for (uint64_t i = 0; i < max_height; i++)
        {
            calc *= k_k*k_k;
            max_bits += calc;
        }

        uint64_t C_t_size = max_bits < t_size_A + t_size_B ? max_bits : t_size_A + t_size_B;
        bit_vector C_t(C_t_size);

        calc *= k_k*k_k;
        max_bits += calc;

        uint64_t C_l_size = max_bits < l_size_A + l_size_B ? max_bits : l_size_A + l_size_B;
        bit_vector C_l(C_l_size);
        ////////

        // Q Initialization
        std::queue<union_node> Q;
        Q.push({0, 1, 1});
        ////////

        union_node next;
        uint8_t bA, bB;
        uint64_t pA, pB, idx_t, idx_l;
        pA = 0;
        pB = 0;
        idx_l = 0;
        idx_t = 0;
        
        int n_total_edges = 0;
        while (!Q.empty()) {
            next = Q.front();
            Q.pop();
            ++next.level;
            for (auto i = 0; i < k_k * k_k; ++i) {
                bA = 0;
                bB = 0;
                if (next.rA == 1) {
                    if (next.level < max_height)
                        bA = k_t[pA];
                    else
                        bA = k_l[pA - t_size_A];
                    ++pA;
                }
                if (next.rB == 1) {
                    if (next.level < max_height)
                        bB = k2_B->k_t[pB];
                    else
                        bB = k2_B->k_l[pB - t_size_B];
                    ++pB;
                }

                if(bA || bB) {
                    if(next.level < max_height) {
                        Q.push({next.level, bA, bB});
                        C_t[idx_t] = 1;
                    } else {
                        C_l[idx_l] = 1;
                        ++n_total_edges;
                    }
                }
                next.level < max_height ? ++idx_t : ++idx_l;
            }
        }
        assert(C_t_size >= idx_t);
        C_t.resize(idx_t);

        assert(C_l_size >= idx_l);
        C_l.resize(idx_l);
        
        k_t.swap(C_t);
        k_l.swap(C_l);
        k_t_rank = t_rank(&k_t);
        k_height = max_height;
        n_marked_edges = 0;
        n_edges = n_total_edges;
        k_t_size = idx_t;
        k_l_size = idx_l;
    }

    //! Move assignment operator
    k2_tree &operator=(k2_tree &&tr)
    {
        if (this != &tr)
        {
            k_t = std::move(tr.k_t);
            k_l = std::move(tr.k_l);
            k_k = std::move(tr.k_k);
            k_height = std::move(tr.k_height);
            k_t_rank = t_rank(&k_t);
            n_vertices = std::move(tr.n_vertices);
            n_marked_edges = std::move(tr.n_marked_edges);
            n_edges = std::move(tr.n_edges);
            max_level = std::move(tr.max_level);
            div_level_table = std::move(tr.div_level_table);
            last_level_rank = std::move(tr.last_level_rank);
            k_t_size = std::move(tr.k_t_size);
            k_l_size = std::move(tr.k_l_size);
        }
        return *this;
    }

    //! Assignment operator
    k2_tree &operator=(k2_tree &tr)
    {
        if (this != &tr)
        {
            k_t = tr.k_t;
            k_l = tr.k_l;
            k_k = tr.k_k;
            k_height = tr.k_height;
            k_t_rank = t_rank(&k_t);
            n_vertices = tr.n_vertices;
            n_marked_edges = tr.n_marked_edges;
            n_edges = tr.n_edges;
            max_level = tr.max_level;
            div_level_table = tr.div_level_table;
            last_level_rank = tr.last_level_rank;
            k_t_size = tr.k_t_size;
            k_l_size = tr.k_l_size;
        }
        return *this;
    }

    bool equal(const k2_tree &tr) const {
        bool result = true;

        result &= n_vertices == tr.n_vertices;
        result &= k_k == tr.k_k;
        result &= k_height == tr.k_height;
        
        if(k_t.size() != tr.k_t.size())
            return false;
        else {
            for(unsigned int i = 0; i < k_t.size(); ++i)
                result &= k_t[i] == tr.k_t[i];
        }

        if(k_l.size() != tr.k_l.size())
            return false;
        else {
            for(unsigned int i = 0; i < k_l.size(); ++i)
                result &= k_l[i] == tr.k_l[i];
        }
        return result;
    }

    //! Swap operator
    void swap(k2_tree &tr)
    {
        if (this != &tr)
        {
            std::swap(k_t, tr.k_t);
            std::swap(k_l, tr.k_l);
            util::swap_support(k_t_rank, tr.k_t_rank, &k_t, &(tr.k_t));
            std::swap(k_k, tr.k_k);
            std::swap(k_height, tr.k_height);
            std::swap(n_vertices, tr.n_vertices);
            std::swap(n_marked_edges, tr.n_marked_edges);
            std::swap(n_edges, tr.n_edges);
            std::swap(max_level, tr.max_level);
            std::swap(div_level_table, tr.div_level_table);
            std::swap(last_level_rank, tr.last_level_rank);
            std::swap(k_t_size, tr.k_t_size);
            std::swap(k_l_size, tr.k_l_size);

        }
    }

    //! Equal operator
    bool operator==(const k2_tree &tr) const
    {
        // TODO check the rank support equality?
        if (k_k != tr.k_k || k_height != tr.k_height)
            return false;
        if (k_t.size() != tr.k_t.size() || k_l.size() != tr.k_l.size())
            return false;
        for (unsigned i = 0; i < k_t.size(); i++)
            if (k_t[i] != tr.k_t[i])
                return false;
        for (unsigned i = 0; i < k_l.size(); i++)
            if (k_l[i] != tr.k_l[i])
                return false;
        return true;
    }

    bool operator!=(const k2_tree &tr) const
    {
        return !(*this == tr);
    }


    //! Indicates whether node j is adjacent to node i or not.
    /*!
         *  \param i Node i.
         *  \param j Node j.
         *  \returns true if there is an edge going from node i to node j,
         *           false otherwise.
         */
    bool adj(idx_type i, idx_type j) const
    {
        uint t_size = k_t.size();
        if (t_size == 0 && k_l.size() == 0)
            return false;
        size_type n = std::pow(k_k, k_height - 1);
        size_type k_2 = std::pow(k_k, 2);
        idx_type col, row;

        // This is duplicated to avoid an extra if at the loop. As idx_type
        // is unsigned and rank has an offset of one, is not possible to run
        // k_t_rank with zero as parameter at the first iteration.
        row = std::floor(i / static_cast<double>(n));
        col = std::floor(j / static_cast<double>(n));
        i = i % n;
        j = j % n;
        idx_type level = k_k * row + col;
        n = n / k_k;

        while (level < t_size)
        {
            if (k_t[level] == 0)
                return false;
            row = std::floor(i / static_cast<double>(n));
            col = std::floor(j / static_cast<double>(n));
            i = i % n;
            j = j % n;
            level = k_t_rank(level + 1) * k_2 + k_k * row + col;
            n = n / k_k;
        }
        return k_l[level - t_size] == 1;
    }


    //! Indicates whether node j is adjacent to node i or not, and makes
    // it adjacent if possible.
    /*!
         *  \param i Node i.
         *  \param j Node j.
         *  \returns true if there is an edge going from node i to node j,
         *           false otherwise.
         */
    bool adj_forced(idx_type i, idx_type j)
    {
        uint t_size = k_t.size();
        if (t_size == 0 && k_l.size() == 0)
            return false;
        size_type n = std::pow(k_k, k_height - 1);
        size_type k_2 = std::pow(k_k, 2);
        idx_type col, row;

        // This is duplicated to avoid an extra if at the loop. As idx_type
        // is unsigned and rank has an offset of one, is not possible to run
        // k_t_rank with zero as parameter at the first iteration.
        row = std::floor(i / static_cast<double>(n));
        col = std::floor(j / static_cast<double>(n));
        i = i % n;
        j = j % n;
        idx_type level = k_k * row + col;
        n = n / k_k;

        while (level < t_size)
        {
            if (k_t[level] == 0)
                return false;
            row = std::floor(i / static_cast<double>(n));
            col = std::floor(j / static_cast<double>(n));
            i = i % n;
            j = j % n;
            level = k_t_rank(level + 1) * k_2 + k_k * row + col;
            n = n / k_k;
        }
        k_l[level - t_size] = 1;
        return true;
    }


    bool contains(uint p, uint q)
    {
        if(k_l.size() > 0) {
            max_level = floor(log(n_vertices)/log(k_k));
            if(max_level != 0 && floor(log(n_vertices)/log(k_k)) == (log(n_vertices)/log(k_k)))
                max_level = max_level-1;

            div_level_table = std::vector<uint64_t>(max_level+1);
            for(int64_t i = 0; i <= max_level; i++)
                div_level_table[i] = exp_pow(k_k, max_level-i);

            std::vector<uint64_t> pL (max_level+1);
            pL[0] = 0;
            pL[1] = k_k*k_k;

            for(int16_t i = 2; i < max_level; i++) {
                pL[i] = (k_t_rank(pL[i-1])+1)*k_k*k_k;
            }
            pL[max_level] = 0;

            if(max_level > 0)
                last_level_rank = pL[max_level-1] == 0? 0 : k_t_rank(pL[max_level-1]);
            else
                last_level_rank = 0;

            return adj_rec(p, q, 0, 0);
        }
        return false;
    }

    bool adj_rec(uint p, uint q, uint64_t node, int level)
    {
        int div_level = div_level_table[level];
        uint64_t newnode = p / div_level * k_k + q / div_level;
        newnode += node;

        if(level == max_level && k_l[newnode]) {
            return true;
        }
        else if (level < max_level - 1 && k_t[newnode])
        {
            return adj_rec(p % div_level, q % div_level, k_t_rank(newnode+1) * k_k * k_k, level + 1);
        }
        else if(level == max_level - 1 && k_t[newnode])
        {
            uint64_t posInf = (k_t_rank(newnode) - last_level_rank) * k_k * k_k;
            uint64_t shift = (q % k_k + (p % k_k) * k_k);
            return k_l[posInf + shift];
        }
        return false;
    }

    //! Returns a list of neighbors of node i.
    /*!
         *  \param i Node to get neighbors from.
         *  \returns A list of neighbors of node i.
         */
    std::vector<idx_type> neigh(idx_type i) const
    {  
        std::vector<idx_type> acc{};
        if (k_l.size() == 0 && k_t.size() == 0)
            return acc;
        // n = k^h / k
        // k^h - dimension n of matrix nxn
        // /k  - to calculate div only once and not for for all parameter again, always (n/k)
        size_type n =
            static_cast<size_type>(std::pow(k_k, k_height)) / k_k;
        // y = k * i/n
        idx_type y = k_k * std::floor(i / static_cast<double>(n));
        for (unsigned j = 0; j < k_k; j++) {
            _neigh(n / k_k, i % n, n * j, y + j, acc);
        }
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
        if (k_l.size() == 0 && k_t.size() == 0)
            return acc;
        // Size of the first square division
        size_type n =
            static_cast<size_type>(std::pow(k_k, k_height)) / k_k;
        idx_type y = std::floor(i / static_cast<double>(n));
        for (unsigned j = 0; j < k_k; j++)
            _reverse_neigh(n / k_k, n * j, i % n, y + j * k_k, acc);

        return acc;
    }

    std::vector<std::pair<idx_type, idx_type>> range(
        idx_type row1, idx_type row2,
        idx_type col1, idx_type col2) const
    {
        std::vector<std::pair<idx_type, idx_type>> res;

        size_type n = static_cast<size_type>(std::pow(k_k, k_height)) / k_k;
        struct state
        {
            idx_type n, row1, row2, col1, col2, dr, dc, z;
            state(idx_type n, idx_type row1, idx_type row2, idx_type col1, idx_type col2,
                  idx_type dr, idx_type dc, idx_type z)
                : n(n), row1(row1), row2(row2), col1(col1), col2(col2), dr(dr), dc(dc), z(z) {}
        };
        std::vector<state> states;
        states.reserve(k_height); // minimum
        states.emplace_back(n, row1, row2, col1, col2, 0, 0, std::numeric_limits<idx_type>::max());

        while (!states.empty())
        {
            state s = states.back();
            states.pop_back();

            //TODO: peel first loop where z==-1 atm
            if (s.z != std::numeric_limits<idx_type>::max() && s.z >= k_t.size())
            { // Last level
                if (k_l[s.z - k_t.size()] == 1)
                {
                    res.emplace_back(s.dr, s.dc);
                }
            }
            else if (s.z == std::numeric_limits<idx_type>::max() || k_t[s.z] == 1)
            {

                auto y = k_t_rank(s.z + 1) * k_k * k_k;

                for (idx_type i = s.row1 / s.n; i <= s.row2 / s.n; ++i)
                {
                    idx_type row1new, row2new;
                    //TODO: loop peeling, first iteration and last iteration special
                    if (i == s.row1 / s.n)
                        row1new = s.row1 % s.n;
                    else
                        row1new = 0;
                    if (i == s.row2 / s.n)
                        row2new = s.row2 % s.n;
                    else
                        row2new = s.n - 1;

                    for (idx_type j = s.col1 / s.n; j <= s.col2 / s.n; ++j)
                    {
                        idx_type col1new, col2new;
                        //TODO: loop peeling, first iteration and last iteration special
                        if (j == s.col1 / s.n)
                            col1new = s.col1 % s.n;
                        else
                            col1new = 0;
                        if (j == s.col2 / s.n)
                            col2new = s.col2 % s.n;
                        else
                            col2new = s.n - 1;

                        states.emplace_back(s.n / k_k, row1new, row2new, col1new, col2new,
                                            s.dr + s.n * i, s.dc + s.n * j, y + k_k * i + j);
                    }
                }
            }
        }

        return res;
    };

    //! Serialize to a stream
    /*! Serialize the k2_tree data structure
         *  \param out Outstream to write the k2_tree.
         *  \param v
         *  \param string_name
         *  \returns The number of written bytes.
         */
    size_type serialize(std::ostream &out, structure_tree_node *v = nullptr,
                        std::string name = "") const
    {
        structure_tree_node *child = structure_tree::add_child(
            v, name, util::class_name(*this));
        size_type written_bytes = 0;

        written_bytes += k_t.serialize(out, child, "t");
        written_bytes += k_l.serialize(out, child, "l");
        written_bytes += k_t_rank.serialize(out, child, "t_rank");
        written_bytes += write_member(k_k, out, child, "k");
        written_bytes += write_member(k_height, out, child, "height");
        written_bytes += write_member(n_vertices, out, child, "n_vertices");
        written_bytes += write_member(n_marked_edges, out, child, "n_marked_edges");
        written_bytes += write_member(n_edges, out, child, "n_edges");

        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    //! Load from istream
    /*! Serialize the k2_tree from the given istream.
         *  \param istream Stream to load the k2_tree from.
         */
    void load(std::istream &in)
    {
        k_t.load(in);
        k_l.load(in);
        k_t_rank.load(in);
        k_t_rank.set_vector(&k_t);
        read_member(k_k, in);
        read_member(k_height, in);
        read_member(n_vertices, in);
        read_member(n_marked_edges, in);
        read_member(n_edges, in);

        max_level = floor(log(n_vertices)/log(k_k));
        if(max_level != 0 && floor(log(n_vertices)/log(k_k)) == (log(n_vertices)/log(k_k)))
            max_level = max_level-1;

        div_level_table = std::vector<uint64_t>(max_level+1);
        for(int64_t i = 0; i <= max_level; i++)
            div_level_table[i] = exp_pow(k_k, max_level-i);
        k_t_size = k_t.size();
        k_l_size = k_l.size();

    }

    bool erase(idx_type i, idx_type j)
    {
        if (k_t.size() == 0 && k_l.size() == 0)
            return false;
        size_type n = std::pow(k_k, k_height - 1);
        size_type k_2 = std::pow(k_k, 2);
        idx_type col, row;

        // This is duplicated to avoid an extra if at the loop. As idx_type
        // is unsigned and rank has an offset of one, is not possible to run
        // k_t_rank with zero as parameter at the first iteration.
        row = std::floor(i / static_cast<double>(n));
        col = std::floor(j / static_cast<double>(n));
        i = i % n;
        j = j % n;
        idx_type level = k_k * row + col;
        n = n / k_k;

        while (level < k_t.size())
        {
            if (k_t[level] == 0)
                return false;
            row = std::floor(i / static_cast<double>(n));
            col = std::floor(j / static_cast<double>(n));
            i = i % n;
            j = j % n;
            level = k_t_rank(level + 1) * k_2 + k_k * row + col;
            n = n / k_k;
        }

        if(k_l[level - k_t.size()] == 1)
        {
            k_l[level - k_t.size()] = 0;
            n_marked_edges++;
            return true;
        }
        return false;
    }

    edg_iterator &edge_begin()
    {
        it_e_begin = edg_iterator(this);
        it_e_end = it_e_begin.end();
        return it_e_begin;
    }

    edg_iterator &edge_end()
    {
        return it_e_end;
    }

    nod_iterator node_begin()
    {
        return nod_iterator(this);
    }

    nod_iterator node_end()
    {
        return nod_iterator(this).end();
    }

    neigh_iterator &neighbour_begin(idx_type node) {
        it_neigh_begin = neigh_iterator(this, node);
        return it_neigh_begin;
    }

    neigh_iterator &neighbour_end() {
        return it_neigh_end;
    }


    uint exp_pow(uint base, uint pow)
    {
        uint i, result = 1;
        for (i = 0; i < pow; i++)
            result *= base;

        return result;
    }

    void edge_it(std::function<void(uint64_t, uint64_t)> func) {
        if(k_l.size() > 0) {
            max_level = floor(log(n_vertices)/log(k_k));
            if(max_level != 0 &&  floor(log(n_vertices)/log(k_k)) == (log(n_vertices)/log(k_k)))
                max_level = max_level-1;

            div_level_table = std::vector<uint64_t>(max_level+1);
            for(int64_t i = 0; i <= max_level; i++)
                div_level_table[i] = exp_pow(k_k, max_level-i);

            pointerL = std::vector<uint64_t>(max_level+1);
            pointerL[0] = 0;
            pointerL[1] = k_k*k_k;

            for(size_t i = 2; i < max_level; i++) {
                pointerL[i] = (k_t_rank(pointerL[i-1])+1)*k_k*k_k;
            }
            pointerL[max_level] = 0;
            _edge_it_rec(0, 0, -1, -1, func);
        }
    }

    void set_edges(uint64_t edges) { n_edges = edges;}
    uint64_t total_edges() { return n_edges;}

    void neigh_it(uint64_t i, std::function<void(uint64_t)> func) {
        if (k_l.size() == 0 && k_t.size() == 0)
            return;
        // n = k^h / k
        // k^h - dimension n of matrix nxn
        // /k  - to calculate div only once and not for for all parameter again, always (n/k)
        size_type n =
            static_cast<size_type>(std::pow(k_k, k_height)) / k_k;
        // y = k * i/n
        idx_type y = k_k * std::floor(i / static_cast<double>(n));
        for (unsigned j = 0; j < k_k; j++) {
            _neigh_rec(n / k_k, i % n, n * j, y + j, func);
        }
    }

    std::vector<idx_type> neigh_test(idx_type i) const
        {
            std::deque<tree_node2> q;
            std::vector<idx_type> acc{};

            if (k_l.size() == 0 && k_t.size() == 0)
                return acc;

            size_type n =
                static_cast<size_type>(std::pow(k_k, k_height)) / k_k;
            // y = k * i/n
            idx_type y = k_k * std::floor(i / static_cast<double>(n));

            for (unsigned j = 0; j < k_k; j++) {
                q.emplace_front(tree_node2(n / k_k, i % n, n * j, y + j));
            }

            ///
            while (!q.empty())
            {
                tree_node2 s = q.front();
                q.pop_front();

                if (s.level >= k_t.size())
                {
                    if (k_l[s.level - k_t.size()] == 1)
                        acc.push_back(s.col);
                    continue; //
                }

                if (k_t[s.level] == 1)
                {
                    idx_type y = k_t_rank(s.level + 1) * std::pow(k_k, 2) +
                                    k_k * std::floor(s.row / static_cast<double>(s.n));
                    for (unsigned j = 0; j < k_k; j++)
                        q.emplace_front(tree_node2(s.n / k_k, s.row % s.n, s.col + s.n * j, y + j));
                }

            }
            return acc;
        }

        class tree_node2
        {
        public:
            size_type n;
            idx_type row, col, level;

            tree_node2() {}
            tree_node2(size_type n, idx_type row, idx_type col, size_type level) : n(n), row(row), col(col), level(level){}
        };
};
} // namespace sdsl

#endif
