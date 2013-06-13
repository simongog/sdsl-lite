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
/*! \file wt_pc.hpp
    \brief wt_pc.hpp contains a class for the wavelet tree of byte sequences.
           The wavelet tree shape is parametrized by a prefix code.
    \author Simon Gog, Timo Beller
*/
#ifndef INCLUDED_SDSL_WT_PC
#define INCLUDED_SDSL_WT_PC

#include "bit_vectors.hpp"
#include "rank_support.hpp"
#include "select_support.hpp"
#include "wt_helper.hpp"
#include <vector>
#include <deque>
#include <queue>
#include <utility>

//! Namespace for the succinct data structure library.
namespace sdsl
{

const int_vector<>::size_type ZoO[2] = {0, (int_vector<>::size_type)-1};

template<class t_shape,
         class t_bitvector   = bit_vector,
         class t_rank        = typename t_bitvector::rank_1_type,
         class t_select      = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type,
         bool  t_dfs_shape   = false>
class wt_pc
{
    public:

        typedef int_vector<>::size_type          size_type;
        typedef unsigned char                    value_type;
        typedef t_bitvector                      bit_vector_type;
        typedef t_rank                           rank_1_type;
        typedef t_select                         select_1_type;
        typedef t_select_zero                    select_0_type;
        typedef wt_tag                           index_category;
        typedef byte_alphabet_tag                alphabet_category;
        typedef typename
        t_shape::template type<wt_pc>::t shape;
        enum { lex_ordered=shape::lex_ordered };

    private:
#ifdef WT_HUFF_CACHE
        mutable value_type m_last_access_answer;
        mutable size_type  m_last_access_i;
        mutable size_type  m_last_access_rl;
#endif

        size_type        m_size  = 0;    // original text size
        size_type        m_sigma = 0;    // alphabet size
        bit_vector_type  m_tree;         // bit vector to store the wavelet tree
        rank_1_type      m_tree_rank;    // rank support for the wavelet tree bit vector
        select_1_type    m_tree_select1; // select support for the wavelet tree bit vector
        select_0_type    m_tree_select0;

        _node<size_type> m_nodes[511];    // nodes for the Huffman tree structure
        uint16_t         m_c_to_leaf[256];// map symbol c to a leaf in the tree structure
        // if m_c_to_leaf[c] == _undef_node the char does
        // not exists in the text
        uint64_t         m_path[256];     // path information for each char; the bits at position
        // 0..55 hold path information; bits 56..63 the length
        // of the path in binary representation

        void copy(const wt_pc& wt) {
            m_size            = wt.m_size;
            m_sigma           = wt.m_sigma;
            m_tree            = wt.m_tree;
            m_tree_rank       = wt.m_tree_rank;
            m_tree_rank.set_vector(&m_tree);
            m_tree_select1    = wt.m_tree_select1;
            m_tree_select1.set_vector(&m_tree);
            m_tree_select0    = wt.m_tree_select0;
            m_tree_select0.set_vector(&m_tree);
            for (size_type i=0; i < 511; ++i)
                m_nodes[i] = wt.m_nodes[i];
            for (size_type i=0; i<256; ++i)
                m_c_to_leaf[i] = wt.m_c_to_leaf[i];
            for (size_type i=0; i<256; ++i) {
                m_path[i] = wt.m_path[i];
            }
        }

        // insert a character into the wavelet tree, see construct method
        void insert_char(uint8_t old_chr, size_type* tree_pos, size_type times, bit_vector& f_tree) {
            uint32_t path_len = (m_path[old_chr]>>56);
            uint64_t p = m_path[old_chr];
            for (uint32_t node=0, l=0; l<path_len; ++l, p >>= 1) {
                if (p&1) {
                    f_tree.set_int(tree_pos[node], 0xFFFFFFFFFFFFFFFFULL,times);
                }
                tree_pos[node] += times;
                node = m_nodes[node].child[p&1];
            }
        }



        // calculates the Huffman tree and returns the size of the WT bit vector
        size_type construct_tree_shape(const size_type* C) {
            std::vector<_node<size_type> > temp_nodes(2*m_sigma-1);  // vector for nodes of the Huffman tree
            size_type node_cnt = shape::construct_tree(C, temp_nodes);
            // Convert Huffman tree into breadth first search order in memory and
            // calculate tree_pos values
            m_nodes[0] = temp_nodes[node_cnt-1];  // insert root at index 0
            size_type tree_size = 0;
            node_cnt = 1;
            uint16_t last_parent = _undef_node;
            std::deque<size_type> q;
            q.push_back(0);
            while (!q.empty()) {
                size_type idx;
                if (!t_dfs_shape) {
                    idx = q.front(); q.pop_front();
                } else {
                    idx = q.back(); q.pop_back();
                }
                size_type frq = m_nodes[idx].tree_pos; // frq_sum was stored in tree_pos
                m_nodes[idx].tree_pos = tree_size;
                if (m_nodes[idx].child[0] != _undef_node)  // if node is not a leaf
                    tree_size += frq;                       // add frequency, as leaves have size 0
                if (idx > 0) { // node is not the root
                    if (last_parent != m_nodes[idx].parent)
                        m_nodes[m_nodes[idx].parent].child[0] = idx;
                    else
                        m_nodes[m_nodes[idx].parent].child[1] = idx;
                    last_parent = m_nodes[idx].parent;
                }
                if (m_nodes[idx].child[0] != _undef_node) { // if node is not a leaf
                    for (size_type k=0; k<2; ++k) {            // add children to tree
                        m_nodes[node_cnt] = temp_nodes[ m_nodes[idx].child[k] ];
                        m_nodes[node_cnt].parent = idx;
                        q.push_back(node_cnt);
                        m_nodes[idx].child[k] = node_cnt++;
                    }
                }
            }

            // initialize m_c_to_leaf
            for (size_type i=0; i<256; ++i)
                m_c_to_leaf[i] = _undef_node; // if c is not in the alphabet m_c_to_leaf[c] = _undef_node
            for (size_type i=0; i < 2*sigma-1; ++i) {
                if (m_nodes[i].child[0] == _undef_node)                 // if node is a leaf
                    m_c_to_leaf[(uint8_t)m_nodes[i].tree_pos_rank] = i; // calculate value
            }
            // initialize path information
            // Note: In the case of a bfs search order,
            // we can classify nodes as right child and left child with an easy criterion:
            //   node is a left child, if node%2==1
            //   node is a right child, if node%2==0
            for (size_type c=0; c<256; ++c) {
                if (m_c_to_leaf[c] != _undef_node) { // if char exists in the alphabet
                    size_type node = m_c_to_leaf[c];
                    uint64_t w = 0; // path
                    uint64_t l = 0; // path len
                    while (node != 0) { // while node is not the root
                        w <<= 1;
                        if (m_nodes[m_nodes[node].parent].child[1] == node) // if the node is a right child
                            w |= 1ULL;
                        ++l;
                        node = m_nodes[node].parent; // go up the tree
                    }
                    if (l > 56) {
                        std::cerr<<"Code tree has max depth > 56!!! ERROR"<<std::endl;
                        throw std::logic_error("Code tree depth is greater than 56!!!");
                    }
                    m_path[c] = w | (l << 56);
                } else {
                    m_path[c] = 0; // i.e. len is also 0, good for special case in rank()
                }
            }
            return tree_size;
        }

        void construct_init_rank_select() {
            util::init_support(m_tree_rank, &m_tree);
            util::init_support(m_tree_select0, &m_tree);
            util::init_support(m_tree_select1, &m_tree);
        }

        void construct_precalc_node_ranks() {
            for (size_type i=0; i<2*m_sigma-1; ++i) {
                if (m_nodes[i].child[0] != _undef_node)  // if node is not a leaf
                    m_nodes[i].tree_pos_rank = m_tree_rank(m_nodes[i].tree_pos);
            }
        }


        // recursive internal version of the method interval_symbols
        void _interval_symbols(size_type i, size_type j, size_type& k,
                               std::vector<value_type>& cs,
                               std::vector<size_type>& rank_c_i,
                               std::vector<size_type>& rank_c_j, uint16_t node) const {
            // invariant: j>i
            // goto right child
            size_type i_new = (m_tree_rank(m_nodes[node].tree_pos + i) - m_nodes[node].tree_pos_rank);
            size_type j_new = (m_tree_rank(m_nodes[node].tree_pos + j) - m_nodes[node].tree_pos_rank);
            // goto left child
            i -= i_new; j -= j_new;
            if (i != j) {
                uint16_t node_new = m_nodes[node].child[0];
                // if node is not a leaf
                if (m_nodes[node_new].child[0] != _undef_node) {
                    _interval_symbols(i, j, k, cs, rank_c_i, rank_c_j, node_new);
                } else {
                    rank_c_i[k] = i;
                    rank_c_j[k] = j;
                    cs[k++] = m_nodes[node_new].tree_pos_rank;
                }
            }
            // goto right child
            if (i_new!=j_new) {
                uint16_t node_new = m_nodes[node].child[1];
                // if node is not a leaf
                if (m_nodes[node_new].child[0] != _undef_node) {
                    _interval_symbols(i_new, j_new, k, cs, rank_c_i, rank_c_j, node_new);
                } else {
                    rank_c_i[k] = i_new;
                    rank_c_j[k] = j_new;
                    cs[k++] = m_nodes[node_new].tree_pos_rank;
                }
            }
        }


    public:

        const size_type&       sigma = m_sigma;
        const bit_vector_type& tree  = m_tree;

        // Default constructor
        wt_pc() {};

        //! Construct the wavelet tree from a file_buffer
        /*! \param input_buf    File buffer of the input.
         *  \param size         The length of the prefix of the random access container, for which the wavelet tree should be build.
         *    \par Time complexity
         *        \f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         */
        wt_pc(int_vector_file_buffer<8>& input_buf, size_type size):m_size(size) {
            if (0 == m_size)
                return;
            // O(n + |\Sigma|\log|\Sigma|) algorithm for calculating node sizes
            size_type C[256] = {0};
            // 1. Count occurrences of characters
            calculate_character_occurences(input_buf, m_size, C);
            // 2. Calculate effective alphabet size
            calculate_effective_alphabet_size(C, m_sigma);
            // 3. Generate tree shape
            size_type tree_size = construct_tree_shape(C);
            // 4. Generate wavelet tree bit sequence m_tree

            bit_vector tmp_tree(tree_size, 0);  // initialize bit_vector for the tree
            //  Calculate starting position of wavelet tree nodes
            size_type tree_pos[511];
            for (size_type i=0; i < 2*sigma-1; ++i) {
                tree_pos[i] = m_nodes[i].tree_pos;
            }
            input_buf.reset();
            if (input_buf.int_vector_size < size) {
                throw std::logic_error("Stream size is smaller than size!");
                return;
            }
            for (size_type i=0, r_sum=0, r = input_buf.load_next_block(); r_sum < m_size;) {
                if (r_sum + r > size) {  // read not more than size chars in the next loop
                    r = size-r_sum;
                }
                uint8_t old_chr = input_buf[i-r_sum], times = 0;
                for (; i < r_sum+r; ++i) {
                    uint8_t chr = input_buf[i-r_sum];
                    if (chr    != old_chr) {
                        insert_char(old_chr, tree_pos, times, tmp_tree);
                        times = 1;
                        old_chr = chr;
                    } else { // chr == old_chr
                        ++times;
                        if (times == 64) {
                            insert_char(old_chr, tree_pos, times, tmp_tree);
                            times = 0;
                        }
                    }
                }
                if (times > 0) {
                    insert_char(old_chr, tree_pos, times, tmp_tree);
                }
                r_sum += r; r = input_buf.load_next_block();
            }
            util::assign(m_tree, tmp_tree);
            // 5. Initialize rank and select data structures for m_tree
            construct_init_rank_select();
            // 6. Finish inner nodes by precalculating the tree_pos_rank values
            construct_precalc_node_ranks();
        }


        //! Copy constructor
        wt_pc(const wt_pc& wt) {
            copy(wt);
        }

        //! Assignment operator
        wt_pc& operator=(const wt_pc& wt) {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Swap operator
        void swap(wt_pc& wt) {
            if (this != &wt) {
                std::swap(m_size, wt.m_size);
                std::swap(m_sigma,  wt.m_sigma);
                m_tree.swap(wt.m_tree);
                util::swap_support(m_tree_rank, wt.m_tree_rank, &m_tree, &(wt.m_tree));

                util::swap_support(m_tree_select1, wt.m_tree_select1, &m_tree, &(wt.m_tree));
                util::swap_support(m_tree_select0, wt.m_tree_select0, &m_tree, &(wt.m_tree));

                for (size_type i=0; i < 511; ++i)
                    std::swap(m_nodes[i], wt.m_nodes[i]);
                for (size_type i=0; i<256; ++i)
                    std::swap(m_c_to_leaf[i], wt.m_c_to_leaf[i]);
                for (size_type i=0; i<256; ++i)
                    std::swap(m_path[i], wt.m_path[i]);
            }
        }

        //! Returns the size of the original vector.
        size_type size()const {
            return m_size;
        }

        //! Returns whether the wavelet tree contains no data.
        bool empty()const {
            return m_size == 0;
        }

        //! Recovers the i-th symbol of the original vector.
        /*! \param i The index of the symbol in the original vector. \f$i \in [0..size()-1]\f$
         *  \return The i-th symbol of the original vector.
         *  \par Time complexity
         *        \f$ \Order{H_0} \f$ on average, where \f$ H_0 \f$ is the zero order entropy
         *      of the sequence.
         */
        value_type operator[](size_type i)const {
            assert(i < size());
            // which stores how many of the next symbols are equal
            // with the current char
            size_type node = 0; // start at root node
            while (m_nodes[node].child[0] != _undef_node) { // while node is not a leaf
                if (m_tree[ m_nodes[node].tree_pos + i]) {  // goto the right child
                    i = m_tree_rank(m_nodes[node].tree_pos + i) - m_nodes[node].tree_pos_rank;
                    node = m_nodes[node].child[1];
                } else { // goto the left child
                    i -= (m_tree_rank(m_nodes[node].tree_pos + i) - m_nodes[node].tree_pos_rank);
                    node = m_nodes[node].child[0];
                }
            }
            return m_nodes[node].tree_pos_rank;
        };

        //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
         *  \param c The symbol to count the occurrences in the prefix.
         *    \return The number of occurrences of symbol c in the prefix [0..i-1] of the supported vector.
         *  \par Time complexity
         *        \f$ \Order{H_0} \f$
         */
        size_type rank(size_type i, value_type c)const {
            assert(i <= size());
            uint64_t p = m_path[c];
            uint32_t path_len = (m_path[c]>>56); // equals zero if char was not present in the original text or m_sigma=1
            if (!path_len and 1 == m_sigma) {    // if m_sigma == 1 return result immediately
                if (m_c_to_leaf[c] == _undef_node) { // if character does not exist return 0
                    return 0;
                }
                return std::min(i, m_size);
            }
            size_type result = i & ZoO[path_len>0]; // important: result has type size_type and ZoO has type size_type
            uint32_t node=0;
            for (uint32_t l=0; l<path_len and result; ++l, p >>= 1) {
                if (p&1) {
                    result     = (m_tree_rank(m_nodes[node].tree_pos+result) -  m_nodes[node].tree_pos_rank);
                } else {
                    result -= (m_tree_rank(m_nodes[node].tree_pos+result) -  m_nodes[node].tree_pos_rank);
                }
                node = m_nodes[node].child[p&1]; // goto child
            }
            return result;
        };

        //! Calculates how many occurrences of symbol wt[i] are in the prefix [0..i-1] of the original sequence.
        /*!
         *  \param i The index of the symbol.
         *  \param c Reference that will contain symbol wt[i].
         *  \return The number of occurrences of symbol wt[i] in the prefix [0..i-1]
         *  \par Time complexity
         *        \f$ \Order{H_0} \f$
         */
        size_type inverse_select(size_type i, value_type& c)const {
            assert(i < size());
            uint32_t node=0;
            while (m_nodes[node].child[0] != _undef_node) { // while node is not a leaf
                if (m_tree[m_nodes[node].tree_pos + i]) { // if bit is set goto right child
                    i     = (m_tree_rank(m_nodes[node].tree_pos + i) -  m_nodes[node].tree_pos_rank);
                    node = m_nodes[node].child[1];
                } else { // goto left child
                    i -= (m_tree_rank(m_nodes[node].tree_pos + i) -  m_nodes[node].tree_pos_rank);
                    node = m_nodes[node].child[0];
                }
            }
            c = m_nodes[node].tree_pos_rank;
            return i;
        }

        //! Calculates the i-th occurrence of the symbol c in the supported vector.
        /*!
         *  \param i The i-th occurrence. \f$i\in [1..rank(size(),c)]\f$.
         *  \param c The symbol c.
         *  \par Time complexity
         *        \f$ \Order{H_0} \f$
         */
        size_type select(size_type i, value_type c)const {
            assert(i > 0);
            assert(i <= rank(size(), c));
            uint16_t node = m_c_to_leaf[c];
            if (node == _undef_node) { // if c was not present in the original text
                return m_size;           // -> return a position right to the end
            }
            if (m_sigma == 1) {
                return std::min(i-1,m_size);
            }
            size_type result = i-1;        // otherwise
            uint64_t p = m_path[c];
            uint32_t path_len = (p>>56);
            p <<= (64-path_len); // Note: path_len > 0, since we have handled m_sigma = 1.
            for (uint32_t l=0; l<path_len; ++l, p <<= 1) {
                if ((p & 0x8000000000000000ULL)==0) { // node was a left child
                    node = m_nodes[node].parent;
                    result = m_tree_select0(m_nodes[node].tree_pos-m_nodes[node].tree_pos_rank + result + 1)
                             - m_nodes[node].tree_pos;
                } else { // node was a right child
                    node = m_nodes[node].parent;
                    result = m_tree_select1(m_nodes[node].tree_pos_rank + result + 1)
                             - m_nodes[node].tree_pos;
                }
            }
            return result;
        };


        //! Calculates for each symbol c in wt[i..j-1], how many times c occurs in wt[0..i-1] and wt[0..j-1].
        /*!
         *  \param i        The start index (inclusive) of the interval.
         *  \param j        The end index (exclusive) of the interval.
         *  \param k        Reference that will contain the number of different
         *                  symbols in wt[i..j-1].
         *  \param cs       Reference to a vector that will contain in
         *                  cs[0..k-1] all symbols that occur in wt[i..j-1] in
         *                  arbitrary order (for Huffman shape) and ascending
         *                  order (for Hu-Tucker shape).
         *  \param rank_c_i Reference to a vector which equals
         *                  rank_c_i[p] = rank(i,cs[p]), for \f$ 0 \leq p < k \f$
         *  \param rank_c_j Reference to a vector which equals
         *                  rank_c_j[p] = rank(j,cs[p]), for \f$ 0 \leq p < k \f$
         *    \par Time complexity
         *        \f$ \Order{\min{\sigma, k \log \sigma}} \f$
         *
         *  \par Precondition
         *       \f$ i \leq j \leq n \f$
         *       \f$ cs.size() \geq \sigma \f$
         *       \f$ rank_{c_i}.size() \geq \sigma \f$
         *       \f$ rank_{c_j}.size() \geq \sigma \f$
         */
        void interval_symbols(size_type i, size_type j, size_type& k,
                              std::vector<value_type>& cs,
                              std::vector<size_type>& rank_c_i,
                              std::vector<size_type>& rank_c_j) const {
            assert(i <= j and j <= size());
            if (i==j) {
                k = 0;
            } else if (1==m_sigma) {
                k = 1;
                cs[0] = m_nodes[0].tree_pos_rank;
                rank_c_i[0] = std::min(i,m_size);
                rank_c_j[0] = std::min(j,m_size);
            } else if ((j-i)==1) {
                k = 1;
                rank_c_i[0] = inverse_select(i, cs[0]);
                rank_c_j[0] = rank_c_i[0]+1;
            } else if ((j-i)==2) {
                rank_c_i[0] = inverse_select(i, cs[0]);
                rank_c_i[1] = inverse_select(i+1, cs[1]);
                if (cs[0]==cs[1]) {
                    k = 1;
                    rank_c_j[0] = rank_c_i[0]+2;
                } else {
                    k = 2;
                    if (lex_ordered and cs[0] > cs[1]) {
                        std::swap(cs[0], cs[1]);
                        std::swap(rank_c_i[0], rank_c_i[1]);
                    }
                    rank_c_j[0] = rank_c_i[0]+1;
                    rank_c_j[1] = rank_c_i[1]+1;
                }
            } else {
                k = 0;
                _interval_symbols(i, j, k, cs, rank_c_i, rank_c_j, 0);
            }
        }


        //! Calculates for symbol c, how many symbols smaller and greater c occur in wt[i..j-1].
        /*!
         *  \param i       The start index (inclusive) of the interval.
         *  \param j       The end index (exclusive) of the interval.
         *  \param c       The symbol to count the occurences in the interval.
         *  \param smaller Reference that will contain the number of symbols smaller than c in wt[i..j-1].
         *  \param greater Reference that will contain the number of symbols greater than c in wt[i..j-1].
         *  \return The number of occurrences of symbol c in wt[0..i-1].
         *
         *  \par Precondition
         *       \f$ i \leq j \leq n \f$
         *       \f$ c must exist in wt \f$
         */
        size_type lex_count(size_type i, size_type j, value_type c, size_type& smaller, size_type& greater)const {
            if (lex_ordered) {
                assert(i <= j and j <= size());
                smaller = 0;
                greater = 0;
                if (1==m_sigma) {
                    return i;
                }
                if (i==j) {
                    return rank(i,c);
                }
                uint64_t p = m_path[c];
                uint32_t path_len = (m_path[c]>>56); // equals zero if char was not present in the original text
                assert(path_len>0);
                size_type res1 = i;
                size_type res2 = j;
                uint32_t node=0;
                for (uint32_t l=0; l<path_len; ++l, p >>= 1) {
                    if (p&1) {
                        size_type r1_1 = (m_tree_rank(m_nodes[node].tree_pos+res1)-m_nodes[node].tree_pos_rank);
                        size_type r1_2 = (m_tree_rank(m_nodes[node].tree_pos+res2)-m_nodes[node].tree_pos_rank);

                        smaller += res2 - r1_2 - res1 + r1_1;

                        res1 = r1_1;
                        res2 = r1_2;
                    } else {
                        size_type r1_1 = (m_tree_rank(m_nodes[node].tree_pos+res1)-m_nodes[node].tree_pos_rank);
                        size_type r1_2 = (m_tree_rank(m_nodes[node].tree_pos+res2)-m_nodes[node].tree_pos_rank);

                        greater += r1_2 - r1_1;

                        res1 -= r1_1;
                        res2 -= r1_2;
                    }
                    node = m_nodes[node].child[p&1];
                }
                return res1;
            } else {
                throw std::logic_error("lex_count is not supported");
            }
        };

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += write_member(m_sigma, out, child, "sigma");
            written_bytes += m_tree.serialize(out, child, "tree");
            written_bytes += m_tree_rank.serialize(out, child, "tree_rank");
            written_bytes += m_tree_select1.serialize(out, child, "tree_select_1");
            written_bytes += m_tree_select0.serialize(out, child, "tree_select_0");
            for (size_type i=0; i < 511; ++i) {               // TODO: use serialize vector
                written_bytes += m_nodes[i].serialize(out);   // is it surely possible to use
            }                                                 // less space
            out.write((char*) m_c_to_leaf, 256*sizeof(m_c_to_leaf[0]));
            written_bytes += 256*sizeof(m_c_to_leaf[0]); // add written bytes from previous loop
            out.write((char*) m_path, 256*sizeof(m_path[0]));
            written_bytes += 256*sizeof(m_path[0]); // add written bytes from previous loop
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            read_member(m_size, in);
            read_member(m_sigma, in);
            m_tree.load(in);
            m_tree_rank.load(in, &m_tree);
            m_tree_select1.load(in, &m_tree);
            m_tree_select0.load(in, &m_tree);
            for (size_type i=0; i < 511; ++i) {
                m_nodes[i].load(in);
            }
            in.read((char*) m_c_to_leaf, 256*sizeof(m_c_to_leaf[0]));
            in.read((char*) m_path, 256*sizeof(m_path[0]));
        }

};

}

#endif
