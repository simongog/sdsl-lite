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
#include <utility>
#include <array>

//! Namespace for the succinct data structure library.
namespace sdsl
{

const int_vector<>::size_type ZoO[2] = {0, (int_vector<>::size_type)-1};

template<class t_shape,
         class t_bitvector   = bit_vector,
         class t_rank        = typename t_bitvector::rank_1_type,
         class t_select      = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type,
         class t_tree_strat  = byte_tree<>
         >
class wt_pc
{
    public:
        typedef typename
        t_tree_strat::template type<wt_pc> tree_strat_type;
        typedef int_vector<>::size_type    size_type;
        typedef typename
        tree_strat_type::value_type        value_type;
        typedef t_bitvector                bit_vector_type;
        typedef t_rank                     rank_1_type;
        typedef t_select                   select_1_type;
        typedef t_select_zero              select_0_type;
        typedef wt_tag                     index_category;
        typedef typename
        tree_strat_type::alphabet_category alphabet_category;
        typedef typename
        t_shape::template type<wt_pc>      shape_type;
        enum { lex_ordered=shape_type::lex_ordered };
        using node_type = typename tree_strat_type::node_type;

    private:

#ifdef WT_HUFF_CACHE
        mutable value_type m_last_access_answer;
        mutable size_type  m_last_access_i;
        mutable size_type  m_last_access_rl;
#endif

        size_type        m_size  = 0;    // original text size
        size_type        m_sigma = 0;    // alphabet size
        bit_vector_type  m_bv;           // bit vector to store the wavelet tree
        rank_1_type      m_bv_rank;      // rank support for the wavelet tree bit vector
        select_1_type    m_bv_select1;   // select support for the wavelet tree bit vector
        select_0_type    m_bv_select0;
        tree_strat_type  m_tree;

        void copy(const wt_pc& wt) {
            m_size            = wt.m_size;
            m_sigma           = wt.m_sigma;
            m_bv            = wt.m_bv;
            m_bv_rank       = wt.m_bv_rank;
            m_bv_rank.set_vector(&m_bv);
            m_bv_select1    = wt.m_bv_select1;
            m_bv_select1.set_vector(&m_bv);
            m_bv_select0    = wt.m_bv_select0;
            m_bv_select0.set_vector(&m_bv);
            m_tree          = wt.m_tree;
        }

        // insert a character into the wavelet tree, see construct method
        void insert_char(value_type old_chr, std::vector<uint64_t>& bv_node_pos,
                         size_type times, bit_vector& bv) {
            uint64_t p = m_tree.bit_path(old_chr);
            uint32_t path_len = p>>56;
            node_type v = m_tree.root();
            for (uint32_t l=0; l<path_len; ++l, p >>= 1) {
                if (p&1) {
                    bv.set_int(bv_node_pos[v], 0xFFFFFFFFFFFFFFFFULL,times);
                }
                bv_node_pos[v] += times;
                v = m_tree.child(v, p&1);
            }
        }



        // calculates the tree shape returns the size of the WT bit vector
        size_type construct_tree_shape(const std::vector<size_type>& C) {
            // vector  for node of the tree
            std::vector<pc_node> temp_nodes; //(2*m_sigma-1);
            shape_type::construct_tree(C, temp_nodes);
            // Convert code tree into BFS order in memory and
            // calculate bv_pos values
            size_type bv_size = 0;
            tree_strat_type temp_tree(temp_nodes, bv_size, this);
            m_tree.swap(temp_tree);
            return bv_size;
        }

        void construct_init_rank_select() {
            util::init_support(m_bv_rank, &m_bv);
            util::init_support(m_bv_select0, &m_bv);
            util::init_support(m_bv_select1, &m_bv);
        }

        // recursive internal version of the method interval_symbols
        void
        _interval_symbols(size_type i, size_type j, size_type& k,
                          std::vector<value_type>& cs,
                          std::vector<size_type>& rank_c_i,
                          std::vector<size_type>& rank_c_j, node_type v) const {
            // invariant: j>i
            // goto right child
            size_type i_new = (m_bv_rank(m_tree.bv_pos(v) + i)
                               - m_tree.bv_pos_rank(v));
            size_type j_new = (m_bv_rank(m_tree.bv_pos(v) + j)
                               - m_tree.bv_pos_rank(v));
            // goto left child
            i -= i_new; j -= j_new;
            if (i != j) {
                node_type v_new = m_tree.child(v, 0);
                if (!m_tree.is_leaf(v_new)) {
                    _interval_symbols(i, j, k, cs, rank_c_i, rank_c_j, v_new);
                } else {
                    rank_c_i[k] = i;
                    rank_c_j[k] = j;
                    cs[k++] = m_tree.bv_pos_rank(v_new);
                }
            }
            // goto right child
            if (i_new!=j_new) {
                node_type v_new = m_tree.child(v, 1);
                if (!m_tree.is_leaf(v_new)) {
                    _interval_symbols(i_new, j_new, k, cs, rank_c_i, rank_c_j,
                                      v_new);
                } else {
                    rank_c_i[k] = i_new;
                    rank_c_j[k] = j_new;
                    cs[k++] = m_tree.bv_pos_rank(v_new);
                }
            }
        }

    public:

        const size_type&       sigma = m_sigma;
        const bit_vector_type& bv  = m_bv;

        // Default constructor
        wt_pc() {};

        //! Construct the wavelet tree from a file_buffer
        /*! \param input_buf    File buffer of the input.
         *  \param size         The length of the prefix.
         *    \par Time complexity
         *        \f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         */
        wt_pc(int_vector_buffer<tree_strat_type::int_width>& input_buf,
              size_type size):m_size(size) {
            if (0 == m_size)
                return;
            // O(n + |\Sigma|\log|\Sigma|) algorithm for calculating node sizes
            // TODO: C should also depend on the tree_strategy. C is just a mapping
            // from a symbol to its frequency. So a map<uint64_t,uint64_t> could be
            // used for integer alphabets...
            std::vector<size_type> C;
            // 1. Count occurrences of characters
            calculate_character_occurences(input_buf, m_size, C);
            // 2. Calculate effective alphabet size
            calculate_effective_alphabet_size(C, m_sigma);
            // 3. Generate tree shape
            size_type tree_size = construct_tree_shape(C);
            // 4. Generate wavelet tree bit sequence m_bv
            bit_vector temp_bv(tree_size, 0);

            // Initializing starting position of wavelet tree nodes
            std::vector<uint64_t> bv_node_pos(m_tree.size(), 0);
            for (size_type v=0; v < m_tree.size(); ++v) {
                bv_node_pos[v] = m_tree.bv_pos(v);
            }
            if (input_buf.size() < size) {
                throw std::logic_error("Stream size is smaller than size!");
                return;
            }
            value_type old_chr = input_buf[0];
            uint32_t times = 0;
            for (size_type i=0; i < m_size; ++i) {
                value_type chr = input_buf[i];
                if (chr != old_chr) {
                    insert_char(old_chr, bv_node_pos, times, temp_bv);
                    times = 1;
                    old_chr = chr;
                } else { // chr == old_chr
                    ++times;
                    if (times == 64) {
                        insert_char(old_chr, bv_node_pos, times, temp_bv);
                        times = 0;
                    }
                }
            }
            if (times > 0) {
                insert_char(old_chr, bv_node_pos, times, temp_bv);
            }
            m_bv = bit_vector_type(std::move(temp_bv));
            // 5. Initialize rank and select data structures for m_bv
            construct_init_rank_select();
            // 6. Finish inner nodes by precalculating the bv_pos_rank values
            m_tree.init_node_ranks(m_bv_rank);
        }


        //! Copy constructor
        wt_pc(const wt_pc& wt) { copy(wt); }

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
                m_bv.swap(wt.m_bv);
                util::swap_support(m_bv_rank, wt.m_bv_rank,
                                   &m_bv, &(wt.m_bv));

                util::swap_support(m_bv_select1, wt.m_bv_select1,
                                   &m_bv, &(wt.m_bv));
                util::swap_support(m_bv_select0, wt.m_bv_select0,
                                   &m_bv, &(wt.m_bv));
                m_tree.swap(wt.m_tree);
            }
        }

        //! Returns the size of the original vector.
        size_type size()const { return m_size; }

        //! Returns whether the wavelet tree contains no data.
        bool empty()const { return m_size == 0; }

        //! Recovers the i-th symbol of the original vector.
        /*! \param i Index in the original vector. \f$i \in [0..size()-1]\f$.
         *  \return The i-th symbol of the original vector.
         *  \par Time complexity
         *        \f$ \Order{H_0} \f$ on average, where \f$ H_0 \f$ is the
         *        zero order entropy of the sequence
         */
        value_type operator[](size_type i)const {
            assert(i < size());
            // which stores how many of the next symbols are equal
            // with the current char
            node_type v = m_tree.root(); // start at root node
            while (!m_tree.is_leaf(v)) {   // while  not a leaf
                if (m_bv[ m_tree.bv_pos(v) + i]) {  // goto right child
                    i = m_bv_rank(m_tree.bv_pos(v) + i)
                        - m_tree.bv_pos_rank(v);
                    v = m_tree.child(v,1);
                } else { // goto the left child
                    i -= (m_bv_rank(m_tree.bv_pos(v) + i)
                          - m_tree.bv_pos_rank(v));
                    v = m_tree.child(v,0);
                }
            }
            // if v is a leaf bv_pos_rank returns symbol itself
            return m_tree.bv_pos_rank(v);
        };

        //! Calculates how many symbols c are in the prefix [0..i-1].
        /*!
         *  \param i Exclusive right bound of the range (\f$i\in[0..size()]\f$).
         *  \param c Symbol c.
         *  \return Number of occurrences of symbol c in the prefix [0..i-1].
         *  \par Time complexity
         *        \f$ \Order{H_0} \f$ on average, where \f$ H_0 \f$ is the
         *        zero order entropy of the sequence
         */
        size_type rank(size_type i, value_type c)const {
            assert(i <= size());
            uint64_t p = m_tree.bit_path(c);
            // path_len == 0, if `c` was not in the text or m_sigma=1
            uint32_t path_len = (p>>56);
            if (!path_len and 1 == m_sigma) {
                if (!m_tree.is_valid(m_tree.c_to_leaf(c))) {   // if `c` was not in the text
                    return 0;
                }
                return std::min(i, m_size); // if m_sigma == 1 answer is trivial
            }
            size_type result = i & ZoO[path_len>0];
            node_type v = m_tree.root();
            for (uint32_t l=0; l<path_len and result; ++l, p >>= 1) {
                if (p&1) {
                    result  = (m_bv_rank(m_tree.bv_pos(v)+result)
                               -  m_tree.bv_pos_rank(v));
                } else {
                    result -= (m_bv_rank(m_tree.bv_pos(v)+result)
                               -  m_tree.bv_pos_rank(v));
                }
                v = m_tree.child(v, p&1); // goto child
            }
            return result;
        };

        //! Calculates how many times symbol wt[i] occurs in the prefix [0..i-1].
        /*!
         *  \param i The index of the symbol.
         *  \param c Reference that will contain the symbol at position i.
         *  \return  Number of occurrences of symbol wt[i] in the prefix [0..i-1].
         *    \par Time complexity
         *        \f$ \Order{H_0} \f$
         */
        size_type inverse_select(size_type i, value_type& c)const {
            assert(i < size());
            node_type v = m_tree.root();
            while (!m_tree.is_leaf(v)) {   // while not a leaf
                if (m_bv[m_tree.bv_pos(v) + i]) {   //  goto right child
                    i = (m_bv_rank(m_tree.bv_pos(v) + i)
                         - m_tree.bv_pos_rank(v));
                    v = m_tree.child(v, 1);
                } else { // goto left child
                    i -= (m_bv_rank(m_tree.bv_pos(v) + i)
                          - m_tree.bv_pos_rank(v));
                    v = m_tree.child(v,0);
                }
            }
            // if v is a leaf bv_pos_rank returns symbol itself
            c = m_tree.bv_pos_rank(v);
            return i;
        }

        //! Calculates the ith occurrence of the symbol c in the supported vector.
        /*!
         *  \param i The ith occurrence. \f$i\in [1..rank(size(),c)]\f$.
         *  \param c The symbol c.
         *  \par Time complexity
         *       \f$ \Order{H_0} \f$ on average, where \f$ H_0 \f$ is the zero order
         *        entropy of the sequence
         */
        size_type select(size_type i, value_type c)const {
            assert(i > 0);
            assert(i <= rank(size(), c));
            node_type v = m_tree.c_to_leaf(c);
            if (!m_tree.is_valid(v)) {   // if c was not in the text
                return m_size;         // -> return a position right to the end
            }
            if (m_sigma == 1) {
                return std::min(i-1,m_size);
            }
            size_type result = i-1;    // otherwise
            uint64_t p = m_tree.bit_path(c);
            uint32_t path_len = (p>>56);
            // path_len > 0, since we have handled m_sigma = 1.
            p <<= (64-path_len);
            for (uint32_t l=0; l<path_len; ++l, p <<= 1) {
                if ((p & 0x8000000000000000ULL)==0) { // node was a left child
                    v  = m_tree.parent(v);
                    result = m_bv_select0(m_tree.bv_pos(v)
                                          - m_tree.bv_pos_rank(v) + result + 1)
                             - m_tree.bv_pos(v);
                } else { // node was a right child
                    v   = m_tree.parent(v);
                    result = m_bv_select1(m_tree.bv_pos_rank(v) + result + 1)
                             - m_tree.bv_pos(v);
                }
            }
            return result;
        };


        //! For each symbol c in wt[i..j-1] get rank(i,c) and rank(j,c).
        /*!
         *  \param i        The start index (inclusive) of the interval.
         *  \param j        The end index (exclusive) of the interval.
         *  \param k        Reference for number of different symbols in [i..j-1].
         *  \param cs       Reference to a vector that will contain in
         *                  cs[0..k-1] all symbols that occur in [i..j-1] in
         *                  arbitrary order (for Huffman shape) and ascending
         *                  order (for Hu-Tucker shape).
         *  \param rank_c_i Reference to a vector which equals
         *                  rank_c_i[p] = rank(i,cs[p]), for \f$ 0 \leq p < k \f$.
         *  \param rank_c_j Reference to a vector which equals
         *                  rank_c_j[p] = rank(j,cs[p]), for \f$ 0 \leq p < k \f$.
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
                cs[0] = m_tree.bv_pos_rank(m_tree.root());
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


        //! How many symbols are lexicographic smaller/greater than c in [i..j-1].
        /*!
         *  \param i       Start index (inclusive) of the interval.
         *  \param j       End index (exclusive) of the interval.
         *  \param c       Symbol c.
         *  \param smaller Reference for symbols smaller than c in [i..j-1].
         *  \param greater Reference for symbols greater than c in [i..j-1].
         *  \return The number of occurrences of symbol c in [0..i-1].
         *
         *  \par Precondition
         *       \f$ i \leq j \leq n \f$
         *       \f$ c must exist in wt \f$
         */
        size_type lex_count(size_type i, size_type j, value_type c,
                            size_type& smaller, size_type& greater) const {
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
                uint64_t p = m_tree.bit_path(c);
                uint32_t path_len = p>>56;
                // path_len equals zero if c was not present
                assert(path_len>0);
                size_type res1 = i;
                size_type res2 = j;
                node_type v = m_tree.root();
                for (uint32_t l=0; l<path_len; ++l, p >>= 1) {
                    if (p&1) {
                        size_type r1_1 = (m_bv_rank(m_tree.bv_pos(v)+res1)
                                          - m_tree.bv_pos_rank(v));
                        size_type r1_2 = (m_bv_rank(m_tree.bv_pos(v)+res2)
                                          - m_tree.bv_pos_rank(v));

                        smaller += res2 - r1_2 - res1 + r1_1;

                        res1 = r1_1;
                        res2 = r1_2;
                    } else {
                        size_type r1_1 = (m_bv_rank(m_tree.bv_pos(v)+res1)
                                          - m_tree.bv_pos_rank(v));
                        size_type r1_2 = (m_bv_rank(m_tree.bv_pos(v)+res2)
                                          - m_tree.bv_pos_rank(v));

                        greater += r1_2 - r1_1;

                        res1 -= r1_1;
                        res2 -= r1_2;
                    }
                    v = m_tree.child(v, p&1);
                }
                return res1;
            } else {
                throw std::logic_error("lex_count is not supported");
            }
        };

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr,
                            std::string name="") const {
            structure_tree_node* child = structure_tree::add_child(
                                             v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size,out,child, "size");
            written_bytes += write_member(m_sigma,out,child, "sigma");
            written_bytes += m_bv.serialize(out,child,"bv");
            written_bytes += m_bv_rank.serialize(out,child,"bv_rank");
            written_bytes += m_bv_select1.serialize(out,child,"bv_select_1");
            written_bytes += m_bv_select0.serialize(out,child,"bv_select_0");
            written_bytes += m_tree.serialize(out,child,"tree");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            read_member(m_size, in);
            read_member(m_sigma, in);
            m_bv.load(in);
            m_bv_rank.load(in, &m_bv);
            m_bv_select1.load(in, &m_bv);
            m_bv_select0.load(in, &m_bv);
            m_tree.load(in);
        }

};

}

#endif
