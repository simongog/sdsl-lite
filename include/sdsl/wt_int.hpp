/* sdsl - succinct data structures library
    Copyright (C) 2009 Simon Gog

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
/*! \file wt_int.hpp
    \brief wt_int.hpp contains a specialized class for a wavelet tree of a
           permutation of the numbers from 0..n. This wavelet tree class takes
           less memory than the generic class for wavelet trees.
    \author Simon Gog, Shanika Kuruppu
*/
#ifndef INCLUDED_SDSL_INT_WAVELET_TREE
#define INCLUDED_SDSL_INT_WAVELET_TREE

#include "sdsl_concepts.hpp"
#include "int_vector.hpp"
#include "rank_support_v.hpp"
#include "select_support_mcl.hpp"
#include "temp_write_read_buffer.hpp"
#include "util.hpp"
#include <set> // for calculating the alphabet size
#include <map> // for mapping a symbol to its lexicographical index
#include <algorithm> // for std::swap
#include <stdexcept>
#include <vector>
#include <queue>
#include <utility>

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A wavelet tree class for integer sequences.
/*!
 *    \par Space complexity
 *        \f$\Order{n\log|\Sigma|}\f$ bits, where \f$n\f$ is the size of the vector the wavelet tree was build for.
 *
 *  \tparam t_bitvector   Type of the bitvector used for representing the wavelet tree.
 *  \tparam t_rank        Type of the support structure for rank on pattern `1`.
 *  \tparam t_select      Type of the support structure for select on pattern `1`.
 *  \tparam t_select_zero Type of the support structure for select on pattern `0`.
 *
 *   @ingroup wt
 */
template<class t_bitvector   = bit_vector,
         class t_rank        = typename t_bitvector::rank_1_type,
         class t_select      = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type>
class wt_int
{
    public:

        typedef int_vector<>::size_type              size_type;
        typedef int_vector<>::value_type             value_type;
        typedef typename t_bitvector::difference_type difference_type;
        typedef random_access_const_iterator<wt_int> const_iterator;
        typedef const_iterator                       iterator;
        typedef t_bitvector                          bit_vector_type;
        typedef t_rank                               rank_1_type;
        typedef t_select                             select_1_type;
        typedef t_select_zero                        select_0_type;
        typedef wt_tag                               index_category;
        typedef int_alphabet_tag                     alphabet_category;
        enum 	{lex_ordered=1};

        typedef std::pair<value_type, size_type>     point_type;
        typedef std::vector<point_type>              point_vec_type;
        typedef std::pair<size_type, point_vec_type> r2d_res_type;

        mutable std::set<size_type> m_access_nodes;
        mutable std::set<size_type> m_greedy_nodes;
        mutable std::set<size_type> m_greedy_sim_nodes;
    protected:

        size_type              m_size  = 0;
        size_type              m_sigma = 0;    //<- \f$ |\Sigma| \f$
        bit_vector_type        m_tree;         // bit vector to store the wavelet tree
        rank_1_type            m_tree_rank;    // rank support for the wavelet tree bit vector
        select_1_type          m_tree_select1; // select support for the wavelet tree bit vector
        select_0_type          m_tree_select0;
        uint32_t               m_max_depth = 0;
        mutable int_vector<64> m_path_off;     // array keeps track of path offset in select-like methods
        mutable int_vector<64> m_path_rank_off;// array keeps track of rank values for the offsets

        void copy(const wt_int& wt) {
            m_size          = wt.m_size;
            m_sigma         = wt.m_sigma;
            m_tree          = wt.m_tree;
            m_tree_rank     = wt.m_tree_rank;
            m_tree_rank.set_vector(&m_tree);
            m_tree_select1  = wt.m_tree_select1;
            m_tree_select1.set_vector(&m_tree);
            m_tree_select0  = wt.m_tree_select0;
            m_tree_select0.set_vector(&m_tree);
            m_max_depth     = wt.m_max_depth;
            m_path_off      = wt.m_path_off;
            m_path_rank_off = wt.m_path_rank_off;
        }

    private:

        void init_buffers(uint32_t max_depth) {
            m_path_off = int_vector<64>(max_depth+1);
            m_path_rank_off = int_vector<64>(max_depth+1);
        }

        // recursive internal version of the method interval_symbols
        void _interval_symbols(size_type i, size_type j, size_type& k,
                               std::vector<value_type>& cs,
                               std::vector<size_type>& rank_c_i,
                               std::vector<size_type>& rank_c_j,
                               size_type depth,
                               size_type path,
                               size_type node_size,
                               size_type offset) const {
            // invariant: j>i

            if (depth >= m_max_depth) {
                rank_c_i[k]= i;
                rank_c_j[k]= j;
                cs[k++]= path;
                return;
            }

            size_type ones_before_o = m_tree_rank(offset);
            size_type ones_before_i = m_tree_rank(offset+i) - ones_before_o;
            size_type ones_before_j = m_tree_rank(offset+j) - ones_before_o;
            size_type ones_before_end = m_tree_rank(offset+ node_size) - ones_before_o;

            // goto left child
            if ((j-i)-(ones_before_j-ones_before_i)>0) {
                size_type new_offset = offset + m_size;
                size_type new_node_size = node_size - ones_before_end;
                size_type new_i = i - ones_before_i;
                size_type new_j = j - ones_before_j;
                _interval_symbols(new_i, new_j, k, cs, rank_c_i, rank_c_j, depth+1, path<<1, new_node_size, new_offset);
            }

            // goto right child
            if ((ones_before_j-ones_before_i)>0) {
                size_type new_offset = offset+(node_size - ones_before_end) + m_size;
                size_type new_node_size = ones_before_end;
                size_type new_i = ones_before_i;
                size_type new_j = ones_before_j;
                _interval_symbols(new_i, new_j, k, cs, rank_c_i, rank_c_j, depth+1, (path<<1)|1, new_node_size, new_offset);
            }
        }

    public:

        const size_type&       sigma = m_sigma; //!< Effective alphabet size of the wavelet tree.
        const bit_vector_type& tree  = m_tree;  //!< A concatenation of all bit vectors of the wavelet tree.

        //! Default constructor
        wt_int() {
            init_buffers(m_max_depth);
        };

        //! Semi-external constructor
        /*! \param buf         File buffer of the int_vector for which the wt_int should be build.
         *  \param size        Size of the prefix of v, which should be indexed.
         *  \param max_depth   Maximal depth of the wavelet tree. If set to 0, determined automatically.
         *    \par Time complexity
         *        \f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         *        I.e. we need \Order{n\log n} if rac is a permutation of 0..n-1.
         *    \par Space complexity
         *        \f$ n\log|\Sigma| + O(1)\f$ bits, where \f$n=size\f$.
         */
        template<uint8_t int_width>
        wt_int(int_vector_buffer<int_width>& buf, size_type size,
               uint32_t max_depth=0) : m_size(size) {
            init_buffers(m_max_depth);
            if (0 == m_size)
                return;
            size_type n = buf.size();  // set n
            if (n < m_size) {
                throw std::logic_error("n="+util::to_string(n)+" < "+util::to_string(m_size)+"=m_size");
                return;
            }
            m_sigma = 0; // init sigma

            std::string dir = util::dirname(buf.filename());

            temp_write_read_buffer<> buf1(5000000, buf.width(), dir);   // buffer for elements in the right node
            int_vector<int_width> rac(m_size, 0, buf.width());          // initialize rac

            value_type x = 1;  // variable for the biggest value in rac
            for (size_type i=0; i < m_size; ++i) { // detect the largest value in rac
                if (buf[i] > x)
                    x = buf[i];
                rac[i] = buf[i];
            }

            if (max_depth == 0) {
                m_max_depth = bits::hi(x)+1; // we need max_depth bits to represent all values in the range [0..x]
            } else {
                m_max_depth = max_depth;
            }
            init_buffers(m_max_depth);

            std::string tree_out_buf_file_name = (dir+"/m_tree"+util::to_string(util::pid())+"_"+util::to_string(util::id()));
            osfstream tree_out_buf(tree_out_buf_file_name, std::ios::binary | std::ios::trunc | std::ios::out);   // open buffer for tree
            size_type bit_size = m_size*m_max_depth;
            tree_out_buf.write((char*) &bit_size, sizeof(bit_size));    // write size of bit_vector

            size_type tree_pos = 0;
            uint64_t tree_word = 0;

            uint64_t        mask_old = 1ULL<<(m_max_depth);
            for (uint32_t k=0; k<m_max_depth; ++k) {
                size_type          start     = 0;
                const uint64_t    mask_new = 1ULL<<(m_max_depth-k-1);
                do {
                    buf1.reset();
                    size_type i           = start;
                    size_type cnt0        =    0;
                    uint64_t  start_value = (rac[i]&mask_old);
                    uint64_t  x;
                    while (i < m_size and((x=rac[i])&mask_old)==start_value) {
                        if (x&mask_new) {
                            tree_word |= (1ULL << (tree_pos&0x3FULL));
                            buf1 << x;
                        } else {
                            rac[start + cnt0++ ] = x;
                        }
                        ++tree_pos;
                        if ((tree_pos & 0x3FULL) == 0) { // if tree_pos % 64 == 0 write old word
                            tree_out_buf.write((char*) &tree_word, sizeof(tree_word));
                            tree_word = 0;
                        }
                        ++i;
                    }
                    buf1.write_close();
                    size_type cnt1 = i-start-cnt0;
                    if (k+1 < m_max_depth) { // inner node
                        for (i=start + cnt0, start = start+cnt0+cnt1; i < start; ++i) {
                            buf1 >> x;
                            rac[ i ] = x;
                        }
                    } else { // leaf node
                        start += cnt0+cnt1;
                        m_sigma += (cnt0>0) + (cnt1>0); // increase sigma for each leaf
                    }

                } while (start < m_size);
                mask_old += mask_new;
            }
            if ((tree_pos & 0x3FULL) != 0) { // if tree_pos % 64 > 0 => there are remaining entries we have to write
                tree_out_buf.write((char*) &tree_word, sizeof(tree_word));
            }
            tree_out_buf.close();
            rac.resize(0);
            bit_vector tree;
            load_from_file(tree, tree_out_buf_file_name);
            sdsl::remove(tree_out_buf_file_name);
            m_tree = bit_vector_type(std::move(tree));
            util::init_support(m_tree_rank, &m_tree);
            util::init_support(m_tree_select0, &m_tree);
            util::init_support(m_tree_select1, &m_tree);
        }

        //! Copy constructor
        wt_int(const wt_int& wt) {
            copy(wt);
        }

        //! Assignment operator
        wt_int& operator=(const wt_int& wt) {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Swap operator
        void swap(wt_int& wt) {
            if (this != &wt) {
                std::swap(m_size, wt.m_size);
                std::swap(m_sigma,  wt.m_sigma);
                m_tree.swap(wt.m_tree);
                util::swap_support(m_tree_rank, wt.m_tree_rank, &m_tree, &(wt.m_tree));
                util::swap_support(m_tree_select1, wt.m_tree_select1, &m_tree, &(wt.m_tree));
                util::swap_support(m_tree_select0, wt.m_tree_select0, &m_tree, &(wt.m_tree));
                std::swap(m_max_depth,  wt.m_max_depth);
                m_path_off.swap(wt.m_path_off);
                m_path_rank_off.swap(wt.m_path_rank_off);
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
        /*! \param i The index of the symbol in the original vector.
         *  \returns The i-th symbol of the original vector.
         *  \par Precondition
         *       \f$ i < size() \f$
         */
        value_type operator[](size_type i)const {
            assert(i < size());
            size_type offset = 0;
            value_type res = 0;
            size_type node_size = m_size;
            size_type id = 1;
            //m_access_nodes.insert(id);
            for (uint32_t k=0; k < m_max_depth; ++k) {
                res <<= 1;
                size_type ones_before_o   = m_tree_rank(offset);
                size_type ones_before_i   = m_tree_rank(offset + i) - ones_before_o;
                size_type ones_before_end = m_tree_rank(offset + node_size) - ones_before_o;
                if (m_tree[offset+i]) { // one at position i => follow right child
                    offset += (node_size - ones_before_end);
                    node_size = ones_before_end;
                    id = id*2+1;
                    i = ones_before_i;
                    res |= 1;
                    //m_access_nodes.insert(id);
                } else { // zero at position i => follow left child
                    node_size = (node_size - ones_before_end);
                    id = id*2;
                    i = (i-ones_before_i);
                    //m_access_nodes.insert(id);
                }
                offset += m_size;
            }
            return res;
        };

        //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
         *  \param c The symbol to count the occurrences in the prefix.
         *    \returns The number of occurrences of symbol c in the prefix [0..i-1] of the supported vector.
         *  \par Time complexity
         *       \f$ \Order{\log |\Sigma|} \f$
         *  \par Precondition
         *       \f$ i \leq size() \f$
         */
        size_type rank(size_type i, value_type c)const {
            assert(i <= size());
            if (((1ULL)<<(m_max_depth))<=c) { // c is greater than any symbol in wt
                return 0;
            }
            size_type offset = 0;
            uint64_t mask = (1ULL) << (m_max_depth-1);
            size_type node_size = m_size;
            for (uint32_t k=0; k < m_max_depth and i; ++k) {
                size_type ones_before_o   = m_tree_rank(offset);
                size_type ones_before_i   = m_tree_rank(offset + i) - ones_before_o;
                size_type ones_before_end = m_tree_rank(offset + node_size) - ones_before_o;
                if (c & mask) { // search for a one at this level
                    offset += (node_size - ones_before_end);
                    node_size = ones_before_end;
                    i = ones_before_i;
                } else { // search for a zero at this level
                    node_size = (node_size - ones_before_end);
                    i = (i-ones_before_i);
                }
                offset += m_size;
                mask >>= 1;
            }
            return i;
        };



        //! Calculates how many occurrences of symbol wt[i] are in the prefix [0..i-1] of the original sequence.
        /*!
         *  \param i The index of the symbol.
         *  \return  Pair (rank(wt[i],i),wt[i])
         *  \par Precondition
         *       \f$ i < size() \f$
         */
        std::pair<size_type, value_type>
        inverse_select(size_type i)const {
            assert(i < size());

            value_type c = 0;
            size_type node_size = m_size, offset = 0;
            for (uint32_t k=0; k < m_max_depth; ++k) {
                size_type ones_before_o   = m_tree_rank(offset);
                size_type ones_before_i   = m_tree_rank(offset + i) - ones_before_o;
                size_type ones_before_end = m_tree_rank(offset + node_size) - ones_before_o;
                c<<=1;
                if (m_tree[offset+i]) { // go to the right child
                    offset += (node_size - ones_before_end);
                    node_size = ones_before_end;
                    i = ones_before_i;
                    c|=1;
                } else { // go to the left child
                    node_size = (node_size - ones_before_end);
                    i = (i-ones_before_i);
                }
                offset += m_size;
            }
            return std::make_pair(i,c);
        }

        //! Calculates the i-th occurrence of the symbol c in the supported vector.
        /*!
         *  \param i The i-th occurrence.
         *  \param c The symbol c.
         *  \par Time complexity
         *       \f$ \Order{\log |\Sigma|} \f$
         *  \par Precondition
         *       \f$ 1 \leq i \leq rank(size(), c) \f$
         */
        size_type select(size_type i, value_type c)const {
            assert(1 <= i and i <= rank(size(), c));
            // possible optimization: if the array is a permutation we can start at the bottom of the tree
            size_type offset = 0;
            uint64_t mask    = (1ULL) << (m_max_depth-1);
            size_type node_size = m_size;
            m_path_off[0] = m_path_rank_off[0] = 0;

            for (uint32_t k=0; k < m_max_depth and node_size; ++k) {
                size_type ones_before_o   = m_tree_rank(offset);
                m_path_rank_off[k] = ones_before_o;
                size_type ones_before_end = m_tree_rank(offset + node_size) - ones_before_o;
                if (c & mask) { // search for a one at this level
                    offset += (node_size - ones_before_end);
                    node_size = ones_before_end;
                } else { // search for a zero at this level
                    node_size = (node_size - ones_before_end);
                }
                offset += m_size;
                m_path_off[k+1] = offset;
                mask >>= 1;
            }
            if (node_size < i) {
                throw std::logic_error("select("+util::to_string(i)+","+util::to_string(c)+"): c does not occur i times in the WT");
                return m_size;
            }
            mask = 1ULL;
            for (uint32_t k=m_max_depth; k>0; --k) {
                offset = m_path_off[k-1];
                size_type ones_before_o = m_path_rank_off[k-1];
                if (c & mask) { // right child => search i'th
                    i = m_tree_select1(ones_before_o + i) - offset + 1;
                } else { // left child => search i'th zero
                    i = m_tree_select0(offset - ones_before_o + i) - offset + 1;
                }
                mask <<= 1;
            }
            return i-1;
        };


        //! For each symbol c in wt[i..j-1] get rank(i,c) and rank(j,c).
        /*!
         * \param i        The start index (inclusive) of the interval.
         * \param j        The end index (exclusive) of the interval.
         * \param k        Reference for number of different symbols in [i..j-1].
         * \param cs       Reference to a vector that will contain in
         *                 cs[0..k-1] all symbols that occur in [i..j-1] in
         *                 ascending order.
         * \param rank_c_i Reference to a vector which equals
         *                 rank_c_i[p] = rank(i,cs[p]), for \f$ 0 \leq p < k \f$.
         * \param rank_c_j Reference to a vector which equals
         *                 rank_c_j[p] = rank(j,cs[p]), for \f$ 0 \leq p < k \f$.
         * \par Time complexity
         *      \f$ \Order{\min{\sigma, k \log \sigma}} \f$
         *
         * \par Precondition
         *      \f$ i \leq j \leq size() \f$
         *      \f$ cs.size() \geq \sigma \f$
         *      \f$ rank_{c_i}.size() \geq \sigma \f$
         *      \f$ rank_{c_j}.size() \geq \sigma \f$
         */
        void interval_symbols(size_type i, size_type j, size_type& k,
                              std::vector<value_type>& cs,
                              std::vector<size_type>& rank_c_i,
                              std::vector<size_type>& rank_c_j) const {
            assert(i <= j and j <= size());
            k=0;
            if (i==j) {
                return;
            }
            if ((i+1)==j) {
                auto res = inverse_select(i);
                cs[0]=res.second;
                rank_c_i[0]=res.first;
                rank_c_j[0]=res.first+1;
                k=1;
                return;
            }

            _interval_symbols(i, j, k, cs, rank_c_i, rank_c_j, 0, 0, m_size, 0);

        }

        //! How many symbols are lexicographic smaller/greater than c in [i..j-1].
        /*!
         * \param i       Start index (inclusive) of the interval.
         * \param j       End index (exclusive) of the interval.
         * \param c       Symbol c.
         * \return A triple containing:
         *         * rank(c,i)
         *         * #symbols smaller than c in [i..j-1]
         *         * #symbols greater than c in [i..j-1]
         *
         * \par Precondition
         *      \f$ i \leq j \leq size() \f$
         */
        template<class t_ret_type = std::tuple<size_type, size_type, size_type>>
        t_ret_type lex_count(size_type i, size_type j, value_type c)const {
            assert(i <= j and j <= size());
            if (((1ULL)<<(m_max_depth))<=c) { // c is greater than any symbol in wt
                return t_ret_type {0, j-i, 0};
            }
            size_type offset  = 0;
            size_type smaller = 0;
            size_type greater = 0;
            uint64_t mask     = (1ULL) << (m_max_depth-1);
            size_type node_size = m_size;
            for (uint32_t k=0; k < m_max_depth; ++k) {
                size_type ones_before_o   = m_tree_rank(offset);
                size_type ones_before_i   = m_tree_rank(offset + i) - ones_before_o;
                size_type ones_before_j   = m_tree_rank(offset + j) - ones_before_o;
                size_type ones_before_end = m_tree_rank(offset + node_size) - ones_before_o;
                if (c & mask) { // search for a one at this level
                    offset += (node_size - ones_before_end);
                    node_size = ones_before_end;
                    smaller += j-i-ones_before_j+ones_before_i;
                    i = ones_before_i;
                    j = ones_before_j;
                } else { // search for a zero at this level
                    node_size -= ones_before_end;
                    greater += ones_before_j-ones_before_i;
                    i -= ones_before_i;
                    j -= ones_before_j;
                }
                offset += m_size;
                mask >>= 1;
            }
            return t_ret_type {i, smaller, greater};
        };

        //! How many symbols are lexicographic smaller than c in [0..i-1].
        /*!
         * \param i Exclusive right bound of the range.
         * \param c Symbol c.
         * \return A tuple containing:
         *         * rank(c,i)
         *         * #symbols smaller than c in [0..i-1]
         * \par Precondition
         *      \f$ i \leq size() \f$
         */
        template<class t_ret_type = std::tuple<size_type, size_type>>
        t_ret_type lex_smaller_count(size_type i, value_type c) const {
            assert(i <= size());
            if (((1ULL)<<(m_max_depth))<=c) { // c is greater than any symbol in wt
                return t_ret_type {0, i};
            }
            size_type offset = 0;
            size_type result = 0;
            uint64_t mask    = (1ULL) << (m_max_depth-1);
            size_type node_size = m_size;
            for (uint32_t k=0; k < m_max_depth and i; ++k) {
                size_type ones_before_o   = m_tree_rank(offset);
                size_type ones_before_i   = m_tree_rank(offset + i) - ones_before_o;
                size_type ones_before_end = m_tree_rank(offset + node_size) - ones_before_o;
                if (c & mask) { // search for a one at this level
                    offset   += (node_size - ones_before_end);
                    node_size = ones_before_end;
                    result   += i - ones_before_i;
                    i         = ones_before_i;
                } else { // search for a zero at this level
                    node_size = (node_size - ones_before_end);
                    i        -= ones_before_i;
                }
                offset += m_size;
                mask >>= 1;
            }
            return t_ret_type {i, result};
        }

        //! range_search_2d searches points in the index interval [lb..rb] and value interval [vlb..vrb].
        /*! \param lb     Left bound of index interval (inclusive)
         *  \param rb     Right bound of index interval (inclusive)
         *  \param vlb    Left bound of value interval (inclusive)
         *  \param vrb    Right bound of value interval (inclusive)
         *  \param report Should the matching points be returned?
         *  \return Pair (#of found points, vector of points), the vector is empty when
         *          report = false.
         */

        std::pair<size_type, std::vector<std::pair<value_type, size_type>>>
        range_search_2d(size_type lb, size_type rb, value_type vlb, value_type vrb,
                        bool report=true) const {
            size_type offsets[m_max_depth+1];
            size_type ones_before_os[m_max_depth+1];
            offsets[0] = 0;
            if (vrb > (1ULL << m_max_depth))
                vrb = (1ULL << m_max_depth);
            if (vlb > vrb)
                return make_pair(0, point_vec_type());
            size_type cnt_answers = 0;
            point_vec_type point_vec;
            _range_search_2d(lb, rb, vlb, vrb, 0, 0, m_size, offsets, ones_before_os, 0, point_vec, report, cnt_answers);
            return make_pair(cnt_answers, point_vec);
        }

        // add parameter path
        // ilb interval left bound
        // irb interval right bound
        void
        _range_search_2d(size_type lb, size_type rb, value_type vlb, value_type vrb, size_type depth,
                         size_type ilb, size_type node_size, size_type offsets[],
                         size_type ones_before_os[], size_type path,
                         point_vec_type& point_vec, bool report, size_type& cnt_answers)
        const {
            if (lb > rb)
                return;
            if (depth == m_max_depth) {
                if (report) {
                    for (size_type j=1; j <= node_size; ++j) {
                        size_type i = j;
                        size_type c = path;
                        for (uint32_t k=m_max_depth; k>0; --k) {
                            size_type offset = offsets[k-1];
                            size_type ones_before_o = ones_before_os[k-1];
                            if (c&1) {
                                i = m_tree_select1(ones_before_o + i) - offset + 1;
                            } else {
                                i = m_tree_select0(offset - ones_before_o + i) - offset + 1;
                            }
                            c >>= 1;
                        }
                        point_vec.emplace_back(i-1, path); // add resulting index; -1 cause of 0 based indexing
                    }
                }
                cnt_answers += node_size;
                return;
            }
            size_type irb = ilb + (1ULL << (m_max_depth-depth));
            size_type mid = (irb + ilb)>>1;

            size_type offset = offsets[depth];

            size_type ones_before_o    = m_tree_rank(offset);
            ones_before_os[depth]      = ones_before_o;
            size_type ones_before_lb   = m_tree_rank(offset + lb);
            size_type ones_before_rb   = m_tree_rank(offset + rb + 1);
            size_type ones_before_end  = m_tree_rank(offset + node_size);
            size_type zeros_before_o   = offset - ones_before_o;
            size_type zeros_before_lb  = offset + lb - ones_before_lb;
            size_type zeros_before_rb  = offset + rb + 1 - ones_before_rb;
            size_type zeros_before_end = offset + node_size - ones_before_end;
            if (vlb < mid and mid) {
                size_type nlb    = zeros_before_lb - zeros_before_o;
                size_type nrb    = zeros_before_rb - zeros_before_o;
                offsets[depth+1] = offset + m_size;
                if (nrb)
                    _range_search_2d(nlb, nrb-1, vlb, std::min(vrb,mid-1), depth+1, ilb, zeros_before_end - zeros_before_o, offsets, ones_before_os, path<<1, point_vec, report, cnt_answers);
            }
            if (vrb >= mid) {
                size_type nlb     = ones_before_lb - ones_before_o;
                size_type nrb     = ones_before_rb - ones_before_o;
                offsets[depth+1]  = offset + m_size + (zeros_before_end - zeros_before_o);
                if (nrb)
                    _range_search_2d(nlb, nrb-1, std::max(mid, vlb), vrb, depth+1, mid, ones_before_end - ones_before_o, offsets, ones_before_os, (path<<1)+1 , point_vec, report, cnt_answers);
            }
        }

        //! Returns a const_iterator to the first element.
        const_iterator begin()const {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        const_iterator end()const {
            return const_iterator(this, size());
        }


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
            written_bytes += write_member(m_max_depth, out, child, "max_depth");
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
            read_member(m_max_depth, in);
            init_buffers(m_max_depth);
        }

        //! Returns the element in T[lb..rb] with rank quantile.
        /*!
         *  \param lb left array bound in T
         *  \param rb right array bound in T
         *  \param quantile 'quantile' smallest symbol (starts with 0)
         */
        std::pair<value_type,size_type>
        quantile_freq(size_type lb, size_type rb,size_type quantile) const {

            size_type  offset = 0;
            value_type sym = 0;
            size_type freq = 0;
            size_type node_size = m_size;

            for (size_t k=0; k < m_max_depth; ++k) {
                sym <<= 1;

                /* number of 1s before the level offset and after the node */
                size_type ones_before_offset   = m_tree_rank(offset);
                size_type ones_before_end = m_tree_rank(offset + node_size) - ones_before_offset;

                /* number of 1s before T[l..r] */
                size_type rank_before_left = 0;
                if (offset+lb>0) rank_before_left = m_tree_rank(offset + lb);

                /* number of 1s before T[r] */
                size_type rank_before_right   = m_tree_rank(offset + rb + 1);

                /* number of 1s in T[l..r] */
                size_type num_ones = rank_before_right - rank_before_left;
                /* number of 0s in T[l..r] */
                size_type num_zeros = (rb-lb+1) - num_ones;

                /* if there are more than q 0s we go right. left otherwise */
                if (quantile >= num_zeros) { /* go right */
                    freq = num_ones; /* calc freq */
                    /* set bit to 1 in sym */
                    sym |= 1;
                    /* number of 1s before T[l..r] within the current node */
                    lb = rank_before_left - ones_before_offset;
                    /* number of 1s in T[l..r] */
                    rb = lb + num_ones - 1;
                    quantile = quantile - num_zeros;

                    /* calc starting pos of right childnode */
                    offset += (node_size - ones_before_end);
                    node_size = ones_before_end;

                } else { /* go left q = q // sym == sym */
                    freq = num_zeros; /* calc freq */
                    /* number of zeros before T[l..r] within the current node */
                    lb = lb - (rank_before_left - ones_before_offset);
                    /* number of zeros in T[l..r] + left bound */
                    rb = lb + num_zeros - 1;

                    /* calc end pos of left childnode */
                    node_size = (node_size - ones_before_end);
                }
                offset += m_size; /* next level */
            }
            return {sym,freq};
        };


        //! Returns the top k most frequent documents in T[lb..rb]
        /*!
         *  \param lb left array bound in T
         *  \param rb right array bound in T
         *  \param k the number of documents to return
         *  \returns the top-k items in descending order.
         *  \par Time complexity
         *      \f$ \Order{\log |\Sigma|} \f$
         */
        class topk_greedy_range_t
        {
            public:
                bool operator<(const topk_greedy_range_t& r) const {
                    if ((rb-lb+1) != (r.rb-r.lb+1))
                        return ((rb-lb+1) < (r.rb-r.lb+1));
                    return sym > r.sym;
                }
            public:
                value_type sym = 0;
                size_type id = 0;
                size_type lb = 0;
                size_type rb = 0;
                size_type offset = 0;
                size_type node_size = 0;
                size_type level = 0;
                size_type freq = 0;
        };

        //! Returns the top k most frequent documents in T[lb..rb]
        /*!
         *  \param lb left array bound in T
         *  \param rb right array bound in T
         *  \param k the number of documents to return
         *  \returns the top-k items in ascending order.
         */
        std::vector< std::pair<value_type,size_type> >
        topk_greedy(size_type lb, size_type rb,size_type k) const {

            std::vector< std::pair<value_type,size_type> > results;
            std::priority_queue<topk_greedy_range_t> heap;

            /* add the initial range */
            topk_greedy_range_t ir;
            ir.id = 1;
            ir.node_size = m_size;
            ir.level = 0;
            ir.lb = lb;
            ir.rb = rb;
            heap.push(ir);

            while (! heap.empty()) {
                topk_greedy_range_t r = heap.top(); heap.pop();
                //m_greedy_nodes.insert(r.id);
                if (r.level == m_max_depth) { /* leaf node */
                    results.emplace_back(r.sym,r.freq);
                    if (results.size()==k) {
                        /* we got the top-k */
                        break;
                    }
                    continue; /* we processed this range */
                }

                /* number of 1s before the level offset and after the node */
                size_type ones_before_offset = m_tree_rank(r.offset);
                size_type ones_before_end = m_tree_rank(r.offset + r.node_size) - ones_before_offset;

                /* number of 1s before T[l..r] */
                size_type rank_before_left = m_tree_rank(r.offset + r.lb);

                /* number of 1s before T[r] */
                size_type rank_before_right   = m_tree_rank(r.offset + r.rb + 1);

                /* number of 1s in T[l..r] */
                size_type num_ones = rank_before_right - rank_before_left;
                /* number of 0s in T[l..r] */
                size_type num_zeros = (r.rb-r.lb+1) - num_ones;

                if (num_ones) { /* map to right child */
                    topk_greedy_range_t nr;
                    nr.sym = (r.sym<<1)|1;
                    nr.id = r.id*2+1;
                    nr.freq = num_ones;
                    nr.offset = m_size + r.offset + (r.node_size - ones_before_end);
                    nr.node_size = ones_before_end;
                    nr.level = r.level + 1;
                    /* number of 1s before T[l..r] within the current node */
                    nr.lb = rank_before_left - ones_before_offset;
                    /* number of 1s in T[l..r] */
                    nr.rb = nr.lb + num_ones - 1;

                    heap.push(nr);
                }
                if (num_zeros) { /* map to left child */
                    topk_greedy_range_t nr;
                    nr.id = r.id*2;
                    nr.sym = r.sym<<1;
                    nr.freq = num_zeros;
                    nr.offset = m_size + r.offset;
                    nr.node_size = (r.node_size - ones_before_end);
                    nr.level = r.level + 1;
                    /* number of 1s before T[l..r] within the current node */
                    nr.lb = r.lb - (rank_before_left - ones_before_offset);
                    /* number of 1s in T[l..r] */
                    nr.rb = nr.lb + num_zeros - 1;
                    heap.push(nr);
                }
            }
            return results;
        };

        //! Returns the top k most frequent documents in T[lb..rb]
        /*!
         *  \param lb left array bound in T
         *  \param rb right array bound in T
         *  \param k the number of documents to return
         *  \returns the top-k items in ascending order.
         */
        std::vector< std::pair<value_type,size_type> >
        topk_qprobing(size_type lb, size_type rb,size_type k) const {
            using p_t = std::pair<value_type,size_type>;
            std::vector<p_t> results;
            auto comp = [](p_t& a,p_t& b) { return a.second > b.second; };
            std::priority_queue<p_t,std::vector<p_t>,decltype(comp)> heap(comp);
            bit_vector seen(1 << m_max_depth); // TODO: better idea?

            /* we start probing using the largest power smaller than len */
            size_type len = rb-lb+1;
            size_type power2greaterlen = 1 << (bits::hi(len)+1);
            size_type probe_interval = power2greaterlen >> 1;

            /* we probe the smallest elem (pos 0 in sorted array) only once */
            auto qf = quantile_freq(lb,rb,0);
            heap.push(qf);
            seen[qf.first] = 1;

            qf = quantile_freq(lb,rb,probe_interval);
            if (!seen[qf.first]) heap.push(qf);
            seen[qf.first] = 1;

            while (probe_interval > 1) {
                size_type probe_pos = probe_interval >> 1;
                while (probe_pos < len) {
                    qf = quantile_freq(lb,rb,probe_pos);
                    if (!seen[qf.first]) { /* not in heap */
                        if (heap.size()<k) {
                            heap.push(qf);
                            seen[qf.first] = 1;
                        } else {
                            /* throw out the smallest and add the new one */
                            if (heap.top().second < qf.second) {
                                heap.pop();
                                heap.push(qf);
                                seen[qf.first] = 1;
                            }
                        }
                    }
                    probe_pos += probe_interval;
                }
                probe_interval >>= 1;
                /* we have enough or can't find anything better */
                if (heap.size() == k && probe_interval-1 <= heap.top().second) break;
            }
            /* populate results */
            while (!heap.empty())  {
                results.emplace(results.begin() , heap.top());
                heap.pop();
            }
            return results;
        };

        //! Returns the intersection of T[lb1..rb1],T[lb2..rb2]...T[lbm..rbm]
        /*!
         *  \param ranges the sp,ep ranges to intersect
         *  \param the threshold t of how many ranges have to be at least
         *         still be present at the leaf level.
         *  \param the results are stored in the results parameter.
         */
        struct intersect_range_t {
            using p_t = std::pair<size_type,size_type>;

            intersect_range_t() {}
            intersect_range_t(size_type off,size_type ns, size_type lvl,
                              value_type _sym, std::vector<p_t>& r)
                :  offset(off) , node_size(ns) , level(lvl), sym(_sym), ranges(r)
            {}
            intersect_range_t(size_type off,size_type ns,size_type lvl,value_type _sym)
                : offset(off) , node_size(ns) , level(lvl), sym(_sym)  {}

            intersect_range_t(const intersect_range_t& r)
                : offset(r.offset), node_size(r.node_size), level(r.level),
                  sym(r.sym), ranges(r.ranges) {}

            intersect_range_t(intersect_range_t&& r)
                : offset(r.offset), node_size(r.node_size), level(r.level),
                  sym(r.sym), ranges(std::move(r.ranges)) {}

            intersect_range_t& operator=(const intersect_range_t& r) {
                offset    = r.offset;
                node_size = r.node_size;
                level     = r.level;
                sym       = r.sym;
                ranges    = r.ranges;
                return *this;
            }

            size_type offset = 0;
            size_type node_size = 0;
            size_type level = 0;
            value_type sym = 0;
            std::vector<p_t> ranges;
        };


        std::vector< std::pair<value_type,size_type> >
        intersect(std::vector< std::pair<size_type,size_type> >& ranges, size_type threshold=0) const {
            using p_t = std::pair<value_type,size_type>;
            std::vector<p_t> results;

            if (threshold==0) { /* default: all ranges must be present */
                threshold = ranges.size();
            }

            std::vector<intersect_range_t> intervals;
            size_type n = m_size;
            intervals.emplace_back(intersect_range_t(0,n,0,0,ranges));

            while (!intervals.empty()) {
                intersect_range_t cr = intervals[intervals.size()-1]; intervals.pop_back();

                if (cr.level == m_max_depth) {
                    if (threshold <= cr.ranges.size()) {
                        /* we found a symbol  */
                        size_type freq = 0;
                        for (auto& r : cr.ranges) freq += (r.second - r.first + 1);
                        results.emplace_back(cr.sym,freq);
                    }
                } else {
                    /* map each range to the corresponding range at level + 1 */

                    /* number of 1s before the level offset and after the node */
                    size_type ones_before_offset = m_tree_rank(cr.offset);
                    size_type ones_before_end = m_tree_rank(cr.offset + cr.node_size) - ones_before_offset;

                    size_type offset_zero = m_size + cr.offset;
                    size_type offset_one = m_size + cr.offset + (cr.node_size - ones_before_end);
                    size_type node_size_zero = m_size + cr.offset;
                    size_type node_size_one = m_size + cr.offset + (cr.node_size - ones_before_end);

                    intersect_range_t range_zero(offset_zero,node_size_zero,cr.level+1,cr.sym<<1);
                    intersect_range_t range_one(offset_one,node_size_one,cr.level+1,(cr.sym<<1)|1);

                    for (size_t i=0; i<cr.ranges.size(); i++) {
                        std::pair<size_t,size_t> r = cr.ranges[i];
                        size_type lb = r.first, rb = r.second;
                        /* number of 1s before T[l..r] */
                        size_type rank_before_left = 0;
                        if (cr.offset+lb>0) rank_before_left = m_tree_rank(cr.offset + lb);
                        /* number of 1s before T[r] */
                        size_type rank_before_right = m_tree_rank(cr.offset + rb + 1);
                        /* number of 1s in T[l..r] */
                        size_type num_ones = rank_before_right - rank_before_left;
                        /* number of 0s in T[l..r] */
                        size_type num_zeros = (rb-lb+1) - num_ones;

                        if (num_ones) { /* map to right child */
                            /* number of 1s before T[l..r] within the current node */
                            size_type lb_one = rank_before_left - ones_before_offset;
                            /* number of 1s in T[l..r] */
                            size_type rb_one = lb_one + num_ones - 1;

                            range_one.ranges.emplace_back(lb_one,rb_one);
                        }
                        if (num_zeros) { /* map to left child */
                            /* number of 1s before T[l..r] within the current node */
                            size_type lb_zero = lb - (rank_before_left - ones_before_offset);
                            /* number of 1s in T[l..r] */
                            size_type rb_zero = lb_zero + num_zeros - 1;

                            range_zero.ranges.emplace_back(lb_zero,rb_zero);
                        }
                    }
                    if (range_zero.ranges.size() >= threshold) intervals.emplace_back(range_zero);
                    if (range_one.ranges.size() >= threshold) intervals.emplace_back(range_one);
                }
            }
            return results;
        }



        //! Returns the top k most relevant documents from ranges T[x_0,y_0]
        //!  ... T[x_1,y_1] ... T[x_m,y_m]
        /*!
         *  \param k the number of documents to return
         *  \returns the top-k items in descending order.
         */
        struct greedy_tfidf_range_t {
            size_type sp;
            size_type ep;
            size_type f_t;
            double w_qt;
            greedy_tfidf_range_t(const greedy_tfidf_range_t& r) : sp(r.sp) , ep(r.ep) , f_t(r.f_t) , w_qt(r.w_qt) {}
            greedy_tfidf_range_t(size_type s,size_type e,size_type nd,size_type N) : sp(s) , ep(e) , f_t(nd) {
                w_qt = log(1+(N/(double)f_t));
            }
            double w_dt() { // estimate w_dt
                return 1 + log((double)(ep-sp+1));
            }
        };

        std::vector< std::pair<value_type,double> >
        topk_greedy_similarity(std::vector<greedy_tfidf_range_t>& input_ranges,size_type k) const {
            std::vector< std::pair<value_type,double> > results;

            struct topk_greedy_tfidf {
                std::vector<greedy_tfidf_range_t> ranges;
                size_type id = 1;
                size_type offset = 0;
                size_type node_size = 0;
                size_type level = 0;
                value_type sym = 0;
                double score = 0.0f;
                topk_greedy_tfidf() {};
                topk_greedy_tfidf(size_type off,size_type ns,size_type lvl,value_type _sym)
                    : offset(off) , node_size(ns) , level(lvl), sym(_sym)  {}
                bool operator < (const topk_greedy_tfidf& r) const { return score < r.score; }
            };

            auto cosine_score = [](double score,greedy_tfidf_range_t& r) -> double {
                return score += (r.w_qt * r.w_dt());
            };

            auto print_range = [](const topk_greedy_tfidf& r) {
                std::cout << "=================================================" << std::endl;
                std::cout << "ranges = ";
                for (const auto& rr : r.ranges) {
                    std::cout << "<" << rr.sp << "," << rr.ep << ">";
                }
                std::cout << std::endl;
                std::cout << "offset = " << r.offset << std::endl;
                std::cout << "node_size = " << r.node_size << std::endl;
                std::cout << "level = " << r.level << std::endl;
                std::cout << "sym = " << r.sym << std::endl;
                std::cout << "score = " << r.score << std::endl;
                std::cout << "=================================================" << std::endl;
            };

            std::priority_queue<topk_greedy_tfidf> heap;

            /* add the initial ranges */
            topk_greedy_tfidf ir;
            ir.node_size = m_size;
            ir.level = 0;
            ir.sym = 0;
            ir.offset = 0;
            ir.id = 1;
            for (const auto& r : input_ranges) {
                ir.ranges.emplace_back(r);
            }
            ir.score = std::accumulate(ir.ranges.begin(),ir.ranges.end(),0.0f, cosine_score);
            heap.push(ir);


            while (! heap.empty()) {
                auto current = heap.top(); heap.pop();
                //m_greedy_sim_nodes.insert(current.id);
                //print_range(current);

                if (current.level == m_max_depth) { /* leaf node */
                    results.emplace_back(current.sym,current.score);
                    if (results.size()==k) {
                        /* we got the top-k */
                        break;
                    }
                    continue; /* we processed this range */
                }
                /* map each range to the corresponding range at level + 1 */

                /* number of 1s before the level offset and after the node */
                size_type ones_before_offset = m_tree_rank(current.offset);
                size_type ones_before_end = m_tree_rank(current.offset + current.node_size) - ones_before_offset;

                size_type offset_zero = m_size + current.offset;
                size_type offset_one = m_size + current.offset + (current.node_size - ones_before_end);
                size_type node_size_zero = current.node_size - ones_before_end;
                size_type node_size_one = ones_before_end;

                topk_greedy_tfidf ranges_zero(offset_zero,node_size_zero,current.level+1,current.sym<<1);
                topk_greedy_tfidf ranges_one(offset_one,node_size_one,current.level+1,(current.sym<<1)|1);

                for (size_t i=0; i<current.ranges.size(); i++) {
                    auto& r = current.ranges[i];
                    size_type lb = r.sp, rb = r.ep;
                    /* number of 1s before T[l..r] */
                    size_type rank_before_left = 0;
                    if (current.offset+lb>0) rank_before_left = m_tree_rank(current.offset + lb);
                    /* number of 1s before T[r] */
                    size_type rank_before_right = m_tree_rank(current.offset + rb + 1);
                    /* number of 1s in T[l..r] */
                    size_type num_ones = rank_before_right - rank_before_left;
                    /* number of 0s in T[l..r] */
                    size_type num_zeros = (rb-lb+1) - num_ones;

                    if (num_ones) { /* map to right child */
                        /* number of 1s before T[l..r] within the current node */
                        greedy_tfidf_range_t r_one = r;
                        r_one.sp = rank_before_left - ones_before_offset;
                        /* number of 1s in T[l..r] */
                        r_one.ep = r_one.sp + num_ones - 1;
                        ranges_one.ranges.push_back(r_one);
                    }
                    if (num_zeros) { /* map to left child */
                        /* number of 1s before T[l..r] within the current node */
                        greedy_tfidf_range_t r_zero = r;
                        r_zero.sp = lb - (rank_before_left - ones_before_offset);
                        /* number of 1s in T[l..r] */
                        r_zero.ep = r_zero.sp + num_zeros - 1;

                        ranges_zero.ranges.push_back(r_zero);
                    }
                }
                if (!ranges_zero.ranges.empty()) {
                    ranges_zero.score = std::accumulate(ranges_zero.ranges.begin(),ranges_zero.ranges.end(),0.0f, cosine_score);
                    ranges_zero.id = current.id*2;
                    heap.push(ranges_zero);
                }
                if (!ranges_one.ranges.empty()) {
                    ranges_one.score = std::accumulate(ranges_one.ranges.begin(),ranges_one.ranges.end(),0.0f, cosine_score);
                    ranges_one.id = current.id*2+1;
                    heap.push(ranges_one);
                }
            }
            return results;
        };
};

}// end namespace sdsl
#endif
