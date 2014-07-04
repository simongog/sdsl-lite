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
           sequence of the numbers. This wavelet tree class takes
           less memory than the wt_pc class for large alphabets.
    \author Simon Gog, Shanika Kuruppu
*/
#ifndef INCLUDED_SDSL_INT_WAVELET_TREE
#define INCLUDED_SDSL_INT_WAVELET_TREE

#include "sdsl_concepts.hpp"
#include "int_vector.hpp"
#include "rank_support_v.hpp"
#include "select_support_mcl.hpp"
#include "wt_helper.hpp"
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


    protected:

        size_type              m_size  = 0;
        size_type              m_sigma = 0;    //<- \f$ |\Sigma| \f$
        bit_vector_type        m_tree;         // bit vector to store the wavelet tree
        rank_1_type            m_tree_rank;    // rank support for the wavelet tree bit vector
        select_1_type          m_tree_select1; // select support for the wavelet tree bit vector
        select_0_type          m_tree_select0;
        uint32_t               m_max_level = 0;
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
            m_max_level     = wt.m_max_level;
            m_path_off      = wt.m_path_off;
            m_path_rank_off = wt.m_path_rank_off;
        }

    private:

        void init_buffers(uint32_t max_level) {
            m_path_off = int_vector<64>(max_level+1);
            m_path_rank_off = int_vector<64>(max_level+1);
        }

        // recursive internal version of the method interval_symbols
        void _interval_symbols(size_type i, size_type j, size_type& k,
                               std::vector<value_type>& cs,
                               std::vector<size_type>& rank_c_i,
                               std::vector<size_type>& rank_c_j,
                               size_type level,
                               size_type path,
                               size_type node_size,
                               size_type offset) const {
            // invariant: j>i

            if (level >= m_max_level) {
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
                _interval_symbols(new_i, new_j, k, cs, rank_c_i, rank_c_j, level+1, path<<1, new_node_size, new_offset);
            }

            // goto right child
            if ((ones_before_j-ones_before_i)>0) {
                size_type new_offset = offset+(node_size - ones_before_end) + m_size;
                size_type new_node_size = ones_before_end;
                size_type new_i = ones_before_i;
                size_type new_j = ones_before_j;
                _interval_symbols(new_i, new_j, k, cs, rank_c_i, rank_c_j, level+1, (path<<1)|1, new_node_size, new_offset);
            }
        }

    public:

        const size_type&       sigma = m_sigma;         //!< Effective alphabet size of the wavelet tree.
        const bit_vector_type& tree  = m_tree;          //!< A concatenation of all bit vectors of the wavelet tree.
        const uint32_t&        max_level = m_max_level; //!< Maximal level of the wavelet tree.

        //! Default constructor
        wt_int() {
            init_buffers(m_max_level);
        };

        //! Semi-external constructor
        /*! \param buf         File buffer of the int_vector for which the wt_int should be build.
         *  \param size        Size of the prefix of v, which should be indexed.
         *  \param max_level   Maximal level of the wavelet tree. If set to 0, determined automatically.
         *    \par Time complexity
         *        \f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         *        I.e. we need \Order{n\log n} if rac is a permutation of 0..n-1.
         *    \par Space complexity
         *        \f$ n\log|\Sigma| + O(1)\f$ bits, where \f$n=size\f$.
         */
        template<uint8_t int_width>
        wt_int(int_vector_buffer<int_width>& buf, size_type size,
               uint32_t max_level=0) : m_size(size) {
            init_buffers(m_max_level);
            if (0 == m_size)
                return;
            size_type n = buf.size();  // set n
            if (n < m_size) {
                throw std::logic_error("n="+util::to_string(n)+" < "+util::to_string(m_size)+"=m_size");
                return;
            }
            m_sigma = 0;
            int_vector<int_width> rac(m_size, 0, buf.width());

            value_type x = 1;  // variable for the biggest value in rac
            for (size_type i=0; i < m_size; ++i) {
                if (buf[i] > x)
                    x = buf[i];
                rac[i] = buf[i];
            }

            if (max_level == 0) {
                m_max_level = bits::hi(x)+1; // max_level bits to represent all values range [0..x]
            } else {
                m_max_level = max_level;
            }
            init_buffers(m_max_level);

            // buffer for elements in the right node
            int_vector_buffer<> buf1(tmp_file(buf.filename(), "_wt_constr_buf"),
                                     std::ios::out, 10*(1<<20), buf.width());
            std::string tree_out_buf_file_name = tmp_file(buf.filename(), "_m_tree");
            osfstream tree_out_buf(tree_out_buf_file_name, std::ios::binary|
                                   std::ios::trunc|std::ios::out);

            size_type bit_size = m_size*m_max_level;
            tree_out_buf.write((char*) &bit_size, sizeof(bit_size));// write size of bit_vector

            size_type tree_pos = 0;
            uint64_t tree_word = 0;

            uint64_t mask_old = 1ULL<<(m_max_level);
            for (uint32_t k=0; k<m_max_level; ++k) {
                size_type          start     = 0;
                const uint64_t    mask_new = 1ULL<<(m_max_level-k-1);
                do {
                    size_type i           = start;
                    size_type cnt0        = 0;
                    size_type cnt1        = 0;
                    uint64_t  start_value = (rac[i]&mask_old);
                    uint64_t  x;
                    while (i < m_size and((x=rac[i])&mask_old)==start_value) {
                        if (x&mask_new) {
                            tree_word |= (1ULL << (tree_pos&0x3FULL));
                            buf1[cnt1++] = x;
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
                    if (k+1 < m_max_level) { // inner node
                        for (size_type j=0; j<cnt1; ++j) {
                            rac[start+cnt0+j] = buf1[j];
                        }
                    } else { // leaf node
                        m_sigma += (cnt0>0) + (cnt1>0); // increase sigma for each leaf
                    }
                    start += cnt0+cnt1;
                } while (start < m_size);
                mask_old += mask_new;
            }
            if ((tree_pos & 0x3FULL) != 0) { // if tree_pos % 64 > 0 => there are remaining entries we have to write
                tree_out_buf.write((char*) &tree_word, sizeof(tree_word));
            }
            buf1.close(true); // remove temporary file
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

        //! Copy constructor
        wt_int(wt_int&& wt) {
            *this = std::move(wt);
        }

        //! Assignment operator
        wt_int& operator=(const wt_int& wt) {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Assignment move operator
        wt_int& operator=(wt_int&& wt) {
            if (this != &wt) {
                m_size          = wt.m_size;
                m_sigma         = wt.m_sigma;
                m_tree          = std::move(wt.m_tree);
                m_tree_rank     = std::move(wt.m_tree_rank);
                m_tree_rank.set_vector(&m_tree);
                m_tree_select1  = std::move(wt.m_tree_select1);
                m_tree_select1.set_vector(&m_tree);
                m_tree_select0  = std::move(wt.m_tree_select0);
                m_tree_select0.set_vector(&m_tree);
                m_max_level     = std::move(wt.m_max_level);
                m_path_off      = std::move(wt.m_path_off);
                m_path_rank_off = std::move(wt.m_path_rank_off);
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
                std::swap(m_max_level,  wt.m_max_level);
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
            for (uint32_t k=0; k < m_max_level; ++k) {
                res <<= 1;
                size_type ones_before_o   = m_tree_rank(offset);
                size_type ones_before_i   = m_tree_rank(offset + i) - ones_before_o;
                size_type ones_before_end = m_tree_rank(offset + node_size) - ones_before_o;
                if (m_tree[offset+i]) { // one at position i => follow right child
                    offset += (node_size - ones_before_end);
                    node_size = ones_before_end;
                    i = ones_before_i;
                    res |= 1;
                } else { // zero at position i => follow left child
                    node_size = (node_size - ones_before_end);
                    i = (i-ones_before_i);
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
            if (((1ULL)<<(m_max_level))<=c) { // c is greater than any symbol in wt
                return 0;
            }
            size_type offset = 0;
            uint64_t mask = (1ULL) << (m_max_level-1);
            size_type node_size = m_size;
            for (uint32_t k=0; k < m_max_level and i; ++k) {
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
            for (uint32_t k=0; k < m_max_level; ++k) {
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
            uint64_t mask    = (1ULL) << (m_max_level-1);
            size_type node_size = m_size;
            m_path_off[0] = m_path_rank_off[0] = 0;

            for (uint32_t k=0; k < m_max_level and node_size; ++k) {
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
            if (0ULL == node_size or node_size < i) {
                throw std::logic_error("select("+util::to_string(i)+","+util::to_string(c)+"): c does not occur i times in the WT");
                return m_size;
            }
            mask = 1ULL;
            for (uint32_t k=m_max_level; k>0; --k) {
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
         *         * rank(i,c)
         *         * #symbols smaller than c in [i..j-1]
         *         * #symbols greater than c in [i..j-1]
         *
         * \par Precondition
         *      \f$ i \leq j \leq size() \f$
         */
        template<class t_ret_type = std::tuple<size_type, size_type, size_type>>
        t_ret_type lex_count(size_type i, size_type j, value_type c)const {
            assert(i <= j and j <= size());
            if (((1ULL)<<(m_max_level))<=c) { // c is greater than any symbol in wt
                return t_ret_type {0, j-i, 0};
            }
            size_type offset  = 0;
            size_type smaller = 0;
            size_type greater = 0;
            uint64_t mask     = (1ULL) << (m_max_level-1);
            size_type node_size = m_size;
            for (uint32_t k=0; k < m_max_level; ++k) {
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
         *         * rank(i,c)
         *         * #symbols smaller than c in [0..i-1]
         * \par Precondition
         *      \f$ i \leq size() \f$
         */
        template<class t_ret_type = std::tuple<size_type, size_type>>
        t_ret_type lex_smaller_count(size_type i, value_type c) const {
            assert(i <= size());
            if (((1ULL)<<(m_max_level))<=c) { // c is greater than any symbol in wt
                return t_ret_type {0, i};
            }
            size_type offset = 0;
            size_type result = 0;
            uint64_t mask    = (1ULL) << (m_max_level-1);
            size_type node_size = m_size;
            for (uint32_t k=0; k < m_max_level and i; ++k) {
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
            size_type offsets[m_max_level+1];
            size_type ones_before_os[m_max_level+1];
            offsets[0] = 0;
            if (vrb > (1ULL << m_max_level))
                vrb = (1ULL << m_max_level);
            if (vlb > vrb)
                return make_pair(0, point_vec_type());
            size_type cnt_answers = 0;
            point_vec_type point_vec;
            _range_search_2d(lb, rb, vlb, vrb, 0, 0, m_size, offsets, ones_before_os, 0, point_vec, report, cnt_answers);
            return make_pair(cnt_answers, point_vec);
        }

        void
        _range_search_2d(size_type lb, size_type rb, value_type vlb, value_type vrb, size_type level,
                         size_type ilb, size_type node_size, size_type offsets[],
                         size_type ones_before_os[], size_type path,
                         point_vec_type& point_vec, bool report, size_type& cnt_answers)
        const {
            if (lb > rb)
                return;
            if (level == m_max_level) {
                if (report) {
                    for (size_type j=lb+1; j <= rb+1; ++j) {
                        size_type i = j;
                        size_type c = path;
                        for (uint32_t k=m_max_level; k>0; --k) {
                            size_type offset = offsets[k-1];
                            size_type ones_before_o = ones_before_os[k-1];
                            if (c&1) {
                                i = m_tree_select1(ones_before_o + i) - offset + 1;
                            } else {
                                i = m_tree_select0(offset - ones_before_o + i) - offset + 1;
                            }
                            c >>= 1;
                        }
                        point_vec.emplace_back(i-1, path);
                    }
                }
                cnt_answers += rb-lb+1;
                return;
            }
            size_type irb = ilb + (1ULL << (m_max_level-level));
            size_type mid = (irb + ilb)>>1;

            size_type offset = offsets[level];

            size_type ones_before_o    = m_tree_rank(offset);
            ones_before_os[level]      = ones_before_o;
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
                offsets[level+1] = offset + m_size;
                if (nrb)
                    _range_search_2d(nlb, nrb-1, vlb, std::min(vrb,mid-1), level+1, ilb, zeros_before_end - zeros_before_o, offsets, ones_before_os, path<<1, point_vec, report, cnt_answers);
            }
            if (vrb >= mid) {
                size_type nlb     = ones_before_lb - ones_before_o;
                size_type nrb     = ones_before_rb - ones_before_o;
                offsets[level+1]  = offset + m_size + (zeros_before_end - zeros_before_o);
                if (nrb)
                    _range_search_2d(nlb, nrb-1, std::max(mid, vlb), vrb, level+1, mid, ones_before_end - ones_before_o, offsets, ones_before_os, (path<<1)+1 , point_vec, report, cnt_answers);
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
            written_bytes += write_member(m_max_level, out, child, "max_level");
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
            read_member(m_max_level, in);
            init_buffers(m_max_level);
        }

        //! Represents a node in the wavelet tree
        struct node_type {
            size_type  offset   = 0;
            size_type  size     = 0;
            size_type  level    = 0;
            value_type sym      = 0;

            // Default constructor
            node_type(size_type o=0, size_type sz=0, size_type l=0,
                      value_type sy=0) :
                offset(o), size(sz), level(l), sym(sy) {}

            // Copy constructor
            node_type(const node_type&) = default;

            // Move copy constructor
            node_type(node_type&&) = default;

            // Assignment operator
            node_type& operator=(const node_type&) = default;

            // Move assignment operator
            node_type& operator=(node_type&&) = default;

            // Comparator operator
            bool operator==(const node_type& v) const {
                return offset == v.offset;
            }

            // Smaller operator
            bool operator<(const node_type& v) const {
                return offset < v.offset;
            }

            // Greater operator
            bool operator>(const node_type& v) const {
                return offset > v.offset;
            }
        };

        //! Checks if the node is a leaf node
        bool is_leaf(const node_type& v) const {
            return v.level == m_max_level;
        }

        value_type sym(const node_type& v) const {
            return v.sym;
        }

        bool empty(const node_type& v) const {
            return v.size == (size_type)0;
        }

        //! Return the root node
        node_type root() const {
            return node_type(0, m_size, 0, 0);
        }

        //! Returns the two child nodes of an inner node
        /*! \param v An inner node of a wavelet tree.
         *  \return Return a pair of nodes (left child, right child).
         *  \pre !is_leaf(v)
         */
        std::pair<node_type, node_type>
        expand(const node_type& v) const {
            node_type v_right = v;
            return expand(std::move(v_right));
        }

        //! Returns the two child nodes of an inner node
        /*! \param v An inner node of a wavelet tree.
         *  \return Return a pair of nodes (left child, right child).
         *  \pre !is_leaf(v)
         */
        std::pair<node_type, node_type>
        expand(node_type&& v) const {
            node_type v_left;
            size_type offset_rank = m_tree_rank(v.offset);
            size_type ones        = m_tree_rank(v.offset + v.size) - offset_rank;

            v_left.offset = v.offset + m_size;
            v_left.size   = v.size - ones;
            v_left.level  = v.level + 1;
            v_left.sym    = v.sym<<1;

            v.offset = v.offset + m_size + v_left.size;
            v.size   = ones;
            v.level  = v.level + 1;
            v.sym    = (v.sym<<1)|1;

            return std::make_pair(std::move(v_left), v);
        }

        //! Returns for each range its left and right child ranges
        /*! \param v      An inner node of an wavelet tree.
         *  \param ranges A vector of ranges. Each range [s,e]
         *                has to be contained in v=[v_s,v_e].
         *  \return A vector a range pairs. The first element of each
         *          range pair correspond to the original range
         *          mapped to the left child of v; the second element to the
         *          range mapped to the right child of v.
         *  \pre !is_leaf(v) and s>=v_s and e<=v_e
         */
        std::pair<range_vec_type, range_vec_type>
        expand(const node_type& v,
               const range_vec_type& ranges) const {
            auto ranges_copy = ranges;
            return expand(v, std::move(ranges_copy));
        }

        //! Returns for each range its left and right child ranges
        /*! \param v      An inner node of an wavelet tree.
         *  \param ranges A vector of ranges. Each range [s,e]
         *                has to be contained in v=[v_s,v_e].
         *  \return A vector a range pairs. The first element of each
         *          range pair correspond to the original range
         *          mapped to the left child of v; the second element to the
         *          range mapped to the right child of v.
         *  \pre !is_leaf(v) and s>=v_s and e<=v_e
         */
        std::pair<range_vec_type, range_vec_type>
        expand(const node_type& v,
               range_vec_type&& ranges) const {
            auto v_sp_rank = m_tree_rank(v.offset);  // this is already calculated in expand(v)
            range_vec_type res(ranges.size());
            size_t i = 0;
            for (auto& r : ranges) {
                auto sp_rank    = m_tree_rank(v.offset + r.first);
                auto right_size = m_tree_rank(v.offset + r.second + 1)
                                  - sp_rank;
                auto left_size  = (r.second-r.first+1)-right_size;

                auto right_sp = sp_rank - v_sp_rank;
                auto left_sp  = r.first - right_sp;

                r = range_type(left_sp, left_sp + left_size - 1);
                res[i++] = range_type(right_sp, right_sp + right_size - 1);
            }
            return make_pair(ranges, std::move(res));
        }

        //! Returns for a range its left and right child ranges
        /*! \param v An inner node of an wavelet tree.
         *  \param r A ranges [s,e], such that [s,e] is
         *           contained in v=[v_s,v_e].
         *  \return A range pair. The first element of the
         *          range pair correspond to the original range
         *          mapped to the left child of v; the second element to the
         *          range mapped to the right child of v.
         *  \pre !is_leaf(v) and s>=v_s and e<=v_e
         */
        std::pair<range_type, range_type>
        expand(const node_type& v, const range_type& r) const {
            auto v_sp_rank = m_tree_rank(v.offset);  // this is already calculated in expand(v)
            auto sp_rank    = m_tree_rank(v.offset + r.first);
            auto right_size = m_tree_rank(v.offset + r.second + 1)
                              - sp_rank;
            auto left_size  = (r.second-r.first+1)-right_size;

            auto right_sp = sp_rank - v_sp_rank;
            auto left_sp  = r.first - right_sp;

            return make_pair(range_type(left_sp, left_sp + left_size - 1),
                             range_type(right_sp, right_sp + right_size - 1));
        }

        //! return the path to the leaf for a given symbol
        std::pair<uint64_t,uint64_t> path(value_type c) const {
            return {m_max_level,c};
        }
};

}// end namespace sdsl
#endif
