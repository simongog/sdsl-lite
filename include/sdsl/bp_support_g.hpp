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
/*! \file bp_support_g.hpp
    \brief bp_support_g.hpp contains an implementation of a balanced parentheses support data structure.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_BP_SUPPORT_G
#define INCLUDED_SDSL_BP_SUPPORT_G

#include "int_vector.hpp"
#include "nearest_neighbour_dictionary.hpp"
#include "rmq_support.hpp"
#include "rank_support.hpp"
#include "select_support.hpp"
#include "bp_support_algorithm.hpp"
#include "util.hpp"
#include <stack>
#include <map>
#include <set>
#include <utility>
#include <stdexcept>

namespace sdsl
{

//! A class that provides support for bit_vectors that represent a BP sequence.
/*! This data structure supports the following operations:
 *   - find_open
 *   - find_close
 *   - enclose
 *   - double_enclose
 *   - rank
 *   - select
 *   - excess
 *   - rr_enclose
 *  An opening parenthesis in the balanced parentheses sequence is represented by a 1 in the bit_vector
 *  and a closing parenthesis by a 0.
 *
 *  \tparam t_nnd     Type which supports rank and select with little space on sparse populated bit_vectors.
 *  \tparam t_rank    Type of rank support structure.
 *  \tparam t_select  Type of select support structure.
 *  \tparam t_rmq     Type which supports range maximum queries on a int_vector<>.
 * \par Reference
 *      Richard F. Geary, Naila Rahman, Rajeev Raman, Venkatesh Raman:
 *      A Simple Optimal Representation for Balanced Parentheses.
 *      CPM 2004: 159-172
 *
 *  @ingroup bps
 */
template<class t_nnd = nearest_neighbour_dictionary<30>,
         class t_rank = rank_support_v5<>,
         class t_select = select_support_mcl<>,
         class t_rmq = range_maximum_support_sparse_table<>,
         uint32_t t_bs=840>
class bp_support_g
{
        static_assert(t_bs > 2, "bp_support_g: block size must be greater than 2.");
    public:
        typedef bit_vector::size_type size_type;
        typedef t_nnd                 nnd_type;
        typedef t_rank                rank_type;
        typedef t_select              select_type;
        typedef t_rmq                 rmq_type;
    private:
        const bit_vector* m_bp;             // the supported BP sequence as bit_vector
        rank_type         m_rank_bp;        // rank support for the BP sequence => see excess() and rank()
        select_type       m_select_bp;      // select support for the BP sequence => see select()

        nnd_type          m_nnd;            // nearest neighbour dictionary for pioneers bit_vector

        bit_vector        m_pioneer_bp;     // first level of recursion: BP sequence of the pioneers
        rank_type         m_rank_pioneer_bp;// rank for the BP sequence of the pioneers
        nnd_type          m_nnd2;           // nearest neighbour dictionary for pioneers of pioneers bit_vector
        int_vector<>      m_match;          //
        int_vector<>      m_enclose;        //
        rmq_type          m_range_max_match;// range maximum support for m_match

        size_type m_size;
        size_type m_blocks; // number of blocks

        void copy(const bp_support_g& bp_support) {
            m_bp = bp_support.m_bp;
            m_rank_bp = bp_support.m_rank_bp;
            m_rank_bp.set_vector(m_bp);
            m_select_bp = bp_support.m_select_bp;
            m_select_bp.set_vector(m_bp);

            m_nnd = bp_support.m_nnd;

            m_pioneer_bp = bp_support.m_pioneer_bp;
            m_rank_pioneer_bp = bp_support.m_rank_pioneer_bp;
            m_rank_pioneer_bp.set_vector(&m_pioneer_bp);
            m_nnd2 = bp_support.m_nnd2;
            m_match = bp_support.m_match;
            m_enclose = bp_support.m_enclose;
            m_range_max_match = bp_support.m_range_max_match;
            m_range_max_match.set_vector(&m_match);

            m_size = bp_support.m_size;
            m_blocks = bp_support.m_blocks;
        }

        /*! Calculates the excess value at index i in the pioneer bitmap.
         * \param i The index of which the excess value should be calculated.
         */
        inline size_type excess_pioneer(size_type i)const {
            return (m_rank_pioneer_bp(i+1)<<1)-i-1;
        }

    public:
        const rank_type&   bp_rank   = m_rank_bp;
        const select_type& bp_select = m_select_bp;

        //! Constructor
        explicit
        bp_support_g(const bit_vector* bp = nullptr) : m_bp(bp),
            m_size(bp==nullptr?0:bp->size()), m_blocks((m_size+t_bs-1)/t_bs) {
            if (bp == nullptr)
                return;
            util::init_support(m_rank_bp, bp);
            util::init_support(m_select_bp, bp);
            bit_vector pioneer = calculate_pioneers_bitmap(*m_bp, t_bs);
            m_nnd = nnd_type(pioneer);
            m_pioneer_bp.resize(m_nnd.ones());
            for (size_type i=1; i<= m_nnd.ones(); ++i)
                m_pioneer_bp[i-1] = (*m_bp)[m_nnd.select(i)];
            util::init_support(m_rank_pioneer_bp, &m_pioneer_bp);
            pioneer = calculate_pioneers_bitmap(m_pioneer_bp, t_bs);
            m_nnd2 = nnd_type(pioneer);

            bit_vector pioneer_bp2 = bit_vector(m_nnd2.ones());
            for (size_type i=1; i<= m_nnd2.ones(); ++i)
                pioneer_bp2[i-1] = m_pioneer_bp[m_nnd2.select(i)];
            calculate_matches(pioneer_bp2, m_match);
            calculate_enclose(pioneer_bp2, m_enclose);
            m_range_max_match = rmq_type(&m_match);
        }

        //! Copy constructor
        bp_support_g(const bp_support_g& bp_support) {
            copy(bp_support);
        }

        //! Move constructor
        bp_support_g(bp_support_g&& bp_support) {
            *this = std::move(bp_support);
        }

        //! Assignment operator
        bp_support_g& operator=(const bp_support_g& bp_support) {
            if (this != &bp_support) {
                copy(bp_support);
            }
            return *this;
        }

        //! Assignment operator
        bp_support_g& operator=(bp_support_g&& bp_support) {
            if (this != &bp_support) {
                m_bp = std::move(bp_support.m_bp);
                m_rank_bp = std::move(bp_support.m_rank_bp);
                m_rank_bp.set_vector(m_bp);
                m_select_bp = std::move(bp_support.m_select_bp);
                m_select_bp.set_vector(m_bp);

                m_nnd = std::move(bp_support.m_nnd);

                m_pioneer_bp = std::move(bp_support.m_pioneer_bp);
                m_rank_pioneer_bp = std::move(bp_support.m_rank_pioneer_bp);
                m_rank_pioneer_bp.set_vector(&m_pioneer_bp);
                m_nnd2 = std::move(bp_support.m_nnd2);
                m_match = std::move(bp_support.m_match);
                m_enclose = std::move(bp_support.m_enclose);
                m_range_max_match = std::move(bp_support.m_range_max_match);
                m_range_max_match.set_vector(&m_match);

                m_size = std::move(bp_support.m_size);
                m_blocks = std::move(bp_support.m_blocks);
            }
            return *this;
        }

        void swap(bp_support_g& bp_support) {
            m_rank_bp.swap(bp_support.m_rank_bp);
            m_select_bp.swap(bp_support.m_select_bp);

            m_nnd.swap(bp_support.m_nnd);

            m_pioneer_bp.swap(bp_support.m_pioneer_bp);
            util::swap_support(m_rank_pioneer_bp, bp_support.m_rank_pioneer_bp,
                               &m_pioneer_bp, &(bp_support.m_pioneer_bp));

            m_nnd2.swap(bp_support.m_nnd2);

            m_match.swap(bp_support.m_match);
            m_enclose.swap(bp_support.m_enclose);
            util::swap_support(m_range_max_match, bp_support.m_range_max_match,
                               &m_match, &(bp_support.m_match));

            std::swap(m_size, bp_support.m_size);
            std::swap(m_blocks, bp_support.m_blocks);
        }

        void set_vector(const bit_vector* bp) {
            m_bp = bp;
            m_rank_bp.set_vector(bp);
            m_select_bp.set_vector(bp);
        }

        /*! Calculates the excess value at index i.
         * \param i The index of which the excess value should be calculated.
         */
        inline size_type excess(size_type i)const {
            return (m_rank_bp(i+1)<<1)-i-1;
        }

        /*! Returns the number of opening parentheses up to and including index i.
         * \pre{ \f$ 0\leq i < size() \f$ }
         */
        size_type rank(size_type i)const {
            return m_rank_bp(i+1);
        }

        /*! Returns the index of the i-th opening parenthesis.
         * \param i Number of the parenthesis to select.
         * \pre{ \f$1\leq i < rank(size())\f$ }
         * \post{ \f$ 0\leq select(i) < size() \f$ }
         */
        size_type select(size_type i)const {
            return m_select_bp(i);
        }

        /*! Calculate the index of the matching closing parenthesis to the parenthesis at index i.
         * \param i Index of an parenthesis. 0 <= i < size().
         * \return * i, if the parenthesis at index i is closing,
         *         * the position j of the matching closing parenthesis, if a matching parenthesis exists,
         *         * size() if no matching closing parenthesis exists.
         */
        size_type find_close(size_type i)const {
            assert(i < m_size);
            if (!(*m_bp)[i]) {// if there is a closing parenthesis at index i return i
                return i;
            }
            size_type mi = 0; // match for i
            if ((mi=near_find_close(*m_bp, i, t_bs))==i) {
                const size_type i2 = m_nnd.rank(i+1)-1; // lemma that this gives us an opening pioneer
                assert(m_pioneer_bp[i2]==1); // assert that i2 is an opening parenthesis
                size_type mi2 = 0; // match for i2
                if ((mi2=near_find_close(m_pioneer_bp, i2, t_bs)) == i2) {
                    const size_type i3 = m_nnd2.rank(i2+1)-1;
                    const size_type mi3 = m_match[i3];     assert(mi3>i3); // assert that i3 is an opening parenthesis
                    mi2 = m_nnd2.select(mi3+1); // matching pioneer position in pioneer_bp
                    mi2 = (mi2/t_bs)*t_bs;
                    size_type epb = excess_pioneer(mi2);// excess of first parenthesis in the pioneer block
                    const size_type ei = excess_pioneer(i2);// excess of pioneer
                    /* invariant: epb >= ei-1 */ assert(epb+1 >= ei);
                    while (epb+1 != ei) {
                        assert(mi2 < m_pioneer_bp.size());
                        if (m_pioneer_bp[++mi2])
                            ++epb;
                        else
                            --epb;
                    }
                }
                mi = m_nnd.select(mi2+1);  /* matching pioneer position in bp */ assert((*m_bp)[mi]==0);
                mi = (mi/t_bs)*t_bs;
                size_type epb = excess(mi); // excess of first parenthesis in the pioneer block
                const size_type ei = excess(i);  // excess at position i
                /* invariant: epb >= ei-1 */ assert(epb+1 >= ei);
                while (epb+1 != ei) {
                    assert(mi < m_size);
                    if ((*m_bp)[++mi])
                        ++epb;
                    else
                        --epb;
                }
            }
            return mi;
        }

        //! Calculate the matching opening parenthesis to the closing parenthesis at position i
        /*! \param i Index of a closing parenthesis.
          * \return * i, if the parenthesis at index i is closing,
          *         * the position j of the matching opening parenthesis, if a matching parenthesis exists,
          *         * size() if no matching closing parenthesis exists.
          */
        size_type find_open(size_type i)const {
            assert(i < m_size);
            if ((*m_bp)[i]) {// if there is a opening parenthesis at index i return i
                return i;
            }
            size_type mi = 0; // match for i
            if ((mi=near_find_open(*m_bp, i, t_bs)) == i) {
                const size_type i2 = m_nnd.rank(i); // lemma that this gives us an closing pioneer
                assert(m_pioneer_bp[i2]==0); // assert that i2 is an opening parenthesis
                const size_type mi2 = find_open_in_pioneers(i2);         assert(m_pioneer_bp[mi2]==1);
                mi = m_nnd.select(mi2+1);  /* matching pioneer position in bp */ assert((*m_bp)[mi]==1);
                mi = (mi/t_bs)*t_bs + t_bs - 1;     assert(mi < i);
                size_type epb = excess(mi); // excess of last parenthesis in the pioneer block
                const size_type ei = excess(i);  // excess at position i
                /*invariant: epb >= ei+1*/      assert(epb >= ei+1);
                while (epb != ei) {
                    assert(mi < m_size);
                    if ((*m_bp)[mi--])
                        --epb;
                    else
                        ++epb;
                }
                ++mi;
            }
            return mi;
        }

        inline size_type find_open_in_pioneers(size_type i)const {
            size_type mi = 0; // match for i
            if ((mi=near_find_open(m_pioneer_bp, i, t_bs))==i) {
                const size_type i3 = m_nnd2.rank(i);
                const size_type mi3 = m_match[i3];      assert(mi3<i3); // assert that i3 is an closing parenthesis
                mi = m_nnd2.select(mi3+1); // matching pioneer position in pioneer_bp
                mi = (mi/t_bs)*t_bs + t_bs - 1;
                size_type epb2 = excess_pioneer(mi);// excess of last parenthesis in the pioneer block
                const size_type ei = excess_pioneer(i);// excess of pioneer
                /* invariant: epb2 >= ei+1 */  assert(epb2 >= ei+1);
                while (epb2 != ei) {
                    assert(mi < m_pioneer_bp.size());
                    if (m_pioneer_bp[mi--])
                        --epb2;
                    else
                        ++epb2;
                }
                ++mi;
            }
            return mi;
        }

        //! Calculate the index of the opening parenthesis corresponding to the closest matching parenthesis pair enclosing i.
        /*! \param i Index of an opening parenthesis.
         *  \return The index of the opening parenthesis corresponding to the closest matching parenthesis pair enclosing i,
         *          or size() if no such pair exists.
         */
        size_type enclose(size_type i)const {
            assert(i < m_size);
            if (!(*m_bp)[i]) { // if there is closing parenthesis at position i
                return find_open(i);
            }
            const size_type exi = excess(i);
            if (exi == 1)  // if i is not enclosed by a parentheses pair..
                return size();
            size_type ei; // enclose  for i
            if ((ei=near_enclose(*m_bp, i, t_bs))==i) {
                const size_type i2 = m_nnd.rank(i); // next parenthesis in the pioneer bitmap
                size_type ei2; // enclose for i2
                if (m_pioneer_bp[i2]) { // search enclose in the pioneer bp
                    if ((ei2=near_enclose(m_pioneer_bp, i2, t_bs))==i2) {
                        const size_type i3  = m_nnd2.rank(i2); // next parenthesis in the pioneer2 bitmap
                        const size_type ei3 = m_enclose[i3];                              assert(ei3<i3);    // assert that enclose answer is valid
                        ei2 = m_nnd2.select(ei3+1);                                       assert(m_pioneer_bp[ei2] == 1);
                        ei2 = (ei2/t_bs)*t_bs + t_bs - 1;         assert(ei2 < i2);
                        size_type epb2 = excess_pioneer(ei2);// excess of the last parenthesis in the pioneer block
                        const size_type exi2 = excess_pioneer(i2);// excess
                        /* invariant epb2+1 >= exi2 */                                    assert(epb2+1 >= exi2);
                        while (epb2+2 != exi2) {
                            if (m_pioneer_bp[ei2--])
                                --epb2;
                            else
                                ++epb2;
                        }
                        ++ei2;
                    }
                } else {
                    // if the next parenthesis in the pioneer bitmap is an closing parenthesis findopen on m_pioneer_bp
                    ei2 = find_open_in_pioneers(i2);
                }
                assert(m_pioneer_bp[ei2]==1);
                ei = m_nnd.select(ei2+1);                                  assert((*m_bp)[ei]==1);
                ei = (ei/t_bs)*t_bs + t_bs - 1;    assert(ei < i);
                size_type epb = excess(ei); // excess of the last parenthesis in the pioneer block
                /* invariant epb+1 >= exi */ assert(epb+1 >= exi);
                while (epb+2 != exi) {
                    if ((*m_bp)[ei--])
                        --epb;
                    else
                        ++epb;
                }
                ++ei;
            }
            return ei;
        }

        //! The range restricted enclose operation
        /*! \param i Index of an opening parenthesis.
         *  \param j Index of an opening parenthesis/ \f$ i<j \wedge findclose(i) < j \f$.
         *  \return The smallest index, say k, of an opening parenthesis such that findclose(i) < k < j and
         *  findclose(j) < findclose(k). If such a k does not exists, restricted_enclose(i,j) returns size().
         * \par Time complexity
         *     \f$ \Order{block\_size} \f$
         */
        size_type rr_enclose(const size_type i, const size_type j)const {
            assert(j > i and j < m_size);
            const size_type mip1 = find_close(i)+1;
            if (mip1 >= j)
                return size();
            return rmq_open(mip1, j);
        }

        /*! Search the interval [l,r-1] for an opening parenthesis, say i, such that find_close(i) >= r.
         * \param l The left end (inclusive) of the interval to search for the result.
         * \param r The right end (exclusive) of the interval to search for the result.
         * \return the minimal opening parenthesis i with \f$ \ell \leq i < r \f$ and \f$ find_close(i) \geq r \f$;
         *         if no such i exists size() is returned.
         * The algorithm consists of 4 steps:
         * 1. scan back from position r to the begin of that block
         * 2. recursively scan back the pioneers of the blocks which lie in between the blocks of l and r
         * 3. scan from position l to the end of the block, which contains l
         * 4. check if there exists a valid solution and return
         * \par Time complexity
         *     \f$ \Order{block\_size} \f$
         */
        size_type rmq_open(const size_type l, const size_type r)const {
            if (l >= r)
                return size();
            size_type        min_ex_pos = r;

            if (l/t_bs == r/t_bs) {
                min_ex_pos = near_rmq_open(*m_bp, l, r);
            } else { // parentheses pair does not start in the same block
//                assert( l>1 ); // mi is at greater or equal than 1
                // note: mi and r are not in the same block
                size_type        k, ex;    // helper variables
                size_type         min_ex = excess(r);// + 2*((*m_bp[r])==0); // minimal excess
                const size_type bl = (l/t_bs+1)*t_bs;  // leftmost position of the leftmost block between the blocks of l and r
                const size_type br = (r/t_bs)*t_bs;       // leftmost position of the block of r


                // 1.2
                size_type l_ = m_nnd.rank(l);   //  l_ inclusive
                size_type r_     = m_nnd.rank(r);   // r_ exclusive

                if (r_ > l_) {
                    size_type min_ex_pos_ = r_;
                    if (l_/t_bs == r_/t_bs) {
                        min_ex_pos_ = near_rmq_open(m_pioneer_bp, l_, r_);
                    } else if (r_ < m_pioneer_bp.size()) {
                        size_type min_ex_      = excess_pioneer(r_)+2*(m_pioneer_bp[r_]==0);
                        const size_type bl_ = (l_/t_bs+1)*t_bs;
                        const size_type br_ = (r_/t_bs)*t_bs;

                        // 2.2
                        size_type l__ = m_nnd2.rank(l_);   // l__ inclusive
                        size_type r__     = m_nnd2.rank(r_);   // r__ exclusive
                        if (r__ > l__) {
                            size_type max_match = 0;
                            k = m_range_max_match(l__, r__-1);
                            max_match = m_match[k];
                            if (max_match >= r__) {
                                k = m_nnd2.select(k+1);
                                if (k < r_ and(ex=excess_pioneer(k)) < min_ex_) {
                                    min_ex_ = ex; min_ex_pos_ = k;
                                }
                            }
                        }
                        if (min_ex_pos_ == r_) {
                            // 2.1
                            k = near_rmq_open(m_pioneer_bp, br_, r_);
                            if (k < r_ and(ex=excess_pioneer(k)) < min_ex_) {
                                min_ex_ = ex; min_ex_pos_ = k;
                            }
                        }
                        // 2.3
                        k = near_rmq_open(m_pioneer_bp, l_, bl_);
                        if (k < bl_ and(ex=excess_pioneer(k)) < min_ex_) {
                            min_ex_ = ex; min_ex_pos_ = k;
                        }
                    }
                    // 2.4
                    if (min_ex_pos_ < r_) {
                        k = m_nnd.select(min_ex_pos_ + 1);
                        if ((ex=excess(k)) < min_ex) {
                            min_ex = ex; min_ex_pos = k;
                        }
                    }
                }
                if (min_ex_pos == r) {
                    // 1.1
                    k = near_rmq_open(*m_bp, br, r);
                    if (k < r and(ex=excess(k)) < min_ex) {
                        min_ex        = ex; min_ex_pos     = k;
                    }
                }
                // 1.3
                k = near_rmq_open(*m_bp, l, bl);
                if (k < bl and(ex=excess(k)) < min_ex) {
                    min_ex = ex; min_ex_pos = k;
                }
            }
            // 1.4
            if (min_ex_pos < r)
                return min_ex_pos;
            else
                return size();
        }

        //! The range restricted enclose operation
        /*! \param i Index of an opening parenthesis.
         *  \param j Index of an opening parenthesis/ \f$ i<j \wedge findclose(i) < j \f$.
         *  \return The smallest index, say k, of an opening parenthesis such that findclose(i) < k < j and
         *  findclose(j) < findclose(k). If such a k does not exists, restricted_enclose(i,j) returns size().
        */
        size_type rr_enclose_naive(size_type i, size_type j)const {
            assert(j > i and j < m_size);
            size_type mi = find_close(i); // matching parenthesis to i
            assert(mi > i and mi < j);
            assert(find_close(j) > j);
            size_type k = enclose(j);
            if (k == m_size or k < i) // there exists no opening parenthesis at position mi<k<j.
                return m_size;
            size_type kk;
            do {
                kk = k;
                k = enclose(k);
            } while (k != m_size and k > mi);
            return kk;
        }

        //! The range minimum query (rmq) returns the index of the parenthesis with minimal excess in the range \f$[l..r]\f$
        /*! \param l The left border of the interval \f$[l..r]\f$ (\f$l\leq r\f$).
         *  \param r The right border of the interval \f$[l..r]\f$ (\f$l \leq r\f$).
         */
        size_type rmq(size_type l, size_type r)const {
            assert(l<=r);
            size_type m = rmq_open(l, r+1);
            if (m==l)
                return l;
            else { // m>l and m<=r
                assert(0 == (*m_bp)[m-1]);
                size_type prev_open = m_rank_bp(m);
                if (prev_open == 0 or m_select_bp(prev_open) < l) { // if there exists no opening parenthesis to the left of m which is greater or equal than l
                    return l;
                } else {
                    return m-1;
                }
            }
        }

        //! The double enclose operation
        /*! \param i Index of an opening parenthesis.
         *  \param j Index of an opening parenthesis \f$ i<j \wedge findclose(i) < j \f$.
         *  \return The maximal opening parenthesis, say k, such that \f$ k<j \wedge k>findclose(j) \f$.
         *          If such a k does not exists, double_enclose(i,j) returns size().
         */
        size_type double_enclose(size_type i, size_type j)const {
            assert(j > i);
            assert((*m_bp)[i]==1 and(*m_bp)[j]==1);
            size_type k = rr_enclose(i, j);
            if (k == size())
                return enclose(j);
            else
                return enclose(k);
        }

        //! Return the number of zeros which procede position i in the balanced parentheses sequence.
        /*! \param i Index of an parenthesis.
         */
        size_type preceding_closing_parentheses(size_type i)const {
            assert(i < m_size);
            if (!i) return 0;
            size_type ones = m_rank_bp(i);
            if (ones) { // ones > 0
                assert(m_select_bp(ones) < i);
                return i - m_select_bp(ones) - 1;
            } else {
                return i;
            }
        }

        /*! The size of the supported balanced parentheses sequence.
         * \return the size of the supported balanced parentheses sequence.
         */
        size_type size() const {
            return m_size;
        }

        //! Serializes the bp_support_g to a stream.
        /*!
         * \param out The outstream to which the data structure is written.
         * \return The number of bytes written to out.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_rank_bp.serialize(out, child, "bp_rank");
            written_bytes += m_select_bp.serialize(out, child, "bp_select");
            written_bytes += m_nnd.serialize(out, child,"nearest_neighbor_dictionary");

            written_bytes += m_pioneer_bp.serialize(out, child, "pioneer_bp");
            written_bytes += m_rank_pioneer_bp.serialize(out, child, "pioneer_bp_rank");
            written_bytes += m_nnd2.serialize(out, child, "nearest_neighbor_dictionary2");
            written_bytes += m_match.serialize(out, child, "match_answers");
            written_bytes += m_enclose.serialize(out, child, "enclose_answers");
            written_bytes += m_range_max_match.serialize(out, child, "rmq_answers");

            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += write_member(m_blocks, out, child, "block_cnt");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load the bp_support_g for a bit_vector v.
        /*!
         * \param in The instream from which the data structure is read.
         * \param bp Bit vector representing a balanced parentheses sequence that is supported by this data structure.
         */
        void load(std::istream& in, const bit_vector* bp) {
            m_bp = bp;
            m_rank_bp.load(in, m_bp);
            m_select_bp.load(in, m_bp);
            m_nnd.load(in);

            m_pioneer_bp.load(in);
            m_rank_pioneer_bp.load(in, &m_pioneer_bp);
            m_nnd2.load(in);
            m_match.load(in);
            m_enclose.load(in);
            m_range_max_match.load(in, &m_match);
            read_member(m_size, in);
            assert(m_size == bp->size());
            read_member(m_blocks, in);
        }
};

}// end namespace




#endif
