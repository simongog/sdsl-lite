/* sdsl - succinct data structures library
    Copyright (C) 2011-2013 Simon Gog

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
/*! \file rrr_vector.hpp
   \brief rrr_vector.hpp contains the sdsl::rrr_vector class, and
          classes which support rank and select for rrr_vector.
   \author Simon Gog, Matthias Petri
*/
#ifndef INCLUDED_SDSL_RRR_VECTOR
#define INCLUDED_SDSL_RRR_VECTOR

#include "int_vector.hpp"
#include "util.hpp"
#include "rrr_helper.hpp" // for binomial helper class
#include "iterators.hpp"
#include <vector>
#include <algorithm> // for next_permutation
#include <iostream>

//! Namespace for the succinct data structure library
namespace sdsl
{

// forward declaration needed for friend declaration
template<uint8_t t_b=1, uint16_t t_bs=15, class t_rac=int_vector<>, uint16_t t_k=32>
class rank_support_rrr;                // in rrr_vector

// forward declaration needed for friend declaration
template<uint8_t t_b=1, uint16_t t_bs=15, class t_rac=int_vector<>, uint16_t t_k=32>
class select_support_rrr;                // in rrr_vector

//! A \f$H_0f$-compressed bitvector representation.
/*!
 *   \tparam t_bs   Size of a basic block.
 *   \tparam t_rac  Random access integer vector. Use to store the block types.
 *                  It is possible to use WTs for t_rac.
 *   \tparam t_k    A rank sample value is stored before every t_k-th basic block.
 *
 *   References:
 *    - Rasmus Pagh
 *      Low redundancy in dictionaries with O(1) worst case lookup time
 *      Technical Report 1998.
 *      ftp://ftp.cs.au.dk/BRICS/Reports/RS/98/28/BRICS-RS-98-28.pdf, Section 2.
 *    - Rajeev Raman, V. Raman and S. Srinivasa Rao
 *      Succinct Indexable Dictionaries with Applications to representations
 *      of k-ary trees and multi-sets.
 *      SODA 2002.
 *    - Francisco Claude, Gonzalo Navarro:
 *      Practical Rank/Select Queries over Arbitrary Sequences.
 *      SPIRE 2008: 176-187
 *    - On the fly-decoding and encoding was discovered in;
 *      Gonzalo Navarro, Eliana Providel:
 *      Fast, Small, Simple Rank/Select on Bitmaps.
 *      SEA 2012
 *
 *    In this version the block size can be adjust by the template parameter t_bs!
 *    \sa sdsl::rrr_vector for a specialized version for block_size=15
 */
template<uint16_t t_bs=63, class t_rac=int_vector<>, uint16_t t_k=32>
class rrr_vector
{
        static_assert(t_bs >= 3 and t_bs <= 256 , "rrr_vector: block size t_bs must be 3 <= t_bs <= 256.");
        static_assert(t_k > 1, "rrr_vector: t_k must be > 0.");
    public:
        typedef bit_vector::size_type                    size_type;
        typedef bit_vector::value_type                   value_type;
        typedef bit_vector::difference_type              difference_type;
        typedef t_rac                                    rac_type;
        typedef random_access_const_iterator<rrr_vector> iterator;
        typedef bv_tag                                   index_category;

        typedef rank_support_rrr<1, t_bs, t_rac, t_k>   rank_1_type;
        typedef rank_support_rrr<0, t_bs, t_rac, t_k>   rank_0_type;
        typedef select_support_rrr<1, t_bs, t_rac, t_k> select_1_type;
        typedef select_support_rrr<0, t_bs, t_rac, t_k> select_0_type;

        friend class rank_support_rrr<0, t_bs, t_rac, t_k>;
        friend class rank_support_rrr<1, t_bs, t_rac, t_k>;
        friend class select_support_rrr<0, t_bs, t_rac, t_k>;
        friend class select_support_rrr<1, t_bs, t_rac, t_k>;

        typedef rrr_helper<t_bs> rrr_helper_type;
        typedef typename rrr_helper_type::number_type number_type;

        enum { block_size = t_bs };
    private:
        size_type    m_size = 0;  // Size of the original bit_vector.
        rac_type     m_bt;     // Vector for the block types (bt). bt equals the
        // number of set bits in the block.
        bit_vector   m_btnr;   // Compressed block type numbers.
        int_vector<> m_btnrp;  // Sample pointers into m_btnr.
        int_vector<> m_rank;   // Sample rank values.
        bit_vector   m_invert; // Specifies if a superblock (i.e. t_k blocks)
        // have to be considered as inverted i.e. 1 and
        // 0 are swapped

        void copy(const rrr_vector& rrr) {
            m_size = rrr.m_size;
            m_bt = rrr.m_bt;
            m_btnr = rrr.m_btnr;
            m_btnrp = rrr.m_btnrp;
            m_rank = rrr.m_rank;
            m_invert = rrr.m_invert;
        }

    public:
        const rac_type& bt     = m_bt;
        const bit_vector& btnr = m_btnr;

        //! Default constructor
        rrr_vector() {};

        //! Copy constructor
        rrr_vector(const rrr_vector& rrr) {
            copy(rrr);
        }

        //! Move constructor
        rrr_vector(rrr_vector&& rrr) : m_size(std::move(rrr.m_size)),
            m_bt(std::move(rrr.m_bt)),
            m_btnr(std::move(rrr.m_btnr)), m_btnrp(std::move(rrr.m_btnrp)),
            m_rank(std::move(rrr.m_rank)), m_invert(std::move(rrr.m_invert)) {}

        //! Constructor
        /*!
        *  \param bv  Uncompressed bitvector.
        *  \param k   Store rank samples and pointers each k-th blocks.
        */
        rrr_vector(const bit_vector& bv) {
            m_size = bv.size();
            int_vector<> bt_array;
            bt_array.width(bits::hi(t_bs)+1);
            bt_array.resize((m_size+t_bs)/((size_type)t_bs)); // blocks for the bt_array + a dummy block at the end,
            // if m_size%t_bs == 0

            // (1) calculate the block types and store them in m_bt
            size_type pos = 0, i = 0, x;
            size_type btnr_pos = 0;
            size_type sum_rank = 0;
            while (pos + t_bs <= m_size) { // handle all blocks full blocks
                bt_array[ i++ ] = x = rrr_helper_type::get_bt(bv, pos, t_bs);
                sum_rank += x;
                btnr_pos += rrr_helper_type::space_for_bt(x);
                pos += t_bs;
            }
            if (pos < m_size) { // handle last not full block
                bt_array[ i++ ] = x = rrr_helper_type::get_bt(bv, pos, m_size - pos);
                sum_rank += x;
                btnr_pos += rrr_helper_type::space_for_bt(x);
            }
            m_btnr  = bit_vector(std::max(btnr_pos, (size_type)64), 0);      // max necessary for case: t_bs == 1
            m_btnrp = int_vector<>((bt_array.size()+t_k-1)/t_k, 0,  bits::hi(btnr_pos)+1);
            m_rank  = int_vector<>((bt_array.size()+t_k-1)/t_k + ((m_size % (t_k*t_bs))>0), 0, bits::hi(sum_rank)+1);
            //                                                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            //   only add a finishing block, if the last block of the superblock is not a dummy block
            m_invert = bit_vector((bt_array.size()+t_k-1)/t_k, 0);

            // (2) calculate block type numbers and pointers into btnr and rank samples
            pos = 0; i = 0;
            btnr_pos= 0, sum_rank = 0;
            bool invert = false;
            while (pos + t_bs <= m_size) {  // handle all full blocks
                if ((i % t_k) == (size_type)0) {
                    m_btnrp[ i/t_k ] = btnr_pos;
                    m_rank[ i/t_k ] = sum_rank;
                    // calculate invert bit for that superblock
                    if (i+t_k <= bt_array.size()) {
                        size_type gt_half_t_bs = 0; // counter for blocks greater than half of the blocksize
                        for (size_type j=i; j < i+t_k; ++j) {
                            if (bt_array[j] > t_bs/2)
                                ++gt_half_t_bs;
                        }
                        if (gt_half_t_bs > (t_k/2)) {
                            m_invert[ i/t_k ] = 1;
                            for (size_type j=i; j < i+t_k; ++j) {
                                bt_array[j] = t_bs - bt_array[j];
                            }
                            invert = true;
                        } else {
                            invert = false;
                        }
                    } else {
                        invert = false;
                    }
                }
                uint16_t space_for_bt = rrr_helper_type::space_for_bt(x=bt_array[i++]);
                sum_rank += (invert ? (t_bs - x) : x);
                if (space_for_bt) {
                    number_type bin = rrr_helper_type::decode_btnr(bv, pos, t_bs);
                    number_type nr = rrr_helper_type::bin_to_nr(bin);
                    rrr_helper_type::set_bt(m_btnr, btnr_pos, nr, space_for_bt);
                }
                btnr_pos += space_for_bt;
                pos += t_bs;
            }
            if (pos < m_size) { // handle last not full block
                if ((i % t_k) == (size_type)0) {
                    m_btnrp[ i/t_k ] = btnr_pos;
                    m_rank[ i/t_k ] = sum_rank;
                    m_invert[ i/t_k ] = 0; // default: set last block to not inverted
                    invert = false;
                }
                uint16_t space_for_bt = rrr_helper_type::space_for_bt(x=bt_array[i++]);
//          no extra dummy block added to bt_array, therefore this condition should hold
                assert(i == bt_array.size());
                sum_rank += invert ? (t_bs - x) : x;
                if (space_for_bt) {
                    number_type bin = rrr_helper_type::decode_btnr(bv, pos, m_size-pos);
                    number_type nr = rrr_helper_type::bin_to_nr(bin);
                    rrr_helper_type::set_bt(m_btnr, btnr_pos, nr, space_for_bt);
                }
                btnr_pos += space_for_bt;
                assert(m_rank.size()-1 == ((i+t_k-1)/t_k));
            } else { // handle last empty full block
                assert(m_rank.size()-1 == ((i+t_k-1)/t_k));
            }
            // for technical reasons we add a last element to m_rank
            m_rank[ m_rank.size()-1 ] = sum_rank; // sum_rank contains the total number of set bits in bv
            m_bt = bt_array;
        }

        //! Swap method
        void swap(rrr_vector& rrr) {
            if (this != &rrr) {
                std::swap(m_size, rrr.m_size);
                m_bt.swap(rrr.m_bt);
                m_btnr.swap(rrr.m_btnr);
                m_btnrp.swap(rrr.m_btnrp);
                m_rank.swap(rrr.m_rank);
                m_invert.swap(rrr.m_invert);
            }
        }

        //! Accessing the i-th element of the original bit_vector
        /*! \param i An index i with \f$ 0 \leq i < size()  \f$.
           \return The i-th bit of the original bit_vector
        */
        value_type operator[](size_type i)const {
            size_type bt_idx = i/t_bs;
            uint16_t bt = m_bt[bt_idx];
            size_type sample_pos = bt_idx/t_k;
            if (m_invert[sample_pos])
                bt = t_bs - bt;
#ifndef RRR_NO_OPT
            if (bt == 0 or bt == t_bs) { // very effective optimization
                return bt>0;
            }
#endif
            uint16_t off = i % t_bs; //i - bt_idx*t_bs;
            size_type btnrp = m_btnrp[ sample_pos ];
            for (size_type j = sample_pos*t_k; j < bt_idx; ++j) {
                btnrp += rrr_helper_type::space_for_bt(m_bt[j]);
            }
            uint16_t btnrlen = rrr_helper_type::space_for_bt(bt);
            number_type btnr = rrr_helper_type::decode_btnr(m_btnr, btnrp, btnrlen);
            return rrr_helper_type::decode_bit(bt, btnr, off);
        }

        //! Get the integer value of the binary string of length len starting at position idx.
        /*! \param idx Starting index of the binary representation of the integer.
         *  \param len Length of the binary representation of the integer. Default value is 64.
         *   \returns The integer value of the binary string of length len starting at position idx.
         *
         *  \pre idx+len-1 in [0..size()-1]
         *  \pre len in [1..64]
         */
        uint64_t get_int(size_type idx, uint8_t len=64)const {
            uint64_t res = 0;
            size_type bb_idx = idx/t_bs; // begin block index
            size_type bb_off = idx%t_bs; // begin block offset
            uint16_t bt = m_bt[bb_idx];
            size_type sample_pos = bb_idx/t_k;
            size_type eb_idx = (idx+len-1)/t_bs; // end block index
            if (bb_idx == eb_idx) {  // extract only in one block
                if (m_invert[sample_pos])
                    bt = t_bs - bt;
                if (bt == 0) {   // all bits are zero
                    res = 0;
                } else if (bt == t_bs and t_bs <= 64) { // all bits are zero
                    res = bits::lo_set[len];
                } else {
                    size_type btnrp = m_btnrp[ sample_pos ];
                    for (size_type j = sample_pos*t_k; j < bb_idx; ++j) {
                        btnrp += rrr_helper_type::space_for_bt(m_bt[j]);
                    }
                    uint16_t btnrlen = rrr_helper_type::space_for_bt(bt);
                    number_type btnr = rrr_helper_type::decode_btnr(m_btnr, btnrp, btnrlen);
                    res =  rrr_helper_type::decode_int(bt, btnr, bb_off, len);
                }
            } else { // solve multiple block case by recursion
                uint16_t b_len = t_bs-bb_off; // remaining bits in first block
                uint16_t b_len_sum = 0;
                do {
                    res |= get_int(idx, b_len) << b_len_sum;
                    idx += b_len;
                    b_len_sum += b_len;
                    len -= b_len;
                    b_len = t_bs;
                    b_len = std::min((uint16_t)len, b_len);
                } while (len > 0);
            }
            return res;
        }

        //! Assignment operator
        rrr_vector& operator=(const rrr_vector& rrr) {
            if (this != &rrr) {
                copy(rrr);
            }
            return *this;
        }

        //! Move assignment operator
        rrr_vector& operator=(rrr_vector&& rrr) {
            swap(rrr);
            return *this;
        }

        //! Returns the size of the original bit vector.
        size_type size()const {
            return m_size;
        }

        //! Answers select queries
        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += m_bt.serialize(out, child, "bt");
            written_bytes += m_btnr.serialize(out, child, "btnr");
            written_bytes += m_btnrp.serialize(out, child, "btnrp");
            written_bytes += m_rank.serialize(out, child, "rank_samples");
            written_bytes += m_invert.serialize(out, child, "invert");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            read_member(m_size, in);
            m_bt.load(in);
            m_btnr.load(in);
            m_btnrp.load(in);
            m_rank.load(in);
            m_invert.load(in);
        }

        iterator begin() const {
            return iterator(this, 0);
        }

        iterator end() const {
            return iterator(this, size());
        }
};

template<uint8_t t_bit_pattern>
struct rank_support_rrr_trait {
    typedef bit_vector::size_type size_type;
    static size_type adjust_rank(size_type r, SDSL_UNUSED size_type n) {
        return r;
    }
};

template<>
struct rank_support_rrr_trait<0> {
    typedef bit_vector::size_type size_type;
    static size_type adjust_rank(size_type r, size_type n) {
        return n - r;
    }
};

//! rank_support for the rrr_vector class
/*!
* \tparam t_b   The bit pattern of size one. (so `0` or `1`)
* \tparam t_bs  The block size of the corresponding rrr_vector
* \tparam t_rac Type used to store the block type in the corresponding rrr_vector.
*  TODO: Test if the binary search can be speed up by
*        saving the (n/2)-th rank value in T[0], the (n/4)-th in T[1],
*        the (3n/4)-th in T[2],... for small number of rank values
*    is this called hinted binary search???
*    or is this called
*/
template<uint8_t t_b, uint16_t t_bs, class t_rac, uint16_t t_k>
class rank_support_rrr
{
        static_assert(t_b == 1u or t_b == 0u , "rank_support_rrr: bit pattern must be `0` or `1`");
    public:
        typedef rrr_vector<t_bs, t_rac, t_k> bit_vector_type;
        typedef typename bit_vector_type::size_type size_type;
        typedef typename bit_vector_type::rrr_helper_type rrr_helper_type;
        typedef typename rrr_helper_type::number_type number_type;
        enum { bit_pat = t_b };

    private:
        const bit_vector_type* m_v; //!< Pointer to the rank supported rrr_vector

    public:
        //! Standard constructor
        /*! \param v Pointer to the rrr_vector, which should be supported
         */
        explicit rank_support_rrr(const bit_vector_type* v=nullptr) {
            set_vector(v);
        }

        //! Answers rank queries
        /*! \param i Argument for the length of the prefix v[0..i-1], with \f$0\leq i \leq size()\f$.
           \returns Number of 1-bits in the prefix [0..i-1] of the original bit_vector.
           \par Time complexity
                \f$ \Order{ sample\_rate of the rrr\_vector} \f$
        */
        const size_type rank(size_type i)const {
            assert(m_v != nullptr);
            assert(i <= m_v->size());
            size_type bt_idx = i/t_bs;
            size_type sample_pos = bt_idx/t_k;
            size_type btnrp = m_v->m_btnrp[ sample_pos ];
            size_type rank  = m_v->m_rank[ sample_pos ];
            if (sample_pos+1 < m_v->m_rank.size()) {
                size_type diff_rank  = m_v->m_rank[ sample_pos+1 ] - rank;
#ifndef RRR_NO_OPT
                if (diff_rank == (size_type)0) {
                    return  rank_support_rrr_trait<t_b>::adjust_rank(rank, i);
                } else if (diff_rank == (size_type)t_bs*t_k) {
                    return  rank_support_rrr_trait<t_b>::adjust_rank(
                                rank + i - sample_pos*t_k*t_bs, i);
                }
#endif
            }
            const bool inv = m_v->m_invert[ sample_pos ];
            for (size_type j = sample_pos*t_k; j < bt_idx; ++j) {
                uint16_t r = m_v->m_bt[j];
                rank  += (inv ? t_bs - r: r);
                btnrp += rrr_helper_type::space_for_bt(r);
            }
            uint16_t off = i % t_bs;
            if (!off) {   // needed for special case: if i=size() is a multiple of t_bs
                // the access to m_bt would cause a invalid memory access
                return rank_support_rrr_trait<t_b>::adjust_rank(rank, i);
            }
            uint16_t bt = inv ? t_bs - m_v->m_bt[ bt_idx ] : m_v->m_bt[ bt_idx ];

            uint16_t btnrlen = rrr_helper_type::space_for_bt(bt);
            number_type btnr = rrr_helper_type::decode_btnr(m_v->m_btnr, btnrp, btnrlen);
            uint16_t popcnt  = rrr_helper_type::decode_popcount(bt, btnr, off);
            return rank_support_rrr_trait<t_b>::adjust_rank(rank + popcnt, i);
        }

        //! Short hand for rank(i)
        const size_type operator()(size_type i)const {
            return rank(i);
        }

        //! Returns the size of the original vector
        const size_type size()const {
            return m_v->size();
        }

        //! Set the supported vector.
        void set_vector(const bit_vector_type* v=nullptr) {
            m_v = v;
        }

        rank_support_rrr& operator=(const rank_support_rrr& rs) {
            if (this != &rs) {
                set_vector(rs.m_v);
            }
            return *this;
        }

        void swap(rank_support_rrr&) { }

        //! Load the data structure from a stream and set the supported vector.
        void load(std::istream&, const bit_vector_type* v=nullptr) {
            set_vector(v);
        }

        //! Serializes the data structure into a stream.
        size_type serialize(std::ostream&, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            structure_tree::add_size(child, 0);
            return 0;
        }
};


//! Select support for the rrr_vector class.
/*
* \tparam t_b   The bit pattern of size one. (so `0` or `1`)
* \tparam t_bs  The block size of the corresponding rrr_vector
* \tparam t_rac Type used to store the block type in the corresponding rrr_vector.
*
* Possible TODO: Add heap which contains the 10 first items of
* each binary search could increase performance.
* Experiments on select_support_interleaved showed about
* 25%.
*/
template<uint8_t t_b, uint16_t t_bs, class t_rac, uint16_t t_k>
class select_support_rrr
{
        static_assert(t_b == 1u or t_b == 0u , "select_support_rrr: bit pattern must be `0` or `1`");
    public:
        typedef rrr_vector<t_bs, t_rac, t_k> bit_vector_type;
        typedef typename bit_vector_type::size_type size_type;
        typedef typename bit_vector_type::rrr_helper_type rrr_helper_type;
        typedef typename rrr_helper_type::number_type number_type;
        enum { bit_pat = t_b };
    private:
        const bit_vector_type* m_v; //!< Pointer to the rank supported rrr_vector

        size_type select1(size_type i)const {
            if (m_v->m_rank[m_v->m_rank.size()-1] < i)
                return size();
            //  (1) binary search for the answer in the rank_samples
            size_type begin=0, end=m_v->m_rank.size()-1; // min included, max excluded
            size_type idx, rank;
            // invariant:  m_rank[end]   >= i
            //             m_rank[begin]  < i
            while (end-begin > 1) {
                idx  = (begin+end) >> 1; // idx in [0..m_rank.size()-1]
                rank = m_v->m_rank[idx];
                if (rank >= i)
                    end = idx;
                else { // rank < i
                    begin = idx;
                }
            }
            //   (2) linear search between the samples
            rank = m_v->m_rank[begin]; // now i>rank
            idx = begin * t_k; // initialize idx for select result
            size_type diff_rank  = m_v->m_rank[end] - rank;
#ifndef RRR_NO_OPT
            if (diff_rank == (size_type)t_bs*t_k) {// optimisation for select<1>
                return idx*t_bs + i-rank -1;
            }
#endif
            const bool inv = m_v->m_invert[ begin ];
            size_type btnrp = m_v->m_btnrp[ begin ];
            uint16_t bt = 0, btnrlen = 0; // temp variables for block_type and space for block type
            while (i > rank) {
                bt = m_v->m_bt[idx++]; bt = inv ? t_bs-bt : bt;
                rank += bt;
                btnrp += (btnrlen=rrr_helper_type::space_for_bt(bt));
            }
            rank -= bt;
            number_type btnr = rrr_helper_type::decode_btnr(m_v->m_btnr, btnrp-btnrlen, btnrlen);
            return (idx-1) * t_bs + rrr_helper_type::decode_select(bt, btnr, i-rank);
        }

        size_type select0(size_type i)const {
            if ((size() - m_v->m_rank[m_v->m_rank.size()-1]) < i) {
                return size();
            }
            //  (1) binary search for the answer in the rank_samples
            size_type begin=0, end=m_v->m_rank.size()-1; // min included, max excluded
            size_type idx, rank;
            // invariant:  m_rank[end] >= i
            //             m_rank[begin] < i
            while (end-begin > 1) {
                idx  = (begin+end) >> 1; // idx in [0..m_rank.size()-1]
                rank = idx*t_bs*t_k - m_v->m_rank[idx];
                if (rank >= i)
                    end = idx;
                else { // rank < i
                    begin = idx;
                }
            }
            //   (2) linear search between the samples
            rank = begin*t_bs*t_k - m_v->m_rank[begin]; // now i>rank
            idx = begin * t_k; // initialize idx for select result
            if (m_v->m_rank[end] == m_v->m_rank[begin]) {      // only for select<0>
                return idx*t_bs +  i-rank -1;
            }
            const bool inv = m_v->m_invert[ begin ];
            size_type btnrp = m_v->m_btnrp[ begin ];
            uint16_t bt = 0, btnrlen = 0; // temp variables for block_type and space for block type
            while (i > rank) {
                bt = m_v->m_bt[idx++]; bt = inv ? t_bs-bt : bt;
                rank += (t_bs-bt);
                btnrp += (btnrlen=rrr_helper_type::space_for_bt(bt));
            }
            rank -= (t_bs-bt);
            number_type btnr = rrr_helper_type::decode_btnr(m_v->m_btnr, btnrp-btnrlen, btnrlen);
            return (idx-1) * t_bs + rrr_helper_type::template decode_select_bitpattern<0, 1>(bt, btnr, i-rank);
        }



    public:
        explicit select_support_rrr(const bit_vector_type* v=nullptr) {
            set_vector(v);
        }

        //! Answers select queries
        size_type select(size_type i)const {
            return  t_b ? select1(i) : select0(i);
        }

        const size_type operator()(size_type i)const {
            return select(i);
        }

        const size_type size()const {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr) {
            m_v = v;
        }

        select_support_rrr& operator=(const select_support_rrr& rs) {
            if (this != &rs) {
                set_vector(rs.m_v);
            }
            return *this;
        }

        void swap(select_support_rrr&) { }

        void load(std::istream&, const bit_vector_type* v=nullptr) {
            set_vector(v);
        }

        size_type serialize(std::ostream&, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            structure_tree::add_size(child, 0);
            return 0;
        }
};

}// end namespace sdsl
#include "rrr_vector_15.hpp" // include specialization

#endif
