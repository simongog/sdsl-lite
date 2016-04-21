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
   \brief rrr_vector.hpp contains a specialisation of the sdsl::rrr_vector class,
          with block size k=15 and lookup table access.
   \author Simon Gog
*/
#ifndef INCLUDED_SDSL_RRR_VECTOR_15
#define INCLUDED_SDSL_RRR_VECTOR_15

#include "int_vector.hpp"
#include "util.hpp"
#include "rrr_helper.hpp" // for binomial helper class
#include "rrr_vector.hpp"
#include "iterators.hpp"
#include <vector>
#include <algorithm> // for next_permutation
#include <iostream>

//! Namespace for the succinct data structure library
namespace sdsl
{

// Helper class for the binomial coefficients \f$ 15 \choose k \f$
/*
 * Size of lookup tables:
 *  * m_nr_to_bin: 64 kB = (2^15 entries x 2 bytes)
 *  * m_bin_to_nr: 64 kB = (2^15 entries x 2 bytes)
 */
class binomial15
{
    public:
        typedef uint32_t number_type;
    private:

        static class impl
        {
            public:
                static const int n = 15;
                static const int MAX_SIZE=32;
                uint8_t m_space_for_bt[16];
                uint8_t m_space_for_bt_pair[256];
                uint64_t m_C[MAX_SIZE];
                int_vector<16> m_nr_to_bin;
                int_vector<16> m_bin_to_nr;

                impl()
                {
                    m_nr_to_bin.resize(1<<n);
                    m_bin_to_nr.resize(1<<n);
                    for (int i=0, cnt=0, class_cnt=0; i<=n; ++i) {
                        m_C[i] = cnt;
                        class_cnt = 0;
                        std::vector<bool> b(n,0);
                        for (int j=0; j<i; ++j) b[n-j-1] = 1;
                        do {
                            uint32_t x=0;
                            for (int k=0; k<n; ++k)
                                x |= ((uint32_t)b[n-k-1])<<(n-1-k);
                            m_nr_to_bin[cnt] = x;
                            m_bin_to_nr[x] = class_cnt;
                            ++cnt;
                            ++class_cnt;
                        } while (next_permutation(b.begin(), b.end()));
                        if (class_cnt == 1)
                            m_space_for_bt[i] = 0;
                        else
                            m_space_for_bt[i] = bits::hi(class_cnt)+1;
                    }
                    if (n == 15) {
                        for (int x=0; x<256; ++x) {
                            m_space_for_bt_pair[x] = m_space_for_bt[x>>4] + m_space_for_bt[x&0x0F];
                        }
                    }
                }
        } iii;

    public:

        static inline uint8_t space_for_bt(uint32_t i)
        {
            return iii.m_space_for_bt[i];
        }

        static inline uint32_t nr_to_bin(uint8_t k, uint32_t nr)
        {
            return iii.m_nr_to_bin[iii.m_C[k]+nr];
        }

        static inline uint32_t bin_to_nr(uint32_t bin)
        {
            return iii.m_bin_to_nr[bin];
        }

        static inline uint8_t space_for_bt_pair(uint8_t x)
        {
            return iii.m_space_for_bt_pair[x];
        }
};


//! A specialization of the rrr_vector class for a block_size of 15.
/*!
 *   \tparam t_rac  Random access integer vector. Use to store the block types.
 *
 *  Several tricks were used to speed-up the operations:
 *  * Whenever possible 2 4-bit blocks are decoded at once.
 *  * When the rank position lies in a block which consists only of zeros or
 *    ones (a uniform block), then we only have to sum up the values of the
 *    block type array between the last sampled position and the
 *    destination block. That can be done by using bit-parallelism on
 *    64-bit words.
*/
template<class t_rac, uint16_t t_k>
class rrr_vector<15, t_rac, t_k>
{
    public:
        typedef bit_vector::size_type                    size_type;
        typedef bit_vector::value_type                   value_type;
        typedef bit_vector::difference_type              difference_type;
        typedef t_rac                                    rac_type;
        typedef random_access_const_iterator<rrr_vector> iterator;
        typedef bv_tag                                   index_category;

        friend class rank_support_rrr<0, 15, t_rac, t_k>;
        friend class rank_support_rrr<1, 15, t_rac, t_k>;
        friend class select_support_rrr<0, 15, t_rac, t_k>;
        friend class select_support_rrr<1, 15, t_rac, t_k>;

        typedef rank_support_rrr<1, 15, t_rac, t_k> rank_1_type;
        typedef rank_support_rrr<0, 15, t_rac, t_k> rank_0_type;
        typedef select_support_rrr<1, 15, t_rac, t_k> select_1_type;
        typedef select_support_rrr<0, 15, t_rac, t_k> select_0_type;

        enum { block_size = 15 };
        typedef binomial15 bi_type;
    private:
        size_type    m_size = 0; // Size of the original bit_vector.
        rac_type     m_bt;       // Vector for block types (bt). bt equals the
        // number of set bits in the block.
        bit_vector   m_btnr;     // Compressed block type numbers.
        int_vector<> m_btnrp;    // Sample pointers into m_btnr.
        int_vector<> m_rank;     // Sample rank values.

        void copy(const rrr_vector& rrr)
        {
            m_size = rrr.m_size;
            m_bt = rrr.m_bt;
            m_btnr = rrr.m_btnr;
            m_btnrp = rrr.m_btnrp;
            m_rank = rrr.m_rank;
        }
    public:
        const rac_type& bt     = m_bt;
        const bit_vector& btnr = m_btnr;

        //! Default constructor
        /*! \param k Store rank samples and pointers each k-th blocks.
         */
        rrr_vector() {};

        //! Copy constructor
        rrr_vector(const rrr_vector& rrr)
        {
            copy(rrr);
        }

        //! Move constructor
        rrr_vector(rrr_vector&& rrr) : m_size(std::move(rrr.m_size)),
            m_bt(std::move(rrr.m_bt)),
            m_btnr(std::move(rrr.m_btnr)), m_btnrp(std::move(rrr.m_btnrp)),
            m_rank(std::move(rrr.m_rank)) {}

        //! Constructor
        /*!
        *  \param bv Uncompressed bitvector.
        *  \param k  Store rank samples and pointers each k-th blocks.
        */
        rrr_vector(const bit_vector& bv)
        {
            m_size = bv.size();
            int_vector<> bt_array;
            bt_array = int_vector<>(m_size/block_size+1, 0, bits::hi(block_size)+1);

            // (1) calculate the block types and store them in m_bt
            size_type pos = 0, i = 0, x;
            size_type btnr_pos = 0;
            size_type sum_rank = 0;
            while (pos + block_size <= m_size) { // handle all full blocks
                bt_array[ i++ ] = x = bits::cnt(bv.get_int(pos, block_size));
                sum_rank += x;
                btnr_pos += bi_type::space_for_bt(x);
                pos += block_size;
            }
            if (pos < m_size) { // handle last full block
                bt_array[ i++ ] = x = bits::cnt(bv.get_int(pos, m_size - pos));
                sum_rank += x;
                btnr_pos += bi_type::space_for_bt(x);
            }
            m_btnr  = bit_vector(std::max(btnr_pos, (size_type)64), 0); // max necessary for case: block_size == 1
            m_btnrp = int_vector<>((bt_array.size()+t_k-1)/t_k, 0,  bits::hi(btnr_pos)+1);

            m_rank  = int_vector<>((bt_array.size()+t_k-1)/t_k + ((m_size % (t_k*block_size))>0), 0, bits::hi(sum_rank)+1);
            //                                                                                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            //                                                                      only add a finishing block, if the last block of the superblock is not a dummy block
            // (2) calculate block type numbers and pointers into btnr and rank samples
            pos = 0; i = 0;
            btnr_pos= 0, sum_rank = 0;
            while (pos + block_size <= m_size) { // handle all full blocks
                if ((i % t_k) == 0) {
                    m_btnrp[ i/t_k ] = btnr_pos;
                    m_rank[ i/t_k ] = sum_rank;
                }
                uint16_t space_for_bt = bi_type::space_for_bt(x=bt_array[i++]);
                sum_rank += x;
                if (space_for_bt) {
                    m_btnr.set_int(btnr_pos, bi_type::bin_to_nr(bv.get_int(pos, block_size)), space_for_bt);
                }
                btnr_pos += space_for_bt;
                pos += block_size;
            }
            if (pos < m_size) { // handle last not full block
                if ((i % t_k) == 0) {
                    m_btnrp[ i/t_k ] = btnr_pos;
                    m_rank[ i/t_k ] = sum_rank;
                }
                uint16_t space_for_bt = bi_type::space_for_bt(x=bt_array[i++]);
                sum_rank += x;
                if (space_for_bt) {
                    m_btnr.set_int(btnr_pos, bi_type::bin_to_nr(bv.get_int(pos, m_size - pos)), space_for_bt);
                }
                btnr_pos += space_for_bt;
                assert(m_rank.size()-1 == ((i+t_k-1)/t_k));
            } else { // handle last empty full block
                assert(m_rank.size()-1 == ((i+t_k-1)/t_k));
            }
            // for technical reasons add an additional element to m_rank
            m_rank[ m_rank.size()-1 ] = sum_rank; // sum_rank contains the total number of set bits in bv
            m_bt = rac_type(std::move(bt_array));
        }

        //! Swap method
        void swap(rrr_vector& rrr)
        {
            if (this != &rrr) {
                std::swap(m_size, rrr.m_size);
                m_bt.swap(rrr.m_bt);
                m_btnr.swap(rrr.m_btnr);
                m_btnrp.swap(rrr.m_btnrp);
                m_rank.swap(rrr.m_rank);
            }
        }

        //! Accessing the i-th element of the original bit_vector
        /*! \param i An index i with \f$ 0 \leq i < size()  \f$.
           \return The i-th bit of the original bit_vector
        */
        value_type operator[](size_type i)const
        {
            size_type bt_idx = i/block_size ;
            uint8_t* bt = (uint8_t*)(m_bt.data());
            uint32_t i_bt = *(bt + (bt_idx/2));
            if (bt_idx%2 == 1) {
                i_bt >>= 4;
            } else {
                i_bt &= 0x0F;
            }
            if (i_bt == 0 or i_bt == block_size) {
                return i_bt > 0;
            }
            size_type sample_pos = bt_idx/t_k;
            size_type btnrp = m_btnrp[ sample_pos ];
            size_type j = (sample_pos*t_k);
            bt += j/2;
            if (j%2 == 1 and j < bt_idx) {
                btnrp += bi_type::space_for_bt((*bt++)>>4);
                ++j;
            }
            while (j+1 < bt_idx) {
                btnrp += bi_type::space_for_bt_pair(*(bt++));   // decode two entries at once
                j+=2;
            }
            if (j < bt_idx) {
                btnrp += bi_type::space_for_bt((*bt)&0x0F);
            }

            uint32_t btnr = m_btnr.get_int(btnrp, bi_type::space_for_bt(i_bt));

            uint8_t off = (uint8_t)(i % block_size); //i - bt_idx*block_size;
            return (bi_type::nr_to_bin(i_bt, btnr) >> off) & (uint32_t)1;
        }

        //! Get the integer value of the binary string of length len starting at position idx.
        /*! \param idx Starting index of the binary representation of the integer.
         *  \param len Length of the binary representation of the integer. Default value is 64.
         *   \returns The integer value of the binary string of length len starting at position idx.
         *
         *  \pre idx+len-1 in [0..size()-1]
         *  \pre len in [1..64]
         */
        uint64_t get_int(size_type idx, uint8_t len=64)const
        {
            uint64_t res = 0;
            size_type bb_idx = idx/block_size; // begin block index
            size_type bb_off = idx%block_size; // begin block offset
            uint16_t bt = m_bt[bb_idx];
            size_type sample_pos = bb_idx/t_k;
            size_type eb_idx = (idx+len-1)/block_size; // end block index
            if (bb_idx == eb_idx) {  // extract only in one block
                if (bt == 0) {   // all bits are zero
                    res = 0;
                } else if (bt == block_size) {  // all bits are zero
                    res = bits::lo_set[len];
                } else {
                    size_type btnrp = m_btnrp[ sample_pos ];
                    for (size_type j = sample_pos*t_k; j < bb_idx; ++j) {
                        btnrp += bi_type::space_for_bt(m_bt[j]);
                    }
                    uint32_t btnr = m_btnr.get_int(btnrp, bi_type::space_for_bt(bt));
                    res = (bi_type::nr_to_bin(bt, btnr) >> bb_off) & bits::lo_set[len];
                }
            } else { // solve multiple block case by recursion
                uint8_t b_len = block_size-bb_off;
                uint8_t b_len_sum = 0;
                do {
                    res |= get_int(idx, b_len) << b_len_sum;
                    idx += b_len;
                    b_len_sum += b_len;
                    len -= b_len;
                    b_len = block_size;
                    b_len = std::min(len, b_len);
                } while (len > 0);
            }
            return res;
        }


        //! Assignment operator
        rrr_vector& operator=(const rrr_vector& rrr)
        {
            if (this != &rrr) {
                copy(rrr);
            }
            return *this;
        }

        //! Move assignment
        rrr_vector& operator=(rrr_vector&& rrr)
        {
            swap(rrr);
            return *this;
        }

        //! Returns the size of the original bit vector.
        size_type size()const
        {
            return m_size;
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            size_type written_bytes = 0;
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += m_bt.serialize(out, child, "bt");
            written_bytes += m_btnr.serialize(out, child, "btnr");
            written_bytes += m_btnrp.serialize(out, child, "btnrp");
            written_bytes += m_rank.serialize(out, child, "rank_samples");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in)
        {
            read_member(m_size, in);
            m_bt.load(in);
            m_btnr.load(in);
            m_btnrp.load(in);
            m_rank.load(in);
        }

        iterator begin() const
        {
            return iterator(this, 0);
        }

        iterator end() const
        {
            return iterator(this, size());
        }
};


//! rank_support for the specialized rrr_vector class of block size 15.
/*! The first template parameter is the bit pattern of size one.
*/
template<uint8_t t_b, class t_rac, uint16_t t_k>
class rank_support_rrr<t_b, 15, t_rac, t_k>
{
        static_assert(t_b == 1u or t_b == 0u , "rank_support_rrr: bit pattern must be `0` or `1`");
    public:
        typedef rrr_vector<15, t_rac, t_k> bit_vector_type;
        typedef typename bit_vector_type::size_type size_type;
        typedef typename bit_vector_type::bi_type bi_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = (uint8_t)1 };

    private:
        const bit_vector_type* m_v; //!< Pointer to the rank supported rrr_vector
        // TODO cache for sequential ranks
//        mutable size_type m_last_bt;
//        mutable size_type m_last_w; // store the last decoded word
//        uint16_t m_space_for_bt[256];
    public:
        //! Standard constructor
        /*! \param v Pointer to the rrr_vector, which should be supported
         */
        explicit rank_support_rrr(const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        //! Answers rank queries
        /*! \param i Argument for the length of the prefix v[0..i-1], with \f$0\leq i \leq size()\f$.
           \returns Number of 1-bits in the prefix [0..i-1] of the original bit_vector.
           \par Time complexity
                \f$ \Order{ sample\_rate of the rrr\_vector} \f$
        */
        const size_type rank(size_type i)const
        {
            size_type bt_idx = i/bit_vector_type::block_size;
            size_type sample_pos = bt_idx/t_k;
            size_type btnrp = m_v->m_btnrp[ sample_pos ];
            size_type rank  = m_v->m_rank[ sample_pos ];
            if (sample_pos+1 < m_v->m_rank.size()) {
                size_type diff_rank  = m_v->m_rank[ sample_pos+1 ] - rank;
                if (diff_rank == 0) {
                    return  rank_support_rrr_trait<t_b>::adjust_rank(rank, i);
                } else if (diff_rank == (size_type)bit_vector_type::block_size*t_k) {
                    return  rank_support_rrr_trait<t_b>::adjust_rank(
                                rank + i - sample_pos*t_k*bit_vector_type::block_size, i);
                }
            }
            uint8_t* bt = (uint8_t*)(m_v->m_bt.data());

            uint8_t last_bt = *(bt + (bt_idx/2));
            if (bt_idx%2 == 1) {
                last_bt >>= 4;
            } else {
                last_bt &= 0x0F;
            }
            // if the final block type consists only of ones or zeros, we don't have to
            // calculate the position pointer and can sum up rank in 64 bit chunks
            if (last_bt == 0 or last_bt == 15) {
                if (last_bt == 15)
                    rank += i % bit_vector_type::block_size;
                size_type j = (sample_pos*t_k) << 2;
                bt_idx = bt_idx << 2;
                if (bt_idx == j)
                    return rank_support_rrr_trait<t_b>::adjust_rank(rank, i);
                // now j < bt_idx
                const uint64_t* bt64 = m_v->m_bt.data() + (j >> 6); // get the word of the start
                uint8_t bt64_off = j & 0x3F; // get the offset in the word of the start
                const uint64_t* bt64_end = m_v->m_bt.data() + (bt_idx >> 6);
                uint8_t bt64_end_off = bt_idx & 0x3F;
                // Case (1)
                if (bt64 == bt64_end) {
                    uint64_t w = ((*bt64) >> bt64_off) & bits::lo_set[bt64_end_off-bt64_off];
                    w = (w & 0x0f0f0f0f0f0f0f0full) + ((w >> 4) & 0x0f0f0f0f0f0f0f0full);
                    rank += ((0x0101010101010101ull*w) >> 56);
                } else { // Case (2)
                    uint64_t w = ((*bt64) >> bt64_off);
                    w = (w & 0x0f0f0f0f0f0f0f0full) + ((w >> 4) & 0x0f0f0f0f0f0f0f0full);
                    rank += ((0x0101010101010101ull*w) >> 56);
                    while ((++bt64) != bt64_end) {
                        w = *bt64;
                        w = (w & 0x0f0f0f0f0f0f0f0full) + ((w >> 4) & 0x0f0f0f0f0f0f0f0full);
                        rank += ((0x0101010101010101ull*w) >> 56);
                    }
                    // now bt64 == bt64_end
                    if (bt64_end_off > 0) {
                        w = *bt64 << (64 - bt64_end_off);
                        w = (w & 0x0f0f0f0f0f0f0f0full) + ((w >> 4) & 0x0f0f0f0f0f0f0f0full);
                        rank += ((0x0101010101010101ull*w) >> 56);
                    }
                }
                return rank_support_rrr_trait<t_b>::adjust_rank(rank, i);  // necessary
            }
            size_type j = sample_pos*t_k;
            bt += j/2;
            if (j%2 == 1 and j < bt_idx) {
                const uint8_t r = (*bt++)>>4;
                rank += r;
                btnrp += bi_type::space_for_bt(r);
                ++j;
            }
            while (j+1 < bt_idx) {
                const uint8_t r = *(bt++);
                rank += (r>>4)+(r&0x0F);
                btnrp += bi_type::space_for_bt_pair(r);   // decode two entries at once
                j+=2;
            }
            if (j < bt_idx) {
                const uint8_t r = (*bt);
                rank += r&0x0F;
                btnrp += bi_type::space_for_bt(r&0x0F);
                ++j;
            }
            uint8_t off = i % bit_vector_type::block_size; //i - bt_idx*bit_vector_type::block_size;
            if (!off) {   // needed for special case: if i=size() is a multiple of block_size
                // the access to m_bt would cause a invalid memory access
                return rank_support_rrr_trait<t_b>::adjust_rank(rank, i);
            }
            uint32_t btnr = m_v->m_btnr.get_int(btnrp, bi_type::space_for_bt(last_bt));
            return rank_support_rrr_trait<t_b>::adjust_rank(rank +
                    bits::cnt(((uint64_t)(bi_type::nr_to_bin(last_bt, btnr))) & bits::lo_set[off]), i);
        }

        //! Short hand for rank(i)
        const size_type operator()(size_type i)const
        {
            return rank(i);
        }

        //! Returns the size of the original vector
        const size_type size()const
        {
            return m_v->size();
        }

        //! Set the supported vector.
        void set_vector(const bit_vector_type* v=nullptr)
        {
            m_v = v;
        }

        rank_support_rrr& operator=(const rank_support_rrr& rs)
        {
            if (this != &rs) {
                set_vector(rs.m_v);
            }
            return *this;
        }

        void swap(rank_support_rrr&) { }

        //! Load the data structure from a stream and set the supported vector.
        void load(std::istream&, const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        //! Serializes the data structure into a stream.
        size_type serialize(std::ostream&, structure_tree_node* v=nullptr, std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            structure_tree::add_size(child, 0);
            return 0;
        }
};


//! Select support for the specialized rrr_vector class of block size 15.
template<uint8_t t_b, class t_rac, uint16_t t_k>
class select_support_rrr<t_b, 15, t_rac, t_k>
{
        static_assert(t_b == 1u or t_b == 0u , "select_support_rrr: bit pattern must be `0` or `1`");
    public:
        typedef rrr_vector<15, t_rac, t_k>          bit_vector_type;
        typedef typename bit_vector_type::size_type size_type;
        typedef typename bit_vector_type::bi_type     bi_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = (uint8_t)1 };
    private:
        const bit_vector_type* m_v; //!< Pointer to the rank supported rrr_vector

        // TODO: hinted binary search
        size_type  select1(size_type i)const
        {
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
            if (diff_rank == (size_type)bit_vector_type::block_size*t_k) {// optimisation for select<1>
                return idx*bit_vector_type::block_size + i-rank -1;
            }
            size_type btnrp = m_v->m_btnrp[ begin ];
            uint8_t bt = 0, s = 0; // temp variables for block_type and space for block type
            while (i > rank) {
                bt = m_v->m_bt[idx++];
                rank += bt;
                btnrp += (s=bi_type::space_for_bt(bt));
            }
            rank -= bt;
            uint32_t btnr = m_v->m_btnr.get_int(btnrp-s, s);
            return (idx-1) * bit_vector_type::block_size + bits::sel(bi_type::nr_to_bin(bt, btnr), i-rank);
        }

        // TODO: hinted binary search
        size_type  select0(size_type i)const
        {
            if ((size()-m_v->m_rank[m_v->m_rank.size()-1]) < i)
                return size();
            //  (1) binary search for the answer in the rank_samples
            size_type begin=0, end=m_v->m_rank.size()-1; // min included, max excluded
            size_type idx, rank;
            // invariant:  m_rank[end] >= i
            //             m_rank[begin] < i
            while (end-begin > 1) {
                idx  = (begin+end) >> 1; // idx in [0..m_rank.size()-1]
                rank = idx*bit_vector_type::block_size*t_k - m_v->m_rank[idx];
                if (rank >= i)
                    end = idx;
                else { // rank < i
                    begin = idx;
                }
            }
            //   (2) linear search between the samples
            rank = begin*bit_vector_type::block_size*t_k - m_v->m_rank[begin]; // now i>rank
            idx = begin * t_k; // initialize idx for select result
            if (m_v->m_rank[end] == m_v->m_rank[begin]) {  // only for select<0>
                return idx*bit_vector_type::block_size +  i-rank -1;
            }
            size_type btnrp = m_v->m_btnrp[ begin ];
            uint8_t bt = 0, s = 0; // temp variables for block_type and space for block type
            while (i > rank) {
                bt = m_v->m_bt[idx++];
                rank += (bit_vector_type::block_size-bt);
                btnrp += (s=bi_type::space_for_bt(bt));
            }
            rank -= (bit_vector_type::block_size-bt);
            uint32_t btnr = m_v->m_btnr.get_int(btnrp-s, s);
            return (idx-1) * bit_vector_type::block_size + bits::sel(~((uint64_t)bi_type::nr_to_bin(bt, btnr)), i-rank);
        }



    public:
        select_support_rrr(const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        //! Answers select queries
        size_type select(size_type i)const
        {
            return  t_b ? select1(i) : select0(i);
        }


        const size_type operator()(size_type i)const
        {
            return select(i);
        }

        const size_type size()const
        {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr)
        {
            m_v = v;
        }

        select_support_rrr& operator=(const select_support_rrr& rs)
        {
            if (this != &rs) {
                set_vector(rs.m_v);
            }
            return *this;
        }

        void swap(select_support_rrr&) { }

        void load(std::istream&, const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        size_type serialize(std::ostream&, structure_tree_node* v=nullptr, std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            structure_tree::add_size(child, 0);
            return 0;
        }
};

}// end namespace sdsl

#endif
