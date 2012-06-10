/* sdsl - succinct data structures library
    Copyright (C) 2011 Simon Gog

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
          with block_size=15 and lookup table access.
   \author Simon Gog
*/
#ifndef INCLUDED_SDSL_RRR_VECTOR_15
#define INCLUDED_SDSL_RRR_VECTOR_15

#include "int_vector.hpp"
#include "bitmagic.hpp"
#include "util.hpp"
#include "rrr_helper.hpp" // for binomial helper class
#include "rrr_vector.hpp"
#include <vector>
#include <algorithm> // for next_permutation
#include <iostream>

//! Namespace for the succinct data structure library
namespace sdsl
{

// Helper class
template<uint8_t n>
class binomial
{
    public:
        typedef uint32_t number_type;
    private:

        static class impl
        {
            public:
                static const int MAX_SIZE=32;
                uint8_t m_space_for_bt[n+1];
                uint8_t m_space_for_bt_pair[256*(n==15)];
                uint64_t m_C[MAX_SIZE];
                int_vector<32> m_nr_to_bin;
                int_vector<32> m_bin_to_nr;

                impl() {
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
                            m_space_for_bt[i] = bit_magic::l1BP(class_cnt)+1;
                        //			cout << cnt << " " << class_cnt << " " << (int)m_space_for_bt[i] << endl;
                    }
                    if (n == 15) {
                        for (int x=0; x<256; ++x) {
                            m_space_for_bt_pair[x] = m_space_for_bt[x>>4] + m_space_for_bt[x&0x0F];
                        }
                    }
                }
        } iii;

    public:

        static inline uint8_t space_for_bt(uint32_t i) {
            return iii.m_space_for_bt[i];
        }

        static inline uint32_t nr_to_bin(uint8_t k, uint32_t nr) {
            return iii.m_nr_to_bin[iii.m_C[k]+nr];
        }

        static inline uint32_t bin_to_nr(uint32_t bin) {
            return iii.m_bin_to_nr[bin];
        }

        static inline uint8_t space_for_bt_pair(uint8_t x) {
            return iii.m_space_for_bt_pair[x];
        }

        static inline uint8_t popcount(number_type x) {
            return bit_magic::b1Cnt(x);
        }

        static inline uint8_t select(number_type x, uint32_t i) {
            return bit_magic::i1BP(x, i);
        }

        template<class bit_vector_type>
        static inline number_type decode_btnr(const bit_vector_type& bv,
                                              typename bit_vector_type::size_type btnrp,
                                              uint8_t btnrlen) {
            return bv.get_int(btnrp, btnrlen);
        }

        template<class bit_vector_type>
        static inline uint8_t get_bt(const bit_vector_type& bv,
                                     typename bit_vector_type::size_type pos,
                                     uint8_t block_size) {
            return bit_magic::b1Cnt(bv.get_int(pos, block_size));
        }

        template<class bit_vector_type>
        static void set_bt(bit_vector_type& bv,
                           typename bit_vector_type::size_type pos,
                           number_type bt,
                           uint8_t space_for_bt) {
            bv.set_int(pos, bt, space_for_bt);
        }

        static inline number_type Li1Mask(uint8_t off) {
            return bit_magic::Li1Mask[off];
        }
};
template<uint8_t n>
typename binomial<n>::impl binomial<n>::iii;



//! A specialization of the rrr_vector class for a block_size of 15.
/*!
 * Several tricks were used to speed-up the operations:
 *  * Whenever possible 2 4-bit blocks are decoded at once.
 *  * When the rank position lies in a block which consists only of zeros or
 *    ones (a uniform block), then we only have to sum up the values of the
 *    block type array between the last sampled position and the
 *    destination block. That can be done by using bit-parallelism on
 *    64 bit words.
*/
template<class wt_type>
class rrr_vector<15, wt_type>
{
    public:
        typedef bit_vector::size_type size_type;
        typedef bit_vector::value_type value_type;

        friend class rrr_rank_support<0, 15, wt_type>;
        friend class rrr_rank_support<1, 15, wt_type>;
        friend class rrr_select_support<0, 15, wt_type>;
        friend class rrr_select_support<1, 15, wt_type>;

        typedef rrr_rank_support<1, 15, wt_type> rank_1_type; // typedef for default types for rank and select
        typedef rrr_rank_support<0, 15, wt_type> rank_0_type;
        typedef rrr_select_support<1, 15, wt_type> select_1_type;
        typedef rrr_select_support<0, 15, wt_type> select_0_type;

        enum { block_size = 15 };
        typedef binomial<block_size> bi_type;
    private:
        size_type      m_size; // length of the original bit_vector
        uint16_t       m_sample_rate;
        wt_type        m_bt; // data structure, which stores the block types (bt). The block type equals the number
        // of ones in a block. Another option for this data structure is wt_huff
        bit_vector     m_btnr; // data structure, which stores the block type numbers of the blocks
        int_vector<>   m_btnrp; // sample pointers into btnr
        int_vector<>   m_rank;  // sample rank values

    public:
        const wt_type& bt;
        const bit_vector& btnr;

        void copy(const rrr_vector& rrr) {
            m_size = rrr.m_size;
            m_sample_rate = rrr.m_sample_rate;
            m_bt = rrr.m_bt;
            m_btnr = rrr.m_btnr;
            m_btnrp = rrr.m_btnrp;
            m_rank = rrr.m_rank;
        }

        //! Default constructor
        rrr_vector(uint16_t sample_rate=32):m_size(0), m_sample_rate(sample_rate), bt(m_bt), btnr(m_btnr) {};

        //! Copy constructor
        rrr_vector(const rrr_vector& rrr): bt(m_bt), btnr(m_btnr)  {
            copy(rrr);
        }

        //! Constructor
        /*!
        *  \param sample_rate Insert a sampled block between sample_rate blocks
        */
        rrr_vector(const bit_vector& bv, uint16_t sample_rate=32): m_sample_rate(sample_rate), bt(m_bt), btnr(m_btnr) {
            m_size = bv.size();
            int_vector<> bt_array;
            util::assign(bt_array, int_vector<>(m_size/block_size+1, 0, bit_magic::l1BP(block_size)+1));
//            bt_array.set_int_width(bit_magic::l1BP(block_size)+1);
//            bt_array.resize((m_size+block_size)/block_size);

            // (1) calculate the block types and store them in m_bt
            size_type pos = 0, i = 0, x;
            size_type btnr_pos = 0;
            size_type sum_rank = 0;
            while (pos + block_size <= m_size) { // handle all full blocks
                bt_array[ i++ ] = x = bit_magic::b1Cnt(bv.get_int(pos, block_size));
                sum_rank += x;
                btnr_pos += bi_type::space_for_bt(x);
                pos += block_size;
            }
            if (pos < m_size) { // handle last full block
                bt_array[ i++ ] = x = bit_magic::b1Cnt(bv.get_int(pos, m_size - pos));
                sum_rank += x;
                btnr_pos += bi_type::space_for_bt(x);
            }
//	 cout << "# bt array initialized "<< endl;
            util::assign(m_btnr, bit_vector(std::max(btnr_pos, (size_type)64), 0));      // max necessary for case: block_size == 1
            util::assign(m_btnrp, int_vector<>((bt_array.size()+m_sample_rate-1)/m_sample_rate, 0,  bit_magic::l1BP(btnr_pos)+1));
            util::assign(m_rank, int_vector<>((bt_array.size()+m_sample_rate-1)/m_sample_rate + 1, 0, bit_magic::l1BP(sum_rank)+1));

            // (2) calculate block type numbers and pointers into btnr and rank samples
            pos = 0; i = 0;
            btnr_pos= 0, sum_rank = 0;
            while (pos + block_size <= m_size) { // handle all full blocks
                if ((i % m_sample_rate) == 0) {
                    m_btnrp[ i/m_sample_rate ] = btnr_pos;
                    m_rank[ i/m_sample_rate ] = sum_rank;
                }
                uint16_t space_for_bt = bi_type::space_for_bt(x=bt_array[i++]);
                sum_rank += x;
                if (space_for_bt) {
                    m_btnr.set_int(btnr_pos, bi_type::bin_to_nr(bv.get_int(pos, block_size)), space_for_bt);
                }
                btnr_pos += space_for_bt;
                pos += block_size;
            }
            if (pos < m_size) { // handle last full block
                if ((i % m_sample_rate) == 0) {
                    m_btnrp[ i/m_sample_rate ] = btnr_pos;
                    m_rank[ i/m_sample_rate ] = sum_rank;
                }
                uint16_t space_for_bt = bi_type::space_for_bt(x=bt_array[i++]);
                sum_rank += x;
                if (space_for_bt) {
                    m_btnr.set_int(btnr_pos, bi_type::bin_to_nr(bv.get_int(pos, m_size - pos)), space_for_bt);
                }
                btnr_pos += space_for_bt;
            }
            // for technical reasons add an additional element to m_rank
            m_rank[ m_rank.size()-1 ] = sum_rank; // sum_rank contains the total number of set bits in bv
            util::assign(m_bt, bt_array);
        }

        //! Accessing the i-th element of the original bit_vector
        /*! \param i An index i with \f$ 0 \leq i < size()  \f$.
           \return The i-th bit of the original bit_vector
        */
        value_type operator[](size_type i)const {
            size_type bt_idx = i/block_size ;
            uint8_t* bt = (uint8_t*)(m_bt.data());
            uint8_t i_bt = *(bt + (bt_idx/2));
            if (bt_idx%2 == 1) {
                i_bt >>= 4;
            } else {
                i_bt &= 0x0F;
            }
            if (i_bt == 0 or i_bt == block_size) {
                return i_bt > 0;
            }
            size_type sample_pos = bt_idx/m_sample_rate;
            size_type btnrp = m_btnrp[ sample_pos ];
            size_type j = (sample_pos*m_sample_rate);
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

            uint8_t off = i % block_size; //i - bt_idx*block_size;
            return (bi_type::nr_to_bin(i_bt, btnr) >> off) & (uint32_t)1;
        }

        //! Assignment operator
        rrr_vector& operator=(const rrr_vector& rrr) {
            if (this != &rrr) {
                copy(rrr);
            }
            return *this;
        }

        //! Returns the size of the original bit vector.
        size_type size()const {
            return m_size;
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out)const {
            size_type written_bytes = 0;
            written_bytes += util::write_member(m_size, out);
            written_bytes += util::write_member(m_sample_rate, out);
            written_bytes += m_bt.serialize(out);
            written_bytes += m_btnr.serialize(out);
            written_bytes += m_btnrp.serialize(out);
            written_bytes += m_rank.serialize(out);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            util::read_member(m_size, in);
            util::read_member(m_sample_rate, in);
            m_bt.load(in);
            m_btnr.load(in);
            m_btnrp.load(in);
            m_rank.load(in);
        }

#ifdef MEM_INFO
        void mem_info(std::string label="")const {
            if (label=="")
                label = "rrr_vector";
            size_type bytes = util::get_size_in_bytes(*this);
            std::cout << "list(label=\""<<label<<"\", size = "<< bytes/(1024.0*1024.0) << "\n,";
            m_bt.mem_info("bt"); std::cout << ",\n";
            m_btnr.mem_info("btnr"); std::cout << ",\n";
            m_btnrp.mem_info("btnrp"); std::cout << ",\n";
            m_rank.mem_info("rank samples"); std::cout << ")\n";
        }
#endif

        // Print information about the object to stdout.
        void print_info()const {
            size_type orig_bv_size = m_size; // size of the original bit vector in bits
            size_type rrr_size 	= util::get_size_in_bytes(*this)*8;
            size_type bt_size 	= util::get_size_in_bytes(m_bt)*8;
            size_type btnr_size 	= util::get_size_in_bytes(m_btnr)*8;
            size_type btnrp_and_rank_size = util::get_size_in_bytes(m_btnrp)*8 + util::get_size_in_bytes(m_rank)*8;
            std::cout << "#block_size\tsample_rate\torig_bv_size\trrr_size\tbt_size\tbtnr_size\tbtnrp_and_rank_size" << std::endl;
            std::cout << (int)block_size << "\t" << m_sample_rate << "\t";
            std::cout << orig_bv_size << "\t" << rrr_size << "\t" << bt_size << "\t" << btnr_size << "\t"
                      << btnrp_and_rank_size << std::endl;
        }
};


//! rank_support for the specialized rrr_vector class of block size 15.
/*! The first template parameter is the bit pattern of size one.
 */
template<uint8_t b, class wt_type>
class rrr_rank_support<b, 15, wt_type>
{
    public:
        typedef rrr_vector<15, wt_type> bit_vector_type;
        typedef typename bit_vector_type::size_type size_type;
        typedef typename bit_vector_type::bi_type bi_type;


    private:
        const bit_vector_type* m_v; //!< Pointer to the rank supported rrr_vector
        uint16_t m_sample_rate;  //!<    "     "   "      "
        // TODO cache for sequential ranks
//		mutable size_type m_last_bt;
//		mutable size_type m_last_w; // store the last decoded word
//		uint16_t m_space_for_bt[256];
    public:
        //! Standard constructor
        /*! \param v Pointer to the rrr_vector, which should be supported
         */
        explicit rrr_rank_support(const bit_vector_type* v=NULL) {
            init(v);
        }

        //! Initialize the data structure with a rrr_vector, which should be supported
        void init(const bit_vector_type* v=NULL) {
            set_vector(v);
        }

        //! Answers rank queries
        /*! \param i Argument for the length of the prefix v[0..i-1], with \f$0\leq i \leq size()\f$.
           \returns Number of 1-bits in the prefix [0..i-1] of the original bit_vector.
           \par Time complexity
        		\f$ \Order{ sample\_rate of the rrr\_vector} \f$
        */
        const size_type rank(size_type i)const {
            size_type bt_idx = i/bit_vector_type::block_size;
            size_type sample_pos = bt_idx/m_sample_rate;
            size_type btnrp = m_v->m_btnrp[ sample_pos ];
            size_type rank  = m_v->m_rank[ sample_pos ];
            size_type diff_rank  = m_v->m_rank[ sample_pos+1 ] - rank;
            if (diff_rank == 0) {
                return  rrr_rank_support_trait<b>::adjust_rank(rank, i);
            } else if (diff_rank == bit_vector_type::block_size*m_sample_rate) {
                return  rrr_rank_support_trait<b>::adjust_rank(
                            rank + i - sample_pos*m_sample_rate*bit_vector_type::block_size, i);
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
                size_type j = (sample_pos*m_sample_rate) << 2;
                bt_idx = bt_idx << 2;
                if (bt_idx == j)
                    return rrr_rank_support_trait<b>::adjust_rank(rank, i);
                // now j < bt_idx
                const uint64_t* bt64 = m_v->m_bt.data() + (j >> 6); // get the word of the start
                uint8_t bt64_off = j & 0x3F; // get the offset in the word of the start
                const uint64_t* bt64_end = m_v->m_bt.data() + (bt_idx >> 6);
                uint8_t bt64_end_off = bt_idx & 0x3F;
                // Case (1)
                if (bt64 == bt64_end) {
                    uint64_t w = ((*bt64) >> bt64_off) & bit_magic::Li1Mask[bt64_end_off-bt64_off];
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
                return rrr_rank_support_trait<b>::adjust_rank(rank, i);  // necessary
            }
            size_type j = sample_pos*m_sample_rate;
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
                return rrr_rank_support_trait<b>::adjust_rank(rank, i);
            }
            uint32_t btnr = m_v->m_btnr.get_int(btnrp, bi_type::space_for_bt(last_bt));
            return rrr_rank_support_trait<b>::adjust_rank(rank +
                    bit_magic::b1Cnt(((uint64_t)(bi_type::nr_to_bin(last_bt, btnr))) & bit_magic::Li1Mask[off]), i);
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
        void set_vector(const bit_vector_type* v=NULL) {
            m_v = v;
            if (v != NULL) {
                m_sample_rate = m_v->m_sample_rate;
            } else {
                m_sample_rate = 0;
            }
        }

        rrr_rank_support& operator=(const rrr_rank_support& rs) {
            if (this != &rs) {
                set_vector(rs.m_v);
                m_sample_rate = rs.m_sample_rate;
            }
            return *this;
        }

        void swap(rrr_rank_support& rs) {
            if (this != &rs) {
                std::swap(m_sample_rate, rs.m_sample_rate);
            }
        }

        bool operator==(const rrr_rank_support& rs)const {
            if (this == &rs)
                return true;
            return m_sample_rate == rs.m_sample_rate;
        }

        bool operator!=(const rrr_rank_support& rs)const {
            return !(*this == rs);
        }

        //! Load the data structure from a stream and set the supported vector.
        void load(std::istream& in, const bit_vector_type* v=NULL) {
            util::read_member(m_sample_rate, in);
            set_vector(v);
        }

        //! Serializes the data structure into a stream.
        size_type serialize(std::ostream& out)const {
            size_type written_bytes = 0;
            written_bytes += util::write_member(m_sample_rate, out);
            return written_bytes;
        }

#ifdef MEM_INFO
        void mem_info(std::string label="")const {
            if (label=="")
                label = "rrr_rank_support";
            size_type bytes = util::get_size_in_bytes(*this);
            std::cout << "list(label=\""<<label<<"\", size = "<< bytes/(1024.0*1024.0) << ")\n";
        }
#endif
};


//! Select support for the specialized rrr_vector class of block size 15.
template<uint8_t b,class wt_type>
class rrr_select_support<b, 15, wt_type>
{
    public:
        typedef rrr_vector<15, wt_type> bit_vector_type;
        typedef typename bit_vector_type::size_type size_type;
        typedef typename bit_vector_type::bi_type bi_type;

    private:
        const bit_vector_type* m_v; //!< Pointer to the rank supported rrr_vector
        uint16_t m_sample_rate;  //!<    "     "   "      "

        // TODO: hinted binary search
        size_type  select1(size_type i)const {
            if (m_v->m_rank[m_v->m_rank.size()-1] < i)
                return size();
            //  (1) binary search for the answer in the rank_samples
            size_type begin=0, end=m_v->m_rank.size()-1; // min included, max excluded
            size_type idx, rank;
            // invariant:  m_rank[end] >= i
            //             m_rank[begin] < i
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
            idx = begin * m_sample_rate; // initialize idx for select result
            size_type diff_rank  = m_v->m_rank[end] - rank;
            if (diff_rank == bit_vector_type::block_size*m_sample_rate) {// optimiziation for select<1>
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
            return (idx-1) * bit_vector_type::block_size + bit_magic::i1BP(bi_type::nr_to_bin(bt, btnr), i-rank);
        }

        // TODO: hinted binary search
        size_type  select0(size_type i)const {
            if ((size()-m_v->m_rank[m_v->m_rank.size()-1]) < i)
                return size();
            //  (1) binary search for the answer in the rank_samples
            size_type begin=0, end=m_v->m_rank.size()-1; // min included, max excluded
            size_type idx, rank;
            // invariant:  m_rank[end] >= i
            //             m_rank[begin] < i
            while (end-begin > 1) {
                idx  = (begin+end) >> 1; // idx in [0..m_rank.size()-1]
                rank = idx*bit_vector_type::block_size*m_sample_rate - m_v->m_rank[idx];
                if (rank >= i)
                    end = idx;
                else { // rank < i
                    begin = idx;
                }
            }
            //   (2) linear search between the samples
            rank = begin*bit_vector_type::block_size*m_sample_rate - m_v->m_rank[begin]; // now i>rank
            idx = begin * m_sample_rate; // initialize idx for select result
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
            return (idx-1) * bit_vector_type::block_size + bit_magic::i1BP(~((uint64_t)bi_type::nr_to_bin(bt, btnr)), i-rank);
        }



    public:
        rrr_select_support(const bit_vector_type* v=NULL) {
            init(v);
        }

        void init(const bit_vector_type* v=NULL) {
            set_vector(v);
        }

        //! Answers select queries
        size_type select(size_type i)const {
            return  b ? select1(i) : select0(i);
        }


        const size_type operator()(size_type i)const {
            return select(i);
        }

        const size_type size()const {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=NULL) {
            m_v = v;
            if (v != NULL) {
                m_sample_rate = m_v->m_sample_rate;
            } else {
                m_sample_rate = 0;
            }
        }

        rrr_select_support& operator=(const rrr_select_support& rs) {
            if (this != &rs) {
                set_vector(rs.m_v);
                m_sample_rate = rs.m_sample_rate;
            }
            return *this;
        }

        void swap(rrr_select_support& rs) {
            if (this != &rs) {
                std::swap(m_sample_rate, rs.m_sample_rate);
            }
        }

        bool operator==(const rrr_select_support& rs)const {
            if (this == &rs)
                return true;
            return m_sample_rate == rs.m_sample_rate;
        }

        bool operator!=(const rrr_select_support& rs)const {
            return !(*this == rs);
        }


        void load(std::istream& in, const bit_vector_type* v=NULL) {
            util::read_member(m_sample_rate, in);
            set_vector(v);
        }

        size_type serialize(std::ostream& out)const {
            size_type written_bytes = 0;
            written_bytes += util::write_member(m_sample_rate, out);
            return written_bytes;
        }

#ifdef MEM_INFO
        void mem_info(std::string label="")const {
            if (label=="")
                label = "rrr_select_support";
            size_type bytes = util::get_size_in_bytes(*this);
            std::cout << "list(label=\""<<label<<"\", size = "<< bytes/(1024.0*1024.0) << ")\n";
        }
#endif
};




}// end namespace sdsl

#endif
