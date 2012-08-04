/* sdsl - succinct data structures library
    Copyright (C) 2012 Simon Gog, Matthias Petri

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
/*!\file bit_vector_interleaved.hpp
   \brief bit_vector_interleaved.hpp contains the sdsl::bit_vector_interleaved class, and
          classes which support rank and select for bit_vector_interleaved.
   \author Matthias Petri, Simon Gog
*/
#ifndef SDSL_BIT_VECTOR_INTERLEAVED
#define SDSL_BIT_VECTOR_INTERLEAVED

#include "int_vector.hpp"
#include "bitmagic.hpp"
#include "util.hpp"

#include <queue>

//! Namespace for the succinct data structure library
namespace sdsl
{

template<uint8_t b=1,uint32_t blockSize=512>// forward declaration needed for friend declaration
class rank_support_interleaved;  // in bit_vector_interleaved

template<uint8_t b=1,uint32_t blockSize=512>// forward declaration needed for friend declaration
class select_support_interleaved;  // in bit_vector_interleaved

//! A bit vector which interleaves the original bit_vector with rank information.
/*!
 * This class is a uncompressed bit vector representation. It copies the original
 * bit_vector and interleaves the data every blockSize bits with a cumulative
 * sum of set bits before the current position. Each cumulative sum is stored
 * in a 64 bit word.
 * \pre blockSize has to be a power of 2 and blockSize >= 64
 */
template<uint32_t blockSize=512>
class bit_vector_interleaved
{
    public:
        typedef bit_vector::size_type   size_type;
        typedef size_type               value_type;

        friend class rank_support_interleaved<1,blockSize>;
        friend class rank_support_interleaved<0,blockSize>;
        friend class select_support_interleaved<1,blockSize>;
        friend class select_support_interleaved<0,blockSize>;

        typedef rank_support_interleaved<1,blockSize>     rank_1_type;
        typedef rank_support_interleaved<0,blockSize>     rank_0_type;
        typedef select_support_interleaved<1,blockSize> select_1_type;
        typedef select_support_interleaved<0,blockSize> select_0_type;
    private:
        size_type m_size;           /* size of the original bit vector */
        size_type m_totalBlocks;    /* total size of m_data in u64s */
        size_type m_superblocks;    /* number of superblocks */
        size_type m_blockShift;
        int_vector<64> m_data;           /* data */
        int_vector<64> m_rank_samples;   /* space for additional rank samples */

        // precondition: m_rank_samples.size() <= m_superblocks
        void init_rank_samples() {
            uint32_t blockSize_U64 = bit_magic::l1BP(blockSize>>6);
            size_type idx = 0;
            std::queue<size_type> lbs, rbs;
            lbs.push(0); rbs.push(m_superblocks);
            while (!lbs.empty()) {
                size_type lb = lbs.front(); lbs.pop();
                size_type rb = rbs.front(); rbs.pop();
                if (/*lb < rb and*/ idx < m_rank_samples.size()) {
                    size_type mid = lb + (rb-lb)/2; // select mid \in [lb..rb)
                    size_type pos = (mid << blockSize_U64) + mid;
                    m_rank_samples[ idx++ ] = m_data[pos];
                    lbs.push(lb); rbs.push(mid);
                    lbs.push(mid+1); rbs.push(rb);
                }
            }
        }

    public:
        bit_vector_interleaved():m_size(0), m_totalBlocks(0), m_superblocks(0), m_blockShift(0) {}

        bit_vector_interleaved(const bit_vector& bv):m_size(0), m_totalBlocks(0), m_superblocks(0), m_blockShift(0) {
            m_size = bv.size();
            /* calculate the number of superblocks */
//          each block of size > 0 gets suberblock in which we store the cumulative sum up to this block
            m_superblocks = (m_size+blockSize) / blockSize;
            m_blockShift = bit_magic::l1BP(blockSize);
            /* allocate new data */
            size_type blocks = (m_size+64)/64;
            size_type mem =  blocks +         m_superblocks + 1;
//                          ^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^   ^
//                          bit vector data | cum. sum data | sum after last block
            util::assign(m_data, int_vector<64>(mem));
            m_totalBlocks = mem;

            /* assign data and calculate super block values */
            const uint64_t* bvp = bv.data();

            size_type j = 0; // 64-bit word counter in the m_data
            size_type cum_sum = 0;
            size_type sample_rate = blockSize/64;
            for (size_type i=0, sample_cnt=sample_rate; i < blocks; ++i, ++sample_cnt) {
                if (sample_cnt == sample_rate) {
                    m_data[j] = cum_sum;
                    sample_cnt = 0;
                    j++;
                }
                m_data[j] = bvp[i];
                cum_sum += bit_magic::b1Cnt(m_data[j]);
                j++;
            }
            m_data[j] = cum_sum; /* last superblock so we can always
                                    get num_ones fast */
            if (m_totalBlocks > 1024*64) {
                // we store at most m_superblocks+1 rank_samples:
                // we do a cache efficient binary search for the select on X=1024
                // or X=the smallest power of two smaller than m_superblock
                m_rank_samples.resize(std::min(1024ULL, 1ULL << bit_magic::l1BP(m_superblocks)));
            }
            init_rank_samples();
        }

        //! Accessing the i-th element of the original bit_vector
        /*! \param i An index i with \f$ 0 \leq i < size()  \f$.
           \return The i-th bit of the original bit_vector
           \par Time complexity
           		\f$ \Order{1} \f$
        */
        value_type operator[](size_type i)const {
            size_type bs = i >> m_blockShift;
            size_type block = bs + (i>>6) + 1;
            return ((m_data[block] >> (i&63)) & 1ULL);
        }

        //! Returns the size of the original bit vector.
        size_type size()const {
            return m_size;
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += util::write_member(m_size, out, child, "size");
            written_bytes += util::write_member(m_totalBlocks, out, child, "totalBlocks");
            written_bytes += util::write_member(m_superblocks, out, child, "superblocks");
            written_bytes += util::write_member(m_blockShift, out, child, "blockShift");
            written_bytes += m_data.serialize(out, child, "data");
            written_bytes += m_rank_samples.serialize(out, child, "rank_samples");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            util::read_member(m_size, in);
            util::read_member(m_totalBlocks, in);
            util::read_member(m_superblocks, in);
            util::read_member(m_blockShift, in);
            m_data.load(in);
            m_rank_samples.load(in);
        }

        void swap(bit_vector_interleaved& bv) {
            if (this != &bv) {
                std::swap(m_size, bv.m_size);
                std::swap(m_totalBlocks, bv.m_totalBlocks);
                std::swap(m_superblocks, bv.m_superblocks);
                std::swap(m_blockShift, bv.m_blockShift);
                m_data.swap(bv.m_data);
                m_rank_samples.swap(bv.m_rank_samples);
            }
        }
};

template<uint8_t b,uint32_t blockSize>
class rank_support_interleaved
{
    public:
        typedef bit_vector::size_type       		size_type;
        typedef bit_vector_interleaved<blockSize> 	bit_vector_type;
    private:
        const bit_vector_type* m_v;
        size_type m_blockShift;
        size_type m_blockMask;
        size_type m_blockSize_U64;  /* blocksize for superblocks in 64 bit words */

        inline size_type rank1(size_type i) const {
            size_type SBlockNum = i >> m_blockShift;
            size_type SBlockPos = (SBlockNum << m_blockSize_U64) + SBlockNum;
            uint64_t resp = m_v->m_data[SBlockPos];
            const uint64_t* B = (m_v->m_data.data() + (SBlockPos+1));
            uint64_t rem = i&63;
            uint64_t bits = (i&m_blockMask) - rem;
            while (bits) {
                resp += bit_magic::b1Cnt(*B++);
                bits -= 64;
            }
            resp += bit_magic::b1Cnt(*B & bit_magic::Li1Mask[rem]);
            return resp;
        }


        inline size_type rank0(size_type i) const {
            size_type SBlockNum = i >> m_blockShift;
            size_type SBlockPos = (SBlockNum << m_blockSize_U64) + SBlockNum;
            uint64_t resp = (SBlockNum << m_blockShift) - m_v->m_data[SBlockPos];
            const uint64_t* B = (m_v->m_data.data() + (SBlockPos+1));
            uint64_t rem = i&63;
            uint64_t bits = (i&m_blockMask) - rem;
            while (bits) {
                resp += bit_magic::b1Cnt(~(*B)); B++;
                bits -= 64;
            }
            resp += bit_magic::b1Cnt((~(*B)) & bit_magic::Li1Mask[rem]);
            return resp;
        }

    public:

        rank_support_interleaved(const bit_vector_type* v=NULL) {
            init(v);
            m_blockShift = bit_magic::l1BP(blockSize);
            m_blockMask = blockSize - 1;
            m_blockSize_U64 = bit_magic::l1BP(blockSize>>6);
        }

        void init(const bit_vector_type* v=NULL) {
            set_vector(v);
        }

        //! Returns the position of the i-th occurrence in the bit vector.
        size_type rank(size_type i) const {
            if (b) return rank1(i);
            return rank0(i);
        }

        const size_type operator()(size_type i)const {
            return rank(i);
        }

        const size_type size()const {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=NULL) {
            m_v = v;
        }

        rank_support_interleaved& operator=(const rank_support_interleaved& rs) {
            if (this != &rs) {
                set_vector(rs.m_v);
            }
            return *this;
        }

        void swap(rank_support_interleaved& rs) { }

        bool operator==(const rank_support_interleaved& ss)const {
            if (this == &ss)
                return true;
            return ss.m_v == m_v;
        }

        bool operator!=(const rank_support_interleaved& rs)const {
            return !(*this == rs);
        }


        void load(std::istream& in, const bit_vector_type* v=NULL) {
            set_vector(v);
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }
};


template<uint8_t b,uint32_t blockSize>
class select_support_interleaved
{
    public:
        typedef bit_vector::size_type size_type;
        typedef bit_vector_interleaved<blockSize> bit_vector_type;
    private:
        const bit_vector_type* m_v;
        size_type m_superblocks;    /* number of superblocks */
        size_type m_blockShift;
        size_type m_blockSize_U64;


        //! Returns the position of the i-th occurrence in the bit vector.
        size_type select1(size_type i) const {
            size_type lb = 0, rb = m_v->m_superblocks; // search interval [lb..rb)
            size_type res = 0;
            size_type idx = 0; // index in m_rank_samples
            /* binary search over super blocks */
            // invariant: lb==0 or m_data[ pos(lb-1) ] < i
            //            m_data[ pos(rb) ] >= i, initial since i < rank(size())
            while (lb < rb) {
                size_type mid = (lb+rb)/2; // select mid \in [lb..rb)
                if (idx < m_v->m_rank_samples.size()) {
                    if (m_v->m_rank_samples[idx] >= i) {
                        idx = (idx<<1) + 1;
                        rb = mid;
                    } else {
                        idx = (idx<<1) + 2;
                        lb = mid + 1;
                    }
                } else {
                    size_type pos = (mid << m_blockSize_U64) + mid;
//                                  ^^^^^^^^^^^^^^^^^^^^^^     ^^^^^^^^^^^^^^^^^^^
//                                    data blocks to jump      superblock position
                    if (m_v->m_data[pos] >= i) {
                        rb = mid;
                    } else {
                        lb = mid + 1;
                    }
                }
            }
            res = (rb-1) << m_blockShift;
            /* iterate in 64 bit steps */
            const uint64_t* w = m_v->m_data.data() + ((rb-1) << m_blockSize_U64) + (rb-1);
            i -= *w;  // subtract the cumulative sum before the superblock
            ++w; /* step into the data */
            size_type ones = bit_magic::b1Cnt(*w);
            while (ones < i) {
                i -= ones; ++w;
                ones = bit_magic::b1Cnt(*w);
                res += 64;
            }
            /* handle last word */
            res += bit_magic::i1BP(*w, i);
            return res;
        }

        //! Returns the position of the i-th occurrence in the bit vector.
        size_type select0(size_type i)const {
            size_type lb = 0, rb = m_v->m_superblocks; // search interval [lb..rb)
            size_type res = 0;
            size_type idx = 0; // index in m_rank_samples
            /* binary search over super blocks */
            // invariant: lb==0 or m_data[ pos(lb-1) ] < i
            //            m_data[ pos(rb) ] >= i, initial since i < rank(size())
            while (lb < rb) {
                size_type mid = (lb+rb)/2; // select mid \in [lb..rb)
                if (idx < m_v->m_rank_samples.size()) {
                    if (((mid << m_blockShift) - m_v->m_rank_samples[idx]) >= i) {
                        idx = (idx<<1) + 1;
                        rb = mid;
                    } else {
                        idx = (idx<<1) + 2;
                        lb = mid + 1;
                    }
                } else {
                    size_type pos = (mid << m_blockSize_U64) + mid;
//                                  ^^^^^^^^^^^^^^^^^^^^^^     ^^^^^^^^^^^^^^^^^^^
//                                    data blocks to jump      superblock position
                    if (((mid << m_blockShift) - m_v->m_data[pos]) >= i) {
                        rb = mid;
                    } else {
                        lb = mid + 1;
                    }
                }
            }
            res = (rb-1) << m_blockShift;

            /* iterate in 64 bit steps */
            const uint64_t* w = m_v->m_data.data() + ((rb-1) << m_blockSize_U64) + (rb-1);
            i = i - (res - *w);  // substract the cumulative sum before the superblock
            ++w; /* step into the data */
            size_type zeros = bit_magic::b1Cnt(~ *w);
            while (zeros < i) {
                i -= zeros; ++w;
                zeros = bit_magic::b1Cnt(~ *w);
                res += 64;
            }
            /* handle last word */
            res += bit_magic::i1BP(~ *w, i);
            return res;
        }

    public:

        select_support_interleaved(const bit_vector_type* v=NULL) {
            init(v);
        }


        void init(const bit_vector_type* v=NULL) {
            set_vector(v);
            m_blockShift = bit_magic::l1BP(blockSize);
            m_blockSize_U64 = bit_magic::l1BP(blockSize>>6);
        }

        //! Returns the position of the i-th occurrence in the bit vector.
        size_type select(size_type i) const {
            if (b) return select1(i);
            return select0(i);
        }

        const size_type operator()(size_type i)const {
            return select(i);
        }

        const size_type size()const {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=NULL) {
            m_v = v;
        }

        select_support_interleaved& operator=(const select_support_interleaved& rs) {
            if (this != &rs) {
                set_vector(rs.m_v);
            }
            return *this;
        }

        void swap(select_support_interleaved& rs) { }

        bool operator==(const select_support_interleaved& ss)const {
            if (this == &ss)
                return true;
            return ss.m_v == m_v;
        }

        bool operator!=(const select_support_interleaved& rs)const {
            return !(*this == rs);
        }


        void load(std::istream& in, const bit_vector_type* v=NULL) {
            set_vector(v);
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }
};

}

#endif
