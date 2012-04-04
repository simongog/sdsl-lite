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

//! Namespace for the succinct data structure library
namespace sdsl
{

template<uint32_t blockSize=512>// forward declaration needed for friend declaration
class rank_support_interleaved;  // in bit_vector_interleaved

template<uint32_t blockSize=512>// forward declaration needed for friend declaration
class select_support_interleaved;  // in bit_vector_interleaved

//! A bit vector which interleaves the original bit_vector with rank information.
/*!
 * This class is a uncompressed bit vector representation. It copies the original
 * bit_vector and interleaves the data every blockSize bits with a cumulative
 * sum of set bits before the current position. Each cumulative sum is stored
 * in a 64 bit word.
 * \pre blockSize has to be a power of 2.
 */
template<uint32_t blockSize=512>
class bit_vector_interleaved
{
    public:
        typedef bit_vector::size_type   size_type;
        typedef size_type               value_type;

        friend class rank_support_interleaved<blockSize>;
        friend class select_support_interleaved<blockSize>;

        typedef rank_support_interleaved<blockSize>   rank_1_type;
        typedef select_support_interleaved<blockSize> select_1_type;
    private:
        size_type m_size;           /* size of the original bitvector */
        size_type m_totalBlocks;    /* total size of m_data in u64s */
        size_type m_blockSize_U64;  /* blocksize for superblocks */
        size_type m_blockMask;      /* blockmask for modulo operation */
        size_type m_superblocks;    /* number of superblocks */
        size_type m_blockShift;
        int_vector<64> m_data;           /* data */

    public:
        bit_vector_interleaved() {}

        bit_vector_interleaved(const bit_vector& bv) {
            m_size = bv.size();
            if (m_size == 0)
                return;

            /* calculate the number of superblocks */
            m_superblocks = m_size / blockSize + 1;

            m_blockShift = bit_magic::l1BP(blockSize);

            /* allocate new data */
            size_type mem = ((m_size+63)>>6)+m_superblocks;
//            m_data = new uint64_t[mem+1];
            util::assign(m_data, int_vector<64>(mem+1));
//			m_data.resize( mem+1 );
            m_totalBlocks = mem+1;

            /* assign data and calculate super block values */
            const uint64_t* bvp = bv.data();

            size_type j = 0;
            size_type popcnt = 0;
            for (size_type i=0; i < (m_size+63)/64; i++) {
                if (i%(blockSize>>6)==0) {
                    m_data[j] = popcnt;
                    j++;
                }
                m_data[j] = bvp[i];
                popcnt += bit_magic::b1Cnt(m_data[j]);
                j++;
            }
            m_data[j] = popcnt; /* last superblock so we can always
                                   get num_ones fast */
        }

        //! Accessing the i-th element of the original bit_vector
        /*! \param i An index i with \f$ 0 \leq i < size()  \f$.
           \return The i-th bit of the original bit_vector
           \par Time complexity
           		\f$ \Order{\log m} \f$, where m equals the number of zeros
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
        size_type serialize(std::ostream& out)const {
            size_type written_bytes = 0;
            written_bytes += util::write_member(m_size, out);
            written_bytes += util::write_member(m_totalBlocks, out);
            written_bytes += util::write_member(m_blockSize_U64, out);
            written_bytes += util::write_member(m_blockMask, out);
            written_bytes += util::write_member(m_superblocks, out);
            written_bytes += util::write_member(m_blockShift, out);
            written_bytes += m_data.serialize(out);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            util::read_member(m_size, in);
            util::read_member(m_totalBlocks, in);
            util::read_member(m_blockSize_U64, in);
            util::read_member(m_blockMask, in);
            util::read_member(m_superblocks, in);
            util::read_member(m_blockShift, in);
            m_data.load(in);
        }

#ifdef MEM_INFO
        void mem_info(std::string label="")const {
            if (label=="")
                label = "bit_vector_interleaved";
            size_type bytes = util::get_size_in_bytes(*this) + m_totalBlocks*sizeof(uint64_t);
            std::cout << "list(label=\""<<label<<"\", size = "<< bytes/(1024.0*1024.0) << "\n,";
        }
#endif

};

template<uint32_t blockSize>
class rank_support_interleaved
{
    public:
        typedef bit_vector::size_type       		size_type;
        typedef bit_vector_interleaved<blockSize> 	bit_vector_type;
    private:
        const bit_vector_type* m_v;
        size_type m_blockShift;
        size_type m_blockMask;
        size_type m_blockSize_U64;  /* blocksize for superblocks */

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

        inline size_type rank(size_type i) const {
            size_type SBlockNum = i >> m_blockShift;
            size_type SBlockPos = (SBlockNum << m_blockSize_U64) + SBlockNum;
            uint64_t resp = m_v->m_data[SBlockPos];
            const uint64_t* B = (m_v->m_data.data() + (SBlockPos+1));
            uint64_t rem = i&63;
            uint64_t bits = (i&m_blockMask) - rem;
// TODO: remove l, and access B by *B and increase B afterwards
            while (bits) {
                resp += bit_magic::b1Cnt( *B++ );
                bits -= 64;
            }
            resp += bit_magic::b1Cnt(*B & bit_magic::Li1Mask[rem]);
            return resp;
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

        size_type serialize(std::ostream& out)const {
            size_type written_bytes = 0;
            return written_bytes;
        }

#ifdef MEM_INFO
        void mem_info(std::string label="")const {
            if (label=="")
                label = "rank_support_interleaved";
            size_type bytes = util::get_size_in_bytes(*this);
            std::cout << "list(label=\""<<label<<"\", size = "<< bytes/(1024.0*1024.0) << ")\n";
        }
#endif

};


template<uint32_t blockSize>
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
        size_type select(size_type i)const {
            size_type l=0, h=m_v->m_superblocks-1;
            size_type mid = l + ((h - l) >> 1);
            size_type pos = (mid << m_blockSize_U64) + mid;
            size_type res = 0;

            /* binary search over super blocks */
            while (l<=h) {
                if (m_v->m_data[pos]<i)
                    l = mid+1;
                else
                    h = mid-1;

                mid = l + ((h - l) >> 1);
                pos = (mid << m_blockSize_U64) + mid;
            }

            res = mid << m_blockShift;

            /* iterate in 64 bit steps */
			const uint64_t *w = m_v->m_data.data() + pos;  // m_v->m_data[pos
			i -= *w; 

            ++w; /* step into the data */
            size_type ones = bit_magic::b1Cnt( *w );
            while (ones < i) {
                i -= ones; ++w;
                ones = bit_magic::b1Cnt( *w );
                res += 64;
            }

            /* last integer */
            res += bit_magic::i1BP( *w, i );
            return res;
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

        size_type serialize(std::ostream& out)const {
            size_type written_bytes = 0;
            return written_bytes;
        }

#ifdef MEM_INFO
        void mem_info(std::string label="")const {
            if (label=="")
                label = "select_support_interleaved";
            size_type bytes = util::get_size_in_bytes(*this);
            std::cout << "list(label=\""<<label<<"\", size = "<< bytes/(1024.0*1024.0) << ")\n";
        }
#endif

};

}

#endif
