/* sdsl - succinct data structures library
    Copyright (C) 2012 Simon Gog 

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
/*! \file gap_vector.hpp
   \brief gap_vector.hpp contains the sdsl::gap_vector class, and 
          classes which support rank and select for gap_vector.
   \author Simon Gog
*/ 
#ifndef SDSL_BIT1L_VECTOR
#define SDSL_BIT1L_VECTOR

#include "int_vector.hpp"
#include "util.hpp"

//! Namespace for the succinct data structure library
namespace sdsl{

template<uint32_t blockSize>// forward declaration needed for friend declaration
class bit1l_rank_support;  // in gap_vector

template<uint32_t blockSize>// forward declaration needed for friend declaration
class bit1l_select_support;  // in gap_vector

//! A bit vector which compresses very sparse populated bit vectors by representing the 1 or 0 by gap encoding
template<uint32_t blockSize=512>	
class bit1l_vector{
	public:
		typedef bit_vector::size_type   size_type;
		typedef size_type               value_type;

		friend class bit1l_rank_support<blockSize>;
		friend class bit1l_select_support<blockSize>;

		typedef bit1l_rank_support<blockSize>   rank_1_type;
		typedef bit1l_select_support<blockSize> select_1_type;
	private:
		size_type m_size;           /* size of the original bitvector */
        size_type m_totalBlocks;    /* total size of m_data in u64s */
        size_type m_blockSize_U64;  /* blocksize for superblocks */
        size_type m_blockMask;      /* blockmask for modulo operation */
        size_type m_superblocks;    /* number of superblocks */
        size_type m_blockShift;
        uint64_t* m_data;           /* data */

	public:
		bit1l_vector(){}

		bit1l_vector(const bit_vector &bv) {
			m_size = bv.size();
			if(m_size == 0)
				return;

            /* calculate the number of superblocks */
            m_superblocks = m_size / blockSize + 1;

            m_blockShift = bit_magic::l1BP(blockSize);

            /* allocate new data */
            size_type mem = ((m_size+63)>>6)+m_superblocks;
            m_data = new uint64_t[mem+1];
            m_totalBlocks = mem+1;

            /* assign data and calculate super block values */
            const uint64_t* bvp = bv.data();

            size_type j = 0;
            size_type popcnt = 0;
            for(size_type i=0;i < (m_size+63)/64;i++) {
                if(i%(blockSize>>6)==0) {
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
		value_type operator[](size_type i)const{
            size_type bs = i >> m_blockShift;
            size_type block = bs + (i>>6) + 1;
            return ((m_data[block] >> (i&63)) & 1ULL);
		}

	   //! Returns the size of the original bit vector.
	   size_type size()const{
		    return m_size;
	    }

		//! Serializes the data structure into the given ostream
	   size_type serialize(std::ostream &out)const{
	       size_type written_bytes = 0;
		   written_bytes += util::write_member(m_size, out);
		   written_bytes += util::write_member(m_totalBlocks, out);

           /* write in blocks */
           uint64_t* p = m_data;
           const static size_type SDSL_BLOCK_SIZE = (1<<20);
           size_type idx = 0;
           while( idx+SDSL_BLOCK_SIZE < m_totalBlocks ) {
               out.write((char*) p, SDSL_BLOCK_SIZE*sizeof(uint64_t) );
               written_bytes += SDSL_BLOCK_SIZE*sizeof(uint64_t);
               p   += SDSL_BLOCK_SIZE;
               idx += SDSL_BLOCK_SIZE;
           }
           out.write((char*) p, (m_totalBlocks-idx)*sizeof(uint64_t));
           written_bytes += (m_totalBlocks-idx)*sizeof(uint64_t);
		   return written_bytes;
	    }

		//! Loads the data structure from the given istream.
		void load(std::istream &in){
			util::read_member(m_size, in);
			util::read_member(m_totalBlocks, in);
            /* calc the rest */
            m_superblocks = m_size / blockSize + 1;
            /* blocksize in U64s */
            m_blockShift = bit_magic::l1BP(blockSize);

            m_data = new uint64_t[m_totalBlocks];

            uint64_t *p = m_data;
            const static size_type SDSL_BLOCK_SIZE = (1<<20);
            size_type idx = 0;

            while( idx+SDSL_BLOCK_SIZE < m_totalBlocks ){
                in.read((char*) p, SDSL_BLOCK_SIZE*sizeof(uint64_t) );
                p   += SDSL_BLOCK_SIZE;
                idx += SDSL_BLOCK_SIZE;
            }
            in.read((char*) p, (m_totalBlocks-idx)*sizeof(uint64_t));
		}	

#ifdef MEM_INFO
		void mem_info(std::string label="")const{
			if(label=="")
				label = "bit1l_vector";
			size_type bytes = util::get_size_in_bytes(*this) + m_totalBlocks*sizeof(uint64_t);
            std::cout << "list(label=\""<<label<<"\", size = "<< bytes/(1024.0*1024.0) << "\n,";
		}
#endif	

};

template<uint32_t blockSize=512>	
class bit1l_rank_support{
	public:
		typedef bit_vector::size_type       size_type;
		typedef bit1l_vector<blockSize>               bit_vector_type;
	private:
		const bit_vector_type* m_v;
        size_type m_blockShift;
        size_type m_blockMask;
        size_type m_blockSize_U64;  /* blocksize for superblocks */

	public:

		bit1l_rank_support(const bit_vector_type *v=NULL){
			init(v);
            m_blockShift = bit_magic::l1BP(blockSize);
            m_blockMask = blockSize - 1;
            m_blockSize_U64 = bit_magic::l1BP(blockSize>>6);
		}

		void init(const bit_vector_type *v=NULL){
			set_vector(v);
		}

		inline size_type rank(size_type i) const{

            /*i++;*/
            size_type SBlockNum = i >> m_blockShift;
            size_type SBlockPos = (SBlockNum << m_blockSize_U64) + SBlockNum;
            uint64_t resp = m_v->m_data[SBlockPos];
            uint64_t* B = (uint64_t*) &(m_v->m_data[SBlockPos+1]);
            uint64_t rem = i&63;
            uint64_t bits = (i&m_blockMask) - rem;
            uint64_t l=0;

            while(bits) {
                resp+=__builtin_popcountll(B[l]);
                bits -= 64;
                l++;
            }
            resp += __builtin_popcountll(B[l]&bit_magic::Li1Mask[rem]);
            return resp;
		} 

		const size_type operator()(size_type i)const{
			return rank(i);
		}

		const size_type size()const{
			return m_v->size();
		}

		void set_vector(const bit_vector_type *v=NULL){
			m_v = v;
		}

		bit1l_rank_support& operator=(const bit1l_rank_support &rs){
			if(this != &rs){
				set_vector(rs.m_v);
			}
			return *this;
		}

		void swap(bit1l_rank_support &rs){ }

		bool operator==(const bit1l_rank_support &ss)const{
			if(this == &ss)
				return true;
			return ss.m_v == m_v;
		}

		bool operator!=(const bit1l_rank_support &rs)const{
			return !(*this == rs);	
		}


		void load(std::istream &in, const bit_vector_type *v=NULL){
			set_vector(v);
		}

		size_type serialize(std::ostream &out)const{
			size_type written_bytes = 0;
			return written_bytes;
		}

#ifdef MEM_INFO
		void mem_info(std::string label="")const{
			if(label=="")
				label = "bit1l_rank_support";
			size_type bytes = util::get_size_in_bytes(*this);
			std::cout << "list(label=\""<<label<<"\", size = "<< bytes/(1024.0*1024.0) << ")\n";
		}
#endif	

};


template<uint32_t blockSize=512>	
class bit1l_select_support{
	public:
		typedef bit_vector::size_type size_type;
		typedef bit1l_vector<blockSize> bit_vector_type;
	private:
		const bit_vector_type *m_v;

	public:

		bit1l_select_support(const bit_vector_type *v=NULL){
			init(v);
		}

		void init(const bit_vector_type *v=NULL){
			set_vector(v);
		}

		//! Returns the position of the i-th occurrence in the bit vector. 
		size_type select(size_type i)const{
            return 2;
		} 

		const size_type operator()(size_type i)const{
			return select(i);
		}

		const size_type size()const{
			return m_v->size();
		}

		void set_vector(const bit_vector_type *v=NULL){
			m_v = v;
		}

		bit1l_select_support& operator=(const bit1l_select_support &rs){
			if(this != &rs){
				set_vector(rs.m_v);
			}
			return *this;
		}

		void swap(bit1l_select_support &rs){ }

		bool operator==(const bit1l_select_support &ss)const{
			if(this == &ss)
				return true;
			return ss.m_v == m_v;
		}

		bool operator!=(const bit1l_select_support &rs)const{
			return !(*this == rs);	
		}


		void load(std::istream &in, const bit_vector_type *v=NULL){
			set_vector(v);
		}

		size_type serialize(std::ostream &out)const{
			size_type written_bytes = 0;
			return written_bytes;
		}

#ifdef MEM_INFO
		void mem_info(std::string label="")const{
			if(label=="")
				label = "bit1l_select_support";
			size_type bytes = util::get_size_in_bytes(*this);
			std::cout << "list(label=\""<<label<<"\", size = "<< bytes/(1024.0*1024.0) << ")\n";
		}
#endif	

};

}

#endif
