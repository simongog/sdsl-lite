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
/*! \file lcp_dac.hpp
    \brief lcp_dac.hpp contains an implementation of a (compressed) lcp array proposed by Brisaboa, Ladra and Navarro
           in the paper "Direct addressable variable-length codes" (SPIRE 2009).
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_LCP_DAC
#define INCLUDED_SDSL_LCP_DAC

#include "lcp.hpp"
#include "int_vector.hpp"
#include "algorithms.hpp"
#include "iterators.hpp"
#include "bitmagic.hpp"
#include "util.hpp"
#include "testutils.hpp"
#include "rank_support_v.hpp"
#include "rank_support_v5.hpp"
#include <iostream>
#include <cassert>
#include <cstring> // for strlen
#include <iomanip>
#include <iterator>
#include <vector>
#include <algorithm> // for max
#include <stdexcept>

namespace sdsl
{

//! A class for the compressed version of lcp information of an suffix array
/*! We use a technique called ,,escaping'' to encode the values.
 *  This is defined as follows (see [1]):
 *  A k-bit integer is split into \f$K=\lceilk/(b-1)\rceil\f$ bits each and
 *  encoded into \f$K\f$ blocks of \f$ b \f$ bits each. All but the last block
 *  are marked with by a 1 in the most significant bit. Escaping with b=8 is
 *  also known as vbyte-coding (see [2]). A experimental study of using escaping
 *  for the LCP array is given in [3].
 *  \par Time complexity
 *		- \f$\Order{\log n/b}\f$ worst case, where b is the number of bits
          in a block
 *  \par References
 *       [1] F. Transier and P. Sanders: ,,Engineering Basic Search Algorithms
 *           of an In-Memory Text Search Engine'', ACM Transactions on
 *           Information Systems, Vol. 29, No.1, Article 2, 2010
 *       [2] H.E. Williams and J. Zobel: ,,Compressing integers for fast file
 *           access'', Computing Journal Vol 43, No.3, 1999
 *       [3] N. Brisboa, S. Ladra, G. Navarro: ,,Directly addressable variable-
 *           length codes'', Proceedings of SPIRE 2009.
 */
template<uint8_t b=4, class rank_support_type=rank_support_v5<> >
class lcp_dac
{
    public:
        typedef typename int_vector<>::value_type		 value_type;	// STL Container requirement
        typedef random_access_const_iterator<lcp_dac>		 const_iterator;// STL Container requirement
        typedef const_iterator 								 iterator;		// STL Container requirement
        typedef const value_type							 const_reference;
        typedef const_reference								 reference;
        typedef const_reference*							 pointer;
        typedef const pointer								 const_pointer;
        typedef int_vector<>::size_type						 size_type;		// STL Container requirement
        typedef ptrdiff_t  									 difference_type; // STL Container requirement

        typedef lcp_plain_tag								 lcp_category; // TODO?

        enum { fast_access = 0,
               text_order  = 0,
               sa_order	  = 1
             }; // as the lcp_dac is not fast for texts with long repetition

        template<class Cst>  // template inner class which is used in CSTs to parametrize lcp classes
        class type           // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef lcp_dac lcp_type;
        };

    private:

        int_vector<b> 		m_data;			  	// vector which holds the block data for every level
        bit_vector	 		m_overflow;			// indicates, if there exists another block for the current number
        rank_support_type 	m_overflow_rank;	// rank data structure for m_overflow
        int_vector<64> 		m_level_pointer_and_rank;
        uint8_t				m_max_level;   // maximal number of levels, at most (log n)/b+1
#ifdef LCP_DAC_CACHING
        int_vector<64>		m_rank_cache;
#endif


        void copy(const lcp_dac& lcp_c) {
            m_data						= lcp_c.m_data;
            m_overflow					= lcp_c.m_overflow;
            m_overflow_rank				= lcp_c.m_overflow_rank;
            m_overflow_rank.set_vector(&m_overflow);
            m_level_pointer_and_rank 	= lcp_c.m_level_pointer_and_rank;
            m_max_level					= lcp_c.m_max_level;
        }

    public:
        //! Default Constructor
        lcp_dac() {
            // has to be initialized for size() methode
            // m_level_pointer_and_rank[2] contains the length of the LCP array
            m_level_pointer_and_rank = int_vector<64>(4,0);
        }

        //! Copy constructor
        lcp_dac(const lcp_dac& lcp_c) {
            copy(lcp_c);
        }

        //! Construct the lcp array from an int_vector_file_buffer
        template<uint8_t int_width, class size_type_class>
        lcp_dac(int_vector_file_buffer<int_width, size_type_class>& lcp_buf);


        //! Number of elements in the instance.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         */
        size_type size()const {
            return m_level_pointer_and_rank[2];
        }

        //! Returns the largest size that lcp_dac can ever have.
        /*! Required for the Container Concept of the STL.
         *  \sa size
         */
        static size_type max_size() {
            return int_vector<>::max_size();
        }

        //! Returns if the data strucutre is empty.
        /*! Required for the Container Concept of the STL.A
         * \sa size
         */
        bool empty()const {
            return m_level_pointer_and_rank[2]==0;
        }

        //! Swap method for lcp_dac
        /*! The swap method can be defined in terms of assignment.
        	This requires three assignments, each of which, for a container type, is linear
        	in the container's size. In a sense, then, a.swap(b) is redundant.
        	This implementation guaranties a run-time complexity that is constant rather than linear.
        	\param lcp_c lcp_dac to swap.

        	Required for the Assignable Conecpt of the STL.
          */
        void swap(lcp_dac& lcp_c);

        //! Returns a const_iterator to the first element.
        /*! Required for the STL Container Concept.
         *  \sa end
         */
        const_iterator begin()const;

        //! Returns a const_iterator to the element after the last element.
        /*! Required for the STL Container Concept.
         *  \sa begin.
         */
        const_iterator end()const;

        //! []-operator
        /*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
         * Time complexity: O(suffix array access)
         * Required for the STL Random Access Container Concept.
         */
        inline value_type operator[](size_type i)const;

        //! Assignment Operator.
        /*!
         *	Required for the Assignable Concept of the STL.
         */
        lcp_dac& operator=(const lcp_dac& lcp_c);

        //! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="") const;

        //! Load from a stream.
        /*! \param in Inputstream to load the data structure from.
         */
        void load(std::istream& in);
};

// == template functions ==


template<uint8_t b, class rank_support_type>
template<uint8_t int_width, class size_type_class>
lcp_dac<b, rank_support_type>::lcp_dac(int_vector_file_buffer<int_width, size_type_class>& lcp_buf)
{
//  (1) Count for each level, how many blocks are needed for the representation
//      Running time: \f$ O(n \times \frac{\log n}{b}  \f$
//      Result is sorted in m_level_pointer_and_rank
    lcp_buf.reset();
    size_type n = lcp_buf.int_vector_size, val=0;
    if (n == 0)
        return;
// 		initialize counter
    m_level_pointer_and_rank.resize(std::max(4*bit_magic::l1BP(2), 2*(((bit_magic::l1BP(n)+1)+b-1) / b)));
    for (size_type i=0; i < m_level_pointer_and_rank.size(); ++i)
        m_level_pointer_and_rank[i] = 0;
    m_level_pointer_and_rank[0] = n; // level 0 has n entries

    uint8_t level_x_2 = 0;
    for (size_type i=0, r_sum=0, r = lcp_buf.load_next_block(); r_sum < n;) {
        for (; i < r_sum+r; ++i) {
            val=lcp_buf[i-r_sum];
            val >>= b; // shift value b bits to the right
            level_x_2 = 2;
            while (val) {
                ++m_level_pointer_and_rank[level_x_2]; // increase counter for current level by 1
                val >>= b; // shift value b bits to the right
                level_x_2 += 2; // increase level by 1
            }
        }
        r_sum += r; r = lcp_buf.load_next_block();
    }

//  (2)	Determine maximum level and prefix sums of level counters
    m_max_level = 0;
    size_type sum_blocks = 0, last_block_size=0;
    for (size_type i=0, t=0; i < m_level_pointer_and_rank.size(); i+=2) {
        t = sum_blocks;
        sum_blocks += m_level_pointer_and_rank[i];
        m_level_pointer_and_rank[i] = t;
//		std::cout<<"i="<<i<<" cnt="<<t<<std::endl;
        if (sum_blocks > t) {
            ++m_max_level;
            last_block_size = sum_blocks - t;
        }
    }
    m_overflow = bit_vector(sum_blocks - last_block_size, 0);
    m_data.resize(sum_blocks);

//	std::cerr<<"m_overflow.size()="<<m_overflow.size()<<" "<<"m_data.size()="<<m_data.size()<<std::endl;
    if (last_block_size == 0) {
        std::cerr<<"ERROR: lcp_dac constructor: last_block_size=0"<<std::endl;
    }

//  (3)	Enter block and overflow data
    int_vector<64> cnt = m_level_pointer_and_rank;
    const uint64_t mask = bit_magic::Li1Mask[b];

    lcp_buf.reset();
    for (size_type i=0,j=0, r_sum=0, r = lcp_buf.load_next_block(); r_sum < n;) {
        for (; i < r_sum+r; ++i) {
            val=lcp_buf[i-r_sum];
            j = cnt[0]++;
            m_data[ j ] =  val & mask;
            val >>= b; // shift value b bits to the right
            level_x_2 = 2;
            while (val) {
//					if(j >= m_overflow.size()){
//						std::cerr<<"i="<<i<<" j="<<j<<" n="<<n<<" m_overflow.size()="<<m_overflow.size()<<" val="<<val<<" level="<<level_x_2/2<<std::endl;
//						return;
//					}
                m_overflow[j] = 1;
                j = cnt[level_x_2]++; // increase counter for current level by 1
//					if(j >= m_data.size()){
//						std::cerr<<"i="<<i<<" j="<<j<<" n="<<n<<" m_data.size()="<<m_data.size()<<" val="<<val<<" level="<<level_x_2/2<<std::endl;
//					}
                m_data[ j ] = val & mask;
                val >>= b; // shift value b bits to the right
                level_x_2 += 2; // increase level by 1
            }
        }
        r_sum += r; r = lcp_buf.load_next_block();
    }

//	std::cerr<<"m_overflow.size()="<<m_overflow.size()<<" "<<"m_data.size()="<<m_data.size()<<std::endl;
//  (4) Initialize rank data structure for m_overflow and precalc rank for
//      pointers
    util::init_support(m_overflow_rank, &m_overflow);
    for (size_type i=0; 2*i < m_level_pointer_and_rank.size() and
         m_level_pointer_and_rank[2*i] < m_overflow.size(); ++i) {
        m_level_pointer_and_rank[2*i+1] = m_overflow_rank(m_level_pointer_and_rank[2*i]);
    }
}

template<uint8_t b, class rank_support_type>
void lcp_dac<b, rank_support_type>::swap(lcp_dac& lcp_c) {
    m_data.swap(lcp_c.m_data);
    m_overflow.swap(lcp_c.m_overflow);
    util::swap_support(m_overflow_rank, lcp_c.m_overflow_rank, &m_overflow, &(lcp_c.m_overflow));

    m_level_pointer_and_rank.swap(lcp_c.m_level_pointer_and_rank);
    std::swap(m_max_level, lcp_c.m_max_level);
}

template<uint8_t b, class rank_support_type>
inline typename lcp_dac<b, rank_support_type>::value_type lcp_dac<b, rank_support_type>::operator[](size_type i)const {
    uint8_t level = 1;
    uint8_t offset = b;
    size_type result = m_data[i];
//	return result;
    const uint64_t* p = m_level_pointer_and_rank.data();
    uint64_t ppi = (*p)+i; // *p plus i
    while (level < m_max_level and m_overflow[ppi]) {
        p += 2;
        ppi = *p + (m_overflow_rank(ppi) - *(p-1));
        result |= (m_data[ppi] << (offset));
        ++level;
        offset += b;
    }
    return result;
}


template<uint8_t b, class rank_support_type>
typename lcp_dac<b, rank_support_type>::size_type lcp_dac<b, rank_support_type>::serialize(std::ostream& out, structure_tree_node* v, std::string name) const {
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_data.serialize(out, child, "data");
    written_bytes += m_overflow.serialize(out, child, "overflow");
    written_bytes += m_overflow_rank.serialize(out, child, "overflow_rank");
    written_bytes += m_level_pointer_and_rank.serialize(out, child, "level_pointer_and_rank");
    written_bytes += util::write_member(m_max_level, out, child, "max_level");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<uint8_t b, class rank_support_type>
void lcp_dac<b, rank_support_type>::load(std::istream& in) {
    m_data.load(in);
    m_overflow.load(in);
    m_overflow_rank.load(in, &m_overflow); // FIXED that at 2012-02-01 15:34
    m_level_pointer_and_rank.load(in);
    util::read_member(m_max_level, in);
}

template<uint8_t b, class rank_support_type>
lcp_dac<b, rank_support_type>& lcp_dac<b, rank_support_type>::operator=(const lcp_dac& lcp_c) {
    if (this != &lcp_c) {
        copy(lcp_c);
    }
    return *this;
}

template<uint8_t b, class rank_support_type>
typename lcp_dac<b, rank_support_type>::const_iterator lcp_dac<b, rank_support_type>::begin()const {
    return const_iterator(this, 0);
}

template<uint8_t b, class rank_support_type>
typename lcp_dac<b, rank_support_type>::const_iterator lcp_dac<b, rank_support_type>::end()const {
    return const_iterator(this, size());
}



} // end namespace sdsl

#endif
