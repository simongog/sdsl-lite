/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

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
/*! \file csa_sada_theo.hpp
    \brief csa_sada_theo.hpp contains an implemenation of the compressed suffix array.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_CSA_SADA_THEO
#define INCLUDED_SDSL_CSA_SADA_THEO

#include "int_vector.hpp"
#include "enc_vector.hpp"
#include "algorithms.hpp"
#include "iterators.hpp"
#include "csa_bitcompressed.hpp"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring> // for strlen
#include <iomanip>
#include <iterator>

namespace sdsl
{

template<class EncVector, class RankSupport>
class csa_sada_theo;

//! A class for the Compressed Suffix Array (CSA) proposed by Sadakane. csa_sada_theo is the one to one implementation of Sadakane's theoretical description with our data structures.
/*! The CSA is parameterized with an enc_vector and a rank_support.
  * The enc_vector holds the \f$\Psi\f$-functions and rank_support helps navigating at the different levels of the data structure.
  * Sadakane proposed to use \f$O(\log\log n)\f$ levels to achive the desired space complexity of the data structure.
  * \sa csa_sada
  * @ingroup csa
 */
template<class EncVector = enc_vector<>, class RankSupport = rank_support_v<> >
class csa_sada_theo
{
    public:
        typedef uint64_t											 value_type;	// STL Container requirement
        typedef random_access_const_iterator<csa_sada_theo> const_iterator;// STL Container requirement
        typedef const_iterator 										 iterator;		// STL Container requirement
        typedef const value_type									 const_reference;
        typedef const_reference										 reference;
        typedef const_reference*									 pointer;
        typedef const pointer										 const_pointer;
        typedef int_vector<>::size_type								 size_type;		// STL Container requirement
        typedef size_type 											 csa_size_type;
        typedef ptrdiff_t  											 difference_type; // STL Container requirement
        typedef const unsigned char*								 pattern_type;

    private:
        uint8_t m_level;   // number of levels of the data structure
//		size_type m_size; // number of eleements of the suffix array --> implicit stored in m_b[0].size()
        EncVector* m_psi;  // psi function at different levels
        int_vector<1>* m_b; // bit vectors supporting navigation at different levels
        RankSupport* m_b_rank; // rank support for m_b
        int_vector<> m_sa; // suffix array at last level (m_level)
        int_vector<> m_last_pos; // isa[m_b[l].size()-1]

        uint16_t 		m_char2comp[256];
        unsigned char 	m_comp2char[256];
        size_type       m_C[257];

        uint8_t m_sigma;

        void copy(const csa_sada_theo<EncVector, RankSupport>& csa);
        // lookup SA at position i in level l
        value_type lookup(size_type i, size_type l=0)const {
//std::cerr<<"i="<<i<<" l="<<l<<std::endl;
            if (l < m_level) {
                if (m_b[l][i]) {
                    return lookup(m_b_rank[l].rank(i), l+1)<<1;
                } else {
                    if (i==m_last_pos[l]) // if i is the index of the greatest entry in SA_l
                        return m_b[l].size()-1;
//						if( i - m_b_rank[l].rank(i) >= m_psi[l].size() ){
//							std::cerr<<"i="<<i<<" l="<<l<<" m_psi.size()="<<m_psi[l].size()<<std::endl;
//						}
                    return lookup(m_psi[l][i - m_b_rank[l].rank(i)], l)-1;
                }
            } else { // at the last level return the stored values
                return m_sa[i];
            }
        }

//		value_type lookup_psi(size_type i, size_type l=0)const{
//
//		}

        // construct the compressed suffix array from a suffix array
        template<typename RandomAccessContainer>
        void construct(RandomAccessContainer& sa, typename RandomAccessContainer::value_type isa_0);


        void setText(const unsigned char* str, size_type len);
    public:
        const uint16_t* char2comp;
        const unsigned char* comp2char;
        const size_type* C;

        //! Default Constructor
        csa_sada_theo():m_level(0),m_psi(NULL),m_b(NULL),m_b_rank(NULL) { }

        //! Copy constructor
        csa_sada_theo(const csa_sada_theo<EncVector, RankSupport>& csa);

        //! Constructor for the CSA taking an already calculated suffix array sa
        /*! \param sa A Container containing an already calculated suffix array.
         *  \param str String for which sa is a suffix array.
         */
        template<typename RandomAccessContainer>
        csa_sada_theo(const RandomAccessContainer& sa, const unsigned char* str);

        //! Constructor for the CSA taking a string for that the CSA should be calculated
        csa_sada_theo(const unsigned char* str);

        //! Default destructor.
        ~csa_sada_theo();

        //! Number of elements in the instance.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size
         */
        size_type size()const {
            return m_level > 0 ? m_b[0].size() : m_sa.size();
        }

        //! Returns the largest size that csa_sada_theo can ever have.
        /*! Required for the Container Concept of the STL.
         *  \sa size
         */
        static size_type max_size() {
            return EncVector::max_size();
        }

        //! Returns if the CSA is empty.
        /*! Required for the Container Concept of the STL.
         */
        bool empty()const {
            return size()==0;
        }

        //! Swap method for csa_sada_theo
        /*! The swap method can be defined in terms of assignment.
        	This requires three assignments, each of which, for a container type, is linear
        	in the container's size. In a sense, then, a.swap(b) is redundant.
        	This implementation guaranties a run-time complexity that is constant rather than linear.
        	\param csa csa_sada_theo to swap.

        	Required for the Assignable Conecpt of the STL.
          */
        void swap(csa_sada_theo<EncVector, RankSupport>& csa);

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
         *
         * Required for the STL Random Access Container Concept.
         */
        inline value_type operator[](size_type i)const;

        //! Assignment Operator
        /*!
         *	Required for the Assignable Concept of the STL.
         */
        csa_sada_theo& operator=(const csa_sada_theo& csa);

        //! Serialize the SDSCompressedSuffixArray to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out) const;

        //! Load SDSCompressedSuffixArray from a stream.
        void load(std::istream& in);
};

// == template functions ==

template<class EncVector, class RankSupport>
void csa_sada_theo<EncVector,RankSupport>::swap(csa_sada_theo<EncVector, RankSupport>& csa)
{
    if (this!=&csa) {
        std::swap(m_level, csa.m_level);
        m_sa.swap(csa.m_sa);
        m_last_pos.swap(csa.m_last_pos);
        std::swap(m_psi 	, csa.m_psi);     // only pointers => swap it with std::swap
        std::swap(m_b	 	, csa.m_b);       // only pointers => swap it with std::swap
        std::swap(m_b_rank	, csa.m_b_rank);  // only pointers => swap it with std::swap
        for (uint16_t i=0; i<256; ++i) {
            std::swap(m_char2comp[i], csa.m_char2comp[i]);
            std::swap(m_comp2char[i], csa.m_comp2char[i]);
            std::swap(m_C[i], csa.m_C[i]);
        }
        std::swap(m_C[256], csa.m_C[256]);
    }
}


template<class EncVector, class RankSupport>
csa_sada_theo<EncVector,RankSupport>::csa_sada_theo(const csa_sada_theo<EncVector,RankSupport>& csa):m_psi(NULL),m_b(NULL),m_b_rank(NULL),char2comp(m_char2comp),comp2char(m_comp2char),C(m_C)
{
    copy(csa);
}

template<class EncVector, class RankSupport>
csa_sada_theo<EncVector,RankSupport>& csa_sada_theo<EncVector,RankSupport>::operator=(const csa_sada_theo<EncVector,RankSupport>& csa)
{
    if (this != &csa) {
        copy(csa);
    }
    return *this;
}


template<class EncVector, class RankSupport>
void csa_sada_theo<EncVector, RankSupport>::setText(const unsigned char* str, size_type len)
{
    for (uint16_t i=0; i<256; ++i)
        m_char2comp[i] = m_comp2char[i] = m_C[i] = 0;
    m_C[256] = 0;
    if (str == NULL or len ==0)
        return;
    const unsigned char* p = str;
    for (size_type i=0; i<len-1; ++i) {
        m_char2comp[*p++] = 1;
    }
    size_type value = 0;
    for (uint16_t i=0; i<256; ++i)
        if (m_char2comp[i])
            m_char2comp[i] = value++;
    m_sigma = value;
    for (uint16_t i=0; i<256; ++i)
        if (m_char2comp[i])
            m_comp2char[m_char2comp[i]] = i;
    for (uint16_t i=0; i<257; ++i) m_C[i]=0;
    for (size_type i=0; i<len; ++i) ++m_C[m_char2comp[str[i]]];
    for (uint16_t i=256; i>0; --i) m_C[i] = m_C[i-1];
    m_C[0]=0;
    for (uint16_t i=1; i<257; ++i) m_C[i] += m_C[i-1];
    assert(m_C[256]==len);
}

template<class EncVector, class RankSupport>
void csa_sada_theo<EncVector,RankSupport>::copy(const csa_sada_theo<EncVector,RankSupport>& csa)
{
    // first free the data structure if it is already in use
    if (m_b_rank != NULL) {
        delete [] m_b_rank;
        m_b_rank = NULL;
    }
    if (m_b		!= NULL) {
        delete [] m_b;
        m_b = NULL;
    }
    if (m_psi	!= NULL) {
        delete [] m_psi;
        m_psi = NULL;
    }
    // second copy the data from csa to this
    m_level 	= csa.m_level;
    m_sa		= csa.m_sa;
    m_last_pos  = csa.m_last_pos;

    if (m_level > 0) {
        // copy arrays
        m_psi	= new EncVector[m_level];
        for (size_type i=0; i < m_level; ++i) {
            m_psi[i]	= csa.m_psi[i];
        }

        m_b		= new int_vector<1>[m_level];
        for (size_type i=0; i < m_level; ++i) {
            m_b[i]	= csa.m_b[i];
        }

        m_b_rank = new RankSupport[m_level];
        for (size_type i=0; i < m_level; ++i) {
            m_b_rank[i] = csa.m_b_rank[i];
            m_b_rank[i].set_vector(&m_b[i]);
        }
    }
    for (int i=0; i<256; ++i) {
        m_char2comp[i] 	= csa.m_char2comp[i];
        m_comp2char[i]	= csa.m_comp2char[i];
        m_C[i]			= csa.m_C[i];
    }
    m_C[256] = csa.m_C[256];
}

template<class EncVector, class RankSupport>
typename csa_sada_theo<EncVector,RankSupport>::const_iterator csa_sada_theo<EncVector,RankSupport>::begin()const
{
    return const_iterator(this,0);
}

template<class EncVector, class RankSupport>
typename csa_sada_theo<EncVector,RankSupport>::const_iterator csa_sada_theo<EncVector,RankSupport>::end()const
{
    return const_iterator(this, size());
}

template<class EncVector, class RankSupport>
csa_sada_theo<EncVector, RankSupport>::~csa_sada_theo()
{
    // free data structure in that order: m_b_rank depends on m_b
    if (m_b_rank != NULL) {
        delete [] m_b_rank;
        m_b_rank = NULL;
    }
    if (m_b		!= NULL) {
        delete [] m_b;
        m_b = NULL;
    }
    if (m_psi	!= NULL) {
        delete [] m_psi;
        m_psi = NULL;
    }
}

template<class EncVector, class RankSupport>
template<typename RandomAccessContainer>
void csa_sada_theo<EncVector,RankSupport>::construct(RandomAccessContainer& psi, typename RandomAccessContainer::value_type isa_0)
{
    // first free the data structure if it is already in use
    if (m_b_rank != NULL) {
        delete [] m_b_rank;
        m_b_rank = NULL;
    }
    if (m_b		!= NULL) {
        delete [] m_b;
        m_b = NULL;
    }
    if (m_psi	!= NULL) {
        delete [] m_psi;
        m_psi = NULL;
    }
    m_sa.resize(0); m_last_pos.resize(0); m_level = 0;
    // if psi is empty there is nothing to do....
    if (psi.empty()) {
        return;
    }
    typename RandomAccessContainer::size_type n = psi.size();
    EncVector psiz(psi);// Store psi compressed
    psi.resize(0);		// free the memory for the uncompressed psi array
    m_level = 0;
    while (bit_magic::l1BP(n)+1 > (1ULL<<m_level))
        ++m_level;
    // special case if m_level equals 0!
    if (m_level==0) {
        algorithm::psi2sa(psiz, isa_0, m_sa);
        m_b =  new int_vector<1>[1]; m_b[0].resize(1); m_b[0][0] = 1;
        return;
    }
    m_psi		= new EncVector[m_level];
    m_b 		= new int_vector<1>[m_level];
    m_b_rank 	= new RankSupport[m_level];
    m_last_pos	= int_vector<>(m_level, 0, bit_magic::l1BP(n)+1);
    util::set_zero_bits(m_last_pos);

    size_type nn = n;
    int_vector<> temp_psi;
    for (size_type l = 0; l < m_level; ++l) { // build data structure for each level
        m_b[l].resize(nn); util::set_zero_bits(m_b[l]); // initialize b to zero values
        typename EncVector::value_type isa_2k = isa_0; // variable for values of isa[0], isa[2],...,isa[2k]
        // walk through the conceptual suffix array and calculate b at level l
        for (size_type k = 0; k < nn-2; k+=2, isa_2k = psiz[psiz[isa_2k]]) { // isa[2(k+1)] = \psi(\psi(isa[2k]))
            m_b[l][ isa_2k ] = 1;
        }
        m_b[l][ isa_2k ] = 1;
        m_last_pos[l] = nn&1 ? isa_2k : psiz[isa_2k]; // save position of entry nn-1 in SA
        util::init_support(m_b_rank[l], &m_b[l]); // build rank structure for b
        temp_psi.set_int_width(bit_magic::l1BP(nn)+1);
        temp_psi.resize(nn/2);
        for (size_type i=0, j=0; i<nn; ++i) {
            // collect the psi values where SA[i] is odd
            if (!m_b[l][i]) temp_psi[j++] = psiz[i];
        }
        m_psi[l] = EncVector(temp_psi); // save compressed psi at level l
        // calculate psi function of the next level
        temp_psi = int_vector<>((nn+1)/2, 0,bit_magic::l1BP((nn+1)/2)+1);
        isa_2k	 = isa_0;//psiz[0];
        for (size_type k=1, t; k+1 <= nn; k+=2) {
            temp_psi[ m_b_rank[l].rank(isa_2k) ] = m_b_rank[l].rank((t=psiz[psiz[isa_2k]]));
            isa_2k = t;
        }
        // if nn is odd add one entry
        if (nn&1) temp_psi[ m_b_rank[l].rank(isa_2k) ] = m_b_rank[l].rank(psiz[isa_2k]);
        nn = (nn+1)/2; // calculate new nn

        isa_0 = m_b_rank[l].rank(isa_0);
        int_vector<> ttt; algorithm::psi2sa(temp_psi, isa_0, ttt);
        psiz = EncVector(temp_psi);
    }
    algorithm::psi2sa(temp_psi, isa_0, m_sa);// calc sa of the last level
}

template<class EncVector, class RankSupport>
template<typename RandomAccessContainer>
csa_sada_theo<EncVector,RankSupport>::csa_sada_theo(const RandomAccessContainer& sa, const unsigned char* str):m_level(0),m_psi(NULL),m_b(NULL),m_b_rank(NULL),char2comp(m_char2comp),comp2char(m_comp2char),C(m_C)
{
    assert(str != NULL);
    size_type n = 0;
    n = strlen((const char*)str)+1;
    setText(str, n);
    assert(n == sa.size());
    typename RandomAccessContainer::value_type isa_0 = 0;
    for (typename RandomAccessContainer::const_iterator it = sa.begin(), end = sa.end(); it!=end; ++it)
        if (*it==0) {
            isa_0 = it-sa.begin();
            break;
        }
//std::cerr<<"verify sa = "<<(SDSAlgorithm::verify_sa(str, n, sa))<<std::endl;
    int_vector<> psi(sa.size(), 0, bit_magic::l1BP(sa.size())+1);
    algorithm::sa2psi(sa, psi); // calculate psi out of sa
    construct(psi, isa_0);
}

template<class EncVector, class RankSupport>
typename csa_sada_theo<EncVector, RankSupport>::value_type csa_sada_theo<EncVector, RankSupport>::operator[](size_type i)const
{
    value_type result= 0;
    uint64_t history = 0;
    uint8_t  history_len = 0;
    uint8_t l 	 	 = 0;

    for (l=0; l < m_level;) {
        history <<= 1;
        ++history_len;
        if (m_b[l][i]) {
            i = m_b_rank[l].rank(i);
            ++history;
            ++l;
        } else {
            if (i==m_last_pos[l]) {
                result = m_b[l].size()-1;
                history >>= 1; --history_len;
                goto ready;
            }
            i = m_psi[l][i-m_b_rank[l].rank(i)];
        }
    }
    result = m_sa[i];
ready:
    while (history_len-- > 0) {
        if (history&1)
            result<<=1;
        else
            result-=1;
        history>>=1;
    }
    return result;
}


template<class EncVector, class RankSupport>
typename csa_sada_theo<EncVector, RankSupport>::size_type csa_sada_theo<EncVector, RankSupport>::serialize(std::ostream& out)const
{
    size_type written_bytes = 0;
    out.write((char*) &m_level, sizeof(m_level));
    written_bytes += sizeof(m_level);
    written_bytes += m_sa.serialize(out);
    if (m_level > 0) {
        for (size_type l=0; l<m_level; ++l) written_bytes += m_psi[l].serialize(out);
        for (size_type l=0; l<m_level; ++l) written_bytes += m_b[l].serialize(out);
        //	for(size_type l=0; l<m_level; ++l) written_bytes += m_b_rank[l].serialize(out);
        written_bytes += m_last_pos.serialize(out);
    }
    size_type wb   = sizeof(m_char2comp[0])*256;
    out.write((char*)m_char2comp, wb);
    written_bytes += wb;
    wb			   = sizeof(m_comp2char[0])*256;
    out.write((char*)m_comp2char, wb);
    written_bytes += wb;
    wb			   = sizeof(m_C[0])*257;
    out.write((char*)C, wb);
    written_bytes += wb;
    return written_bytes;
}

template<class EncVector, class RankSupport>
void csa_sada_theo<EncVector, RankSupport>::load(std::istream& in)
{
    in.read((char*) &m_level, sizeof(m_level));
    m_sa.load(in);
    if (m_psi 	 != NULL) {
        delete [] m_psi;
        m_psi = NULL;
    }
    if (m_b   	 != NULL) {
        delete [] m_b;
        m_b = NULL;
    }
    if (m_b_rank != NULL) {
        delete [] m_b_rank;
        m_b_rank = NULL;
    }
    if (m_level > 0) {
        m_psi 		= new EncVector[m_level];
        m_b   		= new int_vector<1>[m_level];
        m_b_rank 	= new RankSupport[m_level];
        for (size_type l=0; l<m_level; ++l) m_psi[l].load(in);
        for (size_type l=0; l<m_level; ++l) m_b[l].load(in);
        //	for(size_type l=0; l<m_level; ++l) m_b_rank[l].load(in, &m_b[l]);
        for (size_type l=0; l<m_level; ++l) util::init_support(m_b_rank[l], &m_b[l]);
        m_last_pos.load(in);
    }
    in.read((char*)m_char2comp, sizeof(m_char2comp[0])*256);
    in.read((char*)m_comp2char, sizeof(m_comp2char[0])*256);
    in.read((char*)m_C, sizeof(m_C[0])*257);
}

template<class EncVector, class RankSupport>
csa_sada_theo<EncVector,RankSupport>::csa_sada_theo(const unsigned char* str):m_level(0),m_psi(NULL),m_b(NULL),m_b_rank(NULL),char2comp(m_char2comp),comp2char(m_comp2char),C(m_C)
{
//	m_b = NULL; m_psi = NULL; m_b_rank = NULL;
//	assert(str!=NULL);
    csa_bitcompressed<> sa(str);
//	uint64_t n = strlen((const char*)str); // get size of the text
    setText(str, sa.size());
//	int_vector<> sa(n, 0, bit_magic::l1BP(n)+1 );
//	algorithm::calculate_sa(str, n, sa); // calculate the suffix array sa of str

    typename int_vector<>::value_type isa_0 = 0;
    for (csa_bitcompressed<>::size_type i=0; i!=sa.size(); ++i)
        if (sa[i]==0) {
            isa_0 = i;
            break;
        }
//std::cerr<<"verify sa = "<<(SDSAlgorithm::verify_sa(str, n, sa))<<std::endl;
    int_vector<> psi(sa.size(), bit_magic::l1BP(sa.size())+1);
    algorithm::sa2psi(sa, psi); // calculate psi out of sa
    sa = csa_bitcompressed<>();
    construct(psi, isa_0);
}

} // end namespace sdsl

#endif
