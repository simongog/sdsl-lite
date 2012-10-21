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
/*! \file csa_wt.hpp
    \brief csa_wt.hpp contains an implementation of the compressed suffix array based on a wavelet tree.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_CSA_WT
#define INCLUDED_SDSL_CSA_WT

#ifdef SDSL_DEBUG
#define SDSL_DEBUG_CSA_WT
#endif

#include "wt_huff.hpp"
#include "algorithms.hpp"
#include "iterators.hpp"
#include "util.hpp"
#include "suffixarrays.hpp"
#include "bwt_construct.hpp"
#include "fast_cache.hpp"
#include "csa_sampling_strategy.hpp"
#include <iostream>
#include <algorithm> // for std::swap
#include <cassert>
#include <cstring> // for strlen
#include <iomanip>
#include <iterator>

namespace sdsl
{

template<class WaveletTree = wt_huff<>,                       // Wavelet tree type
	     uint32_t SampleDens = 32,                            // Sample density for suffix array (SA) values
		 uint32_t InvSampleDens = 64,                         // Sample density for inverse suffix array (ISA) values
		 class SaSamplingStrategy = sa_order_sa_sampling<>,   // Policy class for the SA sampling. Alternative text_order_sa_sampling.
		 class IsaSampleContainer = int_vector<>              // Container for the ISA samples.
		 > // forward declaration
class csa_wt;


//! A wrapper class for the \f$\Psi\f$ and LF function for (compressed) suffix arrays that are based on a wavelet tree (like sdsl::csa_wt).
template<class CsaWT>
class psi_of_csa_wt
{
    public:
        typedef typename CsaWT::value_type value_type;
        typedef typename CsaWT::size_type size_type;
        typedef typename CsaWT::char_type char_type;
        typedef typename CsaWT::difference_type difference_type;
        typedef random_access_const_iterator<psi_of_csa_wt> const_iterator;// STL Container requirement
    private:
        CsaWT* m_csa_wt; //<- pointer to the (compressed) suffix array that is based on a wavelet tree
    public:

        //! Constructor
        psi_of_csa_wt(CsaWT* csa_wt=NULL) {
            m_csa_wt = csa_wt;
        }

        //! Copy constructor
        psi_of_csa_wt(const psi_of_csa_wt& psi_of_csa) {
            m_csa_wt = psi_of_csa.m_csa_wt;
        }

        //! Calculate the \f$\Psi\f$ value at position i.
        /*!	\param i The index for which the \f$\Psi\f$ value should be calculated, \f$i\in [0..size()-1]\f$.
         *	\par Time complexity
         *		\f$ \Order{\log |\Sigma|} \f$
         */
        value_type operator[](size_type i)const {
            assert(m_csa_wt != NULL);
            char_type c = algorithm::get_ith_character_of_the_first_row(i, *m_csa_wt);
            return m_csa_wt->wavelet_tree.select(i - m_csa_wt->C[m_csa_wt->char2comp[c]] + 1 , c);
        }

        //! Apply \f$\Psi\f$ k times to the value at position i.
        /*!	\param i The index for which \f$\Psi\f$ should be applied k times, \f$i\in [0..size()-1]\f$.
         *  \param k Number of times \f$\Psi\f$ should be applied
         *	\par Time complexity
         *		\f$ \Order{\min\{k\cdot \log |\Sigma|, (s_{\SUF}+s_{\ISA})\log|\Sigma|\}} \f$
         */
        value_type psi_k(size_type i, size_type k)const {
            assert(m_csa_wt != NULL);
            if (k < m_csa_wt->sa_sample_dens + m_csa_wt->isa_sample_dens) {
                for (size_type j=0; j<k; ++j) {
                    i = (*this)[i];
                }
                return i;
            } else {
                size_type x = (*m_csa_wt)[i];
                x += k;
                if (x >= m_csa_wt->size()) {
                    x -= m_csa_wt->size();
                }
                return (*m_csa_wt)(x);
            }
        }

        //! Calculate the LF mapping at position i.
        /*! The LF mapping function is the inverse to the \f$\Psi\f$ function. That is \f$LF[\Psi(i)]=i\f$.
         *  \param i The index for which the LF value should be calculated, \f$i\in [0..size()-1]\f$.
         *	\par Time complexity
         *		\f$ \Order{\log |\Sigma|} \f$
         */
        value_type operator()(size_type i)const {
            assert(m_csa_wt != NULL);
            assert(i < size());
            // (1) get the ith character in the bwt, this takes \f$\log |\Sigma| \f$ time with use of the wavelet tree
//		value_type c = m_csa_wt->m_wavelet_tree[i];
            // (2) next count how many symbols c are in the prefix [0..i-1] of the bwt, this takes \f$\log |\Sigma| \f$ time with use of the wavelet tree.
//		size_type j = m_csa_wt->m_wavelet_tree.rank(i, c);
            // TODO: (1) and (2) can be calculated in only one request to the wavelet tree -> method rank_ith_symbol !!!
            // (3) calculate the position of the j+1 th c in the sorted string of the bwt in constant time.
            typename CsaWT::char_type c;
            size_type j = m_csa_wt->m_wavelet_tree.rank_ith_symbol(i,c); // see documentation of rank_ith_symbol in wt_huff
            return m_csa_wt->C[ m_csa_wt->char2comp[c] ] + j;
        }

        //! Apply LF k times to the value at position i.
        /*!	\param i The index for which LF should be applied k times, \f$i\in [0..size()-1]\f$.
         *  \param k Number of times LF should be applied
         *	\par Time complexity
         *		\f$ \Order{\min\{k\cdot \log |\Sigma|, (s_{\SUF}+s_{\ISA})\log|\Sigma|\}} \f$
         */
        value_type lf_k(size_type i, size_type k)const {
            assert(m_csa_wt != NULL);
            if (k < m_csa_wt->sa_sample_dens + m_csa_wt->isa_sample_dens) {
                for (size_type j=0; j<k; ++j) {
                    i = (*this)(i);
                }
                return i;
            } else {
                size_type x = (*m_csa_wt)[i];
                if (x < k) {
                    x += m_csa_wt->size();
                }
                x -= k;
                return (*m_csa_wt)(x);
            }
        }

        //! Assignment operator
        psi_of_csa_wt& operator=(const psi_of_csa_wt& psi_of_csa) {
            if (this != &psi_of_csa) {
                m_csa_wt = psi_of_csa.m_csa_wt;
            }
            return *this;
        }

        //! Equality operator
        /*! return Always true, since all wrapper objects are equal only the reference to the supported csa differs.
         */
        bool operator==(const psi_of_csa_wt& psi) {
            return true;
        }

        //! Returns the size of the \f$\Psi\f$ function.
        size_type size()const {
            return m_csa_wt->size();
        }

        //! Returns if the \f$\Psi\f$ function is empty.
        size_type empty()const {
            return m_csa_wt->empty();
        }

        //! Swap operation require by the STL.
        void swap(psi_of_csa_wt& psi_of_csa) {
            if (this != &psi_of_csa) {
                ;// std::swap(m_csa_wt, psi_of_csa.m_csa_wt); // do not exchange the supported structure!
            }
        }

        //! Returns a const_iterator to the first element.
        /*! Required for the STL Container Concept.
         *  \sa end
         */
        const_iterator begin()const {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        /*! Required for the STL Container Concept.
         *  \sa begin.
         */
        const_iterator end()const {
            return const_iterator(this, size());
        }

        // Get the number of sampled \f$\Psi\f$ values. In the wavelet tree approach the number of sampled (i.e. explicitly stored) \f$\Psi\f$ values is zero.
        uint32_t get_sample_dens()const {
            return 0;
        }
};

template<class CsaWT>
class bwt_of_csa_wt
{
    public:
        typedef const typename CsaWT::char_type value_type;
        typedef typename CsaWT::size_type size_type;
        typedef typename CsaWT::difference_type difference_type;
        typedef random_access_const_iterator<bwt_of_csa_wt> const_iterator;// STL Container requirement
    private:
        CsaWT* m_csa_wt; //<- pointer to the (compressed) suffix array that is based on a wavelet tree
    public:

        //! Constructor
        bwt_of_csa_wt(CsaWT* csa_wt=NULL) {
            m_csa_wt = csa_wt;
        }

        //! Copy constructor
        bwt_of_csa_wt(const bwt_of_csa_wt& bwt_of_csa) {
            m_csa_wt = bwt_of_csa.m_csa_wt;
        }

        //! Calculate the Burrows Wheeler Transform (BWT) at position i.
        /*!	\param i The index for which the \f$\Psi\f$ value should be calculated, \f$i\in [0..size()-1]\f$.
         *	\par Time complexity
         *		\f$ \Order{\log |\Sigma|} \f$
         */
        value_type operator[](size_type i)const {
            assert(m_csa_wt != NULL);
            assert(i < size());
            return m_csa_wt->m_wavelet_tree[i];
        }

        //! Assignment operator
        bwt_of_csa_wt& operator=(const bwt_of_csa_wt& bwt_of_csa) {
            if (this != &bwt_of_csa) {
                m_csa_wt = bwt_of_csa.m_csa_wt;
            }
            return *this;
        }

        //! Returns the size of the \f$\Psi\f$ function.
        size_type size()const {
            return m_csa_wt->size();
        }

        //! Returns if the bwt function is empty.
        size_type empty()const {
            return m_csa_wt->empty();
        }

        //! Swap operation require by the STL.
        void swap(bwt_of_csa_wt& bwt_of_csa) {
            if (this != &bwt_of_csa) {
                ;// std::swap(m_csa_wt, bwt_of_csa.m_csa_wt); // do not exchange the supported structure!
            }
        }

        //! Returns a const_iterator to the first element.
        /*! Required for the STL Container Concept.
         *  \sa end
         */
        const_iterator begin()const {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        /*! Required for the STL Container Concept.
         *  \sa begin.
         */
        const_iterator end()const {
            return const_iterator(this, size());
        }
};


//! A class for the Compressed Suffix Array (CSA) based on a Wavelet Tree (WT) of the Burrow Wheeler Transform of the orignal text.
/*! The CSA is parameterized with an WavletTree, the sample density SampleDens (\f$s_{SA}\f$), and the sample density for inverse suffix array entries (\f$s_{SA^{-1}}\f$).
  * I.e. every \f$s_{SA}th\f$ value from the original suffix array is explicitly stored with \f$\log n\f$ bits.
  *
  * \todo example, code example
  *  \sa sdsl::csa_sada, sdsl::csa_uncompressed
  * @ingroup csa
 */
template<class WaveletTree, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer>
class csa_wt
{
    public:
        enum { sa_sample_dens = SampleDens,
               isa_sample_dens = InvSampleDens
             };

        typedef uint64_t											                value_type;	// STL Container requirement
        typedef random_access_const_iterator<csa_wt>		 		                const_iterator;// STL Container requirement
        typedef const_iterator 										                iterator;		// STL Container requirement
        typedef const value_type									                const_reference;
        typedef const_reference										                reference;
        typedef const_reference*									                pointer;
        typedef const pointer										                const_pointer;
        typedef int_vector<>::size_type								                size_type;		// STL Container requirement
        typedef size_type 											                csa_size_type;
        typedef ptrdiff_t  											                difference_type; // STL Container requirement
        typedef psi_of_csa_wt<csa_wt>								                psi_type;
        typedef bwt_of_csa_wt<csa_wt>								                bwt_type;
        typedef WaveletTree											                wavelet_tree_type;
        typedef unsigned char 										                char_type;
        typedef const char_type*									                pattern_type;
        typedef typename SaSamplingStrategy::template type<csa_wt>::sample_type     sa_sample_type;
        typedef IsaSampleContainer  												isa_sample_type;

        typedef csa_tag																index_category;

        friend class psi_of_csa_wt<csa_wt>;
        friend class bwt_of_csa_wt<csa_wt>;
//	static const uint32_t sample_dens = SampleDens;
    private:
        WaveletTree		m_wavelet_tree; // the wavelet tree
        psi_type m_psi;  // psi function
        bwt_type m_bwt;  // bwt
        sa_sample_type  m_sa_sample; // suffix array samples
        isa_sample_type m_isa_sample; // inverse suffix array samples
        int_vector<8>	m_char2comp;
        int_vector<8>	m_comp2char;
        int_vector<64>  m_C; // counts for the compact alphabet [0..sigma-1]
        uint16_t		m_sigma;
//		uint32_t m_sample_dens; // additional to SampleDens value
//#define USE_CSA_CACHE
#ifdef USE_CSA_CACHE
        mutable fast_cache csa_cache;
#endif

        void copy(const csa_wt& csa) {
            m_wavelet_tree			= csa.m_wavelet_tree;
            m_sa_sample 			= csa.m_sa_sample;
            m_isa_sample 			= csa.m_isa_sample;
            m_sigma					= csa.m_sigma;
            m_char2comp				= csa.m_char2comp;
            m_comp2char				= csa.m_comp2char;
            m_C = csa.m_C;
            m_psi = psi_type(this);
            m_bwt = bwt_type(this);
        }

    public:
        const int_vector<8>& char2comp;
        const int_vector<8>& comp2char;
        const int_vector<64>& C;
        const uint16_t& sigma;
        const psi_type& psi;
        const bwt_type& bwt;
        const sa_sample_type& sa_sample;
        const isa_sample_type& isa_sample;
        const wavelet_tree_type& wavelet_tree;

        //! Default Constructor
        csa_wt():char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma),
		         psi(m_psi), bwt(m_bwt),sa_sample(m_sa_sample), isa_sample(m_isa_sample), wavelet_tree(m_wavelet_tree) {}
        //! Default Destructor
        ~csa_wt() {}
        //! Copy constructor
        csa_wt(const csa_wt& csa): m_sa_sample(csa.m_sa_sample), m_isa_sample(csa.m_isa_sample),  char2comp(m_char2comp), 
		                           comp2char(m_comp2char), C(m_C), sigma(m_sigma), psi(m_psi), bwt(m_bwt), sa_sample(m_sa_sample), 
								   isa_sample(m_isa_sample), wavelet_tree(m_wavelet_tree) {
            copy(csa);
        }

		//! Constructor taking a file_map of containing the location of BWT and SA on disk
        csa_wt(tMSS& file_map, const std::string& dir, const std::string& id);

        //! Number of elements in the \f$\CSA\f$.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         *  \par Time complexity
         *      \f$ \Order{1} \f$
         */
        size_type size()const {
            return m_wavelet_tree.size();
        }

        //! Returns the largest size that csa_wt can ever have.
        /*! Required for the Container Concept of the STL.
         *  \sa size
         */
        static size_type max_size() {
            return bit_vector::max_size();
        }

        //! Returns if the data strucutre is empty.
        /*! Required for the Container Concept of the STL.
         * \sa size
         */
        bool empty()const {
            return m_wavelet_tree.empty();
        }

        //! Swap method for csa_wt
        /*! The swap method can be defined in terms of assignment.
        	This requires three assignments, each of which, for a container type, is linear
        	in the container's size. In a sense, then, a.swap(b) is redundant.
        	This implementation guaranties a run-time complexity that is constant rather than linear.
        	\param csa csa_wt to swap.

        	Required for the Assignable Conecpt of the STL.
          */
        void swap(csa_wt& csa);

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
         * Required for the STL Random Access Container Concept.
         * \par Time complexity
         *      \f$ \Order{s_{SA}\cdot t_{\Psi}} \f$, where every \f$s_{SA}\f$th suffix array entry is sampled and \f$t_{\Psi}\f$
         *           is the access time for an element in the \f$\Psi\f$-function.
         */
        inline value_type operator[](size_type i)const;

        //! ()-operator return inverse suffix array values
        /*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
         *   \par Time complexity
         *      \f$ \Order{s_{SA^{-1}}\cdot t_{\Psi}} \f$, where every \f$s_{SA^{-1}}\f$th suffix array entry is sampled and \f$t_{\Psi}\f$
         *           is the access time for an element in the \f$\Psi\f$-function.
         */
        inline value_type operator()(size_type i)const;

        //! Assignment Operator.
        /*!
         *	Required for the Assignable Concept of the STL.
         */
        csa_wt& operator=(const csa_wt& csa);

        //! Equality Operator
        /*! Two Instances of csa_wt are equal if
         *  all their members are equal.
         *  \par Required for the Equality Comparable Concept of the STL.
         *  \sa operator!=
         */
        bool operator==(const csa_wt& csa)const;

        //! Unequality Operator
        /*! Two Instances of csa_wt are equal if
         *  not all their members are equal.
         *  \par Required for the Equality Comparable Concept of the STL.
         *  \sa operator==
         */
        bool operator!=(const csa_wt& csa)const;

        //! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const;

        //! Load from a stream.
        /*! \param in Inputstream to load the data structure from.
         */
        void load(std::istream& in);

//		uint32_t get_sample_dens() const;
//		void set_sample_dens(const uint32_t sample_dens);

        uint32_t get_psi_sample_dens() const;
        void set_psi_sample_dens(const uint32_t sample_dens);

        //! Calculates how many symbols c are in the prefix [0..i-1] of the BWT of the original text.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in [0..size()]\f$.
         *  \param c The symbol to count the occurences in the prefix.
         *	\returns The number of occurences of symbol c in the prefix [0..i-1] of the BWT.
         *  \par Time complexity
         *		\f$ \Order{\log |\Sigma|} \f$
         */
        size_type rank_bwt(size_type i, const char_type c)const {
            return m_wavelet_tree.rank(i, c);
        }

        //! Calculates the ith occurrence of symbol c in the BWT of the original text.
        /*!
         *  \param i The ith occurrence. \f$i\in [1..rank(size(),c)]\f$.
         *  \param c The symbol c.
         *	\returns The ith occurrence symbol c in the the BWT or size() if no ith occurence of the symbol exists in the BWT.
         *  \par Time complexity
         *		\f$ \Order{t_{\Psi}} \f$
         */
        size_type select_bwt(size_type i, const char_type c)const {
            assert(i > 0);
            char_type cc = m_char2comp[c];
            if (cc==0 and c!=0)  // character is not in the text => return size()
                return size();
            assert(cc != 255);
            if (m_C[cc]+i-1 <  m_C[cc+1]) {
                return m_wavelet_tree.select(i, c);
            } else
                return size();
        }
};

// == template functions ==

template<class WaveletTree, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer>
csa_wt<WaveletTree, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>::csa_wt(tMSS& file_map, const std::string& dir, const std::string& id):char2comp(m_char2comp), comp2char(m_comp2char),C(m_C), sigma(m_sigma), psi(m_psi), bwt(m_bwt),sa_sample(m_sa_sample),isa_sample(m_isa_sample),wavelet_tree(m_wavelet_tree)
{
    if (file_map.find("bwt") == file_map.end()) { // if bwt is not already stored on disk => construct bwt
        construct_bwt(file_map, dir, id);
//		construct_bwt2(file_map, dir, id);
    }
    int_vector_file_buffer<8> bwt_buf(file_map["bwt"].c_str());
    int_vector_file_buffer<>  sa_buf(file_map["sa"].c_str());
    size_type n = bwt_buf.int_vector_size;
    algorithm::set_text<csa_wt>(bwt_buf, n, m_C, m_char2comp, m_comp2char, m_sigma);
//	m_wavelet_tree = WaveletTree(bwt_buf, n);
    write_R_output("csa", "construct WT", "begin", 1, 0);
	util::assign(m_wavelet_tree, wavelet_tree_type(bwt_buf, n)  );
    write_R_output("csa", "construct WT", "end", 1, 0);

	util::assign(m_sa_sample, sa_sample_type(sa_buf));
    algorithm::set_isa_samples<csa_wt>(sa_buf, m_isa_sample);

    m_psi = psi_type(this);
    m_bwt = bwt_type(this);
}


template<class WaveletTree, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer>
uint32_t csa_wt<WaveletTree, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>::get_psi_sample_dens()const
{
    return m_psi.get_sample_dens();
}

template<class WaveletTree, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer>
void csa_wt<WaveletTree, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>::set_psi_sample_dens(const uint32_t sample_dens)
{
}

template<class WaveletTree, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer>
typename csa_wt<WaveletTree, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>::const_iterator csa_wt<WaveletTree, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>::begin()const
{
    return const_iterator(this, 0);
}

template<class WaveletTree, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer>
typename csa_wt<WaveletTree, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>::const_iterator csa_wt<WaveletTree, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>::end()const
{
    return const_iterator(this, size());
}

template<class WaveletTree, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer>
inline typename csa_wt<WaveletTree, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>::value_type csa_wt<WaveletTree, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>::operator[](size_type i)const
{
    size_type off = 0;
    while ( !m_sa_sample.is_sampled(i) ) {
        i = m_psi(i);
        ++off;
    }
    value_type result = m_sa_sample.sa_value(i); 
    if (result + off < size()) {
        return result + off;
    } else {
        return result + off - size();
    }
    /*
    #ifdef USE_CSA_CACHE
    	size_type r = 0;
    	if( csa_cache.exists(i, r) ){
    		return r;
    	}
    #endif
    	size_type off = 0;
    	while( i % SampleDens ){// while i mod SampleDens != 0 (SA[i] is not sampled)
    		i = m_psi[i];       // go to the position where SA[i]+1 is located
    		++off;              // add 1 to the offset
    	}
    	value_type result = m_sa_sample[i/SampleDens];
    	// TODO: try LF-function for iteration
    	if( result < off ){
    #ifdef USE_CSA_CACHE
    		r = m_psi.size()-(off-result);
    		return r;
    #endif
    		return m_psi.size()-(off-result);
    	}
    	else{
    #ifdef USE_CSA_CACHE
    		r = result-off;
    		return r;
    #endif
    		return result-off;
    	}
    */
}

template<class WaveletTree, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer>
inline typename csa_wt<WaveletTree, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>::value_type csa_wt<WaveletTree, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>::operator()(size_type i)const
{
    size_type ii;
    value_type result = m_isa_sample[ ii = ((i+InvSampleDens-1)/InvSampleDens) ]; // get the leftmost sampled isa value to the right of i
    ii *= InvSampleDens;
    if (ii >= size()) {
        i = size() - 1 - i;
    } else {
        i = ii - i;
    }
    while (i--) {
        result = m_psi(result);
    }
    return result;
}

template<class WaveletTree, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer>
csa_wt<WaveletTree,SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>& csa_wt<WaveletTree, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>::operator=(const csa_wt<WaveletTree,SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>& csa)
{
    if (this != &csa) {
        copy(csa);
    }
    return *this;
}


template<class WaveletTree, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer>
typename csa_wt<WaveletTree, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>::size_type csa_wt<WaveletTree, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_wavelet_tree.serialize(out, child, "wavelet_tree");
    written_bytes += m_sa_sample.serialize(out, child, "sa_samples");
    written_bytes += m_isa_sample.serialize(out, child, "isa_samples");
    written_bytes += m_char2comp.serialize(out, child, "char2comp");
    written_bytes += m_comp2char.serialize(out, child, "comp2char");
    written_bytes += m_C.serialize(out, child, "C");
    written_bytes += util::write_member(m_sigma, out, child, "sigma");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<class WaveletTree, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer>
void csa_wt<WaveletTree, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>::load(std::istream& in)
{
    m_wavelet_tree.load(in);
    m_sa_sample.load(in);
    m_isa_sample.load(in);
    m_char2comp.load(in);
    m_comp2char.load(in);
    m_C.load(in);
    util::read_member(m_sigma, in);
    m_psi = psi_type(this);
    m_bwt = bwt_type(this);
}

template<class WaveletTree, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer>
void csa_wt<WaveletTree, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>::swap(csa_wt<WaveletTree, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer>& csa)
{
    if (this != &csa) {
        m_wavelet_tree.swap(csa.m_wavelet_tree);
        m_sa_sample.swap(csa.m_sa_sample);
        m_isa_sample.swap(csa.m_isa_sample);
        m_char2comp.swap(csa.m_char2comp);
        m_comp2char.swap(csa.m_comp2char);
        m_C.swap(csa.m_C);
        std::swap(m_sigma, csa.m_sigma);
//        m_psi.swap(csa.m_psi);
        m_psi = psi_type(this);
        csa.m_psi = psi_type(&csa);
        m_bwt = bwt_type(this);
        csa.m_bwt = bwt_type(&csa);
    }
}

} // end namespace sdsl

#endif
