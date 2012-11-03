/* sdsl - succinct data structures library
    Copyright (C) 2008-2012 Simon Gog

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
/*! \file csa_sada.hpp
    \brief csa_sada.hpp contains an implementation of the compressed suffix array.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_CSA_SADA
#define INCLUDED_SDSL_CSA_SADA

#include "enc_vector.hpp"
#include "int_vector.hpp"
#include "algorithms.hpp"
#include "iterators.hpp"
#include "suffixarrays.hpp"
#include "suffixarray_helper.hpp"
#include "util.hpp"
#include "testutils.hpp"
#include "bwt_construct.hpp"
#include "csa_sampling_strategy.hpp"
#include "csa_alphabet_strategy.hpp"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring> // for strlen
#include <iomanip>
#include <iterator>

namespace sdsl
{

//! A class for the Compressed Suffix Array (CSA) proposed by Sadakane for practical implementation.
/*! The CSA is parameterized with an EncVector and the sample density SampleDens (\f$s_{SA}\f$).
  * I.e. every \f$s_{SA}th\f$ value from the original suffix array is explicitly stored with \f$\log n\f$ bits.
  *
  * The EncVector (default is sdsl::enc_vector) holds the \f$\Psi\f$-function and can be parametrized with \f$s_{\Psi}\f$.
  * \todo example, code example
  *  \sa csa_sada_theo
  * @ingroup csa
 */
template<class EncVector 			= enc_vector<>, 		 // Vector type used to store the Psi-function 
	     uint32_t SampleDens 		= 32,                    // Sample density for suffix array (SA) values
		 uint32_t InvSampleDens 	= 64,                    // Sample density for inverse suffix array (ISA) values
		 class SaSamplingStrategy 	= sa_order_sa_sampling<>,// Policy class for the SA sampling. Alternative text_order_sa_sampling.
		 class IsaSampleContainer 	= int_vector<>,          // Container for the ISA samples.
		 class AlphabetStrategy   	= byte_alphabet_strategy // Policy class for the representation of the alphabet.
		>
class csa_sada {
    public:
        enum { sa_sample_dens = SampleDens,
               isa_sample_dens = InvSampleDens
             };

        typedef uint64_t											                value_type;	// STL Container requirement
        typedef random_access_const_iterator<csa_sada> 				                const_iterator;// STL Container requirement
        typedef const_iterator 										                iterator;		// STL Container requirement
        typedef const value_type									                const_reference;
        typedef const_reference										                reference;
        typedef const_reference*									                pointer;
        typedef const pointer										                const_pointer;
        typedef int_vector<>::size_type								                size_type;		// STL Container requirement
        typedef size_type 											                csa_size_type;
        typedef ptrdiff_t  											                difference_type; // STL Container requirement
        typedef EncVector											                enc_vector_type;
        typedef psi_of_csa_psi<csa_sada>						 	                psi_type;
        typedef bwt_of_csa_psi<csa_sada>						 	                bwt_type;
        typedef typename SaSamplingStrategy::template type<csa_sada>::sample_type   sa_sample_type;
        typedef IsaSampleContainer  												isa_sample_type;
		typedef AlphabetStrategy													alphabet_type;
		typedef typename alphabet_type::alphabet_category  							alphabet_category;
		typedef typename alphabet_type::comp_char_type								comp_char_type;
        typedef typename alphabet_type::char_type 								    char_type; // Note: This is the char type of the CSA not the WT!
        typedef const char_type*									                pattern_type;

        typedef csa_tag													            index_category;


        friend class psi_of_csa_psi<csa_sada>;
        friend class bwt_of_csa_psi<csa_sada>;

		static const uint32_t linear_decode_limit = 100000;
    private:
        enc_vector_type m_psi;  // psi function
        sa_sample_type 	m_sa_sample; // suffix array samples
        isa_sample_type m_isa_sample; // inverse suffix array samples
		alphabet_type   m_alphabet;   // alphabet component 

        uint64_t *m_psi_buf; //[SampleDens+1]; // buffer for decoded psi values

        void copy(const csa_sada& csa) {
            m_psi         = csa.m_psi;
            m_sa_sample   = csa.m_sa_sample;
            m_isa_sample  = csa.m_isa_sample;
			m_alphabet	  = csa.m_alphabet;
        };

		void create_buffer(){
			if ( enc_vector_type::sample_dens < linear_decode_limit  ){
				m_psi_buf = new uint64_t[enc_vector_type::sample_dens+1];
			}else{
				m_psi_buf = NULL;
			}
		}

		void delete_buffer(){
			if ( m_psi_buf != NULL ){
				delete [] m_psi_buf;
			}
		}

    public:
		const typename alphabet_type::char2comp_type&   char2comp;
		const typename alphabet_type::comp2char_type&  	comp2char;
		const typename alphabet_type::C_type& 			C;
		const typename alphabet_type::sigma_type& 		sigma;
        const psi_type 									psi;
        const bwt_type 									bwt;
        const sa_sample_type& 							sa_sample;
        const isa_sample_type& 							isa_sample;


        //! Default Constructor
        csa_sada(): char2comp(m_alphabet.char2comp), comp2char(m_alphabet.comp2char), C(m_alphabet.C), sigma(m_alphabet.sigma), 
		            psi(this), bwt(this), sa_sample(m_sa_sample), isa_sample(m_isa_sample) {
			create_buffer();
        }
        //! Default Destructor
        ~csa_sada() {
			delete_buffer();
		}

        //! Copy constructor
        csa_sada(const csa_sada& csa): char2comp(m_alphabet.char2comp), comp2char(m_alphabet.comp2char), C(m_alphabet.C), sigma(m_alphabet.sigma),
								       psi(this), bwt(this), sa_sample(m_sa_sample), isa_sample(m_isa_sample) {
			create_buffer();
            copy(csa);
        }

        csa_sada(tMSS& file_map, const std::string& dir, const std::string& id);

        //! Number of elements in the \f$\CSA\f$.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         *  \par Time complexity
         *      \f$ \Order{1} \f$
         */
        size_type size()const {
            return m_psi.size();
        }

        //! Returns the largest size that csa_sada can ever have.
        /*! Required for the Container Concept of the STL.
         *  \sa size
         */
        static size_type max_size() {
            return EncVector::max_size();
        }

        //! Returns if the data strucutre is empty.
        /*! Required for the Container Concept of the STL.A
         * \sa size
         */
        bool empty()const {
            return m_psi.empty();
        }

        //! Swap method for csa_sada
        /*! The swap method can be defined in terms of assignment.
        	This requires three assignments, each of which, for a container type, is linear
        	in the container's size. In a sense, then, a.swap(b) is redundant.
        	This implementation guaranties a run-time complexity that is constant rather than linear.
        	\param csa csa_sada to swap.

        	Required for the Assignable Conecpt of the STL.
          */
        void swap(csa_sada& csa);

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
        csa_sada& operator=(const csa_sada& csa);

        //! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const;

        //! Load from a stream.
        /*! \param in Inputstream to load the data structure from.
         */
        void load(std::istream& in);

        uint32_t get_sample_dens() const;

        uint32_t get_psi_sample_dens() const;

        //! Calculates how many symbols c are in the prefix [0..i-1] of the BWT of the original text.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in [0..size()]\f$.
         *  \param c The symbol to count the occurences in the prefix.
         *	\returns The number of occurences of symbol c in the prefix [0..i-1] of the BWT.
         *  \par Time complexity
         *		\f$ \Order{\log n t_{\Psi}} \f$
         */
        size_type rank_bwt(size_type i, const unsigned char c)const {
            unsigned char cc = char2comp[c];
            if (cc==0 and c!=0)  // character is not in the text => return 0
                return 0;
            if (i == 0)
                return 0;
            assert(i <= size());

            size_type lower_b, upper_b; // lower_b inclusive, upper_b exclusive

            const size_type sd = m_psi.get_sample_dens();
            size_type lower_sb = (C[cc]+sd-1)/sd; // lower_sb inclusive
            size_type upper_sb = (C[cc+1]+sd-1)/sd; // upper_sb exclusive
            while (lower_sb+1 < upper_sb) {
                size_type mid = (lower_sb+upper_sb)/2;
                if (m_psi.sample(mid) >= i)
                    upper_sb = mid;
                else
                    lower_sb = mid;
            }

            if (lower_sb == upper_sb) { // the interval was smaller than sd
                lower_b = C[cc]; upper_b = C[cc+1];
            } else if (lower_sb > (C[cc]+sd-1)/sd) { // main case
// TODO: don't use get_inter_sampled_values if SampleDens is really
//       large				
                lower_b = lower_sb*sd;
				if ( m_psi_buf == NULL ){
					upper_b = std::min(upper_sb*sd, C[cc+1]);
					goto finish;
				}
                uint64_t* p=m_psi_buf;
                // extract the psi values between two samples
                m_psi.get_inter_sampled_values(lower_sb, p);
                p = m_psi_buf;
                uint64_t smpl = m_psi.sample(lower_sb);
                // handle border cases
                if (lower_b + m_psi.get_sample_dens() >= C[cc+1])
                    m_psi_buf[ C[cc+1]-lower_b ] = size()-smpl;
                else
                    m_psi_buf[ m_psi.get_sample_dens() ] = size()-smpl;
                // search the result linear
                while ((*p++)+smpl < i);

                return p-1-m_psi_buf + lower_b - C[cc];
            } else { // lower_b == (m_C[cc]+sd-1)/sd and lower_sb < upper_sb
                if (m_psi.sample(lower_sb) >= i) {
                    lower_b = C[cc];
                    upper_b = lower_sb * sd + 1;
                } else {
                    lower_b = lower_sb * sd;
                    upper_b = std::min(upper_sb*sd, C[cc+1]);
                }
            }
finish:
            // binary search the interval [C[cc]..C[cc+1]-1] for the result
//			size_type lower_b = m_C[cc], upper_b = m_C[cc+1]; // lower_b inclusive, upper_b exclusive
            while (lower_b+1 < upper_b) {
                size_type mid = (lower_b+upper_b)/2;
                if (m_psi[mid] >= i)
                    upper_b = mid;
                else
                    lower_b = mid;
            }
            if (lower_b > C[cc])
                return lower_b - C[cc] + 1;
            else { // lower_b == m_C[cc]
                return m_psi[lower_b] < i;// 1 if m_psi[lower_b]<i, 0 otherwise
            }
        }

        //! Calculates the ith occurence of symbol c in the BWT of the original text.
        /*!
         *"  \param i The ith occurence. \f$i\in [1..rank(size(),c)]\f$.
         *  \param c The symbol c.
         *	\returns The ith occurence symbol c in the the BWT or size() if no ith occurence of the symbol exists in the BWT.
         *  \par Time complexity
         *		\f$ \Order{t_{\Psi}} \f$
         */
        size_type select_bwt(size_type i, const unsigned char c)const {
            assert(i > 0);
            unsigned char cc = char2comp[c];
            if (cc==0 and c!=0)  // character is not in the text => return 0
                return size();
            assert(cc != 255);
            if (C[cc]+i-1 <  C[cc+1]) {
                return m_psi[C[cc]+i-1];
            } else
                return size();
        }
};

// == template functions ==

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer, class AlphabetStrategy>
csa_sada<EncVector, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>::csa_sada(tMSS& file_map, const std::string& dir, const std::string& id):
    char2comp(m_alphabet.char2comp), comp2char(m_alphabet.comp2char),C(m_alphabet.C), sigma(m_alphabet.sigma), psi(this), bwt(this), sa_sample(m_sa_sample), isa_sample(m_isa_sample) {
	create_buffer();
    if (file_map.find(constants::KEY_BWT) == file_map.end()) { // if bwt is not already stored on disk => construct bwt
        construct_bwt(file_map, dir, id);
    }
    int_vector_file_buffer<8> bwt_buf(file_map[constants::KEY_BWT].c_str()); // 8==special case for byte alphabet/int alphabet result in 0 here
    size_type n = bwt_buf.int_vector_size;
    write_R_output("csa", "construct alphabet", "begin", 1, 0);
	util::assign(m_alphabet, alphabet_type(bwt_buf, n));
    write_R_output("csa", "construct alphabet", "end", 1, 0);

	int_vector<> cnt_chr(sigma, 0, bit_magic::l1BP(n)+1);
    for (typename alphabet_type::sigma_type i=0; i < sigma; ++i){
        cnt_chr[i] = C[i];
	}
    write_R_output("csa", "construct PSI","begin",1,0);
    // calculate psi
    {
        bwt_buf.reset();
        int_vector<> psi(n, 0, bit_magic::l1BP(n)+1);
        for (size_type i=0, r_sum=0, r=bwt_buf.load_next_block(); r_sum < n;) {
            for (; i < r_sum+r; ++i) {
                psi[ cnt_chr[ char2comp[bwt_buf[i-r_sum]] ]++ ] = i;
            }
            r_sum += r; r = bwt_buf.load_next_block();
        }
		string psi_file = dir+constants::KEY_PSI+"_"+id;
        if (!util::store_to_file(psi, psi_file.c_str())) {
            throw std::ios_base::failure("#csa_sada: Cannot store PSI to file system!");
        } else {
            file_map[constants::KEY_PSI] = psi_file;
        }
    }
    write_R_output("csa", "construct PSI","end");
    int_vector_file_buffer<> psi_buf(file_map[constants::KEY_PSI].c_str());
	write_R_output("csa", "encoded PSI", "begin");
	util::assign(m_psi, EncVector(psi_buf));
	write_R_output("csa", "encoded PSI", "end");
    int_vector_file_buffer<>  sa_buf(file_map[constants::KEY_SA].c_str());
	util::assign(m_sa_sample, sa_sample_type(sa_buf));
    algorithm::set_isa_samples<csa_sada>(sa_buf, m_isa_sample);
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer, class AlphabetStrategy>
uint32_t csa_sada<EncVector, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>::get_sample_dens()const {
    return SampleDens;
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer, class AlphabetStrategy>
uint32_t csa_sada<EncVector, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>::get_psi_sample_dens()const {
    return m_psi.get_sample_dens();
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer, class AlphabetStrategy>
typename csa_sada<EncVector, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>::const_iterator csa_sada<EncVector, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>::begin()const {
    return const_iterator(this, 0);
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer, class AlphabetStrategy>
typename csa_sada<EncVector, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>::const_iterator csa_sada<EncVector, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>::end()const {
    return const_iterator(this, size());
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer, class AlphabetStrategy>
inline typename csa_sada<EncVector, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>::value_type csa_sada<EncVector, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>::operator[](size_type i)const {
    size_type off = 0;
    while ( !m_sa_sample.is_sampled(i) ) {// while i mod SampleDens != 0 (SA[i] is not sampled)   SG: auf keinen Fall get_sample_dens nehmen, ist total langsam
        i = m_psi[i];       // go to the position where SA[i]+1 is located
        ++off;              // add 1 to the offset
    }
    value_type result = m_sa_sample.sa_value(i); 
    if (result < off) {
        return m_psi.size()-(off-result);
    } else
        return result-off;
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer, class AlphabetStrategy>
inline typename csa_sada<EncVector, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>::value_type csa_sada<EncVector, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>::operator()(size_type i)const {
    value_type result = m_isa_sample[i/InvSampleDens]; // get the rightmost sampled isa value
    i = i % InvSampleDens;
    while (i--) {
        result = m_psi[result];
    }
    return result;
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer, class AlphabetStrategy>
csa_sada<EncVector,SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>& csa_sada<EncVector, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>::operator=(const csa_sada<EncVector,SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>& csa) {
    if (this != &csa) {
        copy(csa);
    }
    return *this;
}


template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer, class AlphabetStrategy>
typename csa_sada<EncVector, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>::size_type csa_sada<EncVector, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const {
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_psi.serialize(out, child, "psi");
    written_bytes += m_sa_sample.serialize(out, child, "sa_samples");
    written_bytes += m_isa_sample.serialize(out, child, "isa_samples");
	written_bytes += m_alphabet.serialize(out, child, "alphabet");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer, class AlphabetStrategy>
void csa_sada<EncVector, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>::load(std::istream& in) {
    m_psi.load(in);
    m_sa_sample.load(in);
    m_isa_sample.load(in);
	m_alphabet.load(in);
}

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, class SaSamplingStrategy, class IsaSampleContainer, class AlphabetStrategy>
void csa_sada<EncVector, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>::swap(csa_sada<EncVector, SampleDens, InvSampleDens, SaSamplingStrategy, IsaSampleContainer, AlphabetStrategy>& csa) {
    if (this != &csa) {
        m_psi.swap(csa.m_psi);
        m_sa_sample.swap(csa.m_sa_sample);
        m_isa_sample.swap(csa.m_isa_sample);
		m_alphabet.swap(csa.m_alphabet);
    }
}

} // end namespace sdsl

#endif
