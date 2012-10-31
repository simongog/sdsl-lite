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
/*! \file csa_bitcompressed.hpp
    \brief csa_bitcompressed.hpp contains an implementation of an bitcompressed suffix array providing information of the suffix array, the inverse suffix array and the psi function.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_CSA_UNCOMPRESSED
#define INCLUDED_SDSL_CSA_UNCOMPRESSED

#include "sdsl_concepts.hpp"
#include "suffixarray_helper.hpp"
#include "int_vector.hpp"
#include "algorithms.hpp"
#include "iterators.hpp"
#include "util.hpp"
#include "testutils.hpp"
#include "csa_sampling_strategy.hpp"
#include "csa_alphabet_strategy.hpp"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring> // for strlen
#include <string>
#include <iomanip>
#include <iterator>


namespace sdsl
{

//! A class for the uncompressed suffix array (SA).
/*!
 * This class stores the information of the suffix array and the inverse suffix array in uncompressed form.
 * In contrast to this class, classes like  sdsl::csa_sada_theo, sdsl::csa_sada, and sdsl::csa_wt store
 * the suffix array and inverse suffix array data  in compressed form.
 *
 * The interface of this class is exactly the same as for the compressed indexes. This is the reason
 * why it is in the group of compressed suffix arrays.
 *
 * \par Space complexity
 *		\f$ 2n\cdot \log n\f$ bits, where \f$n\f$ equals the \f$size()\f$ of the suffix array.
 * @ingroup csa
 */
template<class AlphabetStrategy=byte_alphabet_strategy>	
class csa_bitcompressed
{
    public:
        typedef uint64_t										value_type;	// STL Container requirement
        typedef random_access_const_iterator<csa_bitcompressed> const_iterator;// STL Container requirement
        typedef const_iterator 									iterator;		// STL Container requirement
        typedef const value_type								const_reference;
        typedef const_reference									reference;
        typedef const_reference*								pointer;
        typedef const pointer									const_pointer;
        typedef int_vector<>::size_type							size_type;		// STL Container requirement
        typedef size_type				 						csa_size_type;
        typedef ptrdiff_t  										difference_type; // STL Container requirement
        typedef psi_of_sa_and_isa<csa_bitcompressed>			psi_type;
        typedef bwt_of_csa_psi<csa_bitcompressed>				bwt_type;
        typedef _sa_order_sampling_strategy<1,0>				sa_sample_type;
        typedef int_vector<>									isa_sample_type;
		typedef AlphabetStrategy								alphabet_type;
        typedef typename alphabet_type::char_type 				char_type; // Note: This is the char type of the CSA not the WT!
		typedef typename alphabet_type::comp_char_type			comp_char_type;
        typedef const char_type*								pattern_type;

        typedef csa_tag											index_category;

        enum { sa_sample_dens = 1,
               isa_sample_dens = 1
             };

        friend class psi_of_sa_and_isa<csa_bitcompressed>;
        friend class bwt_of_csa_psi<csa_bitcompressed>;
    private:
        sa_sample_type	m_sa;  // vector for suffix array values
        isa_sample_type	m_isa; // vector for inverse suffix array values
        psi_type		m_psi; // wrapper class for psi function values
		alphabet_type   m_alphabet;

        void copy(const csa_bitcompressed& csa){
		    m_sa 		 = csa.m_sa;
			m_isa 		 = csa.m_isa;
			m_alphabet   = csa.m_alphabet;
			m_psi 		 = psi_type(this);
		}
    public:
		const typename alphabet_type::char2comp_type&   char2comp;
		const typename alphabet_type::comp2char_type&  	comp2char;
		const typename alphabet_type::C_type& 			C;
		const typename alphabet_type::sigma_type& 		sigma;
        const psi_type& 								psi;
        const bwt_type 									bwt;
        const sa_sample_type& 							sa_sample;
        const isa_sample_type& 							isa_sample;

        //! Default Constructor
        csa_bitcompressed() :char2comp(m_alphabet.char2comp), comp2char(m_alphabet.comp2char), C(m_alphabet.C), sigma(m_alphabet.sigma),
		                   psi(m_psi), bwt(this), sa_sample(m_sa), isa_sample(m_isa) {}

        //! Copy constructor
        csa_bitcompressed(const csa_bitcompressed& csa) : char2comp(m_alphabet.char2comp), comp2char(m_alphabet.comp2char), C(m_alphabet.C), sigma(m_alphabet.sigma), 
		                                               psi(m_psi), bwt(this), sa_sample(m_sa), isa_sample(m_isa) {
            copy(csa);
        }

		//! Constructor
        csa_bitcompressed(tMSS& file_map, const std::string& dir, const std::string& id) : char2comp(m_alphabet.char2comp), comp2char(m_alphabet.comp2char), 
																						  C(m_alphabet.C), sigma(m_alphabet.sigma), psi(m_psi), bwt(this), 
																						  sa_sample(m_sa), isa_sample(m_isa) 
		{
			int_vector_file_buffer<8> text_buf(file_map[constants::KEY_TEXT].c_str());
			int_vector_file_buffer<>  sa_buf(file_map[constants::KEY_SA].c_str());
			size_type n = text_buf.int_vector_size;
			util::assign(m_alphabet, alphabet_type(text_buf, n));
			util::assign(m_sa, sa_sample_type(sa_buf));
			algorithm::set_isa_samples<csa_bitcompressed>(sa_buf, m_isa);
			m_psi = psi_type(this);
			write_R_output("csa", "store ISA","begin",1,0);
			std::string isa_file = dir+constants::KEY_ISA+"_"+id;
			if (!util::store_to_file(m_isa, isa_file.c_str(), true)) {
				throw std::ios_base::failure("#csa_bitcompressed: Cannot store ISA to file system!");
			} else {
				file_map[constants::KEY_ISA] = isa_file;
			}
			write_R_output("csa", "store ISA","end",1,0);
		}


        //! Number of elements in the instance.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         */
        size_type size()const {
            return m_sa.size();
        }

        //! Returns the largest size that csa_bitcompressed can ever have.
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
            return m_sa.empty();
        }

        //! Swap method for csa_bitcompressed
        /*! The swap method can be defined in terms of assignment.
        	This requires three assignments, each of which, for a container type, is linear
        	in the container's size. In a sense, then, a.swap(b) is redundant.
        	This implementation guaranties a run-time complexity that is constant rather than linear.
        	\param csa csa_bitcompressed to swap.

        	Required for the Assignable Conecpt of the STL.
          */
        void swap(csa_bitcompressed& csa){
			if (this != &csa) {
				m_sa.swap(csa.m_sa);
				m_isa.swap(csa.m_isa);
				m_alphabet.swap(csa.m_alphabet);
				m_psi = psi_type(this);
				csa.m_psi = psi_type(&csa);
			}	
		}

        //! Returns a const_iterator to the first element.
        /*! Required for the STL Container Concept.
         *  \sa end
         */
        const_iterator begin()const{
    		return const_iterator(this, 0);
		}

        //! Returns a const_iterator to the element after the last element.
        /*! Required for the STL Container Concept.
         *  \sa begin.
         */
        const_iterator end()const{
			return const_iterator(this, size());
		}

        //! []-operator
        /*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
         *
         * Required for the STL Random Access Container Concept.
         */
        inline value_type operator[](size_type i)const {
            return m_sa[i];
        }

        //! ()-operator return inverse suffix array values
        /*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
         */
        inline value_type operator()(size_type i)const {
            return m_isa[i];
        }

        //! Assignment Operator.
        /*!
         *	Required for the Assignable Concept of the STL.
         */
        csa_bitcompressed& operator=(const csa_bitcompressed& csa){
			if (this != &csa) {
        		copy(csa);
			}
			return *this;
		}

        //! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const{
			structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
			size_type written_bytes = 0;
			written_bytes += m_sa.serialize(out, child, "m_sa");
			written_bytes += m_isa.serialize(out, child, "m_isa");
			written_bytes += m_alphabet.serialize(out, child, "m_alphabet");
			structure_tree::add_size(child, written_bytes);
			return written_bytes;	
		}

        //! Load from a stream.
        /*! \param in Inputstream to load the data structure from.
         */
        void load(std::istream& in){
			m_sa.load(in);
			m_isa.load(in);
			m_alphabet.load(in);
			m_psi = psi_type(this);
		}

        size_type get_sample_dens()const {
            return 1;
        }

        //! Calculates how many symbols c are in the prefix [0..i-1] of the BWT of the original text.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in [0..size()]\f$.
         *  \param c The symbol to count the occurences in the prefix.
         *	\returns The number of occurences of symbol c in the prefix [0..i-1] of the BWT.
         *  \par Time complexity
         *		\f$ \Order{\log n} \f$
         */
        size_type rank_bwt(size_type i, const char_type c) const{
			// TODO: special case if c == BWT[i-1] we can use LF to get a constant time answer
			comp_char_type cc = char2comp[c];
			if (cc==0 and c!=0)  // character is not in the text => return 0
				return 0;
			// binary search the interval [C[cc]..C[cc+1]-1] for the result
			size_type lower_b = C[cc], upper_b = C[((size_type)1)+cc]; // lower_b inclusive, upper_b exclusive
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
         *  \param i The ith occurence. \f$i\in [1..rank(size(),c)]\f$.
         *  \param c The symbol c.
         *	\returns The ith occurence symbol c in the the BWT or size() if no ith occurence of the symbol exists in the BWT.
         *  \par Time complexity
         *		\f$ \Order{t_{\Psi}} \f$
         */
        size_type select_bwt(size_type i, const char_type c) const{
			comp_char_type cc = char2comp[c];
			if (cc==0 and c!=0)  // character is not in the text => return size()
				return size();
			if (C[cc]+i-1 <  C[((size_type)1)+cc]) {
				return m_psi[C[cc]+i-1];
			} 
			return size();	
		}
};



} // end namespace sdsl

#endif
