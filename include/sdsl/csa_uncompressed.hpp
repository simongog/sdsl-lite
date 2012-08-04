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
/*! \file csa_uncompressed.hpp
    \brief csa_uncompressed.hpp contains an implementation of an uncompressed suffix array providing information of the suffix array, the inverse suffix array and the psi function.
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
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring> // for strlen
#include <iomanip>
#include <iterator>


namespace sdsl
{

//! A class for the uncmpressed suffix array (SA).
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
class csa_uncompressed
{
    public:
        typedef uint64_t											 value_type;	// STL Container requirement
        typedef random_access_const_iterator<csa_uncompressed> 		 const_iterator;// STL Container requirement
        typedef const_iterator 										 iterator;		// STL Container requirement
        typedef const value_type									 const_reference;
        typedef const_reference										 reference;
        typedef const_reference*									 pointer;
        typedef const pointer										 const_pointer;
        typedef int_vector<>::size_type								 size_type;		// STL Container requirement
        typedef size_type				 							 csa_size_type;
        typedef ptrdiff_t  											 difference_type; // STL Container requirement
        typedef psi_of_sa_and_isa<csa_uncompressed>					 psi_type;
        typedef bwt_of_csa_psi<csa_uncompressed>					 bwt_type;
        typedef const unsigned char*								 pattern_type;
        typedef unsigned char										 char_type;
        typedef int_vector<>										 sa_sample_type;
        typedef int_vector<>										 isa_sample_type;

        typedef csa_tag												index_category;

        enum { sa_sample_dens = 1,
               isa_sample_dens = 1
             };

        friend class psi_of_sa_and_isa<csa_uncompressed>;
        friend class bwt_of_csa_psi<csa_uncompressed>;
    private:
        sa_sample_type	m_sa;  // vector for suffix array values
        isa_sample_type	m_isa; // vector for inverse suffix array values
        psi_type		m_psi; // wrapper class for psi function values
        bwt_type 		m_bwt; // wrapper class for the BWT
        int_vector<8>	m_char2comp;
        int_vector<8> 	m_comp2char;
        int_vector<64>  m_C;
        uint16_t		m_sigma;

        void construct();
        void copy(const csa_uncompressed& csa);
    public:
        const int_vector<8>& char2comp;
        const int_vector<8>& comp2char;
        const int_vector<64>& C;
        const uint16_t& sigma;
        const psi_type& psi;
        const bwt_type& bwt;
        const sa_sample_type& sa_sample;
        const isa_sample_type& isa_sample;

        //! Default Constructor
        csa_uncompressed():char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma) ,psi(m_psi), bwt(m_bwt), sa_sample(m_sa), isa_sample(m_isa) {
            util::assign(m_C, int_vector<64>(257));
            util::assign(m_char2comp, int_vector<8>(256));
            util::assign(m_comp2char, int_vector<8>(256));
        }
        //! Default Destructor
        ~csa_uncompressed() {}
        //! Copy constructor
        csa_uncompressed(const csa_uncompressed& csa):char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma), psi(m_psi), bwt(m_bwt), sa_sample(m_sa), isa_sample(m_isa) {
            copy(csa);
        }

        //! Construct csa_uncompressed from another compressed or uncompressed suffix array
        template<typename RandomAccessContainer>
        csa_uncompressed(const RandomAccessContainer& sa, const unsigned char* str);


        //! Constructor for the CSA taking a string for that the CSA should be calculated
        csa_uncompressed(const unsigned char* str);

        //! Construct the csa_uncompressed form a int_vector_file_buffer for the text and the SA
        template<uint8_t int_width, class size_type_class, class size_type_class_1>
        csa_uncompressed(int_vector_file_buffer<8, size_type_class>& text_buf,
                         int_vector_file_buffer<int_width, size_type_class_1>& sa_buf):char2comp(m_char2comp),
            comp2char(m_comp2char), C(m_C), sigma(m_sigma) ,psi(m_psi), bwt(m_bwt), sa_sample(m_sa), isa_sample(m_isa) {
            text_buf.reset(); sa_buf.reset();
            size_type n = text_buf.int_vector_size;
            algorithm::set_text<csa_uncompressed>(text_buf, n, m_C, m_char2comp, m_comp2char, m_sigma);
            algorithm::set_sa_and_isa_samples<csa_uncompressed>(sa_buf, m_sa, m_isa);
            m_psi = psi_type(this);
            m_bwt = bwt_type(this);
        }

        csa_uncompressed(tMSS& file_map, const std::string& dir, const std::string& id):char2comp(m_char2comp),
            comp2char(m_comp2char), C(m_C), sigma(m_sigma), psi(m_psi), bwt(m_bwt), sa_sample(m_sa), isa_sample(m_isa) {
            construct(file_map, dir, id);
        }


        void construct(tMSS& file_map, const std::string& dir, const std::string& id) {
            int_vector_file_buffer<8> text_buf(file_map["text"].c_str());
            int_vector_file_buffer<>  sa_buf(file_map["sa"].c_str());
            size_type n = text_buf.int_vector_size;
            algorithm::set_text<csa_uncompressed>(text_buf, n, m_C, m_char2comp, m_comp2char, m_sigma);
            algorithm::set_sa_and_isa_samples<csa_uncompressed>(sa_buf, m_sa, m_isa);
            m_psi = psi_type(this);
            m_bwt = bwt_type(this);
            write_R_output("csa", "store ISA","begin",1,0);
            if (!util::store_to_file(m_isa, (dir+"isa_"+id).c_str(), true)) {
                throw std::ios_base::failure("#csa_uncompressed: Cannot store ISA to file system!");
            } else {
                file_map["isa"] = dir+"isa_"+id;
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

        //! Returns the largest size that csa_uncompressed can ever have.
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

        //! Swap method for csa_uncompressed
        /*! The swap method can be defined in terms of assignment.
        	This requires three assignments, each of which, for a container type, is linear
        	in the container's size. In a sense, then, a.swap(b) is redundant.
        	This implementation guaranties a run-time complexity that is constant rather than linear.
        	\param csa csa_uncompressed to swap.

        	Required for the Assignable Conecpt of the STL.
          */
        void swap(csa_uncompressed& csa);

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
        csa_uncompressed& operator=(const csa_uncompressed& csa);

        //! Equality Operator
        /*! Two Instances of csa_uncompressed are equal if
         *  all their members are equal.
         *  \par Required for the Equality Comparable Concept of the STL.
         *  \sa operator!=
         */
        bool operator==(const csa_uncompressed& csa)const;

        //! Unequality Operator
        /*! Two Instances of csa_uncompressed are equal if
         *  not all their members are equal.
         *  \par Required for the Equality Comparable Concept of the STL.
         *  \sa operator==
         */
        bool operator!=(const csa_uncompressed& csa)const;

        //! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out) const;

        //! Load from a stream.
        /*! \param in Inputstream to load the data structure from.
         */
        void load(std::istream& in);

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
        size_type rank_bwt(size_type i, const unsigned char c) const;

        //! Calculates the ith occurence of symbol c in the BWT of the original text.
        /*!
         *  \param i The ith occurence. \f$i\in [1..rank(size(),c)]\f$.
         *  \param c The symbol c.
         *	\returns The ith occurence symbol c in the the BWT or size() if no ith occurence of the symbol exists in the BWT.
         *  \par Time complexity
         *		\f$ \Order{t_{\Psi}} \f$
         */
        size_type select_bwt(size_type i, const unsigned char c) const;
};


template<typename RandomAccessContainer>
csa_uncompressed::csa_uncompressed(const RandomAccessContainer& sa, const unsigned char* str):char2comp(m_char2comp), comp2char(m_comp2char),C(m_C), sigma(m_sigma), psi(m_psi), bwt(m_bwt), sa_sample(m_sa), isa_sample(m_isa)
{
    size_type n = 1;
    if (str != NULL) {
        n = strlen((const char*)str);
    }
    algorithm::set_text<csa_uncompressed>(str, n+1, m_C, m_char2comp, m_comp2char, m_sigma);
//	assert(sa.size() == n+1);
    m_sa = int_vector<>(n+1, 0, bit_magic::l1BP(n+1)+1);
    for (size_type i=0; i<n+1; ++i)
        m_sa[i] = sa[i];
    construct();
    if (n+1 > 0 and n+1 != size())
        throw std::logic_error(util::demangle(typeid(this).name())+": text size differ with sa size!");
}


} // end namespace sdsl

#endif
