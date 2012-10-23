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
/*! \file lcp_support_sada.hpp
    \brief lcp_support_sada.hpp contains an implementation of a compressed lcp array.
	\author Simon Gog
*/
/* Changes:
    - Removed unnecessary rank support
 */
#ifndef INCLUDED_SDSL_LCP_SUPPORT_SADA
#define INCLUDED_SDSL_LCP_SUPPORT_SADA

#include "lcp.hpp"
#include "int_vector.hpp"
#include "algorithms.hpp"
#include "iterators.hpp"
#include "csa_sada.hpp"  // for standard template initialization of lcp_support_sada 
#include "select_support.hpp" // for standard template initialization of lcp_support_sada
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring> // for strlen
#include <iomanip>
#include <iterator>

namespace sdsl
{


//! A class for the compressed version of lcp information of an suffix array of class Csa proposed by Sadakane in the paper "Succinct Representation of lcp Information and Improvements in the Compressed Suffix Arrays".
/*! The class has two template parameters:
 *    - \a Csa is the (compressed) suffix array for which the lcp information should by provided.
 *    - \a SelectSupport is the SelectSupport class which is used in the data structure to calculate the select queries.
 *
 *	\par Space complexity
 *		 	\f$ 2n+o(n) \f$ bits, where 2n is the maximal size of the bitvector for the differences of the PLCP array
 *				and o(n) for the select support data structure.
 *
 *	The representation of the lcp information corresponds to the concept of an immutable random access container of the STL.
 */
template<class Csa = csa_sada<>, class BitVector = bit_vector, class SelectSupport = typename BitVector::select_1_type>
class _lcp_support_sada
{
    public:
        typedef typename Csa::value_type					 value_type;	// STL Container requirement
        typedef random_access_const_iterator<_lcp_support_sada>		 const_iterator;// STL Container requirement
        typedef const_iterator 								 iterator;		// STL Container requirement
        typedef const value_type							 const_reference;
        typedef const_reference								 reference;
        typedef const_reference*							 pointer;
        typedef const pointer								 const_pointer;
        typedef int_vector<>::size_type						 size_type;		// STL Container requirement
        typedef ptrdiff_t  									 difference_type; // STL Container requirement
        typedef BitVector									 bit_vector_type;
        typedef Csa											 csa_type;


        typedef lcp_permuted_tag								 lcp_category;

        enum {	fast_access = 0,
                text_order	= 1,
                sa_order	= 0
             };

        template<class Cst>  // template inner class which is used in CSTs to parametrize lcp classes
        class type           // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef _lcp_support_sada lcp_type;
        };

    private:
        //
        const Csa*		m_csa;
        bit_vector_type m_data;
        SelectSupport  	m_select_support;

        void copy(const _lcp_support_sada& lcp_c) {
            m_csa = lcp_c.m_csa;
            m_data = lcp_c.m_data;
            m_select_support = lcp_c.m_select_support;
            m_select_support.set_vector(&m_data);
        }
    public:
        const Csa*& csa;
        //! Default Constructor
        _lcp_support_sada(): csa(m_csa) {}
        //! Default Destructor
        ~_lcp_support_sada() {}
        //! Copy constructor
        _lcp_support_sada(const _lcp_support_sada& lcp_c):csa(m_csa) {
            copy(lcp_c);
        }

        //! Constructor for the compressed lcp from a suffix array and a text.
        template<class Text, class Sa>
        _lcp_support_sada(const Text& text, const Sa& sa, const Csa* csa);


        //! Construct the lcp array from an lcp array and an int_vector_file_buffer of the inverse suffix array
        template<uint8_t int_width, class size_type_class, uint8_t int_width_1, class size_type_class_1>
        _lcp_support_sada(int_vector_file_buffer<int_width, size_type_class>& lcp_buf,
                          int_vector_file_buffer<int_width_1, size_type_class_1>& isa_buf,
                          const Csa* f_csa);

        void set_csa(const Csa* f_csa) {
            if (f_csa==NULL) {
                std::cerr<<"_lcp_support_sada: Warnung: set m_csa to NULL"<<std::endl;
            }
            m_csa = f_csa;
        }

        //! Number of elements in the instance.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         */
        size_type size()const {
            return m_csa->size();
        }

        //! Returns the largest size that _lcp_support_sada can ever have.
        /*! Required for the Container Concept of the STL.
         *  \sa size
         */
        static size_type max_size() {
            return Csa::max_size();
        }

        //! Returns if the data strucutre is empty.
        /*! Required for the Container Concept of the STL.A
         * \sa size
         */
        bool empty()const {
            return m_csa->empty();
        }

        //! Swap method for _lcp_support_sada
        /*! The swap method can be defined in terms of assignment.
        	This requires three assignments, each of which, for a container type, is linear
        	in the container's size. In a sense, then, a.swap(b) is redundant.
        	This implementation guaranties a run-time complexity that is constant rather than linear.
        	\param lcp_c _lcp_support_sada to swap.

        	Required for the Assignable Conecpt of the STL.
          */
        void swap(_lcp_support_sada& lcp_c);

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
        _lcp_support_sada& operator=(const _lcp_support_sada& lcp_c);

        //! Equality Operator
        /*! Two Instances of _lcp_support_sada are equal if
         *  all their members are equal.
         *  \par Required for the Equality Comparable Concept of the STL.
         *  \sa operator!=
         */
        bool operator==(const _lcp_support_sada& lcp_c)const;

        //! Unequality Operator
        /*! Two Instances of _lcp_support_sada are equal if
         *  not all their members are equal.
         *  \par Required for the Equality Comparable Concept of the STL.
         *  \sa operator==
         */
        bool operator!=(const _lcp_support_sada& lcp_c)const;

        //! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const;

        //! Load from a stream.
        /*! \param in Inputstream to load the data structure from.
         *  \param csa Compressed Suffix Array that is supported.
         */
        void load(std::istream& in, const Csa* csa);
};

// == template functions ==

template<class Csa, class BitVector, class SelectSupport>
template<class Text, class Sa>
_lcp_support_sada<Csa, BitVector, SelectSupport>::_lcp_support_sada(const Text& text, const Sa& sa, const Csa* csa):csa(m_csa)
{
    set_csa(csa);
#ifdef SDSL_DEBUG
    std::cerr<<"start building _lcp_support_sada"<<std::endl;
#endif
    bit_vector data = bit_vector(2*sa.size(), 0);
    size_type data_cnt = 0;
    for (typename Csa::size_type i=0,j=0, sa_1 = sa(0), l=0, oldl=1; i < sa.size(); ++i) {
        if (l) --l;
        if (sa_1) {
            j = sa[sa_1-1];
            while (i+l < sa.size() and j+l < sa.size() and text[i+l]==text[j+l]) ++l;
        } else {
            l = 0;
        }
        data_cnt += l-oldl+1;
        data[data_cnt]=1;
        ++data_cnt;
        sa_1 = sa.psi[sa_1];
        oldl = l;
    }
    data.resize(data_cnt);
    util::assign(m_data, data);
    util::init_support(m_select_support, &m_data);
#ifdef SDSL_DEBUG
    std::cerr<<"finished building _lcp_support_sada"<<std::endl;
#endif
}

template<class Csa, class BitVector, class SelectSupport>
template<uint8_t int_width, class size_type_class, uint8_t int_width_1, class size_type_class_1>
_lcp_support_sada<Csa, BitVector, SelectSupport>::_lcp_support_sada(int_vector_file_buffer<int_width, size_type_class>& lcp_buf,
        int_vector_file_buffer<int_width_1, size_type_class_1>& isa_buf,
        const Csa* f_csa):csa(m_csa)
{
    typedef typename Csa::size_type size_type;
    set_csa(f_csa);
#ifdef SDSL_DEBUG
    std::cerr<<"start building _lcp_support_sada"<<std::endl;
#endif
    int_vector<int_width, size_type_class> lcp;
    util::load_from_file(lcp, lcp_buf.file_name.c_str());
    isa_buf.reset();
    size_type n = lcp.size();
    bit_vector data = bit_vector(2*n, 0);
    size_type data_cnt=0;
    for (size_type i=0, r_sum=0, r = isa_buf.load_next_block(), l=0, old_l=1; r_sum < n;) {
        for (; i < r_sum+r; ++i) {
            l = lcp[isa_buf[i-r_sum]];
            data_cnt += l + 1 - old_l;
            data[data_cnt++] = 1;
            old_l = l;
        }
        r_sum	+=	r;
        r 		=	isa_buf.load_next_block();
    }
    data.resize(data_cnt);
    util::assign(m_data, data);
    util::init_support(m_select_support, &m_data);
#ifdef SDSL_DEBUG
    std::cerr<<"finished building _lcp_support_sada"<<std::endl;
#endif
}

template<class Csa, class BitVector, class SelectSupport>
void _lcp_support_sada<Csa, BitVector, SelectSupport>::swap(_lcp_support_sada& lcp_c)
{
    m_data.swap(lcp_c.m_data);
    util::swap_support(m_select_support, lcp_c.m_select_support, &m_data, &(lcp_c.m_data));
}

template<class Csa, class BitVector, class SelectSupport>
inline typename _lcp_support_sada<Csa, BitVector, SelectSupport>::value_type _lcp_support_sada<Csa, BitVector, SelectSupport>::operator[](size_type i)const
{
    size_type j = (*m_csa)[i];
    size_type s = m_select_support.select(j+1);
    return s-(j<<1);
}


template<class Csa, class BitVector, class SelectSupport>
typename _lcp_support_sada<Csa, BitVector, SelectSupport>::size_type _lcp_support_sada<Csa, BitVector, SelectSupport>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_data.serialize(out, child, "data");
    written_bytes += m_select_support.serialize(out, child, "select_support");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<class Csa, class BitVector, class SelectSupport>
void _lcp_support_sada<Csa, BitVector, SelectSupport>::load(std::istream& in, const Csa* csa)
{
    m_csa = csa;
    m_data.load(in);
    m_select_support.load(in, &m_data);
}


template<class Csa, class BitVector, class SelectSupport>
_lcp_support_sada<Csa, BitVector, SelectSupport>& _lcp_support_sada<Csa, BitVector, SelectSupport>::operator=(const _lcp_support_sada& lcp_c)
{
    if (this != &lcp_c) {
        copy(lcp_c);
    }
    return *this;
}


template<class Csa, class BitVector, class SelectSupport>
bool _lcp_support_sada<Csa, BitVector, SelectSupport>::operator==(const _lcp_support_sada& lcp_c)const
{
    if (this == &lcp_c)
        return true;
    return m_csa == lcp_c.m_csa and m_data == lcp_c.m_data and m_select_support == lcp_c.m_select_support;
}

template<class Csa, class BitVector, class SelectSupport>
bool _lcp_support_sada<Csa, BitVector, SelectSupport>::operator!=(const _lcp_support_sada& lcp_c)const
{
    return !(*this == lcp_c);
}

template<class Csa, class BitVector, class SelectSupport>
typename _lcp_support_sada<Csa, BitVector, SelectSupport>::const_iterator _lcp_support_sada<Csa, BitVector, SelectSupport>::begin()const
{
    return const_iterator(this, 0);
}

template<class Csa, class BitVector, class SelectSupport>
typename _lcp_support_sada<Csa, BitVector, SelectSupport>::const_iterator _lcp_support_sada<Csa, BitVector, SelectSupport>::end()const
{
    return const_iterator(this, size());
}


//! Helper class which provides _lcp_support_sada the context of a CSA.
template<class BitVector = bit_vector, class SelectSupport = typename BitVector::select_1_type>
class lcp_support_sada
{
    public:
        template<class Cst>  // template inner class which is used in CSTs to parametrize lcp classes
        class type           // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef _lcp_support_sada<typename Cst::csa_type, BitVector, SelectSupport> lcp_type;
        };
};

} // end namespace sdsl

#endif
