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
/*! \file lcp_uncompressed.hpp
    \brief lcp_uncompressed.hpp contains an implementation of an uncompressed lcp array.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_LCP_UNCOMPRESSED
#define INCLUDED_SDSL_LCP_UNCOMPRESSED

#include "lcp.hpp"
#include "int_vector.hpp"
#include "algorithms.hpp"
#include "iterators.hpp"
#include "bitmagic.hpp"
#include <iostream>
#include <algorithm> // for lower_bound
#include <cassert>
#include <cstring> // for strlen
#include <iomanip>
#include <iterator>
#include <vector>
#include <utility> // for pair

namespace sdsl
{

//! A class which stores the lcp array uncompressed.
template<uint8_t width=0>
class lcp_uncompressed
{
    public:
        typedef typename int_vector<width>::value_type		 	value_type;	// STL Container requirement
        typedef typename int_vector<width>::size_type		 	size_type;		// STL Container requirement
        typedef random_access_const_iterator<lcp_uncompressed>  const_iterator;// STL Container requirement
        typedef const_iterator 								    iterator;		// STL Container requirement
        typedef const value_type								const_reference;
        typedef const_reference								 	reference;
        typedef const_reference*							 	pointer;
        typedef const pointer								 	const_pointer;
        typedef ptrdiff_t  									 	difference_type; // STL Container requirement

        typedef lcp_plain_tag									lcp_category;

        enum { fast_access = 1,
               text_order  = 0,
               sa_order	  = 1
             };

        template<class Cst>  // template inner class which is used in CSTs to parametrize lcp classes
        class type           // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef lcp_uncompressed lcp_type;
        };

    private:
        int_vector<width>  m_lcp;

        void copy(const lcp_uncompressed& lcp_c) {
            m_lcp	= lcp_c.m_lcp;
        }
    public:
        //! Default Constructor
        lcp_uncompressed() {}
        //! Default Destructor
        ~lcp_uncompressed() {}
        //! Copy constructor
        lcp_uncompressed(const lcp_uncompressed& lcp_c) {
            copy(lcp_c);
        }

        //! Constructor for the compressed lcp from a (compressed) suffix array.
        template<class Text, class Sa>
        lcp_uncompressed(const Text& text, const Sa& sa);

        //! Construct the lcp array from an int_vector_file_buffer
        template<uint8_t int_width, class size_type_class>
        lcp_uncompressed(int_vector_file_buffer<int_width, size_type_class>& lcp_buf);

        template<class Text, class Sa>
        void construct(const Text& text, const Sa& sa);

        template<uint8_t int_width, class size_type_class>
        void construct(int_vector_file_buffer<int_width, size_type_class>& lcp_buf);

        //! Number of elements in the instance.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         */
        size_type size()const {
            return m_lcp.size();
        }

        //! Returns the largest size that lcp_uncompressed can ever have.
        /*! Required for the Container Concept of the STL.
         *  \sa size
         */
        static size_type max_size() {
            return int_vector<width>::max_size();
        }

        //! Returns if the data strucutre is empty.
        /*! Required for the Container Concept of the STL.A
         * \sa size
         */
        bool empty()const {
            return m_lcp.empty();
        }

        //! Swap method for lcp_uncompressed
        /*! The swap method can be defined in terms of assignment.
        	This requires three assignments, each of which, for a container type, is linear
        	in the container's size. In a sense, then, a.swap(b) is redundant.
        	This implementation guaranties a run-time complexity that is constant rather than linear.
        	\param lcp_c lcp_uncompressed to swap.

        	Required for the Assignable Conecpt of the STL.
          */
        void swap(lcp_uncompressed& lcp_c);

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
        lcp_uncompressed& operator=(const lcp_uncompressed& lcp_c);

        //! Equality Operator
        /*! Two Instances of lcp_uncompressed are equal if
         *  all their members are equal.
         *  \par Required for the Equality Comparable Concept of the STL.
         *  \sa operator!=
         */
        bool operator==(const lcp_uncompressed& lcp_c)const;

        //! Inequality Operator
        /*! Two Instances of lcp_uncompressed are equal if
         *  not all their members are equal.
         *  \par Required for the Equality Comparable Concept of the STL.
         *  \sa operator==
         */
        bool operator!=(const lcp_uncompressed& lcp_c)const;

        //! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const;

        //! Load from a stream.
        /*! \param in Inputstream to load the data structure from.
         */
        void load(std::istream& in);
};

// == template functions ==

template<uint8_t width>
template<class Text, class Sa>
lcp_uncompressed<width>::lcp_uncompressed(const Text& text, const Sa& sa)
{
    construct(text, sa);
}

template<uint8_t width>
template<class Text, class Sa>
void lcp_uncompressed<width>::construct(const Text& text, const Sa& sa)
{
    if (sa.size() == 0) {
        return;
    }
    m_lcp = int_vector<width>(sa.size(), 0, width);

    // use Kasai algorithm to compute the lcp values
    typename int_vector<width>::size_type sa_1 = sa(0), l=0;
    for (typename int_vector<width>::size_type i=0,j=0; i < sa.size(); ++i) {
        if (l) --l;
        if (sa_1) {
            j = sa[sa_1-1];
            while (i+l < sa.size() and j+l < sa.size() and text[i+l]==text[j+l]) ++l;
        } else {
            l = 0;
        }
        m_lcp[ sa_1 ] = l;
        sa_1 = sa.psi[sa_1];
    }
}

template<uint8_t width>
template<uint8_t int_width, class size_type_class>
lcp_uncompressed<width>::lcp_uncompressed(int_vector_file_buffer<int_width, size_type_class>& lcp_buf)
{
    construct(lcp_buf);
}

template<uint8_t width>
template<uint8_t int_width, class size_type_class>
void lcp_uncompressed<width>::construct(int_vector_file_buffer<int_width, size_type_class>& lcp_buf)
{
    lcp_buf.reset();
    m_lcp = int_vector<width>(lcp_buf.int_vector_size, 0, lcp_buf.int_width);
    for (size_type i=0, r_sum=0, r = lcp_buf.load_next_block(); r_sum < m_lcp.size();) {
        for (; i < r_sum+r; ++i) {
            m_lcp[i] = lcp_buf[i-r_sum];
        }
        r_sum	+=	r;
        r 		=	lcp_buf.load_next_block();
    }
}

template<uint8_t width>
void lcp_uncompressed<width>::swap(lcp_uncompressed& lcp_c)
{
    m_lcp.swap(lcp_c.m_lcp);
}

template<uint8_t width>
inline typename lcp_uncompressed<width>::value_type lcp_uncompressed<width>::operator[](size_type i)const
{
    return m_lcp[i];
}


template<uint8_t width>
typename lcp_uncompressed<width>::size_type lcp_uncompressed<width>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_lcp.serialize(out, child, "lcp");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<uint8_t width>
void lcp_uncompressed<width>::load(std::istream& in)
{
    m_lcp.load(in);
}


template<uint8_t width>
lcp_uncompressed<width>& lcp_uncompressed<width>::operator=(const lcp_uncompressed& lcp_c)
{
    if (this != &lcp_c) {
        copy(lcp_c);
    }
    return *this;
}


template<uint8_t width>
bool lcp_uncompressed<width>::operator==(const lcp_uncompressed& lcp_c)const
{
    if (this == &lcp_c)
        return true;
    return m_lcp == lcp_c.m_lcp;
}

template<uint8_t width>
bool lcp_uncompressed<width>::operator!=(const lcp_uncompressed& lcp_c)const
{
    return !(*this == lcp_c);
}

template<uint8_t width>
typename lcp_uncompressed<width>::const_iterator lcp_uncompressed<width>::begin()const
{
    return const_iterator(this, 0);
}

template<uint8_t width>
typename lcp_uncompressed<width>::const_iterator lcp_uncompressed<width>::end()const
{
    return const_iterator(this, size());
}



} // end namespace sdsl

#endif
