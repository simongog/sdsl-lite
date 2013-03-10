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
/*! \file lcp_kurtz.hpp
    \brief lcp_kurtz.hpp contains an implementation of a (compressed) lcp array proposed by Stefan Kurtz.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_LCP_KURTZ
#define INCLUDED_SDSL_LCP_KURTZ

#include "lcp.hpp"
#include "int_vector.hpp"
#include "algorithms.hpp"
#include "iterators.hpp"
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

//! A class for the compressed version of lcp information of an suffix array proposed by Stefan Kurtz in the paper "Reducing the Space Requirement of Suffix Trees".
/*! We use 8 bit for each lcp values + \f$ log n \f$ bits for each lcp value which is greater than 254.
 *  \par Time complexity
 *		- \f$\Order{1}\f$ if the value is less than 255 and
		- \f$\Order{\log n}\f$ (\f$n=size()\f$) otherwise.
 */
template<uint8_t width=0>
class lcp_kurtz
{
    public:
        typedef typename int_vector<width>::value_type		 value_type;	// STL Container requirement
        typedef random_access_const_iterator<lcp_kurtz>		 const_iterator;// STL Container requirement
        typedef const_iterator 								 iterator;		// STL Container requirement
        typedef const value_type							 const_reference;
        typedef const_reference								 reference;
        typedef const_reference*							 pointer;
        typedef const pointer								 const_pointer;
        typedef int_vector<>::size_type						 size_type;		// STL Container requirement
        typedef ptrdiff_t  									 difference_type; // STL Container requirement

        typedef lcp_plain_tag								 lcp_category;

        enum { fast_access = 0,
               text_order  = 0,
               sa_order	  = 1
             }; // as the lcp_kurtz is not fast for texts with long repetition

        template<class Cst>  // template inner class which is used in CSTs to parametrize lcp classes
        class type           // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef lcp_kurtz lcp_type;
        };

    private:

        int_vector<8>      m_small_lcp;		// vector for lcp values < 255
        int_vector<width>   m_big_lcp;		// vector for lcp values > 254
        int_vector<width>   m_big_lcp_idx;   // index of the lcp entries in the lcp array

        typedef std::pair<size_type, size_type> tPII;
        typedef std::vector<tPII> tVPII;

        void copy(const lcp_kurtz& lcp_c) {
            m_small_lcp 	= lcp_c.m_small_lcp;
            m_big_lcp		= lcp_c.m_big_lcp;
            m_big_lcp_idx	= lcp_c.m_big_lcp_idx;
        }

    public:
        //! Default Constructor
        lcp_kurtz() {}

        //! Copy constructor
        lcp_kurtz(const lcp_kurtz& lcp_c) {
            copy(lcp_c);
        }

        //! Constructor for the compressed lcp from a compressed suffix array.
        template<class Text, class Sa>
        lcp_kurtz(const Text& text, const Sa& sa);

        //! Construct the lcp array from an int_vector_file_buffer
        template<uint8_t int_width>
        lcp_kurtz(int_vector_file_buffer<int_width>& lcp_buf);

        //! Number of elements in the instance.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         */
        size_type size()const {
            return m_small_lcp.size();
        }

        //! Returns the largest size that lcp_kurtz can ever have.
        /*! Required for the Container Concept of the STL.
         *  \sa size
         */
        static size_type max_size() {
            return int_vector<8>::max_size();
        }

        //! Returns if the data strucutre is empty.
        /*! Required for the Container Concept of the STL.A
         * \sa size
         */
        bool empty()const {
            return m_small_lcp.empty();
        }

        //! Swap method for lcp_kurtz
        /*! The swap method can be defined in terms of assignment.
        	This requires three assignments, each of which, for a container type, is linear
        	in the container's size. In a sense, then, a.swap(b) is redundant.
        	This implementation guaranties a run-time complexity that is constant rather than linear.
        	\param lcp_c lcp_kurtz to swap.

        	Required for the Assignable Conecpt of the STL.
          */
        void swap(lcp_kurtz& lcp_c);

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
        lcp_kurtz& operator=(const lcp_kurtz& lcp_c);

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
lcp_kurtz<width>::lcp_kurtz(const Text& text, const Sa& sa) {
    if (sa.size() == 0) {
        return;
    }
    m_small_lcp = int_vector<8>(sa.size());
    tVPII big_values;

    // use Kasai algorithm to compute the lcp values
    typename Sa::size_type sa_1 = sa(0), l=0, max_l=0, max_big_idx=0;
    for (typename Sa::size_type i=0,j=0; i < sa.size(); ++i) {
        if (l) --l;
        if (sa_1) {
            j = sa[sa_1-1];
            while (i+l < sa.size() and j+l < sa.size() and text[i+l]==text[j+l]) ++l;
        } else {
            l = 0;
        }
        if (l < 255) {
            m_small_lcp[ sa_1 ] = l;
        } else {
            m_small_lcp[ sa_1 ] = 255;
            big_values.push_back(tPII(sa_1, l));
            if (l > max_l) max_l = l;
            if (sa_1 > max_big_idx) max_big_idx = sa_1;
        }
        sa_1 = sa.psi[sa_1];
    }

    sort(big_values.begin(), big_values.end());
    m_big_lcp 		= int_vector<width>(big_values.size(), 0, bit_magic::l1BP(max_l)+1);
    m_big_lcp_idx 	= int_vector<width>(big_values.size(), 0, bit_magic::l1BP(max_big_idx)+1);
    for (tVPII::size_type i=0; i<big_values.size(); ++i) {
        m_big_lcp[i] 	 = big_values[i].second;
        m_big_lcp_idx[i] = big_values[i].first;
    }
}


template<uint8_t width>
template<uint8_t int_width>
lcp_kurtz<width>::lcp_kurtz(int_vector_file_buffer<int_width>& lcp_buf) {
    m_small_lcp = int_vector<8>(lcp_buf.int_vector_size);
    typename int_vector<>::size_type l=0, max_l=0, max_big_idx=0, big_sum=0;
    lcp_buf.reset();

    for (size_type i=0, r_sum=0, r = lcp_buf.load_next_block(); r_sum < m_small_lcp.size();) {
        for (; i < r_sum+r; ++i) {
            if ((l=lcp_buf[i-r_sum]) < 255) {
                m_small_lcp[i] = l;
            } else {
                m_small_lcp[i] = 255;
                if (l > max_l) max_l = l;
                max_big_idx = i;
                ++big_sum;
            }
        }
        r_sum	+=	r;
        r 		=	lcp_buf.load_next_block();
    }
    m_big_lcp 		= int_vector<>(big_sum, 0, bit_magic::l1BP(max_l)+1);
    m_big_lcp_idx 	= int_vector<>(big_sum, 0, bit_magic::l1BP(max_big_idx)+1);

    lcp_buf.reset();
    for (size_type i=0, r_sum=0, r = lcp_buf.load_next_block(),ii=0; r_sum<m_small_lcp.size();) {
        for (; i < r_sum+r; ++i) {
            if ((l=lcp_buf[i-r_sum]) >= 255) {
                m_big_lcp[ii] = l;
                m_big_lcp_idx[ii] = i;
                ++ii;
            }
        }
        r_sum	+=	r;
        r 		=	lcp_buf.load_next_block();
    }
}

template<uint8_t width>
void lcp_kurtz<width>::swap(lcp_kurtz& lcp_c) {
    m_small_lcp.swap(lcp_c.m_small_lcp);
    m_big_lcp.swap(lcp_c.m_big_lcp);
    m_big_lcp_idx.swap(lcp_c.m_big_lcp_idx);
}

template<uint8_t width>
inline typename lcp_kurtz<width>::value_type lcp_kurtz<width>::operator[](size_type i)const {
    if (m_small_lcp[i]!=255) {
        return m_small_lcp[i];
    } else {
        size_type idx = lower_bound(m_big_lcp_idx.begin(), m_big_lcp_idx.end(), i) - m_big_lcp_idx.begin();
        return m_big_lcp[idx];
    }
}


template<uint8_t width>
typename lcp_kurtz<width>::size_type lcp_kurtz<width>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const {
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_small_lcp.serialize(out, child, "small_lcp");
    written_bytes += m_big_lcp.serialize(out, child, "large_lcp");
    written_bytes += m_big_lcp_idx.serialize(out, child, "large_lcp_idx");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<uint8_t width>
void lcp_kurtz<width>::load(std::istream& in) {
    m_small_lcp.load(in);
    m_big_lcp.load(in);
    m_big_lcp_idx.load(in);
}


template<uint8_t width>
lcp_kurtz<width>& lcp_kurtz<width>::operator=(const lcp_kurtz& lcp_c) {
    if (this != &lcp_c) {
        copy(lcp_c);
    }
    return *this;
}

template<uint8_t width>
typename lcp_kurtz<width>::const_iterator lcp_kurtz<width>::begin()const {
    return const_iterator(this, 0);
}

template<uint8_t width>
typename lcp_kurtz<width>::const_iterator lcp_kurtz<width>::end()const {
    return const_iterator(this, size());
}



} // end namespace sdsl

#endif
