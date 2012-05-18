/* sdsl - succinct data structures library
    Copyright (C) 2010 Simon Gog

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
/* \file lcp_vlc.hpp
    \brief lcp_vlc.hpp contains an implementation of a (compressed) lcp array.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_LCP_VLC
#define INCLUDED_SDSL_LCP_VLC

#include "lcp.hpp"
#include "vlc_vector.hpp"
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

// A class for the compressed version of lcp information of an suffix array
/* We use variable length code for the lcp entries.
 *  \par Time complexity
 *		- \f$\Order{d}\f$, where d is the sample density
 */
template<class vlc_vec_type = vlc_vector<> >
class lcp_vlc
{
    public:
        typedef typename vlc_vec_type::value_type		 value_type;	// STL Container requirement
        typedef random_access_const_iterator<lcp_vlc>		 const_iterator;// STL Container requirement
        typedef const_iterator 								 iterator;		// STL Container requirement
        typedef const value_type							 const_reference;
        typedef const_reference								 reference;
        typedef const_reference*							 pointer;
        typedef const pointer								 const_pointer;
        typedef typename vlc_vec_type::size_type			 size_type;		// STL Container requirement
        typedef typename vlc_vec_type::difference_type		 difference_type; // STL Container requirement

        typedef lcp_plain_tag								 lcp_category;

        enum { fast_access = 0,
               text_order  = 0,
               sa_order	  = 1
             };

    private:

        vlc_vec_type		m_vec;

        void copy(const lcp_vlc& lcp_c) {
            m_vec = lcp_c.m_vec;
        }

    public:
        //! Default Constructor
        lcp_vlc() {}
        //! Default Destructor
        ~lcp_vlc() {}
        //! Copy constructor
        lcp_vlc(const lcp_vlc& lcp_c) {
            copy(lcp_c);
        }

        //! Constructor for the compressed lcp from a compressed suffix array.
        template<class Text, class Sa>
        lcp_vlc(const Text& text, const Sa& sa);

        //! Construct the lcp array from an int_vector_file_buffer
        template<uint8_t int_width, class size_type_class>
        lcp_vlc(int_vector_file_buffer<int_width, size_type_class>& lcp_buf);

        template<class Text, class Sa>
        void construct(const Text& text, const Sa& sa);

        template<uint8_t int_width, class size_type_class>
        void construct(int_vector_file_buffer<int_width, size_type_class>& lcp_buf);

        //! Number of elements in the instance.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         */
        size_type size()const {
            return m_vec.size();
        }

        //! Returns the largest size that lcp_vlc can ever have.
        /*! Required for the Container Concept of the STL.
         *  \sa size
         */
        static size_type max_size() {
            return vlc_vec_type::max_size();
        }

        //! Returns if the data strucutre is empty.
        /*! Required for the Container Concept of the STL.A
         * \sa size
         */
        bool empty()const {
            return m_vec.empty();
        }

        //! Swap method for lcp_vlc
        /*! The swap method can be defined in terms of assignment.
        	This requires three assignments, each of which, for a container type, is linear
        	in the container's size. In a sense, then, a.swap(b) is redundant.
        	This implementation guaranties a run-time complexity that is constant rather than linear.
        	\param lcp_c lcp_vlc to swap.

        	Required for the Assignable Conecpt of the STL.
          */
        void swap(lcp_vlc& lcp_c);

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
        lcp_vlc& operator=(const lcp_vlc& lcp_c);

        //! Equality Operator
        /*! Two Instances of lcp_vlc are equal if
         *  all their members are equal.
         *  \par Required for the Equality Comparable Concept of the STL.
         *  \sa operator!=
         */
        bool operator==(const lcp_vlc& lcp_c)const;

        //! Unequality Operator
        /*! Two Instances of lcp_vlc are equal if
         *  not all their members are equal.
         *  \par Required for the Equality Comparable Concept of the STL.
         *  \sa operator==
         */
        bool operator!=(const lcp_vlc& lcp_c)const;

        //! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out) const;

        //! Load from a stream.
        /*! \param in Inputstream to load the data structure from.
         */
        void load(std::istream& in);
#ifdef MEM_INFO
        void mem_info(std::string label="")const {
            if (label=="")
                label = "lcp";
            size_type bytes = util::get_size_in_bytes(*this);
            std::cout << "list(label = \""<<label<<"\", size = "<< bytes/(1024.0*1024.0) <<")\n";
        }
#endif
};

// == template functions ==

template<class vlc_vec_type>
template<class Text, class Sa>
lcp_vlc<vlc_vec_type>::lcp_vlc(const Text& text, const Sa& sa)
{
    construct(text, sa);
}

template<class vlc_vec_type>
template<class Text, class Sa>
void lcp_vlc<vlc_vec_type>::construct(const Text& text, const Sa& sa)
{
    if (sa.size() == 0) {
        return;
    }
    int_vector<> lcp;
    lcp.set_int_width(bit_magic::l1BP(sa.size()) + 1);
    lcp.resize(sa.size());

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
        lcp[ sa_1 ] = l;
        sa_1 = sa.psi[sa_1];
    }
    m_vec = vlc_vec_type(lcp);
}


template<class vlc_vec_type>
template<uint8_t int_width, class size_type_class>
lcp_vlc<vlc_vec_type>::lcp_vlc(int_vector_file_buffer<int_width, size_type_class>& lcp_buf)
{
    construct(lcp_buf);
}

template<class vlc_vec_type>
template<uint8_t int_width, class size_type_class>
void lcp_vlc<vlc_vec_type>::construct(int_vector_file_buffer<int_width, size_type_class>& lcp_buf)
{
    m_vec.init(lcp_buf);
}

template<class vlc_vec_type>
void lcp_vlc<vlc_vec_type>::swap(lcp_vlc& lcp_c)
{
    m_vec.swap(lcp_c.m_vec);
}

template<class vlc_vec_type>
inline typename lcp_vlc<vlc_vec_type>::value_type lcp_vlc<vlc_vec_type>::operator[](size_type i)const
{
    return m_vec[i];
}


template<class vlc_vec_type>
typename lcp_vlc<vlc_vec_type>::size_type lcp_vlc<vlc_vec_type>::serialize(std::ostream& out) const
{
    size_type written_bytes = 0;
    written_bytes += m_vec.serialize(out);
    return written_bytes;
}

template<class vlc_vec_type>
void lcp_vlc<vlc_vec_type>::load(std::istream& in)
{
    m_vec.load(in);
}


template<class vlc_vec_type>
lcp_vlc<vlc_vec_type>& lcp_vlc<vlc_vec_type>::operator=(const lcp_vlc& lcp_c)
{
    if (this != &lcp_c) {
        copy(lcp_c);
    }
    return *this;
}


template<class vlc_vec_type>
bool lcp_vlc<vlc_vec_type>::operator==(const lcp_vlc& lcp_c)const
{
    if (this == &lcp_c)
        return true;
    return m_vec == lcp_c.m_vec;
}

template<class vlc_vec_type>
bool lcp_vlc<vlc_vec_type>::operator!=(const lcp_vlc& lcp_c)const
{
    return !(*this == lcp_c);
}

template<class vlc_vec_type>
typename lcp_vlc<vlc_vec_type>::const_iterator lcp_vlc<vlc_vec_type>::begin()const
{
    return const_iterator(this, 0);
}

template<class vlc_vec_type>
typename lcp_vlc<vlc_vec_type>::const_iterator lcp_vlc<vlc_vec_type>::end()const
{
    return const_iterator(this, size());
}



} // end namespace sdsl

#endif
