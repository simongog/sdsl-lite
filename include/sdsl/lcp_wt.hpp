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
/*! \file lcp_wt.hpp
    \brief lcp_wt.hpp contains an implementation of a (compressed) lcp array proposed by Simon Gog based on the Wavelet Tree.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_LCP_WT
#define INCLUDED_SDSL_LCP_WT

#include "lcp.hpp"
#include "int_vector.hpp"
#include "algorithms.hpp"
#include "iterators.hpp"
#include "bitmagic.hpp"
#include "wt_huff.hpp"
#include "select_support_bs.hpp" // dummy select support for wavelet tree, as we don't use select in this application
#include "util.hpp"
#include "testutils.hpp"
#include <iostream>
#include <algorithm> // for lower_bound
#include <cassert>
#include <cstring> // for strlen
#include <iomanip>
#include <iterator>
#include <vector>
#include <utility> // for pair
#include <stdexcept>

namespace sdsl
{

//! A class for the compressed version of lcp information of an suffix array
/*! We use \f$H_0\f$ bit for each lcp values < 255 and \f$ log n \f$ bits for each lcp value which is greater than 254.
 *  \par Time complexity
 *		- \f$\Order{1}\f$ if the value is less than 255 and
		- \f$\Order{\log n}\f$ (\f$n=size()\f$) otherwise.
 */
template<uint8_t width=0>
class lcp_wt
{
    public:
        typedef typename int_vector<width>::value_type		 value_type;	// STL Container requirement
        typedef random_access_const_iterator<lcp_wt>		 const_iterator;// STL Container requirement
        typedef const_iterator 								 iterator;		// STL Container requirement
        typedef const value_type							 const_reference;
        typedef const_reference								 reference;
        typedef const_reference*							 pointer;
        typedef const pointer								 const_pointer;
        typedef int_vector<>::size_type						 size_type;		// STL Container requirement
        typedef ptrdiff_t  									 difference_type; // STL Container requirement
        typedef select_support_bs< rank_support_v<> > 		 tDummySS;
		typedef wt_huff<bit_vector, rank_support_v<>,
				        tDummySS, tDummySS>					 small_lcp_type;

        typedef lcp_plain_tag								 lcp_category;

        enum { fast_access = 0,
               text_order  = 0,
               sa_order	  = 1
             }; // as the lcp_wt is not fast for texts with long repetition

        template<class Cst>  // template inner class which is used in CSTs to parametrize lcp classes
        class type           // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef lcp_wt lcp_type;
        };

    private:
        small_lcp_type 		m_small_lcp; // vector for lcp values < 255
        int_vector<width>   m_big_lcp;	 // vector for lcp values > 254

        typedef std::pair<size_type, size_type> tPII;
        typedef std::vector<tPII> tVPII;

        void copy(const lcp_wt& lcp_c) {
            m_small_lcp 	= lcp_c.m_small_lcp;
            m_big_lcp		= lcp_c.m_big_lcp;
        }

    public:
        //! Default Constructor
        lcp_wt() {}
        //! Default Destructor
        ~lcp_wt() {}
        //! Copy constructor
        lcp_wt(const lcp_wt& lcp_c) {
            copy(lcp_c);
        }

        //! Constructor for the compressed lcp from a compressed suffix array.
        template<class Text, class Sa>
        lcp_wt(const Text& text, const Sa& sa);

        //! Construct the lcp array from an int_vector_file_buffer
        template<uint8_t int_width, class size_type_class>
        lcp_wt(int_vector_file_buffer<int_width, size_type_class>& lcp_buf);

        //! Number of elements in the instance.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         */
        size_type size()const {
            return m_small_lcp.size();
        }

        //! Returns the largest size that lcp_wt can ever have.
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
            return m_small_lcp.size()==0;
        }

        //! Swap method for lcp_wt
        /*! The swap method can be defined in terms of assignment.
        	This requires three assignments, each of which, for a container type, is linear
        	in the container's size. In a sense, then, a.swap(b) is redundant.
        	This implementation guaranties a run-time complexity that is constant rather than linear.
        	\param lcp_c lcp_wt to swap.

        	Required for the Assignable Conecpt of the STL.
          */
        void swap(lcp_wt& lcp_c);

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
        lcp_wt& operator=(const lcp_wt& lcp_c);

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
lcp_wt<width>::lcp_wt(const Text& text, const Sa& sa)
{
    if (sa.size() == 0) {
        return;
    }
    throw std::logic_error("This constructor of lcp_wt is not yet implemented!");
}


template<uint8_t width>
template<uint8_t int_width, class size_type_class>
lcp_wt<width>::lcp_wt(int_vector_file_buffer<int_width, size_type_class>& lcp_buf)
{
    std::string temp_file = "/tmp/lcp_sml" + util::to_string(util::get_pid()) + "_" + util::to_string(util::get_id()) ;// TODO: remove absolute file name
//	write_R_output("lcp","construct sml","begin");
    typename int_vector<>::size_type l=0, max_l=0, big_sum=0, n = lcp_buf.int_vector_size;
    {
        int_vector<8> small_lcp = int_vector<8>(n);
        lcp_buf.reset();
        for (size_type i=0, r_sum=0, r = lcp_buf.load_next_block(); r_sum < n;) {
            for (; i < r_sum+r; ++i) {
                if ((l=lcp_buf[i-r_sum]) < 255) {
                    small_lcp[i] = l;
                } else {
                    small_lcp[i] = 255;
                    if (l > max_l) max_l = l;
                    ++big_sum;
                }
            }
            r_sum += r; r = lcp_buf.load_next_block();
        }
        util::store_to_file(small_lcp, temp_file.c_str());
    }
//	write_R_output("lcp","construct sml","end");
    int_vector_file_buffer<8> lcp_sml_buf(temp_file.c_str());

	util::assign( m_small_lcp, small_lcp_type(lcp_sml_buf, lcp_sml_buf.int_vector_size) );

    std::remove(temp_file.c_str());
//	write_R_output("lcp","construct big","begin");
    m_big_lcp 		= int_vector<>(big_sum, 0, bit_magic::l1BP(max_l)+1);
    {
        lcp_buf.reset();
        for (size_type i=0, ii=0, r_sum=0, r = lcp_buf.load_next_block(); r_sum < n;) {
            for (; i < r_sum+r; ++i) {
                if (lcp_buf[i-r_sum] >= 255) {
                    m_big_lcp[ ii++ ] = lcp_buf[i-r_sum];
                }
            }
            r_sum += r; r = lcp_buf.load_next_block();
        }
    }
//	write_R_output("lcp","construct big","end");
}

template<uint8_t width>
void lcp_wt<width>::swap(lcp_wt& lcp_c)
{
    m_small_lcp.swap(lcp_c.m_small_lcp);
    m_big_lcp.swap(lcp_c.m_big_lcp);
}

template<uint8_t width>
inline typename lcp_wt<width>::value_type lcp_wt<width>::operator[](size_type i)const
{
    if (m_small_lcp[i]!=255) {
        return m_small_lcp[i];
    } else {
        return m_big_lcp[ m_small_lcp.rank(i, 255) ];
    }
}


template<uint8_t width>
typename lcp_wt<width>::size_type lcp_wt<width>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_small_lcp.serialize(out, child,  "small_lcp");
    written_bytes += m_big_lcp.serialize(out, child, "large_lcp");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<uint8_t width>
void lcp_wt<width>::load(std::istream& in)
{
    m_small_lcp.load(in);
    m_big_lcp.load(in);
}


template<uint8_t width>
lcp_wt<width>& lcp_wt<width>::operator=(const lcp_wt& lcp_c)
{
    if (this != &lcp_c) {
        copy(lcp_c);
    }
    return *this;
}

template<uint8_t width>
typename lcp_wt<width>::const_iterator lcp_wt<width>::begin()const
{
    return const_iterator(this, 0);
}

template<uint8_t width>
typename lcp_wt<width>::const_iterator lcp_wt<width>::end()const
{
    return const_iterator(this, size());
}



} // end namespace sdsl

#endif
