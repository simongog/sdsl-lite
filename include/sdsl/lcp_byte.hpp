/* sdsl - succinct data structures library
    Copyright (C) 2009-2013 Simon Gog

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
/*! \file lcp_byte.hpp
    \brief lcp_byte.hpp contains an implementation of a (compressed) lcp array.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_LCP_BYTE
#define INCLUDED_SDSL_LCP_BYTE

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

//! A class for a simple compressed version of LCP information
/*! Each small LCP value x=LCP[i] (\f$\leq 254\f$) is represented in a byte.
 *  For x=LCP[i] \f$\geq 255\f$ a pair (i, x) is store an list of word  pairs.
 *  This list is binary search to access LCP[i].
 *  \par Time complexity
 *		- \f$\Order{1}\f$ if the value is less than 255 and
 *		- \f$\Order{\log n}\f$ (\f$n=size()\f$) otherwise.
 * \par Reference
 *   Mohamed Ibrahim Abouelhoda, Stefan Kurtz, Enno Ohlebusch:
 *   Replacing suffix trees with enhanced suffix arrays.
 *   J. Discrete Algorithms 2(1): 53-86 (2004)
 */
template<uint8_t t_width=0>
class lcp_byte
{
    public:
        typedef typename int_vector<t_width>::value_type		 value_type;	// STL Container requirement
        typedef random_access_const_iterator<lcp_byte>		 const_iterator;// STL Container requirement
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
             }; // as the lcp_byte is not fast for texts with long repetition

        template<class Cst>  // template inner class which is used in CSTs to parametrize lcp classes
        class type           // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef lcp_byte lcp_type;
        };

    private:

        int_vector<8>      m_small_lcp;		// vector for lcp values < 255
        int_vector<t_width>   m_big_lcp;		// vector for lcp values > 254
        int_vector<t_width>   m_big_lcp_idx;   // index of the lcp entries in the lcp array

        typedef std::pair<size_type, size_type> tPII;
        typedef std::vector<tPII> tVPII;

        void copy(const lcp_byte& lcp_c) {
            m_small_lcp 	= lcp_c.m_small_lcp;
            m_big_lcp		= lcp_c.m_big_lcp;
            m_big_lcp_idx	= lcp_c.m_big_lcp_idx;
        }

    public:
        //! Default Constructor
        lcp_byte() {}

        //! Copy constructor
        lcp_byte(const lcp_byte& lcp_c) {
            copy(lcp_c);
        }

        //! Constructor
        lcp_byte(cache_config& config);

        //! Number of elements in the instance.
        size_type size()const {
            return m_small_lcp.size();
        }

        //! Returns the largest size that lcp_byte can ever have.
        static size_type max_size() {
            return int_vector<8>::max_size();
        }

        //! Returns if the data strucutre is empty.
        bool empty()const {
            return m_small_lcp.empty();
        }

        //! Swap method for lcp_byte
        void swap(lcp_byte& lcp_c);

        //! Returns a const_iterator to the first element.
        const_iterator begin()const {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        const_iterator end()const {
            return const_iterator(this, size());
        }



        //! []-operator
        /*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
         * Time complexity: O(1) for small and O(log n) for large values
         */
        inline value_type operator[](size_type i)const;

        //! Assignment Operator.
        lcp_byte& operator=(const lcp_byte& lcp_c);

        //! Serialize to a stream.
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const;

        //! Load from a stream.
        void load(std::istream& in);
};

// == template functions ==

template<uint8_t t_width>
lcp_byte<t_width>::lcp_byte(cache_config& config)
{
    int_vector_file_buffer<> lcp_buf(util::cache_file_name(constants::KEY_LCP, config));
    m_small_lcp = int_vector<8>(lcp_buf.int_vector_size);
    typename int_vector<>::size_type l=0, max_l=0, max_big_idx=0, big_sum=0;

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

template<uint8_t t_width>
void lcp_byte<t_width>::swap(lcp_byte& lcp_c)
{
    m_small_lcp.swap(lcp_c.m_small_lcp);
    m_big_lcp.swap(lcp_c.m_big_lcp);
    m_big_lcp_idx.swap(lcp_c.m_big_lcp_idx);
}

template<uint8_t t_width>
inline typename lcp_byte<t_width>::value_type lcp_byte<t_width>::operator[](size_type i)const
{
    if (m_small_lcp[i]!=255) {
        return m_small_lcp[i];
    } else {
        size_type idx = lower_bound(m_big_lcp_idx.begin(), m_big_lcp_idx.end(), i) - m_big_lcp_idx.begin();
        return m_big_lcp[idx];
    }
}

template<uint8_t t_width>
typename lcp_byte<t_width>::size_type lcp_byte<t_width>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_small_lcp.serialize(out, child, "small_lcp");
    written_bytes += m_big_lcp.serialize(out, child, "large_lcp");
    written_bytes += m_big_lcp_idx.serialize(out, child, "large_lcp_idx");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<uint8_t t_width>
void lcp_byte<t_width>::load(std::istream& in)
{
    m_small_lcp.load(in);
    m_big_lcp.load(in);
    m_big_lcp_idx.load(in);
}

template<uint8_t t_width>
lcp_byte<t_width>& lcp_byte<t_width>::operator=(const lcp_byte& lcp_c)
{
    if (this != &lcp_c) {
        copy(lcp_c);
    }
    return *this;
}

} // end namespace sdsl
#endif
