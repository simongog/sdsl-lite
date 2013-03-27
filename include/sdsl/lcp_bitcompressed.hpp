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
/*! \file lcp_bitcompressed.hpp
    \brief lcp_bitcompressed.hpp contains an implementation of an uncompressed lcp array.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_LCP_BITCOMPRESSED
#define INCLUDED_SDSL_LCP_BITCOMPRESSED

#include "lcp.hpp"
#include "int_vector.hpp"
#include "algorithms.hpp"
#include "iterators.hpp"
#include <cassert>
#include <utility> // for pair

namespace sdsl
{

//! A class which stores the lcp array uncompressed.
template<uint8_t t_width=0>
class lcp_bitcompressed
{
    public:
        typedef typename int_vector<t_width>::value_type        value_type;    // STL Container requirement
        typedef typename int_vector<t_width>::size_type         size_type;        // STL Container requirement
        typedef random_access_const_iterator<lcp_bitcompressed> const_iterator;// STL Container requirement
        typedef const_iterator                                  iterator;        // STL Container requirement
        typedef const value_type                                const_reference;
        typedef const_reference                                 reference;
        typedef const_reference*                                pointer;
        typedef const pointer                                   const_pointer;
        typedef ptrdiff_t                                       difference_type; // STL Container requirement

        typedef lcp_plain_tag                                   lcp_category;

        enum { fast_access = 1,
               text_order  = 0,
               sa_order    = 1
             };

        template<class Cst>  // template inner class which is used in CSTs to parametrize lcp classes
        class type           // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef lcp_bitcompressed lcp_type;
        };

    private:
        int_vector<t_width>  m_lcp;

        void copy(const lcp_bitcompressed& lcp_c) {
            m_lcp    = lcp_c.m_lcp;
        }
    public:
        //! Default Constructor
        lcp_bitcompressed() {}

        //! Copy constructor
        lcp_bitcompressed(const lcp_bitcompressed& lcp_c) {
            copy(lcp_c);
        }

        //! Constructor taking a cache_config
        lcp_bitcompressed(cache_config& config);

        //! Number of elements in the instance.
        size_type size()const {
            return m_lcp.size();
        }

        //! Returns the largest size that lcp_bitcompressed can ever have.
        static size_type max_size() {
            return int_vector<t_width>::max_size();
        }

        //! Returns if the data strucutre is empty.
        bool empty()const {
            return m_lcp.empty();
        }

        //! Swap method for lcp_bitcompressed
        void swap(lcp_bitcompressed& lcp_c);

        //! Returns a const_iterator to the first element.
        const_iterator begin()const {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        const_iterator end()const {
            return const_iterator(this, size());
        }

        //! Access operator
        /*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
         */
        inline value_type operator[](size_type i)const;

        //! Assignment Operator.
        lcp_bitcompressed& operator=(const lcp_bitcompressed& lcp_c);

        //! Serialize to a stream.
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const;

        //! Load from a stream.
        void load(std::istream& in);
};

// == template functions ==

template<uint8_t t_width>
lcp_bitcompressed<t_width>::lcp_bitcompressed(cache_config& config)
{
    int_vector_file_buffer<> lcp_buf(cache_file_name(constants::KEY_LCP, config));
    m_lcp = int_vector<t_width>(lcp_buf.int_vector_size, 0, lcp_buf.width);
    for (size_type i=0, r_sum=0, r = lcp_buf.load_next_block(); r_sum < m_lcp.size();) {
        for (; i < r_sum+r; ++i) {
            m_lcp[i] = lcp_buf[i-r_sum];
        }
        r_sum += r;
        r      = lcp_buf.load_next_block();
    }
}

template<uint8_t t_width>
void lcp_bitcompressed<t_width>::swap(lcp_bitcompressed& lcp_c)
{
    m_lcp.swap(lcp_c.m_lcp);
}

template<uint8_t t_width>
inline typename lcp_bitcompressed<t_width>::value_type lcp_bitcompressed<t_width>::operator[](size_type i)const
{
    return m_lcp[i];
}


template<uint8_t t_width>
typename lcp_bitcompressed<t_width>::size_type lcp_bitcompressed<t_width>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_lcp.serialize(out, child, "lcp");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<uint8_t t_width>
void lcp_bitcompressed<t_width>::load(std::istream& in)
{
    m_lcp.load(in);
}


template<uint8_t t_width>
lcp_bitcompressed<t_width>& lcp_bitcompressed<t_width>::operator=(const lcp_bitcompressed& lcp_c)
{
    if (this != &lcp_c) {
        copy(lcp_c);
    }
    return *this;
}


} // end namespace sdsl

#endif
