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
    \brief lcp_byte.hpp contains a (compressed) lcp array.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_LCP_BYTE
#define INCLUDED_SDSL_LCP_BYTE

#include "lcp.hpp"
#include "int_vector.hpp"
#include "iterators.hpp"
#include <iostream>
#include <algorithm> // for lower_bound
#include <cassert>
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
 *        - \f$\Order{1}\f$ if the value is less than 255 and
 *        - \f$\Order{\log n}\f$ (\f$n=size()\f$) otherwise.
 * \par Reference
 *   Mohamed Ibrahim Abouelhoda, Stefan Kurtz, Enno Ohlebusch:
 *   Replacing suffix trees with enhanced suffix arrays.
 *   J. Discrete Algorithms 2(1): 53-86 (2004)
 */
template<uint8_t t_width=0>
class lcp_byte
{
    public:

        typedef typename int_vector<t_width>::value_type value_type;
        typedef random_access_const_iterator<lcp_byte>   const_iterator;
        typedef const_iterator                           iterator;
        typedef const value_type                         const_reference;
        typedef const_reference                          reference;
        typedef const_reference*                         pointer;
        typedef const pointer                            const_pointer;
        typedef int_vector<>::size_type                  size_type;
        typedef ptrdiff_t                                difference_type;

        typedef lcp_plain_tag                            lcp_category;
        typedef lcp_tag                                  index_tag;

        enum { fast_access = 0,
               text_order  = 0,
               sa_order    = 1
             }; // as the lcp_byte is not fast for texts with long repetition

        template<class Cst>
        using type = lcp_byte;

    private:

        int_vector<8>       m_small_lcp;   // vector for LCP values < 255
        int_vector<t_width> m_big_lcp;     // vector for LCP values > 254
        int_vector<t_width> m_big_lcp_idx; // index of LCP entries in the LCP array

        typedef std::pair<size_type, size_type> tPII;
        typedef std::vector<tPII> tVPII;

    public:

        //! Default Constructor
        lcp_byte() = default;
        lcp_byte(const lcp_byte&) = default;
        lcp_byte(lcp_byte&&) = default;
        lcp_byte& operator=(const lcp_byte&) = default;
        lcp_byte& operator=(lcp_byte&&) = default;

        //! Constructor
        lcp_byte(cache_config& config)
        {
            std::string lcp_file = cache_file_name(conf::KEY_LCP, config);
            int_vector_buffer<> lcp_buf(lcp_file);
            m_small_lcp = int_vector<8>(lcp_buf.size());
            size_type l=0, max_l=0, max_big_idx=0, big_sum=0;

            for (size_type i=0; i < m_small_lcp.size(); ++i) {
                if ((l=lcp_buf[i]) < 255) {
                    m_small_lcp[i] = l;
                } else {
                    m_small_lcp[i] = 255;
                    if (l > max_l) max_l = l;
                    max_big_idx = i;
                    ++big_sum;
                }
            }
            m_big_lcp     = int_vector<>(big_sum, 0, bits::hi(max_l)+1);
            m_big_lcp_idx = int_vector<>(big_sum, 0, bits::hi(max_big_idx)+1);

            for (size_type i=0,ii=0; i<m_small_lcp.size(); ++i) {
                if ((l=lcp_buf[i]) >= 255) {
                    m_big_lcp[ii] = l;
                    m_big_lcp_idx[ii] = i;
                    ++ii;
                }
            }
        }

        //! Number of elements in the instance.
        size_type size()const
        {
            return m_small_lcp.size();
        }

        //! Returns the largest size that lcp_byte can ever have.
        static size_type max_size()
        {
            return int_vector<8>::max_size();
        }

        //! Returns if the data structure is empty.
        bool empty()const
        {
            return m_small_lcp.empty();
        }

        //! Swap method for lcp_byte
        void swap(lcp_byte& lcp_c)
        {
            m_small_lcp.swap(lcp_c.m_small_lcp);
            m_big_lcp.swap(lcp_c.m_big_lcp);
            m_big_lcp_idx.swap(lcp_c.m_big_lcp_idx);
        }


        //! Returns a const_iterator to the first element.
        const_iterator begin()const
        {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        const_iterator end()const
        {
            return const_iterator(this, size());
        }

        //! []-operator
        /*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
         * Time complexity: O(1) for small and O(log n) for large values
         */
        inline value_type operator[](size_type i)const
        {
            if (m_small_lcp[i]!=255) {
                return m_small_lcp[i];
            } else {
                size_type idx = lower_bound(m_big_lcp_idx.begin(),
                                            m_big_lcp_idx.end(),i)
                                - m_big_lcp_idx.begin();
                return m_big_lcp[idx];
            }
        }

        //! Serialize to a stream.
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr,
                            std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(v, name,
                                         util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_small_lcp.serialize(out, child, "small_lcp");
            written_bytes += m_big_lcp.serialize(out, child, "large_lcp");
            written_bytes += m_big_lcp_idx.serialize(out, child, "large_lcp_idx");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load from a stream.
        void load(std::istream& in)
        {
            m_small_lcp.load(in);
            m_big_lcp.load(in);
            m_big_lcp_idx.load(in);
        }
};

} // end namespace sdsl
#endif
