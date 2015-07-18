/* sdsl - succinct data structures library
    Copyright (C) 2010-2013 Simon Gog

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
    \brief lcp_wt.hpp contains a (compressed) LCP array based on a WT.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_LCP_WT
#define INCLUDED_SDSL_LCP_WT

#include "lcp.hpp"
#include "wt_huff.hpp"
#include "int_vector.hpp"
#include "iterators.hpp"
#include "util.hpp"
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
/*! We use \f$H_0\f$ bit for each lcp values < 255 and \f$ log n \f$ bits for
 *  each LCP value which is greater than 254.
 *  \tparam t_width   Width of int_vector storing the large LCP values.
 *  \par Time complexity
 *        - \f$\Order{1}\f$ if the value is less than 255 and
 *        - \f$\Order{\log n}\f$ (\f$n=size()\f$) otherwise.
 */
template<uint8_t t_width=0>
class lcp_wt
{
    public:

        typedef typename int_vector<t_width>::value_type value_type;
        typedef random_access_const_iterator<lcp_wt>     const_iterator;
        typedef const_iterator                           iterator;
        typedef const value_type                         const_reference;
        typedef const_reference                          reference;
        typedef const_reference*                         pointer;
        typedef const pointer                            const_pointer;
        typedef int_vector<>::size_type                  size_type;
        typedef ptrdiff_t                                difference_type;
        typedef wt_huff<bit_vector, rank_support_v<>,
                select_support_scan<1>,
                select_support_scan<0>>                  small_lcp_type;

        typedef lcp_plain_tag                            lcp_category;
        typedef lcp_tag                                  index_category;

        enum { fast_access = 0,
               text_order  = 0,
               sa_order    = 1
             }; // as the lcp_wt is not fast for texts with long repetition

        template<class Cst>
        using type = lcp_wt;

    private:

        small_lcp_type      m_small_lcp; // vector for lcp values < 255
        int_vector<t_width> m_big_lcp;     // vector for lcp values > 254

        typedef std::pair<size_type, size_type> tPII;
        typedef std::vector<tPII> tVPII;

    public:

        //! Default Constructor
        lcp_wt() = default;
        //! Copy / Move constructor
        lcp_wt(const lcp_wt&) = default;
        lcp_wt(lcp_wt&&) = default;
        lcp_wt& operator=(const lcp_wt&) = default;
        lcp_wt& operator=(lcp_wt&&) = default;

        //! Constructor
        lcp_wt(cache_config& config, std::string other_key="")
        {
            std::string temp_file = tmp_file(config, "_lcp_sml");
            std::string lcp_key  = conf::KEY_LCP;
            if ("" != other_key) {
                lcp_key = other_key;
            }
            int_vector_buffer<> lcp_buf(cache_file_name(lcp_key, config));
            size_type l=0, max_l=0, big_sum=0, n = lcp_buf.size();
            {
                int_vector<8> small_lcp = int_vector<8>(n);
                for (size_type i=0; i < n; ++i) {
                    if ((l=lcp_buf[i]) < 255) {
                        small_lcp[i] = l;
                    } else {
                        small_lcp[i] = 255;
                        if (l > max_l) max_l = l;
                        ++big_sum;
                    }
                }
                store_to_file(small_lcp, temp_file);
            }
            {
                int_vector_buffer<8> lcp_sml_buf(temp_file);
                small_lcp_type tmp(lcp_sml_buf, lcp_sml_buf.size());
                m_small_lcp.swap(tmp);
            }
            sdsl::remove(temp_file);
            m_big_lcp = int_vector<>(big_sum, 0, bits::hi(max_l)+1);
            {
                for (size_type i=0, ii=0; i < n; ++i) {
                    if (lcp_buf[i] >= 255) {
                        m_big_lcp[ ii++ ] = lcp_buf[i];
                    }
                }
            }
        }

        //! Number of elements in the instance.
        size_type size()const
        {
            return m_small_lcp.size();
        }

        //! Returns the largest size that lcp_wt can ever have.
        static size_type max_size()
        {
            return int_vector<8>::max_size();
        }

        //! Returns if the data structure is empty.
        bool empty()const
        {
            return 0==m_small_lcp.size();
        }

        //! Swap method for lcp_wt
        void swap(lcp_wt& lcp_c)
        {
            m_small_lcp.swap(lcp_c.m_small_lcp);
            m_big_lcp.swap(lcp_c.m_big_lcp);
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
         */
        inline value_type operator[](size_type i)const
        {
            if (m_small_lcp[i]!=255) {
                return m_small_lcp[i];
            } else {
                return m_big_lcp[ m_small_lcp.rank(i, 255) ];
            }
        }

        //! Serialize to a stream.
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr,
                            std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(
                                             v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_small_lcp.serialize(out, child, "small_lcp");
            written_bytes += m_big_lcp.serialize(out, child, "large_lcp");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load from a stream.
        void load(std::istream& in)
        {
            m_small_lcp.load(in);
            m_big_lcp.load(in);
        }


};

} // end namespace sdsl
#endif
