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
    \brief lcp_wt.hpp contains an implementation of a (compressed) LCP array based on a WT.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_LCP_WT
#define INCLUDED_SDSL_LCP_WT

#include "lcp.hpp"
#include "wt_huff.hpp"
#include "int_vector.hpp"
#include "algorithms.hpp"
#include "iterators.hpp"
#include "select_support_bs.hpp" // dummy select support for wavelet tree, as we don't use select in this application
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
/*! We use \f$H_0\f$ bit for each lcp values < 255 and \f$ log n \f$ bits for each lcp value which is greater than 254.
 *  \tparam t_width   Width of int_vector storing the large LCP values.
 *  \par Time complexity
 *        - \f$\Order{1}\f$ if the value is less than 255 and
 *        - \f$\Order{\log n}\f$ (\f$n=size()\f$) otherwise.
 */
template<uint8_t t_width=0>
class lcp_wt
{
    public:
        typedef typename int_vector<t_width>::value_type value_type;    // STL Container requirement
        typedef random_access_const_iterator<lcp_wt>     const_iterator;// STL Container requirement
        typedef const_iterator                           iterator;        // STL Container requirement
        typedef const value_type                         const_reference;
        typedef const_reference                          reference;
        typedef const_reference*                         pointer;
        typedef const pointer                            const_pointer;
        typedef int_vector<>::size_type                  size_type;        // STL Container requirement
        typedef ptrdiff_t                                difference_type; // STL Container requirement
        typedef select_support_bs< rank_support_v<> >    tDummySS;
        typedef wt_huff<bit_vector, rank_support_v<>,
                tDummySS, tDummySS>                      small_lcp_type;

        typedef lcp_plain_tag                            lcp_category;

        enum { fast_access = 0,
               text_order  = 0,
               sa_order      = 1
             }; // as the lcp_wt is not fast for texts with long repetition

        template<class Cst>  // template inner class which is used in CSTs to parametrize lcp classes
        class type           // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef lcp_wt lcp_type;
        };

    private:
        small_lcp_type         m_small_lcp; // vector for lcp values < 255
        int_vector<t_width>   m_big_lcp;     // vector for lcp values > 254

        typedef std::pair<size_type, size_type> tPII;
        typedef std::vector<tPII> tVPII;

        void copy(const lcp_wt& lcp_c) {
            m_small_lcp     = lcp_c.m_small_lcp;
            m_big_lcp        = lcp_c.m_big_lcp;
        }

    public:
        //! Default Constructor
        lcp_wt() {}
        //! Copy constructor
        lcp_wt(const lcp_wt& lcp_c) {
            copy(lcp_c);
        }

        //! Constructor
        lcp_wt(cache_config& config, std::string other_key="");

        //! Number of elements in the instance.
        size_type size()const {
            return m_small_lcp.size();
        }

        //! Returns the largest size that lcp_wt can ever have.
        static size_type max_size() {
            return int_vector<8>::max_size();
        }

        //! Returns if the data structure is empty.
        bool empty()const {
            return 0==m_small_lcp.size();
        }

        //! Swap method for lcp_wt
        void swap(lcp_wt& lcp_c);

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
         */
        inline value_type operator[](size_type i)const;

        //! Assignment Operator.
        lcp_wt& operator=(const lcp_wt& lcp_c);

        //! Serialize to a stream.
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const;

        //! Load from a stream.
        void load(std::istream& in);
};

// == template functions ==


template<uint8_t t_width>
lcp_wt<t_width>::lcp_wt(cache_config& config, std::string other_key)
{
    std::string tmp_file = util::tmp_file(config, "_lcp_sml");
    std::string lcp_key  = constants::KEY_LCP;
    if ("" != other_key) {
        lcp_key = other_key;
    }
    int_vector_file_buffer<> lcp_buf(util::cache_file_name(lcp_key, config));
    typename int_vector<>::size_type l=0, max_l=0, big_sum=0, n = lcp_buf.int_vector_size;
    {
        int_vector<8> small_lcp = int_vector<8>(n);
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
        util::store_to_file(small_lcp, tmp_file);
    }
    int_vector_file_buffer<8> lcp_sml_buf(tmp_file);
    {
        small_lcp_type tmp_small_lcp(lcp_sml_buf, lcp_sml_buf.int_vector_size);
        m_small_lcp.swap(tmp_small_lcp);
    }
    sdsl::remove(tmp_file);
    m_big_lcp         = int_vector<>(big_sum, 0, bit_magic::l1BP(max_l)+1);
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
}

template<uint8_t t_width>
void lcp_wt<t_width>::swap(lcp_wt& lcp_c)
{
    m_small_lcp.swap(lcp_c.m_small_lcp);
    m_big_lcp.swap(lcp_c.m_big_lcp);
}

template<uint8_t t_width>
inline typename lcp_wt<t_width>::value_type lcp_wt<t_width>::operator[](size_type i)const
{
    if (m_small_lcp[i]!=255) {
        return m_small_lcp[i];
    } else {
        return m_big_lcp[ m_small_lcp.rank(i, 255) ];
    }
}


template<uint8_t t_width>
typename lcp_wt<t_width>::size_type lcp_wt<t_width>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_small_lcp.serialize(out, child,  "small_lcp");
    written_bytes += m_big_lcp.serialize(out, child, "large_lcp");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<uint8_t t_width>
void lcp_wt<t_width>::load(std::istream& in)
{
    m_small_lcp.load(in);
    m_big_lcp.load(in);
}


template<uint8_t t_width>
lcp_wt<t_width>& lcp_wt<t_width>::operator=(const lcp_wt& lcp_c)
{
    if (this != &lcp_c) {
        copy(lcp_c);
    }
    return *this;
}
} // end namespace sdsl

#endif
