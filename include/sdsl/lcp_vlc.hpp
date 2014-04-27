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
    \brief lcp_vlc.hpp contains an implementation of a (compressed) LCP array.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_LCP_VLC
#define INCLUDED_SDSL_LCP_VLC

#include "lcp.hpp"
#include "vlc_vector.hpp"
#include "int_vector.hpp"
#include "iterators.hpp"
#include <iostream>
#include <cassert>
#include <vector>

namespace sdsl
{

// A class for a compressed LCP array based on variable-length coding.
/*
 * \tparam t_vlc_vec Type of the underlying variable-length coder.
 */
template<class t_vlc_vec = vlc_vector<> >
class lcp_vlc
{
    public:

        typedef typename t_vlc_vec::value_type        value_type;
        typedef random_access_const_iterator<lcp_vlc> const_iterator;
        typedef const_iterator                        iterator;
        typedef const value_type                      const_reference;
        typedef const_reference                       reference;
        typedef const_reference*                      pointer;
        typedef const pointer                         const_pointer;
        typedef typename t_vlc_vec::size_type         size_type;
        typedef typename t_vlc_vec::difference_type   difference_type;
        typedef t_vlc_vec                             vlc_vec_type;

        typedef lcp_plain_tag                         lcp_category;

        enum { fast_access = 0,
               text_order  = 0,
               sa_order    = 1
             };

        template<class Cst>
        using type = lcp_vlc;

    private:

        vlc_vec_type        m_vec;

    public:

        //! Default Constructor
        lcp_vlc() = default;

        //! Copy / Move constructor
        lcp_vlc(const lcp_vlc&) = default;
        lcp_vlc(lcp_vlc&&) = default;
        lcp_vlc& operator=(const lcp_vlc&) = default;
        lcp_vlc& operator=(lcp_vlc&&) = default;

        //! Construct
        lcp_vlc(cache_config& config, std::string other_key="") {
            std::string lcp_key  = conf::KEY_LCP;
            if ("" != other_key) {
                lcp_key = other_key;
            }
            int_vector_buffer<> lcp_buf(cache_file_name(lcp_key, config));
            vlc_vec_type tmp_vec(lcp_buf);
            m_vec.swap(tmp_vec);
        }

        //! Number of elements in the instance.
        size_type size()const {
            return m_vec.size();
        }

        //! Returns the largest size that lcp_vlc can ever have.
        static size_type max_size() {
            return vlc_vec_type::max_size();
        }

        //! Returns if the data strucutre is empty.
        bool empty()const {
            return m_vec.empty();
        }

        //! Swap method for lcp_vlc
        void swap(lcp_vlc& lcp_c) {
            m_vec.swap(lcp_c.m_vec);
        }

        //! Returns a const_iterator to the first element.
        const_iterator begin()const {
            return const_iterator(this, 0);
        }


        //! Returns a const_iterator to the element after the last element.
        const_iterator end()const {
            return const_iterator(this, size());
        }

        //! []-operator
        inline value_type operator[](size_type i)const {
            return m_vec[i];
        }

        //! Serialize to a stream.
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_vec.serialize(out, child, "vec");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load from a stream.
        void load(std::istream& in) {
            m_vec.load(in);
        }
};

} // end namespace sdsl
#endif
