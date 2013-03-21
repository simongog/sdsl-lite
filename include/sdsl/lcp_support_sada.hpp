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
/*! \file lcp_support_sada.hpp
    \brief lcp_support_sada.hpp contains an implementation of a compressed lcp array.
    \author Simon Gog
*/
/* Changes:
    - Removed unnecessary rank support
 */
#ifndef INCLUDED_SDSL_LCP_SUPPORT_SADA
#define INCLUDED_SDSL_LCP_SUPPORT_SADA

#include "lcp.hpp"
#include "int_vector.hpp"
#include "algorithms.hpp"
#include "iterators.hpp"
#include "csa_sada.hpp"  // for standard template initialization of lcp_support_sada 
#include "select_support.hpp" // for standard template initialization of lcp_support_sada
#include <cassert>

namespace sdsl
{


//! A class to represent the LCP array in compressed form.
/*!
 * \tparam t_csa    Type of the Underlying CSA.
 * \tparam t_bitvec Type of the bitvector used to store the unary
 *                  representation of the deltas of the permuted LCP array.
 * \tparam t_select Type of the select structure use to select on the
 *                  bitvector of the unary representation of the PLCP array.
 *
 *    \par Space complexity
 *             \f$ 2n+o(n) \f$ bits, where 2n is the maximal size of the bitvector for the differences of the PLCP array
 *                and o(n) for the select support data structure.
 * \par Reference
 *   Kunihiko Sadakane:
 *   Succinct representations of lcp information and improvements in the compressed suffix arrays.
 *   SODA 2002: 225-232
 */
template<class t_csa = csa_sada<>, class t_bitvec = bit_vector, class t_select = typename t_bitvec::select_1_type>
class _lcp_support_sada
{
    public:
        typedef typename t_csa::value_type                      value_type;    // STL Container requirement
        typedef random_access_const_iterator<_lcp_support_sada> const_iterator;// STL Container requirement
        typedef const_iterator                                  iterator;        // STL Container requirement
        typedef const value_type                                const_reference;
        typedef const_reference                                 reference;
        typedef const_reference*                                pointer;
        typedef const pointer                                   const_pointer;
        typedef int_vector<>::size_type                         size_type;        // STL Container requirement
        typedef ptrdiff_t                                       difference_type; // STL Container requirement
        typedef t_bitvec                                        bit_vector_type;
        typedef t_csa                                           csa_type;
        typedef t_select                                        select_type;


        typedef lcp_permuted_tag                                lcp_category;

        enum { fast_access = 0,
               text_order  = 1,
               sa_order    = 0
             };

        template<class Cst>  // template inner class which is used in CSTs to parametrize lcp classes
        class type           // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef _lcp_support_sada lcp_type;
        };

    private:
        const csa_type* m_csa;
        bit_vector_type m_data;
        select_type     m_select_support;

        void copy(const _lcp_support_sada& lcp_c) {
            m_csa            = lcp_c.m_csa;
            m_data           = lcp_c.m_data;
            m_select_support = lcp_c.m_select_support;
            m_select_support.set_vector(&m_data);
        }
    public:
        const t_csa*& csa;
        //! Default Constructor
        _lcp_support_sada(): csa(m_csa) {}
        //! Default Destructor
        ~_lcp_support_sada() {}
        //! Copy constructor
        _lcp_support_sada(const _lcp_support_sada& lcp_c):csa(m_csa) {
            copy(lcp_c);
        }

        //! Constructor
        _lcp_support_sada(cache_config& config, const t_csa* f_csa);

        void set_csa(const t_csa* f_csa) {
            m_csa = f_csa;
        }

        //! Number of elements in the instance.
        size_type size()const {
            return m_csa->size();
        }

        //! Returns the largest size that _lcp_support_sada can ever have.
        static size_type max_size() {
            return t_csa::max_size();
        }

        //! Returns if the data structure is empty.
        bool empty()const {
            return m_csa->empty();
        }

        //! Swap method for _lcp_support_sada
        void swap(_lcp_support_sada& lcp_c);

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
         * Time complexity: O(suffix array access)
         */
        inline value_type operator[](size_type i)const;

        //! Assignment Operator.
        _lcp_support_sada& operator=(const _lcp_support_sada& lcp_c);

        //! Serialize to a stream.
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const;

        //! Load from a stream.
        void load(std::istream& in, const t_csa* csa);
};

// == template functions ==

template<class t_csa, class t_bitvec, class t_select>
_lcp_support_sada<t_csa, t_bitvec, t_select>::_lcp_support_sada(cache_config& config, const t_csa* f_csa):csa(m_csa)
{
    typedef typename t_csa::size_type size_type;
    set_csa(f_csa);
    int_vector<> lcp;
    util::load_from_file(lcp, util::cache_file_name(constants::KEY_LCP, config));
    if (!util::cache_file_exists(constants::KEY_ISA, config)) {
        construct_isa(config);
    }
    int_vector_file_buffer<> isa_buf(util::cache_file_name(constants::KEY_ISA, config));
    size_type n = lcp.size();
    bit_vector data = bit_vector(2*n, 0);
    size_type data_cnt=0;
    for (size_type i=0, r_sum=0, r = isa_buf.load_next_block(), l=0, old_l=1; r_sum < n;) {
        for (; i < r_sum+r; ++i) {
            l = lcp[isa_buf[i-r_sum]];
            data_cnt += l + 1 - old_l;
            data[data_cnt++] = 1;
            old_l = l;
        }
        r_sum += r;
        r      = isa_buf.load_next_block();
    }
    data.resize(data_cnt);
    util::assign(m_data, data);
    util::init_support(m_select_support, &m_data);
}

template<class t_csa, class t_bitvec, class t_select>
void _lcp_support_sada<t_csa, t_bitvec, t_select>::swap(_lcp_support_sada& lcp_c)
{
    m_data.swap(lcp_c.m_data);
    util::swap_support(m_select_support, lcp_c.m_select_support, &m_data, &(lcp_c.m_data));
}

template<class t_csa, class t_bitvec, class t_select>
inline typename _lcp_support_sada<t_csa, t_bitvec, t_select>::value_type _lcp_support_sada<t_csa, t_bitvec, t_select>::operator[](size_type i)const
{
    size_type j = (*m_csa)[i];
    size_type s = m_select_support.select(j+1);
    return s-(j<<1);
}


template<class t_csa, class t_bitvec, class t_select>
typename _lcp_support_sada<t_csa, t_bitvec, t_select>::size_type _lcp_support_sada<t_csa, t_bitvec, t_select>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_data.serialize(out, child, "data");
    written_bytes += m_select_support.serialize(out, child, "select_support");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<class t_csa, class t_bitvec, class t_select>
void _lcp_support_sada<t_csa, t_bitvec, t_select>::load(std::istream& in, const t_csa* csa)
{
    m_csa = csa;
    m_data.load(in);
    m_select_support.load(in, &m_data);
}


template<class t_csa, class t_bitvec, class t_select>
_lcp_support_sada<t_csa, t_bitvec, t_select>& _lcp_support_sada<t_csa, t_bitvec, t_select>::operator=(const _lcp_support_sada& lcp_c)
{
    if (this != &lcp_c) {
        copy(lcp_c);
    }
    return *this;
}

//! Helper class which provides _lcp_support_sada the context of a CSA.
template<class t_bitvec = bit_vector, class t_select = typename t_bitvec::select_1_type>
class lcp_support_sada
{
    public:
        template<class Cst>  // template inner class which is used in CSTs to parametrize lcp classes
        class type           // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef _lcp_support_sada<typename Cst::csa_type, t_bitvec, t_select> lcp_type;
        };
};

} // end namespace sdsl

#endif
