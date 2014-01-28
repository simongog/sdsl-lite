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
    \brief lcp_support_sada.hpp contains a compressed lcp array.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_LCP_SUPPORT_SADA
#define INCLUDED_SDSL_LCP_SUPPORT_SADA

#include "lcp.hpp"
#include "int_vector.hpp"
#include "iterators.hpp"
#include "csa_sada.hpp"  // for default template initialization
#include "select_support.hpp" // for default template initialization
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
 *         \f$ 2n+o(n) \f$ bits, where 2n is the worst case size of the
 *         unary encoding of the deltas of the PLCP array and o(n) for
 *         the select support data structure.
 * \par Reference
 *   Kunihiko Sadakane:
 *   Succinct representations of lcp information and improvements in the
 *   compressed suffix arrays.
 *   SODA 2002: 225-232
 */
template<class t_csa = csa_sada<>,
         class t_bitvec = bit_vector,
         class t_select = typename t_bitvec::select_1_type>
class _lcp_support_sada
{
    public:
        typedef typename t_csa::value_type     value_type;
        typedef random_access_const_iterator<
        _lcp_support_sada> const_iterator;
        typedef const_iterator                 iterator;
        typedef const value_type               const_reference;
        typedef const_reference                reference;
        typedef const_reference*               pointer;
        typedef const pointer                  const_pointer;
        typedef int_vector<>::size_type        size_type;
        typedef ptrdiff_t                      difference_type;
        typedef t_bitvec                       bit_vector_type;
        typedef t_csa                          csa_type;
        typedef t_select                       select_type;


        typedef lcp_permuted_tag               lcp_category;

        enum { fast_access = 0,
               text_order  = 1,
               sa_order    = 0
             };
        // inner class which is used in CSTs to parametrize lcp classes
        // with information about the CST.
        template<class Cst>
        struct type {
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
        const t_csa*& csa = m_csa;
        //! Default Constructor
        _lcp_support_sada() {}

        //! Copy constructor
        _lcp_support_sada(const _lcp_support_sada& lcp_c) {
            copy(lcp_c);
        }

        //! Move constructor
        _lcp_support_sada(_lcp_support_sada&& lcp_c) {
            *this = std::move(lcp_c);
        }

        //! Constructor
        _lcp_support_sada(cache_config& config, const t_csa* f_csa) {
            typedef typename t_csa::size_type size_type;
            set_csa(f_csa);
            if (!cache_file_exists(conf::KEY_ISA, config)) {
                construct_isa(config);
            }
            int_vector<> lcp;
            load_from_file(lcp, cache_file_name(conf::KEY_LCP, config));
            std::string isa_file = cache_file_name(conf::KEY_ISA, config);
            int_vector_buffer<> isa_buf(isa_file);
            size_type n = lcp.size();
            bit_vector data = bit_vector(2*n, 0);
            size_type data_cnt=0;
            for (size_type i=0, l=0, old_l=1; i < n; ++i) {
                l = lcp[isa_buf[i]];
                data_cnt += l + 1 - old_l;
                data[data_cnt++] = 1;
                old_l = l;
            }
            data.resize(data_cnt);
            m_data = bit_vector_type(data);
            util::init_support(m_select_support, &m_data);
        }

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
        void swap(_lcp_support_sada& lcp_c) {
            m_data.swap(lcp_c.m_data);
            util::swap_support(m_select_support, lcp_c.m_select_support,
                               &m_data, &(lcp_c.m_data));
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
        /*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
         * Time complexity: O(suffix array access)
         */
        inline value_type operator[](size_type i)const {
            size_type j = (*m_csa)[i];
            size_type s = m_select_support.select(j+1);
            return s-(j<<1);
        }

        //! Assignment Operator.
        _lcp_support_sada& operator=(const _lcp_support_sada& lcp_c) {
            if (this != &lcp_c) {
                copy(lcp_c);
            }
            return *this;
        }

        //! Assignment Move Operator.
        _lcp_support_sada& operator=(_lcp_support_sada&& lcp_c) {
            if (this != &lcp_c) {
                m_csa            = std::move(lcp_c.m_csa);
                m_data           = std::move(lcp_c.m_data);
                m_select_support = std::move(lcp_c.m_select_support);
                m_select_support.set_vector(&m_data);
            }
            return *this;
        }

        //! Serialize to a stream.
        size_type
        serialize(std::ostream& out, structure_tree_node* v=nullptr,
                  std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name,
                                         util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_data.serialize(out, child, "data");
            written_bytes += m_select_support.serialize(out, child,
                             "select_support");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load from a stream.
        void load(std::istream& in, const t_csa* csa) {
            m_csa = csa;
            m_data.load(in);
            m_select_support.load(in, &m_data);
        }
};

//! Helper class which provides _lcp_support_sada the context of a CSA.
template<class t_bitvec = bit_vector,
         class t_select = typename t_bitvec::select_1_type>
struct lcp_support_sada
{
    template<class t_cst>
    using type = _lcp_support_sada<typename t_cst::csa_type,
          t_bitvec,
          t_select>;
};

} // end namespace sdsl
#endif
