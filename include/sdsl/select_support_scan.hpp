/* sdsl - succinct data structures library
    Copyright (C) 2012 Simon Gog

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
/*! \file select_support_scan.hpp
    \brief select_support_scan.hpp contains classes that support a sdsl::bit_vector with linear time select.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_SELECT_SUPPORT_SCAN
#define INCLUDED_SDSL_SELECT_SUPPORT_SCAN

#include "int_vector.hpp"
#include "select_support.hpp"

//#define SDSL_DEBUG_SELECT_SUPPORT_JMC

#ifdef SDSL_DEBUG_SELECT_SUPPORT_JMC
#include "testutils.hpp"
#endif

//! Namespace for the succinct data structure library.
namespace sdsl
{


//! A class supporting linear time select queries.
/*! \par Space complexity
 *       Constant.
 * @ingroup select_support_group
 */
template<uint8_t b=1, uint8_t pattern_len=1>
class select_support_scan : public select_support
{
    public:
        typedef bit_vector bit_vector_type;
    public:
        explicit select_support_scan(const bit_vector* v=NULL){set_vector(v);}
        select_support_scan(const select_support_scan<b,pattern_len>& ss){set_vector(ss.m_v);}
        inline const size_type select(size_type i) const;
        //! Alias for select(i).
        inline const size_type operator()(size_type i)const{return select(i);}
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const{return 0;}
        void load(std::istream& in, const bit_vector* v=NULL){}
        void set_vector(const bit_vector* v=NULL){m_v = v;}
        select_support_scan<b, pattern_len>& operator=(const select_support_scan& ss){set_vector(ss.m_v); return *this;}
        void swap(select_support_scan<b, pattern_len>& ss){}
};

template<uint8_t b, uint8_t pattern_len>
inline const typename select_support_scan<b,pattern_len>::size_type select_support_scan<b,pattern_len>::select(size_type i)const {
	const uint64_t* data = m_v->data();
	size_type word_pos = 0;
	size_type word_off = 0;
	uint64_t carry = select_support_trait<b,pattern_len>::init_carry(data, word_pos);
	size_type args = select_support_trait<b,pattern_len>::args_in_the_first_word(*data, word_off, carry);
	if (args >= i) {
		return (word_pos<<6)+select_support_trait<b,pattern_len>::ith_arg_pos_in_the_first_word(*data, i, word_off, carry);
	}
	word_pos+=1;
	size_type sum_args = args;
	carry = select_support_trait<b,pattern_len>::get_carry(*data);
	uint64_t old_carry = carry;
	args = select_support_trait<b,pattern_len>::args_in_the_word(*(++data), carry);
	while (sum_args + args < i) {
		sum_args += args;
		assert(data+1 < m_v->data() + (m_v->capacity()>>6));
		old_carry = carry;
		args = select_support_trait<b,pattern_len>::args_in_the_word(*(++data), carry);
		word_pos+=1;
	}
	return (word_pos<<6) + select_support_trait<b,pattern_len>::ith_arg_pos_in_the_word(*data, i-sum_args, old_carry);
}

} // end namespace
#endif
