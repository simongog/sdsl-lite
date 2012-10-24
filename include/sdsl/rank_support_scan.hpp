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
/*! \file rank_support_scan.hpp
    \brief rank_support_scan.hpp contains rank_support_scan that support a sdsl::bit_vector with linear time rank information.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_RANK_SUPPORT_SCAN
#define INCLUDED_SDSL_RANK_SUPPORT_SCAN

#include "rank_support.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A class supporting rank queries in linear time. 
/*! \par Space complexity
 *       Constant.
 * @ingroup rank_support_group
 */
template<uint8_t b=1, uint8_t pattern_len=1>
class rank_support_scan : public rank_support
{
    public:
        typedef bit_vector bit_vector_type;
    public:
        explicit rank_support_scan(const bit_vector* v = NULL){set_vector(v);}
        rank_support_scan(const rank_support_scan& rs){}
        const size_type rank(size_type idx) const;
        const size_type operator()(size_type idx)const;
        const size_type size()const;
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const{return 0;}
        void load(std::istream& in, const int_vector<1>* v=NULL){ set_vector(v); }
        void set_vector(const bit_vector* v=NULL){m_v=v;}

        //! Assign Operator
        rank_support_scan& operator=(const rank_support_scan& rs);

        //! swap Operator
        void swap(rank_support_scan& rs){}

        //! Equality Operator
        bool operator==(const rank_support_scan& rs)const;
        //! Inequality Operator
        bool operator!=(const rank_support_scan& rs)const;
};

template<uint8_t b, uint8_t pattern_len>
inline const typename rank_support_scan<b, pattern_len>::size_type rank_support_scan<b, pattern_len>::size()const
{
    return m_v->size();
}

template<uint8_t b, uint8_t pattern_len>
inline const typename rank_support_scan<b, pattern_len>::size_type rank_support_scan<b, pattern_len>::rank(size_type idx)const {
    const uint64_t* p 	= m_v->data();
	size_type 	i		= 0;
	size_type   result  = 0;
	while ( i+64 <= idx ){
		result += rank_support_trait<b, pattern_len>::full_word_rank(p, i);
		i += 64;
	}
    return  result+rank_support_trait<b, pattern_len>::word_rank(p, idx);
}


template<uint8_t b, uint8_t pattern_len>
inline const typename rank_support_scan<b, pattern_len>::size_type rank_support_scan<b, pattern_len>::operator()(size_type idx)const {
    return rank(idx);
}

template<uint8_t b, uint8_t pattern_len>
inline rank_support_scan<b, pattern_len>& rank_support_scan<b, pattern_len>::operator=(const rank_support_scan& rs) {
    if (this != &rs) {
        set_vector(rs.m_v);
    }
    return *this;
}

// TODO: == operator remove pointer comparison
template<uint8_t b, uint8_t pattern_len>
inline bool rank_support_scan<b, pattern_len>::operator==(const rank_support_scan& rs)const {
    if (this == &rs)
        return true;
    return *(rs.m_v) == *m_v;
}

template<uint8_t b, uint8_t pattern_len>
inline bool rank_support_scan<b, pattern_len>::operator!=(const rank_support_scan& rs)const {
    return !(*this == rs);
}

}// end namespace sds

#endif // end file 
