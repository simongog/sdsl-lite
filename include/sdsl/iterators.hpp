/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

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
/*! \file iterators.hpp
    \brief iterators.hpp contains an generic iterator for random access containers.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_ITERATORS
#define INCLUDED_SDSL_ITERATORS

#include <iterator>

namespace sdsl
{

//! Generic iterator for a random access container
/*! \tparam t_rac Type of random access container.
 */
template<class t_rac>
class random_access_const_iterator: public std::iterator<std::random_access_iterator_tag, typename t_rac::value_type, typename t_rac::difference_type>
{
    public:
        typedef const typename t_rac::value_type  const_reference;
        typedef typename t_rac::size_type size_type;
        typedef random_access_const_iterator<t_rac> iterator;
        typedef typename t_rac::difference_type difference_type;

    private:
        const t_rac* m_rac;// pointer to the random access container
        typename t_rac::size_type m_idx;

        template<class t_RAC>
        friend typename random_access_const_iterator<t_RAC>::difference_type operator-(const random_access_const_iterator<t_RAC>& x,
                const random_access_const_iterator<t_RAC>& y);


    public:
        //! Constructor
        random_access_const_iterator(const t_rac* rac, size_type idx = 0) : m_rac(rac), m_idx(idx) { }

        //! Dereference operator for the Iterator.
        const_reference operator*()const {
            return (*m_rac)[m_idx];
        }

        //! Prefix increment of the Iterator.
        iterator& operator++() {
            ++m_idx;
            return *this;
        }

        //! Postfix increment of the Iterator.
        iterator operator++(int) {
            random_access_const_iterator it = *this;
            ++(*this);
            return it;
        }

        //! Prefix decrement of the Iterator.
        iterator& operator--() {
            --m_idx;
            return *this;
        }

        //! Postfix decrement of the Iterator.
        iterator operator--(int) {
            random_access_const_iterator it = *this;
            ++(*this);
            return it;
        }

        iterator& operator+=(difference_type i) {
            if (i<0)
                return *this -= (-i);
            m_idx += i;
            return *this;
        }

        iterator& operator-=(difference_type i) {
            if (i<0)
                return *this += (-i);
            m_idx -= i;
            return *this;
        }

        iterator operator+(difference_type i) const {
            iterator it = *this;
            return it += i;
        }

        iterator operator-(difference_type i) const {
            iterator it = *this;
            return it -= i;
        }

        const_reference operator[](difference_type i) const {
            return *(*this + i);
        }

        bool operator==(const iterator& it)const {
            return it.m_rac == m_rac && it.m_idx == m_idx;
        }

        bool operator!=(const iterator& it)const {
            return !(*this==it);
        }

        bool operator<(const iterator& it)const {
            return m_idx < it.m_idx;
        }

        bool operator>(const iterator& it)const {
            return m_idx > it.m_idx;
        }

        bool operator>=(const iterator& it)const {
            return !(*this < it);
        }

        bool operator<=(const iterator& it)const {
            return !(*this > it);
        }

};

template<class t_rac>
inline typename random_access_const_iterator<t_rac>::difference_type operator-(const random_access_const_iterator<t_rac>& x, const random_access_const_iterator<t_rac>& y)
{
    return (typename random_access_const_iterator<t_rac>::difference_type)x.m_idx
           - (typename random_access_const_iterator<t_rac>::difference_type)y.m_idx;
}

template<class t_rac>
inline random_access_const_iterator<t_rac> operator+(typename random_access_const_iterator<t_rac>::difference_type n, const random_access_const_iterator<t_rac>& it)
{
    return it+n;
}


//! Generic uint64_t iterator for bitvectors.
/*! \tparam t_nv Bitvector type.
 */
template<class t_bv>
class uint64_bv_const_iterator: public std::iterator<std::random_access_iterator_tag, uint64_t, typename t_bv::difference_type>
{
        // TODO: static assert for index_category bv_tag
    public:
        typedef const uint64_t const_reference;
        typedef typename t_bv::size_type size_type;
        typedef uint64_bv_const_iterator<t_bv> iterator;
        typedef typename t_bv::difference_type difference_type;

    private:
        const t_bv* m_bv;// pointer to the random access container
        size_type m_idx;

        template<class t_RAC>
        friend typename uint64_bv_const_iterator<t_RAC>::difference_type operator-(const uint64_bv_const_iterator<t_RAC>& x,
                const uint64_bv_const_iterator<t_RAC>& y);


    public:
        //! Constructor
        uint64_bv_const_iterator(const t_bv* bv, size_type idx = 0) : m_bv(bv),
            m_idx(idx) { }

        //! Dereference operator for the Iterator.
        const_reference operator*()const {
            size_type idx = m_idx << 6;
            if (idx+64 < m_bv->size()) {
                return m_bv->get_int(idx, 64);
            } else if (idx < m_bv->size()) {
                return m_bv->get_int(idx, m_bv->size()-idx);
            } else {
                throw std::out_of_range("m_bv->size()="+std::to_string(m_bv->size())+
                                        " idx="+std::to_string(idx));
                return 0;
            }
        }

        //! Prefix increment of the Iterator.
        iterator& operator++() {
            ++m_idx;
            return *this;
        }

        //! Postfix increment of the Iterator.
        iterator operator++(int) {
            uint64_bv_const_iterator it = *this;
            ++(*this);
            return it;
        }

        //! Prefix decrement of the Iterator.
        iterator& operator--() {
            --m_idx;
            return *this;
        }

        //! Postfix decrement of the Iterator.
        iterator operator--(int) {
            uint64_bv_const_iterator it = *this;
            ++(*this);
            return it;
        }

        iterator& operator+=(difference_type i) {
            if (i<0)
                return *this -= (-i);
            m_idx += i;
            return *this;
        }

        iterator& operator-=(difference_type i) {
            if (i<0)
                return *this += (-i);
            m_idx -= i;
            return *this;
        }

        iterator operator+(difference_type i) const {
            iterator it = *this;
            return it += i;
        }

        iterator operator-(difference_type i) const {
            iterator it = *this;
            return it -= i;
        }

        const_reference operator[](difference_type i) const {
            return *(*this + i);
        }

        bool operator==(const iterator& it)const {
            return it.m_bv == m_bv && it.m_idx == m_idx;
        }

        bool operator!=(const iterator& it)const {
            return !(*this==it);
        }

        bool operator<(const iterator& it)const {
            return m_idx < it.m_idx;
        }

        bool operator>(const iterator& it)const {
            return m_idx > it.m_idx;
        }

        bool operator>=(const iterator& it)const {
            return !(*this < it);
        }

        bool operator<=(const iterator& it)const {
            return !(*this > it);
        }

};

template<class t_bv>
inline typename uint64_bv_const_iterator<t_bv>::difference_type operator-(const uint64_bv_const_iterator<t_bv>& x, const uint64_bv_const_iterator<t_bv>& y)
{
    return (typename uint64_bv_const_iterator<t_bv>::difference_type)x.m_idx
           - (typename uint64_bv_const_iterator<t_bv>::difference_type)y.m_idx;
}

template<class t_bv>
inline uint64_bv_const_iterator<t_bv> operator+(typename uint64_bv_const_iterator<t_bv>::difference_type n, const uint64_bv_const_iterator<t_bv>& it)
{
    return it+n;
}


} // end namespace sdsl
#endif
