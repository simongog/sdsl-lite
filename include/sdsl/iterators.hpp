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

        //! Default Constructor
        random_access_const_iterator() : m_rac(nullptr), m_idx(0) { }

        //! Constructor
        random_access_const_iterator(const t_rac* rac, size_type idx = 0) : m_rac(rac), m_idx(idx) { }

        //! Dereference operator for the Iterator.
        const_reference operator*()const
        {
            return (*m_rac)[m_idx];
        }

        //! Prefix increment of the Iterator.
        iterator& operator++()
        {
            ++m_idx;
            return *this;
        }

        //! Postfix increment of the Iterator.
        iterator operator++(int)
        {
            random_access_const_iterator it = *this;
            ++(*this);
            return it;
        }

        //! Prefix decrement of the Iterator.
        iterator& operator--()
        {
            --m_idx;
            return *this;
        }

        //! Postfix decrement of the Iterator.
        iterator operator--(int)
        {
            random_access_const_iterator it = *this;
            --(*this);
            return it;
        }

        iterator& operator+=(difference_type i)
        {
            if (i<0)
                return *this -= (-i);
            m_idx += i;
            return *this;
        }

        iterator& operator-=(difference_type i)
        {
            if (i<0)
                return *this += (-i);
            m_idx -= i;
            return *this;
        }

        iterator operator+(difference_type i) const
        {
            iterator it = *this;
            return it += i;
        }

        iterator operator-(difference_type i) const
        {
            iterator it = *this;
            return it -= i;
        }

        const_reference operator[](difference_type i) const
        {
            return *(*this + i);
        }

        bool operator==(const iterator& it)const
        {
            return it.m_rac == m_rac && it.m_idx == m_idx;
        }

        bool operator!=(const iterator& it)const
        {
            return !(*this==it);
        }

        bool operator<(const iterator& it)const
        {
            return m_idx < it.m_idx;
        }

        bool operator>(const iterator& it)const
        {
            return m_idx > it.m_idx;
        }

        bool operator>=(const iterator& it)const
        {
            return !(*this < it);
        }

        bool operator<=(const iterator& it)const
        {
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


template<typename t_F>
struct random_access_container {
    typedef int_vector<>::size_type                               size_type;
    typedef int_vector<>::difference_type                         difference_type;
    typedef typename std::result_of<t_F(size_type)>::type         value_type;
    typedef random_access_const_iterator<random_access_container> iterator_type;

    t_F f;
    size_type m_size;

    random_access_container() {};
    random_access_container(t_F ff, size_type size) : f(ff), m_size(size) { }

    value_type operator[](size_type i) const { return f(i); }

    size_type size() const { return m_size; }

    iterator_type begin() const
    {
        return iterator_type(this, 0);
    }

    iterator_type end() const
    {
        return iterator_type(this, size());
    }
};

} // end namespace sdsl
#endif
