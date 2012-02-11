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

namespace sdsl{

template<class RandomAccessContainer>
class random_access_const_iterator: public std::iterator<std::random_access_iterator_tag, typename RandomAccessContainer::value_type, typename RandomAccessContainer::difference_type>{
	public:
	typedef const typename RandomAccessContainer::value_type  const_reference;
	typedef typename RandomAccessContainer::size_type size_type;
	typedef random_access_const_iterator<RandomAccessContainer> iterator;
	typedef typename RandomAccessContainer::difference_type difference_type;

	private:
	const RandomAccessContainer *m_rac;// pointer to the random access container
	typename RandomAccessContainer::size_type m_idx;

	template<class RAC>
	friend typename random_access_const_iterator<RAC>::difference_type operator-(const random_access_const_iterator<RAC> &x, 
																				 const random_access_const_iterator<RAC> &y);
	

	public:
	//! Constructor
	random_access_const_iterator(const RandomAccessContainer *rac, size_type idx = 0){
		m_rac = rac;
		m_idx = idx;
	}

	//! Dereference operator for the Iterator.
	const_reference operator*()const{
		return (*m_rac)[m_idx];
	}

	//! Prefix increment of the Iterator.
	iterator& operator++(){
		++m_idx;
		return *this;
	}

	//! Postfix increment of the Iterator.
	iterator operator++(int x){
		random_access_const_iterator it = *this;
		++(*this);
		return it;
	}

	//! Prefix decrement of the Iterator.
	iterator& operator--(){
		--m_idx;
		return *this;
	}

	//! Postfix decrement of the Iterator.
	iterator operator--(int x){
		random_access_const_iterator it = *this;
		++(*this);
		return it;
	}

	iterator& operator+=(difference_type i){
		if(i<0)
			return *this -= (-i);
		m_idx += i;
		return *this;
	}

	iterator& operator-=(difference_type i){
		if(i<0)
			return *this += (-i);
		m_idx -= i;
		return *this;
	}

	iterator operator+(difference_type i) const{
		iterator it = *this;
		return it += i;
	}

	iterator operator-(difference_type i) const{
		iterator it = *this;
		return it -= i;
	}

	const_reference operator[](difference_type i) const{
		return *(*this + i);
	}

	bool operator==(const iterator &it)const{
		return it.m_rac == m_rac && it.m_idx == m_idx;
	}

	bool operator!=(const iterator &it)const{
		return !(*this==it);
	}

	bool operator<(const iterator &it)const{
		return m_idx < it.m_idx;
	}

	bool operator>(const iterator &it)const{
		return m_idx > it.m_idx;
	}

	bool operator>=(const iterator &it)const{
		return !(*this < it);
	}

	bool operator<=(const iterator &it)const{
		return !(*this > it);
	}

};

template<class RandomAccessContainer>
inline typename random_access_const_iterator<RandomAccessContainer>::difference_type operator-(const random_access_const_iterator<RandomAccessContainer> &x, const random_access_const_iterator<RandomAccessContainer> &y){
	return  (typename random_access_const_iterator<RandomAccessContainer>::difference_type)x.m_idx
			- (typename random_access_const_iterator<RandomAccessContainer>::difference_type)y.m_idx;
};

template<class RandomAccessContainer>
inline random_access_const_iterator<RandomAccessContainer> operator+(typename random_access_const_iterator<RandomAccessContainer>::difference_type n, const random_access_const_iterator<RandomAccessContainer>& it){
	return it+n;
}

} // end namespace sdsl

#endif // end of include guard
