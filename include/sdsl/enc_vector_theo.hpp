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
/*! \file enc_vector_theo.hpp
   \brief enc_vector_theo.hpp contains the sdsl::enc_vector_theo class and the iterator class for enc_vector_theo.
   \author Simon Gog
*/  
#include "int_vector.hpp"
#include "elias_delta_coder.hpp"
#include "rank_support.hpp"
#include "select_support.hpp"

#ifndef SDSL_ENC_VECTOR_THEO
#define SDSL_ENC_VECTOR_THEO

//! Namespace for the succinct data structure library.
namespace sdsl{

template<class EncVector>	
class enc_vector_theo_const_iterator; // forward declaration

template<uint8_t>
struct enc_vector_theo_trait{
	typedef int_vector<0> int_vector_type;
};	

template<>
struct enc_vector_theo_trait<32>{
	typedef int_vector<32> int_vector_type;
};

//! A generic immutable space-saving vector class for unsigned integers. It encodes each integer with its self-delimiting code and still provides constant time access.
/*! The values of a enc_vector_theo are immutable after the constructor call. The class
    could be parametrized with a self-delimiting codes Coder and Rank/Select-Support
	for constant time access to the elements.
 */
template<class Coder=coder::elias_delta,
		 uint32_t SampleDens = 1,
		 class RankSupport=rank_support_v<>,
		 class SelectSupport = select_support_mcl<1>,
		 uint8_t fixedIntWidth = 0 >
class enc_vector_theo{
	public:
	typedef uint64_t 							value_type;  	// STL Container requirement
//	typedef enc_vector_theo_const_iterator<Coder, SampleDens, RankSupport, SelectSupport, fixedIntWidth>	iterator;    	// STL Container requirement
	typedef enc_vector_theo_const_iterator<enc_vector_theo<Coder, SampleDens, RankSupport, SelectSupport, fixedIntWidth> >	iterator;    	// STL Container requirement
	typedef iterator							const_iterator; // STL Container requirement
	typedef const value_type		 			reference;
	typedef const value_type 					const_reference;
//	typedef SDSBitVectorReference*				pointer;
	typedef const value_type*					const_pointer;
	typedef ptrdiff_t 							difference_type;// STL Container requirement
	typedef int_vector<>::size_type				size_type;		// STL Container requirement
	typedef Coder								coder;
	typedef RankSupport							rank_support;
	typedef SelectSupport						select_support;
	static  const uint32_t 						sample_dens	= SampleDens;

	friend class enc_vector_theo_const_iterator<enc_vector_theo<Coder, SampleDens, RankSupport, SelectSupport, fixedIntWidth> >;    	// STL Container requirement
	private:
	int_vector<0> 	m_z; 		// compressed bit stream
	int_vector<1> 	m_sample;   // indicator for a sample

	typename enc_vector_theo_trait<fixedIntWidth>::int_vector_type m_sample_pointer;
	typename enc_vector_theo_trait<fixedIntWidth>::int_vector_type m_sample_vals;

	RankSupport		m_sample_rank; // rank support for m_sample
	size_type		m_elements;    // number of elements

	// workaround function for the constructor
	void construct(){
		m_elements = 0;
	}
	void copy(const enc_vector_theo &v);
	public:
		//! Default Constuctor
		enc_vector_theo(){
			construct();
		};
		//! Copy constructor
		/*! \param v The enc_vector_theo to copy.
		  	Required for the Assignable Concept of the STL
		 */
		enc_vector_theo(const enc_vector_theo &v);

		//! Constructor for a Container of positiv integers.
		/*! \param c A container of positive integers.
		    \par The container is used to build the EncVector of the
			     integer sequence.
		  */
		template<class Container>
		enc_vector_theo(const Container &c){
			construct();
			init(c);
		};

		template<class Container>
		void init(const Container &c);

		//! Default Destructor
		~enc_vector_theo(){
		};

		//! The number of elements in the enc_vector_theo.
		/*!
		 	Required for the Container Concept of the STL.
			\sa max_size
		 */
		size_type size()const;

		//! Return the largest size that this container can ever have.
		/*! Required for the Container Concept of the STL.
		 */
		static size_type max_size();

		//!	Returns if the enc_vector_theo is empty.
		/*! Equivalent to size() == 0.
		 *
		 * 	Required for the STL Container Concept.
		 *  \sa size()
		 */
		bool empty() const;

		//! Swap method for enc_vector_theo
		/*! The swap method can be defined in terms of assignment.
			This requires three assignments, each of which, for a container type, is linear
			in the container's size. In a sense, then, a.swap(b) is redundant.
			This implementation guaranties a run-time complexity that is constant rather than linear.
			\param v enc_vector_theo to swap.

			Required for the Assignable Conecpt of the STL.
		  */
		void swap(enc_vector_theo &v);

		//! Iterator that points to the first element of the enc_vector_theo.
		/*!
		 * 	Required for the Container Concept of the STL.
		 *  \sa end()
		 */
		const const_iterator begin()const;

		//! Iterator that points to the position after the last element of the enc_vector_theo.
		/*!
		 *	Required for the Container Concept of the STL
		 *  \sa begin()
	     */
		const const_iterator end()const;

		// Iterator that points to the last element of the enc_vector_theo.
		/*
		 * 	Required for the Container Concept of the STL.
		 *  \sa rend()
		 */
//		reverse_iterator rbegin()const;

		// Iterator that points to the position before the first element of the enc_vector_theo.
		/*
		 *	Required for the Container Concept of the STL
		 *  \sa rbegin()
	     */
//		reverse_iterator rend()const;

		//! []-operator
		/*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
		 *
		 *  Required for the STL Random Access Container Concept.
		 */
		value_type operator[](size_type i)const;

		//! Assignment Operator
		/*!
		 *	Required for the Assignable Concept of the STL.
		 */
		enc_vector_theo& operator=(const enc_vector_theo &v);

		//! Equality Operator
		/*! Two enc_vector_theos are equal if all member variables are equal
		 *  (including the sample density of the enc_vector_theos).
		 *  \note If the sample density is not equal you should use
		 *  SDSAlgorithm::equal_container_values to compare two enc_vector_theos.
		 *
		 * 	Required for the Equality Comparable Concept of the STL.
		 *  \sa operator!=
		 */
		bool operator==(const enc_vector_theo &v)const;

		//! Unequality Operator
		/*! Two enc_vector_theos are unuequal if not all member variables are equal
		 *  (including the sample density of the enc_vector_theos).
		 *  \note If the sample density is not equal you should use
		 *  SDSAlgorithm::equal_container_values to compare two enc_vector_theos.
		 *
		 * 	Required for the Equality Comparable Concept of the STL.
		 *  \sa operator==
		 */
		bool operator!=(const enc_vector_theo &v)const;

		//! Serialzes the enc_vector_theo to a stream.
		/*! \param out Outstream to write the data structure.
		    \return The number of written bytes.
		 */
		size_type serialize(std::ostream &out) const;

		//! Load the enc_vector_theo from a stream.
		void load(std::istream &in);

		//! Returns the ith sample of enc_vector_theo_prac
		/*! \param i The index of the sample. 0 <= i < size()/get_sample_dens()
		 *  \return The value of the ith sample.
		 */		
		value_type sample(const size_type i) const;
};


// TODO: uint64_t als value type ersetzen durch EncVector::value_type geht nicht??? mit using....
template<class EncVector>
struct enc_vector_theo_const_iterator: public std::iterator<std::random_access_iterator_tag, uint64_t,  ptrdiff_t>{
	typedef const reference	const_reference;
	typedef enc_vector_theo_const_iterator<EncVector> const_iterator;
	typedef typename EncVector::size_type size_type;
	typedef typename EncVector::value_type value_type;
	typedef typename EncVector::difference_type difference_type;

	const EncVector *m_v;		// enc_vector_theo the interator belongs to
	value_type m_decoded_val[((size_type(EncVector::sample_dens))<<6)/EncVector::coder::min_codeword_length+1];// buffer for decoded values
	static const size_type m_decoded_val_size = ((size_type(EncVector::sample_dens))<<6)/EncVector::coder::min_codeword_length;// buffer for decoded values
	size_type m_decoded_val_start_idx;		// Index of the first element that is buffered
	size_type m_decoded_val_end_idx;		// Index of the first element after the buffered elements
	size_type m_idx;  // Index of the current element

	//! Constructor
	/*!
	 *	\param v Pointer to the enc_vector_theo that is supported.
	 *	\param idx Index of the i-th element. \f$idx\in [0..v->size()]\f$
	 */
	enc_vector_theo_const_iterator(const EncVector *v, size_type idx=0){
		m_v = v;
		m_idx = idx;
		m_decoded_val_start_idx = 0;
		m_decoded_val_end_idx = 0;
	}

	enc_vector_theo_const_iterator(const enc_vector_theo_const_iterator &it){
		enc_vector_theo_const_iterator(it.m_v, it.m_idx);
	}

	~enc_vector_theo_const_iterator(){ }

	inline const_reference operator*();

	//! Prefix increment of the Iterator
	const_iterator& operator++(){
		++m_idx;
		return *this;
	}

	//! Postfix increment of the Iterator
	const_iterator operator++(int x){
		enc_vector_theo_const_iterator it = *this;
		++(*this);
		return it;
	}

	//! Prefix decrement of the Iterator
	const_iterator& operator--(){
		--m_idx;
		return *this;
	}

	//! Postfix decrement of the Iterator
	const_iterator operator--(int x){
		enc_vector_theo_const_iterator it = *this;
		++(*this);
		return it;
	}

	const_iterator& operator+=(difference_type i){
		if(i<0)
			return *this -= (-i);
		m_idx += i;
		return *this;
	}

	const_iterator& operator-=(difference_type i){
		if(i<0)
			return *this += (-i);
		m_idx -= i;
		return *this += -i;
	}

	const_iterator operator+(difference_type i) const{
		const_iterator it = *this;
		return it += i;
	}

	const_iterator operator-(difference_type i) const{
		const_iterator it = *this;
		return it -= i;
	}

	const_reference operator[](difference_type i) const{
		return *(*this + i);
	}

	bool operator==(const enc_vector_theo_const_iterator &it)const{
		// supported vectors and index have to be equal
		return m_idx == it.m_idx and m_v == it.m_v;
	}

	bool operator!=(const enc_vector_theo_const_iterator &it)const{
		return m_idx != it.m_idx or m_v != it.m_v;
	}

	bool operator<(const enc_vector_theo_const_iterator &it)const{
		return m_idx < it.m_idx;
	}

	bool operator>(const enc_vector_theo_const_iterator &it)const{
		return m_idx > it.m_idx;
	}

	bool operator>=(const enc_vector_theo_const_iterator &it)const{
		return !(*this < it);
	}

	bool operator<=(const enc_vector_theo_const_iterator &it)const{
		return !(*this > it);
	}

};

template<class Coder,
		 uint32_t SampleDens,
		 class RankSupport,
		 class SelectSupport, 
		 uint8_t fixedIntWidth>
inline typename enc_vector_theo<Coder, SampleDens, RankSupport, SelectSupport, fixedIntWidth>::value_type enc_vector_theo<Coder, SampleDens, RankSupport, SelectSupport, fixedIntWidth>::operator[](const size_type i)const {
		if( i+1 == 0 || i >= m_elements  ){
			throw std::out_of_range("OUT_OF_RANGE_ERROR: enc_vector_theo::operator[](size_type); idx >= size()!");
			return 0;
		}
		const size_type sample_rank	= m_sample_rank.rank(i+1);
		if(m_sample[i]){// -> sample_idx == i
			return m_sample_vals[sample_rank-1];
		}
		size_type sample_idx	= bit_magic::prev(m_sample.data(), i);
		size_type idx_of_sample_in_z = m_sample_pointer[sample_rank-1];
		return m_sample_vals[sample_rank-1] + Coder::decode_prefix_sum(m_z.data(), idx_of_sample_in_z, i-sample_idx);
}

template<class Coder,
		 uint32_t SampleDens,
		 class RankSupport,
		 class SelectSupport, 
		 uint8_t fixedIntWidth>
inline typename enc_vector_theo<Coder, SampleDens, RankSupport, SelectSupport, fixedIntWidth>::value_type enc_vector_theo<Coder, SampleDens, RankSupport, SelectSupport, fixedIntWidth>::sample(const size_type i)const {
		if( i >= m_sample_vals.size()  ){
			throw std::out_of_range("OUT_OF_RANGE_ERROR: enc_vector_theo::sample(size_type); i*SampleDens >= size()!");
			return 0;
		}
		return m_sample_vals[i];
}
template<class Coder,
		 uint32_t SampleDens,
		 class RankSupport,
		 class SelectSupport,
		 uint8_t fixedIntWidth>
inline enc_vector_theo<>::size_type enc_vector_theo<Coder, SampleDens, RankSupport, SelectSupport, fixedIntWidth>::size()const
{
	return m_elements;
}

template<class Coder,
		 uint32_t SampleDens,
		 class RankSupport,
		 class SelectSupport,
		 uint8_t fixedIntWidth>
inline enc_vector_theo<>::size_type enc_vector_theo<Coder, SampleDens, RankSupport, SelectSupport, fixedIntWidth>::max_size()
{
	return int_vector<>::max_size()/2; // each element could possible occupy double space with selfdelimiting codes
}

template<class Coder,
		 uint32_t SampleDens,
		 class RankSupport,
		 class SelectSupport,
		 uint8_t fixedIntWidth>
inline bool enc_vector_theo<Coder, SampleDens, RankSupport, SelectSupport, fixedIntWidth>::empty()const
{
	return 0==m_elements;
}


template<class Coder,
		 uint32_t SampleDens,
		 class RankSupport,
		 class SelectSupport,
		 uint8_t fixedIntWidth>
void enc_vector_theo<Coder, SampleDens,RankSupport, SelectSupport, fixedIntWidth>::copy(const enc_vector_theo<Coder, SampleDens,RankSupport, SelectSupport, fixedIntWidth> &v){
	m_z					= v.m_z;				// copy compressed bit stream
	m_sample			= v.m_sample;			// copy indicator for a sample
	m_sample_pointer	= v.m_sample_pointer;	// copy sample pointer vector
	m_sample_vals		= v.m_sample_vals;      // copy sample values
	m_sample_rank		= v.m_sample_rank;		// copy rank support for sample indicator
	m_sample_rank.set_vector(&m_sample);         // set rank supported bit_vector
	m_elements			= v.m_elements;			// copy number of stored elements
}

template<class Coder,
		 uint32_t SampleDens,
		 class RankSupport,
		 class SelectSupport,
		 uint8_t fixedIntWidth>
enc_vector_theo<Coder, SampleDens,RankSupport, SelectSupport, fixedIntWidth>::enc_vector_theo(const enc_vector_theo &v){
	copy(v);
}

template<class Coder,
		 uint32_t SampleDens,
		 class RankSupport,
		 class SelectSupport,
		 uint8_t fixedIntWidth>
enc_vector_theo<Coder, SampleDens,RankSupport, SelectSupport, fixedIntWidth>& enc_vector_theo<Coder, SampleDens,RankSupport, SelectSupport, fixedIntWidth>::operator=(const enc_vector_theo<Coder, SampleDens,RankSupport, SelectSupport, fixedIntWidth> &v){
	if( this != &v ){// if v and _this_ are not the same object
		copy(v);
	}
	return *this;
}

template<class Coder,
		 uint32_t SampleDens,
		 class RankSupport,
		 class SelectSupport,
		 uint8_t fixedIntWidth>
bool enc_vector_theo<Coder, SampleDens,RankSupport, SelectSupport, fixedIntWidth>::operator==(const enc_vector_theo<Coder, SampleDens,RankSupport, SelectSupport, fixedIntWidth> &v)const{
	if(this == &v)
		return true;
	return	 	m_elements == v.m_elements
		 and	m_z == v.m_z
		 and	m_sample == v.m_sample
		 and	m_sample_pointer == v.m_sample_pointer
		 and	m_sample_vals == v.m_sample_vals
		 and	m_sample_rank == v.m_sample_rank;
}

template<class Coder,
		 uint32_t SampleDens,
		 class RankSupport,
		 class SelectSupport,
		 uint8_t fixedIntWidth>
bool enc_vector_theo<Coder, SampleDens,RankSupport, SelectSupport, fixedIntWidth>::operator!=(const enc_vector_theo<Coder, SampleDens,RankSupport, SelectSupport, fixedIntWidth> &v)const{
	return !(*this == v);
}

template<class Coder,
		 uint32_t SampleDens,
		 class RankSupport,
		 class SelectSupport,
		 uint8_t fixedIntWidth>
void enc_vector_theo<Coder, SampleDens,RankSupport, SelectSupport, fixedIntWidth>::swap(enc_vector_theo<Coder, SampleDens,RankSupport, SelectSupport, fixedIntWidth> &v){
	if(this != &v){// if v and _this_ are not the same object
		m_z.swap(v.m_z);					// swap compressed bit streams
		m_sample_pointer.swap(v.m_sample_pointer); // swap the sample pointer vector
		m_sample_vals.swap(v.m_sample_vals);
		m_sample_rank.swap(v.m_sample_rank);	// swap rank support for sample indicator
//		m_sample_select.swap(v.m_sample_select); // swap select support for sample indicator
		m_sample.swap(v.m_sample);			// swap the sample vector
		std::swap(m_elements, v.m_elements);// swap the number of elements
	}
}


template<class Coder,
		 uint32_t SampleDens,
		 class RankSupport,
		 class SelectSupport,
		 uint8_t fixedIntWidth>
template<class Container>
void enc_vector_theo<Coder, SampleDens, RankSupport , SelectSupport, fixedIntWidth>::init(const Container &c){
	// clear BitVectors
	m_z.resize(0);
	m_elements = 0;
//	m_inc_start.resize(0);
	m_sample.resize(0);
	m_sample_pointer.resize(0);
	m_sample_vals.resize(0);
	if( c.empty() ) // if c is empty there is nothing to do...
		return;
//	m_inc_start.resize( c.size() ); util::setZeroBits(m_inc_start);
	m_sample.resize( c.size() );	util::set_zero_bits(m_sample);
	typename Container::const_iterator	it		 	= c.begin(), end = c.end();
	// determine the greatest value in the difference encoding
	typename Container::value_type 		max_value	= 1;//*(it++) + 1;
	// determine the greatest value of the samples
	typename Container::value_type		max_sample_value = (*it++) + 1;
	typename Container::value_type 		v1			= max_sample_value, v2;
	size_type z_size = 0;//Coder::encoding_length(v1);
	size_type sample_cnt = 1, max_sample_pointer = 0;
	m_sample[0]		= 1;
	const size_type threshold = (SampleDens<<6);
//std::cerr<<"calculate max_sample_pointer"<<std::endl;
	for(size_type i=1, z_partial_size = 0, elen=0; it != end; ++it,++i){
		if( ((v2=(*it)+1) <= v1) or (elen = Coder::encoding_length(v2-v1))+z_partial_size > threshold  ){// start of an increasing sequence or force encoding of sample
//			if( max_value < v2 ) max_value = v2;
			if( max_sample_value < v2 ) max_sample_value = v2;
			m_sample[i]		= 1;
			z_size += z_partial_size;
			max_sample_pointer = z_size;
			++sample_cnt;
			z_partial_size = 0;
		}else{// part of an increasing sequence
			z_partial_size += elen;
			if( max_value < v2-v1 ) max_value = v2-v1;
		}
		v1 = v2;
	}
//std::cerr<<"Calculate delta"<<std::endl;
	{
		int_vector<0> delta_c( c.size()-sample_cnt, 0, bit_magic::l1BP(max_value)+1 ); // Vector for difference encoding of c
		m_sample_pointer.set_int_width( bit_magic::l1BP(max_sample_pointer)+1 );
		m_sample_pointer.resize(sample_cnt);
		m_sample_vals.set_int_width(bit_magic::l1BP(max_sample_value)+1);
		m_sample_vals.resize(sample_cnt);
		int_vector<0>::iterator d_it = delta_c.begin();
		typename enc_vector_theo_trait<fixedIntWidth>::int_vector_type::iterator sp_it = m_sample_pointer.begin();
		typename enc_vector_theo_trait<fixedIntWidth>::int_vector_type::iterator sv_it = m_sample_vals.begin();
		it = c.begin();
		z_size = 0;
		for(int_vector<1>::const_iterator s_it = m_sample.begin(); it != c.end(); ++it, ++s_it){
			v2 = *it+1;
			if( *s_it ){// if sample
				*(sv_it++) = v2-1;
				*(sp_it++) = z_size;
//				std::cerr<<"sample  "<<v2<<std::endl;
			}else{
				*d_it = v2-v1;
				z_size += Coder::encoding_length( *d_it );
				++d_it;
			}
			v1 = v2;
		}
//std::cerr<<"Encode delta "<<delta_c.size()<<" integers"<<std::endl;
		Coder::encode(delta_c, m_z); // encode delta_c to m_z
	}
//	delta_c.resize(0);
//std::cerr<<"Calc rank"<<std::endl;
	m_sample_rank.init(&m_sample);  // init rank for m_sample
//std::cerr<<"Calc select"<<std::endl;
//	m_sample_select.init(&m_sample); // init select for m_sample
//std::cerr<<"Finished "<<std::endl;,
	m_elements = c.size();
}

template<class Coder,
		 uint32_t SampleDens,
		 class RankSupport,
		 class SelectSupport,
		 uint8_t fixedIntWidth>
enc_vector_theo<>::size_type enc_vector_theo<Coder, SampleDens,RankSupport, SelectSupport, fixedIntWidth>::serialize(std::ostream &out) const{
	size_type written_bytes = 0;
	out.write((char *) &m_elements, sizeof(m_elements));
	written_bytes += sizeof(m_elements);
	written_bytes += m_z.serialize(out);
	written_bytes += m_sample.serialize(out);
	written_bytes += m_sample_pointer.serialize(out);
	written_bytes += m_sample_vals.serialize(out);
	written_bytes += m_sample_rank.serialize(out);
	return written_bytes;
}

template<class Coder,
		 uint32_t SampleDens,
		 class RankSupport,
		 class SelectSupport,
		 uint8_t fixedIntWidth>
void enc_vector_theo<Coder, SampleDens,RankSupport, SelectSupport, fixedIntWidth>::load(std::istream &in){
	in.read((char *) &m_elements, sizeof(m_elements));
	m_z.load(in);
	m_sample.load(in);
	m_sample_pointer.load(in);
	m_sample_vals.load(in);
	m_sample_rank.load(in, &m_sample);
}

template<class Coder,
		 uint32_t SampleDens,
		 class RankSupport,
		 class SelectSupport,
		 uint8_t fixedIntWidth>
const typename enc_vector_theo<Coder,SampleDens, RankSupport, SelectSupport, fixedIntWidth>::const_iterator enc_vector_theo<Coder, SampleDens,RankSupport, SelectSupport, fixedIntWidth>::begin()const{
	return const_iterator(this, 0);
}

template<class Coder,
		 uint32_t SampleDens,
		 class RankSupport,
		 class SelectSupport,
		 uint8_t fixedIntWidth>
const typename enc_vector_theo<Coder,SampleDens,RankSupport, SelectSupport, fixedIntWidth>::const_iterator enc_vector_theo<Coder, SampleDens,RankSupport, SelectSupport, fixedIntWidth>::end()const{
	return const_iterator(this, this->m_elements);
}

template<class EncVector>
typename enc_vector_theo_const_iterator<EncVector>::const_reference enc_vector_theo_const_iterator<EncVector>::operator*(){
	// if requested element is buffered
	if( m_idx >= m_decoded_val_start_idx and m_idx < m_decoded_val_end_idx ){
		return m_decoded_val[m_idx-m_decoded_val_start_idx];
	}
	else{
		if(m_idx >= m_v->m_elements ){
			throw std::out_of_range("enc_vector_theo_const_iterator: dereferencing failed!");
		}
		size_type sample_rank	= m_v->m_sample_rank.rank(m_idx+1);// get number of samples before m_idx
		size_type sample_idx = bit_magic::prev(m_v->m_sample.data(), m_idx);

		size_type idx_of_sample_in_z = m_v->m_sample_pointer[sample_rank-1];// get idx in z
		size_type n = 0;//m_decoded_val_size;
		if( sample_rank == m_v->m_sample_vals.size() ){// last one
			n = m_v->m_elements - sample_idx;
		}
		else{
			n = bit_magic::next(m_v->m_sample.data(), m_idx) - sample_idx;
		}
		if(n>1){
			EncVector::coder::template decode<false, true, uint64_t* >(m_v->m_z.data(), idx_of_sample_in_z, n-1, m_decoded_val+1);
		}
		m_decoded_val_start_idx = sample_idx;
		m_decoded_val_end_idx	= m_decoded_val_start_idx + n;
		value_type *vals = m_decoded_val;
		*(vals++) 	   = m_v->m_sample_vals[sample_rank-1];
		for(size_type i=1;i<n;++i, ++vals){
				*vals += *(vals-1);
		}
		return m_decoded_val[m_idx-m_decoded_val_start_idx];
	}
}

} // end namespace sdsl

#endif
