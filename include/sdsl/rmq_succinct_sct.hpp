/* sdsl - succinct data structures library
    Copyright (C) 2009 Simon Gog 

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
/*! \file rmq_support_sct.hpp
    \brief rmq_support_sct.hpp contains the class rmq_support_sct which supports range minimum or range maximum queries on a random access container in constant time and \f$2 n+o(n) bits\f$ space. 
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_RMQ_SUPPORT_SCT
#define INCLUDED_SDSL_RMQ_SUPPORT_SCT

#include "rmq_support.hpp"
#include "int_vector.hpp"
#include "algorithms_for_compressed_suffix_trees.hpp"
#include "bp_support_sada.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl{

	
template<class RandomAccessContainer = int_vector<>, bool Minimum = true, class Bp_support = bp_support_sada<> >	
class rmq_support_sct;

template<class RandomAccessContainer = int_vector<>, class Bp_support = bp_support_sada<> >
struct range_maximum_support_sct{
	typedef rmq_support_sct<RandomAccessContainer, false, Bp_support> type;
};

//! A class to support range minimum or range maximum queries on a random access container.
/*!
 * This class takes three template parameters:
 *  - RandomAccessContainer is the type of random access container, for which the range minimum/maximum support should be build,
 *  - minumum 				specifies whether the data structure should answer range minimum queries (mimumum=true) of range maximum queries (maximum=false), and
 *  - Bp_support			is the support structure for the balanced parentheses sequence of the Super-Cartesian tree used internal in the class.
 * \par Time complexity
 *		\f$ \Order{1} \f$ for the range minimum/maximum queries if the balanced parentheses support structure supports constant time operations.
 * \par Space complexity:
 *		\f$ \Order{2n}+o(n) \f$ bits for the data structure ( \f$ n=size() \f$ ).
 *
 */	
template<class RandomAccessContainer, bool Minimum, class Bp_support>	
class rmq_support_sct{
	const RandomAccessContainer *m_v;
	bit_vector					m_sct_bp; 		//!< A bit vector which contains the balanced parentheses sequence of the Super-Cartesian tree of the input container.
	Bp_support					m_sct_bp_support; 	//!< Support structure for the balanced parentheses of the Super-Cartesian tree.

	void construct(){
		if( m_v == NULL ){
			m_sct_bp = bit_vector(0); m_sct_bp_support = Bp_support();
		}else{
#ifdef RMQ_SCT_BUILD_BP_NOT_SUCCINCT			
			// this method takes \f$n\log n\f$ bits extra space in the worst case
			algorithm::construct_supercartesian_tree_bp(*m_v, m_sct_bp);
#else			
			// this method takes only \f$n\f$ bits extra space in all cases 
			algorithm::construct_supercartesian_tree_bp_succinct(*m_v, m_sct_bp); 
			// TODO: falls alle werte im bereich von 0..n liegen sind nur 2n bits noetig
		   //  TODO: constructor mit int_vector_file_buffer	
#endif			
			m_sct_bp_support = Bp_support(&m_sct_bp);
		}
	}

	void copy(const rmq_support_sct &rm){
		m_v = rm.m_v;
		m_sct_bp = rm.m_sct_bp;
		m_sct_bp_support = rm.m_sct_bp_support;
		m_sct_bp_support.set_vector(&m_sct_bp);
	}

	public:
	typedef typename RandomAccessContainer::size_type size_type;
	typedef typename RandomAccessContainer::size_type value_type;

	const bit_vector &sct_bp;
	const Bp_support &sct_bp_support;

	//! Constructor
	rmq_support_sct(const RandomAccessContainer *v=NULL):m_v(v), sct_bp(m_sct_bp), sct_bp_support(m_sct_bp_support){
		construct();
	}

	//! Copy constructor
	rmq_support_sct(const rmq_support_sct &rm){
		if( this != &rm ){ // if v is not the same object
			copy(rm);
		}
	}

	//! Destructor
	~rmq_support_sct(){ }

	rmq_support_sct& operator=(const rmq_support_sct &rm){
		if( this != &rm ){
			copy(rm);
		}
		return *this;
	}

	void set_vector(const RandomAccessContainer *v){
		m_v = v;
	}

	//! Range minimum/maximum query for the supported random access container v.
	/*! 
	 * \param l Leftmost position of the interval \f$[\ell..r]\f$.  
	 * \param r Rightmost position of the interval \f$[\ell..r]\f$.
	 * \return The minimal index i with \f$\ell \leq i \leq r\f$ for which \f$ v[i] \f$ is minimal/maximal. 
	 * \pre 
	 *   - r < size()
	 *   - \f$ \ell \leq r \f$
	 * \par Time complexity
	 *      \f$ \Order{1} \f$
	 */
	size_type operator()(const size_type l, const size_type r)const{
		assert( l <= r ); assert( r < size() );
		if( l==r )
			return l;
		size_type i		= m_sct_bp_support.select(l+1);
		size_type j		= m_sct_bp_support.select(r+1);
		size_type fc_i	= m_sct_bp_support.find_close(i);
		if( j < fc_i ){ // i < j < find_close(j) < find_close(i)
			return l;
		}else{ // if i < find_close(i) < j < find_close(j)
			size_type ec = m_sct_bp_support.rr_enclose(i,j);
			if( ec == m_sct_bp_support.size() ){// no restricted enclosing pair found
				return r;
			}else{// found range restriced enclosing pair
				return m_sct_bp_support.rank(ec)-1; // subtract 1, as the index is 0 based
			}
		}
	}

	size_type size()const{
		return m_sct_bp.size()/2;
	}

	size_type serialize(std::ostream &out)const{
		size_type written_bytes = 0;
		written_bytes += m_sct_bp.serialize(out);
		written_bytes += m_sct_bp_support.serialize(out);
		return written_bytes;
	}

	void load(std::istream &in, const RandomAccessContainer *v){
		m_sct_bp.load(in);
		m_sct_bp_support.load(in, &m_sct_bp);
	}
};


}


#endif
