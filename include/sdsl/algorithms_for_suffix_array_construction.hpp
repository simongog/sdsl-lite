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
/*! \file algorithms_for_suffix_array_construction.hpp
    \brief algorithms_for_suffix_array_construction.hpp contains an interface to access suffix array construction algorithms
	\author Simon Gog
*/ 

#ifndef INCLUDED_SDSL_ALGORITHMS_FOR_SUFFIX_ARRAY_CONSTRUCTION
#define INCLUDED_SDSL_ALGORITHMS_FOR_SUFFIX_ARRAY_CONSTRUCTION

#include "config.hpp"
#include "int_vector.hpp"

#include "divsufsort.h"
#include "divsufsort64.h"

namespace sdsl{

namespace algorithm{	

//	
// Forward declarations
//----------------------------------------------------------

//! Calculates the Suffix Array for a text.
/*!
 * \param c Text (c-string) to calculate the suffix array. The lex. order is given by the ascii-codes of the characters.
 * \param len Length of the text. *(c+len)=0 and for i<len *(c+len)!=0
 * \param sa Reference to a RandomAccessContainer which will contain the result of the calculation. 
 * \pre sa.size() has to be equal to len.
 */
template<class RandomAccessContainer>
static void calculate_sa(const unsigned char *c, typename RandomAccessContainer::size_type len, RandomAccessContainer &sa);

template<uint8_t fixedIntWidth>
static void calculate_sa(const unsigned char *c, typename int_vector<fixedIntWidth>::size_type len, int_vector<fixedIntWidth> &sa);

template<class RandomAccessContainer>
void calculate_sa(const unsigned char *c, typename RandomAccessContainer::size_type len, RandomAccessContainer &sa){
	typedef typename RandomAccessContainer::size_type size_type;
	if(len<=1){ // handle special case 
		sa = RandomAccessContainer(len,0);
		return;
	}
	bool small_file = (sizeof(len) <= 4 or len < 0x7FFFFFFFULL );
	if( small_file ){
		int32_t *sufarray = new int32_t[len];
		divsufsort(c, sufarray, len);		
		for(size_type i=0; i<len; ++i) { sa[i] = sufarray[i]; }
		delete [] sufarray;
	}else{
		int64_t *sufarray = new int64_t[len];
		divsufsort64(c, sufarray, len);		
		for(size_type i=0; i<len; ++i) { sa[i] = sufarray[i]; }
		delete [] sufarray;
	}
}

template<uint8_t fixedIntWidth>
void calculate_sa(const unsigned char *c, typename int_vector<fixedIntWidth>::size_type len, int_vector<fixedIntWidth> &sa){
	typedef typename int_vector<fixedIntWidth>::size_type size_type;
	if(len <= 1){ // handle special case 
		sa = int_vector<fixedIntWidth>(len,0);
		return;
	}
	bool small_file = (sizeof(len) <= 4 or len < 0x7FFFFFFFULL );
	if( small_file ){
		uint8_t oldIntWidth = sa.get_int_width();
		if( 32 == fixedIntWidth or (0==fixedIntWidth and 32 >= oldIntWidth) ){
			sa.set_int_width(32);
			sa.resize( len );
			divsufsort(c, (int32_t*)sa.m_data, len);
			// copy integers back to the right positions
			if(oldIntWidth!=32){
				for(size_type i=0; i<len; ++i) { sa.set_int(i*oldIntWidth, sa.get_int(i<<5, 32), oldIntWidth);  }
				sa.set_int_width(oldIntWidth);
				sa.resize(len);
			}
		}
		else{
			if( sa.get_int_width() < bit_magic::l1BP(len)+1 ){
				throw std::logic_error( "width of int_vector is to small for the text!!!" ); 
			}
			int32_t *sufarray = new int32_t[len];
			divsufsort(c, sufarray, len);		
			for(size_type i=0; i<len; ++i) { sa[i] = sufarray[i]; }
			delete [] sufarray;
		}
	}else{
		uint8_t oldIntWidth = sa.get_int_width();
		sa.set_int_width(64);
		sa.resize( len );
		divsufsort64(c, (int64_t*)sa.m_data, len);
		// copy integers back to the right positions
		if(oldIntWidth!=64){
			for(size_type i=0; i<len; ++i) { sa.set_int(i*oldIntWidth, sa.get_int(i<<6, 64), oldIntWidth);  }
			sa.set_int_width(oldIntWidth);
			sa.resize(len);
		}
	}
}

// \param c Char array pointing to the text
// \param n Length of the text
bool shift_text(char *c, uint64_t n, bool shift=true);

} // end namespace algorithm

} // end namespace sdsl

#endif
