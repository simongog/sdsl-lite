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
    \brief algorithms_for_suffix_array_construction.hpp contains algorithms for the construction of suffixa rrays.
	\author Simon Gog
*/ 
// TODO: force use of DC3 in some situtations

#ifndef INCLUDED_SDSL_ALGORITHMS_FOR_SUFFIX_ARRAY_CONSTRUCTION
#define INCLUDED_SDSL_ALGORITHMS_FOR_SUFFIX_ARRAY_CONSTRUCTION

#include "int_vector.hpp"

#cmakedefine divsufsort_FOUND
#cmakedefine divsufsort64_FOUND

#ifdef divsufsort_FOUND
	#include "divsufsort.h"
#endif

#ifdef divsufsort64_FOUND
	#include "divsufsort64.h"
#endif

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

/* _radixPass stable sorts the sequence a of length n  according 
   to the ordering given by r. The result is written in
   b. K is the size of the alphabet of r. 
 */
template<typename I>
static void _radixPass(I *a, I *b, I *r, I n, I K);	

template<typename I>
static void _calcSuffixArrayDC3(I *s, I *sa, I n, I K);	

//
// Implementation
//----------------------------------------------------------

template<class RandomAccessContainer>
void calculate_sa(const unsigned char *c, typename RandomAccessContainer::size_type len, RandomAccessContainer &sa){
	typedef typename RandomAccessContainer::size_type size_type;
	if(len<=1){ // handle special case 
		sa = RandomAccessContainer(len,0);
		return;
	}
	bool small_file = (sizeof(len) <= 4 or len < 0x7FFFFFFFULL );
	#if !defined(divsufsort_FOUND) || !defined(divsufsort64_FOUND)
		bool divss=false, divss64=false;
		#ifdef divsufsort_FOUND
			divss = true;
		#endif
		#ifdef divsufsort_FOUND64
			divss64 = true;
		#endif
		uint8_t K=0; // size of the alphabet without zero
		uint8_t occ[256] = {0}; // indicator if a char occurs in c
		const unsigned char *cp = c; // char pointer to c
		if( (small_file and !divss) or (!small_file and !divss64) ){
			for(size_type i=0; i<len; ++i){ occ[*(cp++)] = 1; }
			for(uint16_t i=0;i<256; ++i){
				if(occ[i]) ++K;
				if(i) occ[i] += occ[i-1];
			}
			cp = c;
		}
	#endif
	if( small_file ){
		#if defined(divsufsort_FOUND)
			int32_t *sufarray = new int32_t[len];
			divsufsort(c, sufarray, len);		
			for(size_type i=0; i<len; ++i) { sa[i] = sufarray[i]; }
			delete [] sufarray;
		#else	
			uint32_t *s = new uint32_t[len+3], *suffarray = new uint32_t[len];
			for(size_type i=0; i<len; ++i) s[i] = occ[*(cp++)];
			s[len] = s[len+1] = s[len+2] = 0; 
			_calcSuffixArrayDC3(s, suffarray, (uint32_t)len, (uint32_t)K);
			for(size_type i=0; i<len; ++i) sa[i] = suffarray[i];
			delete [] s; delete [] suffarray;
		#endif	
	}else{
		#if defined(divsufsort64_FOUND)
			int64_t *sufarray = new int64_t[len];
			divsufsort64(c, sufarray, len);		
			for(size_type i=0; i<len; ++i) { sa[i] = sufarray[i]; }
			delete [] sufarray;
		#else	
			uint64_t *s = new uint64_t[len+3], *suffarray = new uint64_t[len];
			for(size_type i=0; i<len; ++i) s[i] = occ[(uint8_t)*(cp++)];
			s[len] = s[len+1] = s[len+2]= 0;	
			_calcSuffixArrayDC3(s, suffarray, (uint64_t)len, (uint64_t)K);
			for(size_type i=0; i<len; ++i) sa[i] = suffarray[i];
			delete [] s; delete [] suffarray;
		#endif	
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
		#if defined(divsufsort_FOUND)
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
		#else
			calculate_sa< int_vector<fixedIntWidth> >(c, len, sa);
		#endif	
	}else{
		#if defined(divsufsort64_FOUND)
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
		#else
			calculate_sa< int_vector<fixedIntWidth> >(c, len, sa);
		#endif	
	}
}


// The following code (function _radixPass and _calcSuffixArrayDC3) was originally written by Juha Kärkkainen.
// 
//

template<typename I>
void _radixPass(I *a, I *b, I *r, I n, I K){
	I *c = new I[K+1]; // counter array
	for(I i = 0; i <= K; ++i) c[i] = 0;
	for(I i = 0; i < n; ++i) c[r[a[i]]]++;
	for(I i=0, sum=0; i<=K; ++i){
		I t = c[i]; c[i] = sum; sum += t;
	}
	for(I i=0; i<n; ++i) b[c[r[a[i]]]++] = a[i];
	delete [] c;
}

template<typename I>
void _calcSuffixArrayDC3(I *s, I *sa, I n, I K){
  I n0=(n+2)/3, n1=(n+1)/3, n2=n/3;// number of suffixs mod {0,1,2}==0; n0>=n1>=n2
  I n02=n0+n2; 
  I *s12  = new I[n02 + 3];  s12[n02]= s12[n02+1]= s12[n02+2]=0; 
  I *SA12 = new I[n02 + 3]; SA12[n02]=SA12[n02+1]=SA12[n02+2]=0;
  I* s0   = new I[n0];
  I* SA0  = new I[n0];
  
  I j=0;
  for (I i=0/*, j=0*/;  i < n+(n0-n1);  ++i) if (i%3 != 0) s12[j++] = i;
  // lsb radix sort the mod 1 and mod 2 triples
  _radixPass(s12 , SA12, s+2, n02, K);
  _radixPass(SA12, s12 , s+1, n02, K);  
  _radixPass(s12 , SA12, s  , n02, K);

  // find lexicographic names of triples
  I name = 0, c0 = 0, c1 = 0, c2 = 0;
  if( n02 > 0 ) c0 = s[SA12[0]]+1;
  for (I i = 0;  i < n02;  ++i) {
    if (s[SA12[i]] != c0 || s[SA12[i]+1] != c1 || s[SA12[i]+2] != c2) { 
      name++;  c0 = s[SA12[i]];  c1 = s[SA12[i]+1];  c2 = s[SA12[i]+2];
    }
    if (SA12[i] % 3 == 1) { s12[SA12[i]/3]      = name; } // left half
    else                  { s12[SA12[i]/3 + n0] = name; } // right half
  }
//  cerr<<"found recursive names"<<endl;
  // recurse if names are not yet unique
  if (name < n02) {
    _calcSuffixArrayDC3(s12, SA12, n02, name);
    // store unique names in s12 using the suffix array 
    for (I i = 0;  i < n02;  ++i) s12[SA12[i]] = i + 1;
  } else // generate the suffix array of s12 directly
    for (I i = 0;  i < n02;  ++i) SA12[s12[i] - 1] = i; 

  // stably sort the mod 0 suffixes from SA12 by their first character
  for (I i=0, j=0;  i < n02;  ++i) if (SA12[i] < n0) s0[j++] = 3*SA12[i];
  _radixPass(s0, SA0, s, n0, K);

  // merge sorted SA0 suffixes and sorted SA12 suffixes
  for (I p=0,  t=n0-n1,  k=0;  k < n;  ++k) {
#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
    I i = GetI(); // pos of current offset 12 suffix
    I j = SA0[p]; // pos of current offset 0  suffix
    if (SA12[t] < n0 ? 
        leq(s[i],       s12[SA12[t] + n0], s[j],       s12[j/3]) :
        leq(s[i],s[i+1],s12[SA12[t]-n0+1], s[j],s[j+1],s12[j/3+n0]))
    { // suffix from SA12 is smaller
      sa[k] = i;  t++;
      if (t == n02) { // done --- only SA0 suffixes left
        for (k++;  p < n0;  p++, k++) sa[k] = SA0[p];
      }
    } else { 
      sa[k] = j;  p++; 
      if (p == n0)  { // done --- only SA12 suffixes left
        for (k++;  t < n02;  t++, k++) sa[k] = GetI(); 
      }
    }  
  } 
  delete [] s12; delete [] SA12; delete [] SA0; delete [] s0; 
} 

// \param c Char array pointing to the text
// \param n Length of the text
bool shift_text(char *c, uint64_t n, bool shift=true);

} // end namespace algorithm

} // end namespace sdsl

#endif
