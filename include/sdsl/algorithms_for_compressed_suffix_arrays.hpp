/* sdsl - succinct data structures library
    Copyright (C) 20010 Simon Gog 

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
/*! \file algorithms_for_compressed_suffix_arrays.hpp
    \brief algorithms_for_compressed_suffix_arrays.hpp contains algorithms for compressed suffix arrays.
	\author Simon Gog
*/ 
#ifndef INCLUDED_SDSL_ALGORITHMS_FOR_COMPRESSED_SUFFIX_ARRAYS
#define INCLUDED_SDSL_ALGORITHMS_FOR_COMPRESSED_SUFFIX_ARRAYS

#include "int_vector.hpp" // for bit_vector
#include <stack> // for calculate_supercartesian_tree_bp

#ifdef SDSL_DEBUG_ALGORITHMS_FOR_COMPRESSED_SUFFIX_ARRAYS
#include "testutils.hpp"
#endif

namespace sdsl{

namespace algorithm{	

template<class Csa>	
void set_text(const unsigned char *str, typename Csa::size_type len, typename Csa::size_type *C, unsigned char *char2comp, unsigned char *comp2char, uint8_t &sigma){
#ifdef SDSL_DEBUG_ALGORITHMS_FOR_COMPRESSED_SUFFIX_ARRAYS
	stop_watch sw; sw.start();
#endif	
	typedef typename Csa::size_type size_type;
	for(uint16_t i=0; i<256; ++i)
		char2comp[i] = comp2char[i] = C[i] = 0;
	C[256] = 0;
	if( str == NULL or len ==0 )
		return;
	const unsigned char *p = str;
//	assert(str[len-1]==0); // muss bei bwt z.B. nicht so sein
	for(size_type i=0;i<len;++i){
		++C[*p++];
	}
	assert(1 == C[0]);
	size_type value = 0;
	for(uint16_t i=0; i<256; ++i)
		if(C[i]){
			char2comp[i] 	= value;
			C[value]		= C[i];
			if( i > value )
				C[i]=0;
			++value;	
		}
	sigma = value;
	for(uint16_t i=0; i<256; ++i)
		if(char2comp[i])
			comp2char[char2comp[i]] = i;
	for(uint16_t i=256;i>0; --i) C[i] = C[i-1];
	C[0] = 0;
	for(uint16_t i=1; i<257; ++i) C[i] += C[i-1];	
	if(C[256]!=len){
		std::cerr<<"C[256]="<<C[256]<<" "<<len<<std::endl;
	}
	assert(C[256]==len);
	assert(C[sigma+1]==len);
#ifdef SDSL_DEBUG_ALGORITHMS_FOR_COMPRESSED_SUFFIX_ARRAYS
	sw.stop();
	std::cerr<<"set_text takes "<<sw.get_real_time()<<" ms real time";
	std::cerr<<"set_text takes "<<sw.get_user_time()<<" ms user time";
#endif	
}

template<class Csa, class size_type_class>	
void set_text(int_vector_file_buffer<8, size_type_class> &str_buf, typename Csa::size_type len, typename Csa::size_type *C, unsigned char *char2comp, unsigned char *comp2char, uint8_t &sigma){
#ifdef SDSL_DEBUG_ALGORITHMS_FOR_COMPRESSED_SUFFIX_ARRAYS
	stop_watch sw; sw.start();
#endif	
	typedef typename Csa::size_type size_type;
	for(uint16_t i=0; i<256; ++i)
		char2comp[i] = comp2char[i] = C[i] = 0;
	C[256] = 0;
	str_buf.reset();
	if( 0 == str_buf.int_vector_size )
		return;
	for(size_type i=0, r_sum=0, r = str_buf.load_next_block(); i < len; ){
		for(; i < r_sum+r; ++i){
			++C[str_buf[i-r_sum]];
		}
		r_sum += r; r = str_buf.load_next_block();
	}
	assert(1 == C[0]);
	size_type value = 0;
	for(uint16_t i=0; i<256; ++i)
		if(C[i]){
			char2comp[i] 	= value;
			C[value]		= C[i];
			if( i > value )
				C[i]=0;
			++value;	
		}
	sigma = value;
	for(uint16_t i=0; i<256; ++i)
		if(char2comp[i])
			comp2char[char2comp[i]] = i;
	for(uint16_t i=256;i>0; --i) C[i] = C[i-1];
	C[0] = 0;
	for(uint16_t i=1; i<257; ++i) C[i] += C[i-1];	
	if(C[256]!=len){
		std::cerr<<"C[256]="<<C[256]<<" "<<len<<std::endl;
	}
	assert(C[256]==len);
	assert(C[sigma+1]==len);
#ifdef SDSL_DEBUG_ALGORITHMS_FOR_COMPRESSED_SUFFIX_ARRAYS
	sw.stop();
	std::cerr<<"set_text takes "<<sw.get_real_time()<<" ms real time";
	std::cerr<<"set_text takes "<<sw.get_user_time()<<" ms user time";
#endif	
}

//! Calculates if a character c occures in the text of the compressed suffix array c.
/*!
 * \param csa 	The csa in which we search for the occurence of c.
 * \param c		The character c for which we search in the text of the compressed suffix array.
 * \return 		True, if c occures in the text of the compressed suffix array and False otherwise.
 * \pre \f$ csa.size()>0 \f$ 
 */ 
template<class Csa>
bool char_occures_in_text_of_csa(const Csa &csa, typename Csa::char_type c){
	return (csa.char2comp[c] > 0) or (csa.char2comp[c]==c); 
}

template<class Csa, uint8_t int_width, class size_type_class>
void set_sa_and_isa_samples(int_vector_file_buffer<int_width, size_type_class> &sa_buf, typename Csa::sa_sample_type &sa_sample, typename Csa::isa_sample_type &isa_sample){
	typedef typename Csa::size_type size_type;
	size_type  n = sa_buf.int_vector_size;

	sa_sample.set_int_width( bit_magic::l1BP(n)+1 );
	sa_sample.resize( (n+Csa::sa_sample_dens-1)/Csa::sa_sample_dens );

	isa_sample.set_int_width( bit_magic::l1BP(n)+1 );
	if( n >= 1 ){ // so n+Csa::isa_sample_dens >= 2
		isa_sample.resize( (n-1+Csa::isa_sample_dens-1)/Csa::isa_sample_dens + 1);
	}

	util::set_one_bits(sa_sample);
	util::set_one_bits(isa_sample);

	sa_buf.reset();
	for(size_type i=0, r_sum = 0, r = sa_buf.load_next_block(), cnt_mod=Csa::sa_sample_dens, cnt_sum=0; r_sum < n; ){
		for(; i < r_sum+r; ++i, ++cnt_mod){
			size_type sa = sa_buf[i-r_sum];
			if( (sa % Csa::isa_sample_dens) == 0 ){
				isa_sample[sa/Csa::isa_sample_dens] = i;
			}else if( sa+1 == n ){
				isa_sample[(sa+Csa::isa_sample_dens-1)/Csa::isa_sample_dens] = i;
			}
			if( Csa::sa_sample_dens == cnt_mod  ){
				cnt_mod = 0;
				sa_sample[cnt_sum++] = sa;
			}
		}
		r_sum += r; r = sa_buf.load_next_block();
	}
/*	if(isa_sample.size() < 20 ){
		std::cerr<<"isa_samples = ";
		for(int_vector<>::size_type i=0;i<isa_sample.size(); ++i)
			std::cerr<<isa_sample[i]<<" ";
		std::cerr<<std::endl;
		std::cerr<<"Csa::isa_sample_dens = "<<Csa::isa_sample_dens <<std::endl;
		std::cerr<<"n = "<<n <<std::endl;
		std::cerr<<"isa_sample.size() = "<< isa_sample.size() <<std::endl;
	}
*/	
}

//template<class Csa>
//bool char_at_char_pos_equals_char(const Csa &csa, typename Csa::size_type char_pos, typename Csa::char_type c){
//	typename Csa::
//}

/*
 * \par Time complexity
 *    \f$ \Order{\log \sigma} \f$
 *  TODO: add hinted binary search? Two way binary search?
*/ 
template<class Csa>	
const unsigned char get_ith_character_of_the_first_row(const typename Csa::size_type i, const Csa &csa){
		assert( i < csa.size() );
		if( csa.sigma < 16 ){ //<- if sigma is small search linear
			typename Csa::size_type res=1;
			while( csa.C[res] <= i )
				++res;
			return csa.comp2char[res-1];
		}
		else
		{	// binary search the character with C 
			typename Csa::size_type upper_c = csa.sigma, lower_c = 0; // lower_c inclusive, upper_c exclusive
			typename Csa::size_type res=0;
//std::cout<<" binary search: csa.sigma= "<<upper_c<<" i="<<i<<std::endl; 
			do{
				res = (upper_c+lower_c)/2;
				if( i < csa.C[res] ){
					upper_c = res;	
				}else if( i >= csa.C[res+1] ){
					lower_c = res+1;
				}
			}
			while(  i < csa.C[res] or i >= csa.C[res+1]  ); // i is not in the interval
			return csa.comp2char[res];
		}
}

}// end namespace algorithm

}// end namespace sdsl

#endif

