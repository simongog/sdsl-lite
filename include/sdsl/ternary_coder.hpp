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
/*! \file ternary_coder.hpp
    \brief ternary_coder.hpp contains the class sdsl::coder::ternary
	\author Simon Gog
 */
#ifndef SDSL_TERNARY_CODER
#define SDSL_TERNARY_CODER

#include "int_vector.hpp"

namespace sdsl{
	
namespace coder{

//! A class to encode and decode between ternary and binary code.
class ternary{
public: 
		typedef uint64_t size_type;

		//! Array contains precomputed values for the decoding of a number in the ternary code.
		/*! The 3 most significant bits (call it hi) contain information about how far 
		 *  to shift to get to the next next encoded 
		 *  integer. If the value of the 13 least significant bits is greater than 2187 or \f$ 2^{12}=4096 \f$
		 *  there is no complete complete decoding of an integer. So shift by 16 bits and decode the next word...
		 *  otherwise shift by 2*(hi+1) to decode the next integer.
		 */ 
		static const uint16_t Trny2bin_0_16[1<<16];
		 
		static const uint16_t Trny2bin_0_16_greedy[1<<16];



		static const uint8_t min_codeword_length = 4; // 11 01 represents 1 and is the code word with minumum length

		static uint8_t encoding_length(uint64_t w);

		template<bool sumup, bool increment, class Iterator>
		static uint64_t decode(const uint64_t *data, const size_type start_idx, size_type n, Iterator it=(Iterator)NULL);

		static uint64_t decode_prefix_sum(const uint64_t *data, const size_type start_idx, size_type n);
		static uint64_t decode_prefix_sum(const uint64_t *data, const size_type start_idx, const size_type end_idx, size_type n);

		template<class int_vector>
		static bool encode(const int_vector &v, int_vector &z);
		template<class int_vector>
		static bool decode(const int_vector &z, int_vector &v);

		//! Encode one positive integer x to an int_vector at bit position start_idx.
		/* \param x Positive integer to encode.
		   \param z Raw data of vector to write the encoded form of x.
		   \param start_idx Beginning bit index to write the encoded form ox x in z.
		*/
		static void encode(uint64_t x, uint64_t *&z, uint8_t &offset);

		template<class int_vector>
		static uint64_t* raw_data(int_vector &v){
			return v.m_data;
		};	
};

inline uint8_t ternary::encoding_length(uint64_t w){
	if(w==0) return 2;
	uint8_t res = 2;   // add two for the terminating "11" bit pair
	do{
		w/=3;
		res+=2;
	}while(w!=0);
	return res;
}

template<class int_vector>
bool ternary::encode(const int_vector &v, int_vector &z){
	z.set_int_width( v.get_int_width() );
	typename int_vector::size_type z_bit_size = 0;
	for(typename int_vector::const_iterator it = v.begin(), end = v.end(); it != end; ++it){
		z_bit_size += encoding_length(*it);
	}
	z.bit_resize( z_bit_size ); // Initial size of z
	if( z_bit_size & 0x3F ){ // if z_bit_size % 64 != 0
		*(z.m_data + (z_bit_size>>6)) = 0; // initialize last word
	}
	z_bit_size = 0;
	uint64_t *z_data = z.m_data, *z_temp;
	uint8_t offset=0, off_temp;
	uint64_t w=0;
	typename int_vector::size_type len=0;
	for(typename int_vector::const_iterator it = v.begin(), end=v.end(); it != end; ++it){
		w = *it;
		if( w == 0 ){
			throw std::logic_error("ternary::encode(const SDSBitVector &v, SDSBitVector &z); entry of v equals 0 that cannot be encoded!");
		}
		// (length of bits to represent w) 
		len = encoding_length(*it)-2;
		if( (offset=(offset + len)) > 63){
			z_data += (offset>>6);
		}
		offset &= 0x3F;
		z_temp = z_data;
		off_temp = offset;
		if( (offset=(offset+2)) > 63 ){
			++z_data;
		}
		offset &= 0x3F;
//		bool not_zero = false;
		// write ternary number
		bit_magic::write_int(z_temp, 3, off_temp, 2);
		for(int32_t i=len; i>0; i-=2){
			if( off_temp != 0 )
				off_temp -= 2;
			else{
				off_temp = 62;
				--z_temp;
			}
//			if( (w%3)>0 )
//				not_zero = true;
			bit_magic::write_int(z_temp, w%3, off_temp, 2);
			w/=3;
		}
//		assert(not_zero);
	}
	return true;
}

template<class int_vector>
bool ternary::decode(const int_vector &z, int_vector &v){
	uint64_t n = 0; // n = number of values to be decoded
	const uint64_t *data = z.data();
	// Determine size of v
	if( z.empty() ){// if z is empty we are ready with decoding
		v.set_int_width(z.get_int_width());
		v.resize(0);
		return true;
	}
	for(typename int_vector::size_type i=0; i < (z.capacity()>>6)-1; ++i, ++data){
		n += bit_magic::eB11Cnt( *data );
	}
	if( z.capacity() != z.bit_size() ){
		n += bit_magic::eB11Cnt( (*data) & bit_magic::Li1Mask[z.bit_size()&0x3F] );
	}
	else{
		n += bit_magic::eB11Cnt( *data );
	}
	v.set_int_width( z.get_int_width() ); v.resize( n );
	return decode<false, true>(z.data(), 0, n, v.begin());
}

template<bool sumup, bool increment, class Iterator>
inline uint64_t ternary::decode(const uint64_t *data, const size_type start_idx, size_type n, Iterator it){
	data += (start_idx >> 6);
	size_type i = 0;
	uint64_t w = 0, value = 0, tempv=0;
	uint32_t v;
	int8_t buffered = 0; // bits buffered in w, in 0..64
	int8_t read = start_idx & 0x3F; // read bits in current *data 0..63
	int8_t shift = 0;
	while( i < n ){// while not all values are decoded
		while(buffered < 16 and i + bit_magic::eB11Cnt(w) < n ){
			w |= (((*data)>>read)<<buffered);
			if( read >= buffered ){
				++data;
				buffered += 64-read;
				read = 0;
			}else{ // read < buffered
				read += 64-buffered;
				buffered = 64;
			}
		}
		v = Trny2bin_0_16[w&0xFFFF];
		//if shift is 16 and the most significant bit pair is not 11
		if( ((v>>13)==7) && (w&0xC000)!=0xC000 ){// not end of decoding
			tempv *= 6561; // tempv *= 3**8
			tempv += v&0x1FFF;
			w >>= 16;
			buffered -= 16;
		}else{// end of decoding
			shift = ((v>>13)+1)<<1;
			w >>= shift;
			buffered -= shift;
			tempv *= bit_magic::powerOf3[v>>13];
			tempv += v&0xFFF;
			if( sumup ) 
				value += tempv;
			else
				value = tempv;
			if( increment ) *(it++) = value;
			tempv = 0;
			++i;
		}
	}
	return value;
}

} // end namespace coder
} // end namespace sdsl 

#endif
