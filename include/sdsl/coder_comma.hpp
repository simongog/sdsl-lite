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
/*! \file coder_comma.hpp
    \brief coder_comma.hpp contains the class sdsl::coder::comma
    \author Uwe Baier
 */
#ifndef SDSL_CODER_COMMA_INCLUDED
#define SDSL_CODER_COMMA_INCLUDED

#include <sdsl/bits.hpp>
#include <array>
#include <math.h>

namespace sdsl {

namespace coder {

//! A class to encode and decode between comma code and binary code.
/*! \author Uwe Baier
 *
 *    Comma coding works as the following:
 *    First of all, comma coding needs a parameter t_width which indicates 
 *    how big a encoded digit will be (in bits), let's say t_width = 2.
 *    By use of t_width one can calculate a base for encoding, in detail 
 *    this means 
 *    base = 2^t_width - 1
 *    now, given any number it is encoded as the follows: The number gets displayed
 *    in the calculated base, and each digit of the number is saved with t_width bits.
 *    To indicate the end of the number, a termination digit is used (namely this is 
 *    the value base).
 *    Example:
 *    t_width = 2 => base = 2^2 - 1 = 3
 *    Value to be encoded: 15
 *    15 (base 10) = 120 (base 3)
 *    Encoded value: 120 (plus termination digit) = 01 10 00 11 in binary 
 *    (last digit is termination digit)
 *
 *  \tparam t_width Width of one digit used in comma code
 */
template<uint8_t t_width = 2>
class comma {
	private:
		static_assert(t_width > 1 && t_width <= 32, 
			"comma coder: Width must be in interval [2,32]");

		//base in which numbers are coded
		static const uint32_t base = (1 << t_width) - 1;

		//table needed for computation of encoding lengths.
		//table contains entries of the kind (index, base^index)
		//to know how much digits a number needs to be encoded.
		static const size_t codelentbllen = ceil(64 / log2(base));
		static std::array<uint64_t, codelentbllen> codelentbl;

		//utility function to set up codelen table
		static std::array<uint64_t, codelentbllen> createCodeLenTbl();

		//helper function to encode a single number without
		//termination digit
		static void encode_in_base(uint64_t x, uint64_t *& z, 
					uint8_t& offset); 
	public:
		typedef uint64_t size_type;
		static const uint8_t min_codeword_length = 
			t_width; //0 needs t_width bits as termination

		//// ENCODING /////////////////////////////////////////////////

		//! Get the number of bits that are necessary to encode 
		//  the value w in comma code.
		/*! \param w 64bit int to get the length of its comma encoding.
		*/
		static uint8_t encoding_length(uint64_t w);

		//! Encode one positive integer x to an int_vector 
		//  at bit position start_idx.
		/* \param x Positive integer to encode.
		   \param z Raw data of vector to write the encoded form of x.
		   \param start_idx Beginning bit index to write the encoded form ox x in z.
		*/
		static void encode(uint64_t x, uint64_t*& z, uint8_t& offset);

		//! Encode integers contained in vector v into vector z
		/* \param v vector containing positive integer values
		   \param z vector to put the encoded values
		*/
		template<class int_vector>
		static bool encode(const int_vector& v, int_vector& z);

		//// DECODING /////////////////////////////////////////////////

		//! Decode n comma encoded values beginning at start_idx 
		//  in the bitstring "data"
		/* \param data Bitstring
		   \param start_idx Starting index of the decoding.
		   \param n Number of values to decode from the bitstring.
		   \param it Iterator to store the values.
		*/
		template<bool t_sumup, bool t_inc,class t_iter>
		static uint64_t decode(const uint64_t* data, 
			const size_type start_idx, size_type n, 
			t_iter it=(t_iter)nullptr);

		//! Decode n comma gamma encoded integers 
		//  beginning at start_idx in the bitstring "data" 
		//  and return the sum of these values.
		/*! \param data Pointer to the beginning 
			of the comma encoded bitstring.
		    \param start_idx Index of the first bit 
			to encode the values from.
		    \param n Number of values to decode from the bitstring.
			Attention: There have to be at least n encoded 
			values in the bitstring.
		*/
		static uint64_t decode_prefix_sum(const uint64_t* data, 
			const size_type start_idx, size_type n);

		//! Decode n comma gamma encoded integers 
		//  beginning at start_idx ending at end_idx (exclusive)
		//  in the bitstring "data" 
		//  and return the sum of these values.
		/*! \param data Pointer to the beginning 
			of the comma encoded bitstring.
		    \param start_idx Index of the first bit 
			to encode the values from.
		    \param end_idx Index of the last bit
			to encode the values from.
		    \param n Number of values to decode from the bitstring.
			Attention: There have to be at least n encoded 
			values in the bitstring.
		*/
		static uint64_t decode_prefix_sum(const uint64_t* data, 
			const size_type start_idx, const size_type end_idx, 
			size_type n);

		//! Decode vector z containing comma encoded integers
		//  and store them in vector v.
		/*! \param z vector that contains encoded integers.
		    \param v vector to store the decoded integers
		*/
		template<class int_vector>
		static bool decode(const int_vector& z, int_vector& v);

		//interface needs this function for whatever :>
		template<class int_vector>
		static uint64_t* raw_data(int_vector& v) {
			return v.m_data;
		}
};

//// IMPLEMENTATION ///////////////////////////////////////////////////////////

//// CODELENGTH TABLE SETUP ///////////////////////////////

template<uint8_t t_width>
std::array<uint64_t, comma<t_width>::codelentbllen> comma<t_width>::codelentbl = 
		createCodeLenTbl();

template<uint8_t t_width>
std::array<uint64_t, comma<t_width>::codelentbllen> comma<t_width>::createCodeLenTbl() {
	std::array<uint64_t, codelentbllen> tbl;
	uint64_t n = 1;
	for (size_t i = 0; i < codelentbllen; i++) {
		tbl[i] = n;
		n = (n << t_width) - n; //n = n * base
	}
	return tbl;
}

//// Encoding /////////////////////////////////////////////

template<uint8_t t_width>
inline uint8_t comma<t_width>::encoding_length(uint64_t w) {
	//use function table and binary search to determine the number of digits
	//needed to encode w in given base.
	uint8_t numdigits = 
		std::upper_bound(codelentbl.begin(), codelentbl.end(), w)
		- codelentbl.begin();
	//finally calculate length. 
	//Don't forget termination character on calculations ;)
	return (numdigits + 1) * t_width;
}

template<uint8_t t_width>
void comma<t_width>::encode_in_base(uint64_t x, uint64_t *& z, 
					uint8_t& offset) {
	if (x) {
		uint32_t digit = x % base; //get next digit
		//encode digits with higher order
		encode_in_base(x / base, z, offset);
		//and write own digit
		bits::write_int_and_move(z, digit, offset, t_width);
	}
}

template<uint8_t t_width>
inline void comma<t_width>::encode(uint64_t x, uint64_t*& z, uint8_t& offset) {
	//encode x itself
	encode_in_base(x, z, offset);
	//and append the termination digit	
	bits::write_int_and_move(z, base, offset, t_width);
}

template<uint8_t t_width>
template<class int_vector>
bool comma<t_width>::encode(const int_vector& v, int_vector& z) {
	//first, find out how much bits vector z needs to save values
	typedef typename int_vector::size_type size_type;
	size_type z_bit_size = 0;
	for (typename int_vector::const_iterator it = v.begin(), end = v.end();
	     it != end; ++it) {
		z_bit_size += encoding_length(*it);
	}

	//trim vector z to correct size
	z.width(v.width());
	z.bit_resize(z_bit_size); //for future may check if resizing works

	//iterate again and save values in z
	uint64_t* z_data = z.m_data;
	uint8_t offset = 0;
	for (typename int_vector::const_iterator it = v.begin(), end = v.end();
	     it != end; ++it) {
		encode(*it, z_data, offset);
	}
	return true;
}

//// DECODING /////////////////////////////////////////////

template<uint8_t t_width>
template<bool t_sumup, bool t_inc, class t_iter>
inline uint64_t comma<t_width>::decode(const uint64_t* data, 
		const size_type start_idx, size_type n, t_iter it)	{
	data += (start_idx >> 6); //jump to byte offset
	uint8_t offset = start_idx & 0x3F; //and calculate bit offset
	uint64_t value = 0;
	for (size_type i = 0; i < n; i++) {
		//read next value
		uint64_t v = 0;
		for (uint32_t digit = (uint32_t)bits::read_int_and_move(data, offset, t_width); //read first digit
		     digit != base;                  //while digit is not the terminating digit
                     v = (v << t_width) - v + digit, //v = v * base + digit
 		     digit = (uint32_t)bits::read_int_and_move(data, offset, t_width)); //and read next digit
		//now decide how to handle value
		value = (t_sumup) ? value + v : v;
		if (t_inc)	*(it++) = value;
	}
	return value;
}

template<uint8_t t_width>
uint64_t comma<t_width>::decode_prefix_sum(const uint64_t* data, 
			const size_type start_idx, size_type n)	{
	//easiest seems to be to use already build function decode...
	return decode<true,false,int *>(data, start_idx, n);
	//Note for above: 3rd template parameter ca be any pntr except void *
}

template<uint8_t t_width>
uint64_t comma<t_width>::decode_prefix_sum(const uint64_t* data, 
		const size_type start_idx, 
		SDSL_UNUSED const size_type end_idx, size_type n)	{
	//end index does not change anything here...
	return decode_prefix_sum(data, start_idx, n);
}

template<uint8_t t_width>
template<class int_vector>
bool comma<t_width>::decode(const int_vector& z, int_vector& v)	{
	//check if bit size is dividable through t_width.
	if (z.bit_size() % t_width != 0)	return false;

	//calculate num of overall digits in z (including terminating digit)
	uint64_t numOfDigits = z.bit_size() / t_width;
	//iteration vars for z vector
	const uint64_t *z_data = z.data();
	uint8_t z_offset = 0;
	//utility to count number of entries in z, and last read digit
	uint32_t digit = base;
	typename int_vector::size_type n = 0;	

	//iterate over all digits. each time a termination digit is 
	// detected, a encoded number in vector ends.
	while (numOfDigits--) {
		digit = (uint32_t)bits::read_int_and_move(z_data, z_offset, t_width);
		if (digit == base)	n++;	//termination digit detected
	}

	//also, ensure last read digit was a termination digit
	if (digit != base)	return false;

	//resize vector v
	v.width(z.width());
	v.resize(n);

	//and finally decode and save result in v
	decode<false, true>(z.data(), 0, n, v.begin());
	return true;
}

} //end of namespace coder
} //end of namespace sdsl

#endif
