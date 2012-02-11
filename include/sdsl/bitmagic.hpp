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
/*! \file bitmagic.hpp
    \brief bitmagic.hpp contains the sdsl::bit_magic class. 
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_BITMAGIC
#define INCLUDED_SDSL_BITMAGIC

#include <stdint.h> // for uint64_t uint32_t declaration
#include <iostream>// for cerr
#include <cassert>

//! Namespace for the succinct data structure library.
namespace sdsl{

//! A helper class for bitwise tricks on 64 bit words.
/*!
	bit_magic is a helper class for bitwise tricks and
	techniques. For the basic tricks and techiques we refer to Donald E. Knuth's
	"The Art of Computer Programming", Volume 4A, Chapter 7.1.3 and
	the informative website of Sean E. Anderson about the topic:
	http://www-graphics.stanford.edu/~seander/bithacks.html .

	We have added new functions like: b11Cnt and i11BP.

	All members of this class are static variables or methods.
	This class cannot be instantiated.

	\author Simon Gog
 */
class bit_magic{
	private:
		bit_magic(); // This helper class can not be instantiated
	public:
		//! 64bit mask with all bits set to 1.
		static const int64_t  All1Mask = -1LL;

		//! This constant represents a de Bruijn sequence B(k,n) for k=2 and n=6.
		/*! Details for de Bruijn sequences see
           http://en.wikipedia.org/wiki/De_bruijn_sequence
		   DeBruijn64 is used in combination with the
		   array DeBruijn64ToIndex.
		*/
		static const uint64_t DeBruijn64 = 0x0218A392CD3D5DBFULL;

		//! This table maps a 6-bit subsequence S[idx...idx+5] of constant DeBruijn64 to idx.
		/*! \sa DeBruijn64
		*/
		static const uint32_t DeBruijn64ToIndex[64];

		//! Array containing fibonacci numbers less than \f$2^64\f$.
		static const uint64_t Fib[92];

		//! Array containing precomputed values for the leftmost 1-bit in a 8bit integer.
	    /*! \sa l1BP
	     */
		static const uint32_t L1BP[256];

		//! An array with entry i containing a 64bit integer with the last i bits set.
		static const uint64_t Li1Mask[65];

		//! An array with entry i containing a 64bit integer with the last i bits not set.
		static const uint64_t Li0Mask[65];

		static const uint8_t lookupr1BP[256];

		//! An array with precomputed select queries on a 8-bit integer.
		/*! Entry at idx = 256*j + i equals the position of the
		    (j+1)-th leftmost 1 bit in the integer i. Positions lie in the range \f$[0..7]\f$.
		 */
		static const uint8_t Select256[256*8];

		static const uint64_t PsOverflow[65];

		static const uint32_t powerOf3[9];

		static const uint8_t max_excess_8bit[256]; 

		//! Contains 64 bit constants. Constant i is equal to 8 bytes each set to i. 
		static const uint64_t _8_x_the_byte[65];

		static const uint8_t cover0[1];
		static const uint8_t cover1[2];
		static const uint8_t cover2[3];
		static const uint8_t cover3[4];
		static const uint8_t cover4[5];
		static const uint8_t cover5[7];
		static const uint8_t cover6[9];
		static const uint8_t cover7[13];
		static const uint8_t cover8[20];

		static const uint8_t *covers[9];

		static const uint32_t cover_sizes[9];


		static void generate_first_pos_of_excess_val(){
			for(int i=0;i<9;++i){
				for(int j=0;j<256;++j){
					int x = j;
					int excess = 0;
					int found=0;
					for(int k=0;k<8 and !found;++k, x>>=1){
						if(x&1)
							excess++;
						else
							excess--;
						if(excess == i){
							std::cout<<k<<",";
							found=1;
						}
					}
					if(!found)
						std::cout<<0<<",";
					if((j+1)%16==0)
						std::cout<<std::endl;
				}
				std::cout<<std::endl;	
			}
		}

		static const uint8_t first_pos_of_excess_val[256*9];

		static void generate_last_pos_of_excess_val(){
			std::cout<<"const uint8_t bit_magic::last_pos_of_excess_val[]={"<<std::endl;
			for(int i=-8;i<9;++i){
				for(int j=0;j<256;++j){
					int x = j;
					int excess = 0;
					int idx=0;
					for(int k=0;k<8;++k, x>>=1){
						if(x&1)
							excess++;
						else
							excess--;
						if(excess == i){
							idx=k;
						}
					}
					std::cout<<idx;
					if(i*j!=8*255)
						std::cout<<",";
					if((j+1)%16==0)
						std::cout<<std::endl;
				}
				std::cout<<std::endl;	
			}
			std::cout<<"};"<<std::endl;
		}

		static const uint8_t last_pos_of_excess_val[256*17];

		static void generate_very_near_find_open(){
			std::cout<<"const uint8_t bit_magic::very_near_find_open[]={"<<std::endl;
			for(int w=0;w<256;++w){
				int excess_val=0;
				for(int i=0;i<8;++i){
					if( (w>>i)&1 )
						++excess_val;
					else
						--excess_val;
				}
				int res=0;
				int cur_excess_val=0;
				for(int i=0;i<8;++i){
					if( (w>>i)&1 )
						++cur_excess_val;
					else
						--cur_excess_val;
					if(cur_excess_val==excess_val-1)
						res = 7-i;
				}
				std::cout<<res;
				if(w!=255)
					std::cout<<",";
				if( (w+1)%16==0 )
					std::cout<<std::endl;
			}
			std::cout<<"};"<<std::endl;
		}

		//! Precomputed answers for find_open for blocks of 8 bit preceding the closing parenthesis.
		/*! If there is no matching opening parenthesis within these 8 bits 0 is returned otherwise
		 *  the distance i \f$i\in \{1,3,5,7\}\f$ from the index of the closing parenthesis to the 
		 *  matching opening parenthesis is returned. 
		 */
		static const uint8_t very_near_find_open[256];

		static void generate_very_near_enclose(){
			std::cout<<"const uint8_t bit_magic::very_near_enclose[]={"<<std::endl;
			for(int w=0;w<256;++w){
				int excess_val = 0, res=0;
				for(int i=7;i>=0 and !res;--i){
					if( w&(1<<i) )
						++excess_val;
					else
						--excess_val;
					if(excess_val==1)
						res = 8-i;
				}
				std::cout<<res;
				if(w!=255)
					std::cout<<",";
				if( (w+1)%16==0 )
					std::cout<<std::endl;
			}
			std::cout<<"};"<<std::endl;
		}

		static const uint8_t very_near_enclose[256];
		//! Counts the number of 1-bits in the 64bit integer x.
		/*! This function is also known as sideway addition in literature.
		    E.g. see Donald E. Knuth, "The Art of Computer Programming", Volume 4A.
		    \param x 64bit integer to count the bits.
		    \return The number of 1-bits in x.
		 */
		static uint64_t b1Cnt(uint64_t x);

		//! Counts the number of 1-bits in the 32bit integer x.
		/*! This function is a variant of the method b1Cnt. If
			32bit multiplication is fast, this method beats the b1Cnt.
			for 32bit integers.
			\param x 64bit integer to count the bits.
			\return The number of 1-bits in x.
		 */
		static uint32_t b1Cnt32(uint32_t x);

		//! Naive implementation of the b1Cnt function.
		static uint32_t b1CntNaive(uint64_t x);

		//! Count the number of consecutive and distinct 11 in the 64bit integer x.
		/*!
		  	\param x 64bit integer to count the terminating sequence 11 of a fibonacci code.
			\param c Carry equals msb of the previous 64bit integer.
		 */
		static uint32_t b11Cnt(uint64_t x, uint64_t &c);

		static uint32_t b11CntS(uint64_t x, uint64_t &c);
		static uint32_t b11CntS(uint64_t x);

		//! Count the number of consecutive and distinct 11 in the 64bit integer x.
		/*!
		  	\param x 64bit integer to count the terminating sequence 11 of a fibonacci code.
		 */
		static uint32_t b11Cnt(uint64_t x);

		//! Naive implementation of the b11Cnt function.
		static uint32_t b11CntNaive(uint64_t x);

		//! Count 11 bit pairs starting at even positions in the 64bit integer x.
		static uint32_t eB11Cnt(uint64_t x);

		//! Count 10 bit pairs in the word x.
		/*!
		 * \param x 64bit integer to count the 10 bit pairs.
		 * \param c Carry equals msb of the previous 64bit integer.
		 */
		static uint32_t b10Cnt(uint64_t x, uint64_t &c);

		//! Count 01 bit pairs in the word x.
		/*!
		 * \param x 64bit integer to count the 01 bit pairs.
		 * \param c Carry equals msb of the previous 64bit integer.
		 */
		static uint32_t b01Cnt(uint64_t x, uint64_t &c);

		static uint32_t b10CntNaive(uint64_t x, uint64_t &c);

		//! Map all 10 bit pairs to 01 or 1 if c=1 and the lsb=0. All other pairs are mapped to 00.
		static uint64_t b10Map(uint64_t x, uint64_t c=0);

		//! Map all 01 bit pairs to 01 of 1 if c=1 and the lsb=0. All other pairs are mapped to 00.
		static uint64_t b01Map(uint64_t x, uint64_t c=1);

		//! Calculate the position of the i-th rightmost 1 bit in the 64bit integer x
		/*!
		  	\param x 64bit integer.
			\param i Argument i must be in the range \f$[1..b1Cnt(x)]\f$.
			\pre Argument i must be in the range \f$[1..b1Cnt(x)]\f$.
		  	\sa l1BP, r1BP
		 */
		static uint32_t i1BP(uint64_t x, uint32_t i);

		//! i1BP implementation proposed by Vigna.
		/*! \sa i1BP
		 */
		static uint32_t j1BP(uint64_t x, uint32_t j);

		//! Naive implementation of i1BP.
		static uint32_t i1BPNaive(uint64_t x, uint32_t i);

		//! Calculates the position of the leftmost 1-bit in the 64bit integer x if it exists
		/*! \param x 64 bit integer.
		    \return The position (in 0..63) of the leftmost 1-bit in the 64bit integer x if
			        x>0 and 0 if x equals 0.
			\sa i1BP, r1BP
		*/
		static const uint32_t l1BP(uint64_t x);

		//! Naive implementation of l1BP.
		static const uint32_t l1BPNaive(uint64_t x);

		//! Calculates the position of the rightmost 1-bit in the 64bit integer x if it exists
		/*! This method is e.g. used in the getUnaryBit method.
			\param x 64 bit integer.
			\return The position (in 0..63) of the rightmost 1-bit in the 64bit integer x if
			        x>0 and 0 if x equals 0.
			\sa i1BP, l1BP
		*/
		static const uint32_t r1BP(uint64_t x);

		//! Naive implementation of r1BP.
		static const uint32_t r1BPNaive(uint64_t x);

		//! Calcluates the position of the i-th rightmost 11-bit-pattern which terminates a fibonacci coded integer in x.
		/*!	\param x 64 bit integer.
		    \param i Index of 11-bit-pattern. \f$i \in [1..b11Cnt(x)]\f$
			\param c Carry bit from word before
		 	\return The position (in 1..63) of the i-th 11-bit-pattern which terminates a fibonacci coded integer in x if
			        x contains at least i 11-bit-patterns and a undefined value otherwise.
		    \sa b11Cnt, l11BP, i1BP

		 */
		static const uint32_t i11BP(uint64_t x, uint32_t i, uint32_t c=0);

		//! Naive implementation of i11BP
		static const uint32_t i11BPNaive(uint64_t x, uint32_t i){
			for(uint32_t j=0; j<63; ++j){
				if( (x&3)==3  ){
					i--;
					if(!i) return j+1;
					x>>=1;
					++j;
				}	
				x>>=1;
			}
			return 63;
		}

		static const uint64_t all11BPs(uint64_t x, bool &c);

		//! Calculates the position of the i-th rightmost 11-bit-pattern which terminates a Ternary coded integer in x.
		static const uint32_t eI11BP(uint64_t x, uint32_t i);

		//! Calculates the position of the leftmost 11-bit-pattern which terminates a fibonacci coded integer in x.
	    /*! \param x 64 bit integer.
		    \return The position (in 1..63) of the leftmost 1 of the leftmost 11-bit-pattern which
			        terminates a fibonacci coded integer in x if x contains a 11-bit-pattern
					and 0 otherwise.
			\sa b11Cnt, i11BP
	    */
		static const uint32_t l11BP(uint64_t x);

		//! TODO: Documentation for this function
		static void writeInt(uint64_t *word, uint64_t x, const uint8_t offset=0, const uint8_t len=64);
		static void writeInt2(uint64_t *word, uint64_t x, const uint8_t offset=0, const uint8_t len=64);

		//! TODO: Documentation for this function
		static void writeIntAndMove(uint64_t* &word, uint64_t x, uint8_t &offset, const uint8_t len);

		//! TODO: Documentation for this function
		static uint64_t readInt(const uint64_t *word, uint8_t offset=0, const uint8_t len=64);

		//! TODO: Documentation for this function
		static uint64_t readIntAndMove(const uint64_t* &word, uint8_t &offset, const uint8_t len=64);

		//! TODO: Documentation for this function
		static uint64_t readUnaryInt(const uint64_t* word, uint8_t offset=0);

		//! TODO: Documentation for this function
		static uint64_t readUnaryIntAndMove(const uint64_t* &word, uint8_t &offset);

		//! Move the bit pointer consisting of a uint64_t pointer and an offset len positions to the right.
		/*!
		 * \param word   64-bit word part of the bit pointer
		 * \param offset Offset part of the bit pointer
		 * \param len Position to move to the right. \f$ len \in [0..64] \f$
		 * \sa move_left
		 */
		static void move_right(const uint64_t* &word, uint8_t &offset, const uint8_t len);
		//! Move the bit pointer consisting of a uint64_t pointer and an offset len positions to the left.
		/*!
		 * \param word   64-bit word part of the bit pointer
		 * \param offset Offset part of the bit pointer
		 * \param len Position to move to the left. \f$ len \in [0..64] \f$
		 * \sa move_right
		 */
		static void move_left(const uint64_t* &word, uint8_t &offset, const uint8_t len);

		//! Get the first one bit in the intervall \f$[idx..\infty )\f$
		static uint64_t next(const uint64_t* word, uint64_t idx);

		//! Get the one bit with the greatest position in the intervall \f$[0..idx]\f$
		static uint64_t prev(const uint64_t* word, uint64_t idx);

/*		// precalc the deBruijn64ToIndex table
		static const void calcDeBruijnHash(uint64_t h[]){
			h[0] = 0;
			for(int i=1; i<64; ++i){
				h[((1ULL<<i)*DeBruijn64)>>58] = i;
			}
		}
*/
/*		// precalc the lookup for select256
		static const void select256_table(uint64_t table[]){
			for(int ones=1;ones<=8;++ones){
				for(uint16_t i=0;i<256;++i){
					uint16_t tmp = i, res=0, one_occs=0;
					while(tmp){
						if(tmp&1)
							if( ++one_occs==ones)
								break;
						tmp>>=1;
						++res;
					}
					if(one_occs!=ones)
						res = 0;
					table[((ones-1)<<8) + i] = res;
				}
			}
		}
*/
		// precalac for PSoverflow
/*		static const void psoverflow_table(uint64_t table[]){
			for(int i=0; i<=64; ++i){
				table[i]= 128-i;
				for(int j=1;j<8;++j){
					table[i]<<=8;
					table[i] += (128-i);
				}
			}
		}
*/
		static uint16_t max_excess(uint64_t x, uint16_t &b1Cnt);
		static uint16_t max_excess2(uint64_t x, uint16_t &b1Cnt);
		static uint16_t max_excess3(uint64_t x, uint16_t &b1Cnt);
		/*! Calculates the maximal excess values for the eight byte prefixes of a 64 bit word.
		 *  \param w The 64 bit word representing a parentheses sequence.
		 *  \param max_byte_excesses The \e i th byte contains the maximal excess value that  
		 *                           ..... TODO: description
		 *  \param byte_prefix_sums_x_2 64bit word that contains the doubled prefix sums of the 8 bytes.
		 */
		static void max_byte_excesses(uint64_t w, uint64_t &max_byte_excesses, uint64_t &byte_prefix_sums_x_2);
		static void max_byte_excesses2(uint64_t w, uint64_t &max_byte_excesses, uint64_t &byte_prefix_sums_x_2);

		static void min_max_byte_excesses(uint64_t max2, uint64_t &min_byte_excesses, uint64_t &max_byte_excesses, uint64_t &byte_prefix_sums_x_2);

		/*! Calculates the minimal position where a specific excess_val is reached.
		 * \param w The 64 bit word representing a parentheses sequence.
		 * \param excess_val The excess value which should be reached. \f$0<excess\_val \leq 64 \f$. 
		 * \param byte_prefix_sums_x_2 Contains for each byte i of w the number of ones the the byte + the
		 *                             number of ones in the previous bytes in w. 
		 * \return The minimal index where a the excess value equals excess_val or 64 if there is
		 *         no such index.
		 * \sa first_excess_position_naive
		 */ 
		static uint8_t first_excess_position(uint64_t w, uint8_t excess_val, uint64_t &byte_prefix_sums_x_2);

		/*! Naive implementation of first_excess_position
		 * \sa first_excess_position 
		 */
		static uint8_t first_excess_position_naive(uint64_t w, uint8_t excess_val, uint64_t &byte_prefix_sums_x_2);

		static uint8_t last_excess_position(uint64_t w, uint8_t excess_val, uint64_t &byte_prefix_sums_x_2);
		static uint8_t last_excess_position_naive(uint64_t w, uint8_t excess_val, uint64_t &byte_prefix_sums_x_2);
		/*
		 * \param close_parenthesis_index Index of the closing parenthesis. \f$0 \leq close_parenthesis_index \leq 64\f$
		 */
		static uint8_t find_open(uint64_t w, uint8_t close_parenthesis_index);
		static uint8_t find_open_naive(uint64_t w, uint8_t close_parenthesis_index);

		static uint8_t find_enclose(uint64_t w, uint8_t open_parenthesis_index);
		static uint8_t find_enclose_naive(uint64_t w, uint8_t open_parenthesis_index);
};


// ============= inline - implementations ================

inline uint16_t bit_magic::max_excess(uint64_t max2, uint16_t &b1Cnt){
	uint64_t sum = max2-((max2>>1)&0x5555555555555555ULL);
	uint64_t max = sum&max2;//(~(x&0x5555555555555555ULL&sum))&sum;
	// Maxima fuer 4bit Bloecke
	max2 = ((max>>2)&0x3333333333333333ULL) + 0x6666666666666666ULL + ((sum&0x3333333333333333ULL)<<1);
//	cout<<"max2 ";
//	printDebug(max2);
	uint64_t mask = (max2 - (max = max&0x3333333333333333ULL))&0x8888888888888888ULL;
	mask = (mask | (((mask>>3))^0x1111111111111111ULL))-0x1111111111111111ULL;
//	max2 = ((max2 & mask) | (max & ~mask));
	max2 = max ^ ((max^max2) & mask);
	// Summen fuer 4bit Bloecke
	sum = ((sum & 0x3333333333333333ULL) + ((sum>>2) & 0x3333333333333333ULL))<<1;
	// Maxima fuer 8bit Bloecke
   	max  = ((max2>>4)&0x0707070707070707ULL) + 0x0C0C0C0C0C0C0C0CULL + ((sum&0x0F0F0F0F0F0F0F0FULL));	
	mask = (max - (max2 = max2&0x0707070707070707ULL)) & 0x1010101010101010ULL;
	mask = (mask | ((mask>>4)^0x0101010101010101ULL))-0x0101010101010101ULL;
//	max  = ((max & mask) | (max2 & ~mask));
	max  = max2 ^ ((max2^max) & mask);  
	// Summe fuer 8 bit Bloecke
	sum = ((sum + (sum>>4)) & 0x1E1E1E1E1E1E1E1EULL);
	// Maxima fuer 16bit Bloecke
	max2 = ((max>>8)&0x000F000F000F000FULL) + 0x0018001800180018ULL + ((sum&0x00FF00FF00FF00FFULL));
	mask = (max2 - (max = max&0x000F000F000F000FULL)) & 0x0020002000200020ULL;
	mask = (mask | ((mask>>5)^0x0001000100010001ULL))-0x0001000100010001ULL;
//	max2 = ((max2 & mask) | (max & ~mask));
	max2 = max ^ ((max^max2) & mask);
	// Summe fuer 16 bit Bloecke
	sum = (sum + (sum>>8)) & 0x00FF00FF00FF00FFULL;
	// Maxima fuer 32bit Bloecke
	max = ((max2>>16)&0x0000001F0000001FULL) + 0x0000003000000030ULL + ((sum&0x0000FFFF0000FFFFULL));
	mask = (max - (max2 = max2&0x0000001F0000001FULL)) & 0x0000004000000040ULL;
	mask = (mask | ((mask>>6)^0x0000000100000001ULL))-0x0000000100000001ULL;
//	max = ((max & mask) | (max2 & ~mask)) & 0x0000003F0000003FULL;
	max  = (max2 ^ ((max2^max) & mask)) & 0x0000003F0000003FULL;  
	// Summe fuer 32 bit Bloecke
	sum = (sum + (sum>>16)) & 0x0000FFFF0000FFFFULL;
	b1Cnt = (sum+(sum>>32))>>1;
	// Maxima fuer 64bit Bloecke
	max2 = (max>>32) + 0x0000000000000060ULL + (sum&0x00000000FFFFFFFFULL);
	return ((max2 - (max&0x00000000FFFFFFFFULL)) & 0x0000000000000080ULL) ? max2 & 0x000000000000007FULL : max & 0x00000000FFFFFFFFULL;
/*	mask = (max2 - (max = max&0x00000000FFFFFFFFULL)) & 0x0000000000000080ULL;
	return max2;*/
}

inline uint16_t bit_magic::max_excess2(uint64_t x, uint16_t &b1Cnt){
	uint64_t _max_byte_excesses, byte_prefix_sums_x_2;	
	max_byte_excesses(x, _max_byte_excesses, byte_prefix_sums_x_2);
	b1Cnt = byte_prefix_sums_x_2 >> 57;
	uint8_t maxi1 = _max_byte_excesses, maxi2 = _max_byte_excesses>>32;
	for(uint8_t i=1, m;i<4;++i){
		if( (m = (_max_byte_excesses >>= 8)) > maxi1 ) maxi1 = m;
	}
	_max_byte_excesses >>= 8;
	for(uint8_t i=5, m;i<8;++i){
		if( (m = (_max_byte_excesses >>= 8)) > maxi2 ) maxi2 = m;
	}
	return (maxi1>maxi2)?maxi1&0x7F:maxi2&0x7F;

/*	
	uint8_t maxi[8] = {sum, sum>>8, sum>>16, sum>>24, sum>>32, sum>>40, sum>>48,sum>>56};
	if(maxi[0]<maxi[1])maxi[0]=maxi[1];
	if(maxi[2]<maxi[3])maxi[2]=maxi[3];
	if(maxi[4]<maxi[5])maxi[4]=maxi[5];
	if(maxi[6]<maxi[7])maxi[6]=maxi[7];
	if(maxi[0]<maxi[2])maxi[0]=maxi[2];
	if(maxi[4]<maxi[6])maxi[4]=maxi[6];
	return (maxi[0]<maxi[4])?maxi[4]&0x7F:maxi[0]&0x7F;
*/	
	
/*	max = (sum&0x00FF00FF00FF00FFULL)|0x0100010001000100ULL;
	mask = max-(max2 = ((sum>>4)&0x00FF00FF00FF00FFULL));	
	mask = (mask | ((mask>>5)^0x0001000100010001ULL))-0x0001000100010001ULL;
	max  = max2 ^ (max2^max) & mask;  
	uint8_t maxi1 = max, maxi2 = max>>32;
	if( ((uint8_t)(max>>16)) > maxi1 ) maxi1 = max>>16;
	if( ((uint8_t)(max>>48)) > maxi2 ) maxi2 = max>>48;
	return maxi2 > maxi1 ? maxi2 : maxi1;
*/	
	
/*	for(uint8_t i=1, m;i<4;++i){
		if( (m = (max >>= 16)) > maxi ) maxi = m;
	}
	return maxi&0x7F;
*/	
}

inline void bit_magic::max_byte_excesses(uint64_t max2, uint64_t &max_byte_excesses, uint64_t &byte_prefix_sums_x_2){
	byte_prefix_sums_x_2 = max2-((max2>>1)&0x5555555555555555ULL);
	uint64_t max = byte_prefix_sums_x_2&max2;
	// Maxima fuer 4bit Bloecke
	max2 = ((max>>2)&0x3333333333333333ULL) + 0x6666666666666666ULL + ((byte_prefix_sums_x_2&0x3333333333333333ULL)<<1);
	uint64_t mask = (max2 - (max = max&0x3333333333333333ULL))&0x8888888888888888ULL;
	mask = (mask | ((mask>>3)^0x1111111111111111ULL))-0x1111111111111111ULL;
	max2 = max ^ ((max^max2) & mask);
	// Summen fuer 4bit Bloecke
	byte_prefix_sums_x_2 = ((byte_prefix_sums_x_2 & 0x3333333333333333ULL) + ((byte_prefix_sums_x_2>>2) & 0x3333333333333333ULL))<<1;
	// Maxima fuer 8bit Bloecke
   	max  = ((max2>>4)&0x0707070707070707ULL) + 0x0C0C0C0C0C0C0C0CULL + ((byte_prefix_sums_x_2&0x0F0F0F0F0F0F0F0FULL));	
	mask = (max - (max2 = max2&0x0707070707070707ULL)) & 0x1010101010101010ULL;
	mask = (mask | ((mask>>4)^0x0101010101010101ULL))-0x0101010101010101ULL;
	max  = max2 ^ ((max2^max) & mask);  
	// Summe fuer 8 bit Bloecke
	byte_prefix_sums_x_2 = ((byte_prefix_sums_x_2 + (byte_prefix_sums_x_2>>4)) & 0x1E1E1E1E1E1E1E1EULL);
	byte_prefix_sums_x_2 = 0x0101010101010101ULL*byte_prefix_sums_x_2;// partialsummen der 2*summe der einsen im praefix
	max_byte_excesses = (( (byte_prefix_sums_x_2<<8) | 0x8080808080808080ULL) - 0x3830282018100800ULL) + max; 
}

inline void bit_magic::max_byte_excesses2(uint64_t max2, uint64_t &max_byte_excesses, uint64_t &byte_prefix_sums_x_2){
	byte_prefix_sums_x_2 = max2-((max2>>1)&0x5555555555555555ULL);
//	uint64_t max = (byte_prefix_sums_x_2&max2)+((max2|(max2>>1))&0x5555555555555555ULL);// add one to each max if max != 0
	uint64_t max = ((max2&0x5555555555555555ULL)<<1)|((max2&0x5555555555555555ULL)>>1);
	// Maxima fuer 4bit Bloecke
	max2 = ((max>>2)&0x3333333333333333ULL) + 0x6666666666666666ULL + ((byte_prefix_sums_x_2&0x3333333333333333ULL)<<1);
	uint64_t mask = (max2 - (max = max&0x3333333333333333ULL))&0x8888888888888888ULL;
	mask = (mask | ((mask>>3)^0x1111111111111111ULL))-0x1111111111111111ULL;
	max2 = max ^ ((max^max2) & mask);
	// Summen fuer 4bit Bloecke
	byte_prefix_sums_x_2 = ((byte_prefix_sums_x_2 & 0x3333333333333333ULL) + ((byte_prefix_sums_x_2>>2) & 0x3333333333333333ULL))<<1;
	// Maxima fuer 8bit Bloecke
   	max  = ((max2>>4)&0x0707070707070707ULL) + 0x0C0C0C0C0C0C0C0CULL + ((byte_prefix_sums_x_2&0x0F0F0F0F0F0F0F0FULL));	
	mask = (max - (max2 = max2&0x0707070707070707ULL)) & 0x1010101010101010ULL;
	mask = (mask | ((mask>>4)^0x0101010101010101ULL))-0x0101010101010101ULL;
	max  = max2 ^ ((max2^max) & mask);  
	// Summe fuer 8 bit Bloecke
	byte_prefix_sums_x_2 = ((byte_prefix_sums_x_2 + (byte_prefix_sums_x_2>>4)) & 0x1E1E1E1E1E1E1E1EULL);
	byte_prefix_sums_x_2 = 0x0101010101010101ULL*byte_prefix_sums_x_2;// partialsummen der 2*summe der einsen im praefix
	max_byte_excesses =  (byte_prefix_sums_x_2<<8) + 0x474f575f676F777FULL + max; 
}

inline void bit_magic::min_max_byte_excesses(uint64_t max2, uint64_t &min_byte_excesses, uint64_t &max_byte_excesses, uint64_t &byte_prefix_sums_x_2){
	byte_prefix_sums_x_2 = max2-((max2>>1)&0x5555555555555555ULL);
//	uint64_t max = (byte_prefix_sums_x_2&max2)+((max2|(max2>>1))&0x5555555555555555ULL);// add one to each max if max != 0
	uint64_t max = ((max2&0x5555555555555555ULL)<<1)|((max2>>1)&0x5555555555555555ULL);
	uint64_t min = ~max, min2,mask,mask2;
//std::cerr<<"min = "<<min<<std::endl;	
	// Maxima and minima fuer 4bit Bloecke
	max2 = ((max>>2)&0x3333333333333333ULL) + 0x6666666666666666ULL + ((byte_prefix_sums_x_2&0x3333333333333333ULL)<<1);
	min2 = ((min>>2)&0x3333333333333333ULL) + 0x6666666666666666ULL + (0x4444444444444444ULL-(((byte_prefix_sums_x_2)&0x3333333333333333ULL)<<1));
//	min2 = ((min>>2)&0x3333333333333333ULL) + mask;
	mask = (max2 - (max = max&0x3333333333333333ULL))&0x8888888888888888ULL;
	mask2 = (min2 - (min = min&0x3333333333333333ULL))&0x8888888888888888ULL;
	mask = (mask | ((mask>>3)^0x1111111111111111ULL))-0x1111111111111111ULL;
	mask2 = (mask2 | ((mask2>>3)^0x1111111111111111ULL))-0x1111111111111111ULL;
	max2 = max ^ ((max^max2) & mask);
	min2 = min ^ ((min^min2) & mask2);
//std::cerr<<"min2 = "<<min2<<std::endl;	
	// Minima fuer 4bit Bloecke
	// Summen fuer 4bit Bloecke
	byte_prefix_sums_x_2 = ((byte_prefix_sums_x_2 & 0x3333333333333333ULL) + ((byte_prefix_sums_x_2>>2) & 0x3333333333333333ULL))<<1;
	// Maxima and minima fuer 8bit Bloecke
   	max  = ((max2>>4)&0x0707070707070707ULL) + 0x0C0C0C0C0C0C0C0CULL + ((byte_prefix_sums_x_2&0x0F0F0F0F0F0F0F0FULL));	
   	min  = ((min2>>4)&0x0707070707070707ULL) + 0x0C0C0C0C0C0C0C0CULL + (0x0808080808080808ULL-((byte_prefix_sums_x_2)&0x0F0F0F0F0F0F0F0FULL));	
  // 	min  = ((min2>>4)&0x0707070707070707ULL) + mask;	
	mask = (max - (max2 = max2&0x0707070707070707ULL)) & 0x1010101010101010ULL;
	mask2 = (min - (min2 = min2&0x0707070707070707ULL)) & 0x1010101010101010ULL;
	mask = (mask | ((mask>>4)^0x0101010101010101ULL))-0x0101010101010101ULL;
	mask2 = (mask2 | ((mask2>>4)^0x0101010101010101ULL))-0x0101010101010101ULL;
	max  = max2 ^ ((max2^max) & mask);  
	min  = min2 ^ ((min2^min) & mask2);  
//std::cerr<<"min = "<<min<<std::endl;	
	// Summe fuer 8 bit Bloecke
	byte_prefix_sums_x_2 = ((byte_prefix_sums_x_2 + (byte_prefix_sums_x_2>>4)) & 0x1E1E1E1E1E1E1E1EULL);
	byte_prefix_sums_x_2 = 0x0101010101010101ULL*byte_prefix_sums_x_2;// partialsummen der 2*summe der einsen im praefix
//std::cerr<<"max = "<<max<<std::endl;	
//std::cerr<<"min = "<<min<<std::endl;	
//std::cerr<<"byte_prefix_sums_x = "<<(byte_prefix_sums_x_2) << std::endl; 
	max_byte_excesses =  (byte_prefix_sums_x_2<<8) + 0x474f575f676F777FULL + max; 
	min_byte_excesses =  (byte_prefix_sums_x_2<<8) + 0x4951596169717981ULL - min; 
}


inline uint8_t bit_magic::first_excess_position(uint64_t w, uint8_t excess_val, uint64_t &byte_prefix_sums_x_2){
	uint64_t max_byte_excesses;
	bit_magic::max_byte_excesses(w, max_byte_excesses, byte_prefix_sums_x_2);
	assert(excess_val <= 64);
	max_byte_excesses = (max_byte_excesses-_8_x_the_byte[excess_val])&0x8080808080808080ULL;
	if(max_byte_excesses){
		uint8_t block_of_occurence_x_8 = DeBruijn64ToIndex[((max_byte_excesses&-max_byte_excesses)*DeBruijn64)>>58]&0xF8;
		assert(block_of_occurence_x_8<64);
		uint64_t sums = byte_prefix_sums_x_2<<8;
		excess_val = (excess_val + (block_of_occurence_x_8)) - ((sums >> (block_of_occurence_x_8))&0xFF);
		assert(excess_val <= 8);
		return first_pos_of_excess_val[((w>>(block_of_occurence_x_8))&0xFF)+(excess_val<<8)]+block_of_occurence_x_8;
	}
	else{
		return 64;
	}
}

inline uint8_t bit_magic::first_excess_position_naive(uint64_t w, uint8_t excess_val, uint64_t &byte_prefix_sums_x_2){
	byte_prefix_sums_x_2 = w-( (w>>1) & 0x5555555555555555ull);
	byte_prefix_sums_x_2 = (byte_prefix_sums_x_2 & 0x3333333333333333ull) + ((byte_prefix_sums_x_2 >> 2) & 0x3333333333333333ull);
	byte_prefix_sums_x_2 = (byte_prefix_sums_x_2 + (byte_prefix_sums_x_2 >> 4)) & 0x0f0f0f0f0f0f0f0full;
	byte_prefix_sums_x_2 *= 0x0101010101010101ull;
	byte_prefix_sums_x_2 <<= 1;
	for(uint8_t i=0;i<64;++i, w>>=1){
		excess_val += (1-((w&1)<<1));
		if(excess_val==0)
			return i;
	}
	return 64;
}

inline uint8_t bit_magic::find_open(uint64_t w, uint8_t close_parenthesis_index){
	if(!close_parenthesis_index)
		return 64;
	assert( close_parenthesis_index <= 64 );
	w <<= (64-close_parenthesis_index); 
	uint8_t near_match = very_near_find_open[(w>>56)&0xFF];// find near matches
//	std::cerr<<"near_match="<<(int)near_match<<std::endl;
    if( near_match ){ // the matching opening parenthesis is at most 7 positions before the closing paranthesis
		return close_parenthesis_index-near_match;
	}else{ // the matching opening parenthesis is not among the first 8 preceding bits/parentheses
		// calculate 
		uint64_t min_byte_excesses, max_byte_excesses, byte_prefix_sums_x_2;
//std::cerr<<"w="<<w<<" excess_val="<<(int)excess_val<<std::endl;
		bit_magic::min_max_byte_excesses(w, min_byte_excesses, max_byte_excesses, byte_prefix_sums_x_2);
		int8_t desired_excess = (byte_prefix_sums_x_2 >> 56) - 65;
//		std::cerr<<"desired_excess="<<(int)desired_excess<<"  "<<(int)(64-close_parenthesis_index)<<std::endl;
		if( desired_excess >= 0 ){
			max_byte_excesses = (max_byte_excesses-_8_x_the_byte[desired_excess])&0x8080808080808080ULL;
			max_byte_excesses &= ~(min_byte_excesses-_8_x_the_byte[desired_excess]-_8_x_the_byte[1]);
		}else{ // desire_excess <0
//			std::cerr<<"max_byte_excesses="<<max_byte_excesses<<std::endl;
//			std::cerr<<"min_byte_excesses="<<min_byte_excesses<<std::endl;
			max_byte_excesses = (max_byte_excesses+_8_x_the_byte[-desired_excess])&0x8080808080808080ULL;
			max_byte_excesses &= ~(min_byte_excesses+_8_x_the_byte[-desired_excess]-_8_x_the_byte[1]);
		}
//		std::cerr<<"max_byte_excesses="<<max_byte_excesses<<std::endl;
		if(max_byte_excesses){
			max_byte_excesses |= (max_byte_excesses>>8);
			max_byte_excesses |= (max_byte_excesses>>16);
			max_byte_excesses |= (max_byte_excesses>>32);
			max_byte_excesses -= (max_byte_excesses>>8);
			assert( max_byte_excesses = (max_byte_excesses&-max_byte_excesses) );
			uint8_t block_of_occurence_x_8 = DeBruijn64ToIndex[(max_byte_excesses*DeBruijn64)>>58]-7;
//std::cerr<<"block_of_occurence_x_8="<<(int)block_of_occurence_x_8<<std::endl;
//std::cerr<<"result+(64-close_parenthesis_index)="<<result+(64-close_parenthesis_index)<<std::endl;
			assert(block_of_occurence_x_8 <= 56);
			byte_prefix_sums_x_2 <<= 8;
			int excess_val = (int)(desired_excess + block_of_occurence_x_8 + 8) - ((byte_prefix_sums_x_2 >> block_of_occurence_x_8)&0xFF);
//std::cerr<<"desired_excess+block_of_occurence_x_8="<< (int)(desired_excess + block_of_occurence_x_8) <<std::endl;
//std::cerr<<"excess_so_far ="<< ((int)((byte_prefix_sums_x_2 >> block_of_occurence_x_8)&0xFF))-(block_of_occurence_x_8) <<std::endl;
//std::cerr<<"byte_prefix_sums ="<< byte_prefix_sums_x_2 <<std::endl;
//std::cerr<<"excess_val="<<(int)excess_val<<std::endl;		
			assert(excess_val<=16);
			return last_pos_of_excess_val[((w>>(block_of_occurence_x_8))&0xFF)+(excess_val<<8)]+block_of_occurence_x_8-63+close_parenthesis_index;
		}
		else // no near parenthesis found 
			return 64; 
	}
}

inline uint8_t bit_magic::find_open_naive(uint64_t w, uint8_t close_parenthesis_index){
	for(int i=close_parenthesis_index-1, excess_val=-1; i>=0; --i){
		if( (w>>i)&1 )
			++excess_val;
		else
			--excess_val;
		if(excess_val==0)
			return i;
	}
	return 64;
}

inline uint8_t bit_magic::find_enclose_naive(uint64_t w, uint8_t open_parenthesis_index){
	for(int i=open_parenthesis_index-1, excess_val=0; i>=0; --i){
		if( (w>>i)&1 )
			++excess_val;
		else
			--excess_val;
		if(excess_val==1)
			return i;	
	}
	return 64;
}

inline uint8_t bit_magic::find_enclose(uint64_t w, uint8_t open_parenthesis_index){
	return find_open(w, open_parenthesis_index);
}

inline uint8_t bit_magic::last_excess_position(uint64_t w, uint8_t excess_val, uint64_t &byte_prefix_sums_x_2){
	uint64_t min_byte_excesses, max_byte_excesses;
//std::cerr<<"w="<<w<<" excess_val="<<(int)excess_val<<std::endl;
	bit_magic::min_max_byte_excesses(w, min_byte_excesses, max_byte_excesses, byte_prefix_sums_x_2);
//std::cerr<<"min_byte_accesses="<<min_byte_excesses<<std::endl;
//std::cerr<<"max_byte_accesses="<<max_byte_excesses<<std::endl;
	assert(excess_val <= 64);
	max_byte_excesses = (max_byte_excesses-_8_x_the_byte[excess_val])&0x8080808080808080ULL;
	max_byte_excesses &= ~(min_byte_excesses-_8_x_the_byte[excess_val]-_8_x_the_byte[1]);
//std::cerr<<"max_byte_excesses="<<max_byte_excesses<<std::endl;
	if(max_byte_excesses){
//		std::cerr<<"max_byte_accesses="<<max_byte_excesses<<std::endl;
		uint8_t block_of_occurence_x_8;
		if(max_byte_excesses&0x8080808000000000ULL){// 4..7
			if(max_byte_excesses&0x8080000000000000ULL){// 6..7
				block_of_occurence_x_8 = 48 + ((max_byte_excesses>>60)&0x8);
			}else{// 4..5
				block_of_occurence_x_8 = 32 + ((max_byte_excesses>>44)&0x8);
			}
		}else{// 0..3
			if(max_byte_excesses&0x8080808080800000ULL){//2..3
				block_of_occurence_x_8 = 16 + ((max_byte_excesses>>28)&0x8);	
			}else{//0..1
				block_of_occurence_x_8 = (max_byte_excesses>>12)&0x8;
			}
		}
//			std::cerr<<" block_of_oc...="<<(int)block_of_occurence_x_8<<std::endl;
#ifndef NDEBUG		
#endif		
		assert(block_of_occurence_x_8<64);
		uint64_t sums = byte_prefix_sums_x_2<<8;
		excess_val = (excess_val + (block_of_occurence_x_8)+8) - ((sums >> (block_of_occurence_x_8))&0xFF);
//		if(excess_val>8+8){
//			std::cerr<<"block_of_occurence_x_8="<< (int)block_of_occurence_x_8  <<" excess_val="<<(int)excess_val<<std::endl;
//			std::cerr<<"(sums >> (block_of_occurence_x_8))&0xFF = "<< ((sums >> (block_of_occurence_x_8))&0xFF) << std::endl;
//		}
		assert(excess_val <= 16);
//		std::cerr<<"excess_val="<<(int)((int)excess_val)-8<<std::endl;
//		std::cerr<<"pos_inside_8_bit_block "<<(int)last_pos_of_excess_val[((w>>(block_of_occurence_x_8))&0xFF)+(excess_val<<8)]<<std::endl;
//		std::cerr<<"pos_inside_8_bit_block+block_of... "<<(int)last_pos_of_excess_val[((w>>(block_of_occurence_x_8))&0xFF)+(excess_val<<8)]+block_of_occurence_x_8<<std::endl;
		return last_pos_of_excess_val[((w>>(block_of_occurence_x_8))&0xFF)+(excess_val<<8)]+block_of_occurence_x_8;
	}
	else{
		return 64;
	}
}


inline uint8_t bit_magic::last_excess_position_naive(uint64_t w, uint8_t excess_val, uint64_t &byte_prefix_sums_x_2){
	byte_prefix_sums_x_2 = w-( (w>>1) & 0x5555555555555555ull);
	byte_prefix_sums_x_2 = (byte_prefix_sums_x_2 & 0x3333333333333333ull) + ((byte_prefix_sums_x_2 >> 2) & 0x3333333333333333ull);
	byte_prefix_sums_x_2 = (byte_prefix_sums_x_2 + (byte_prefix_sums_x_2 >> 4)) & 0x0f0f0f0f0f0f0f0full;
	byte_prefix_sums_x_2 *= 0x0101010101010101ull;
	byte_prefix_sums_x_2 <<= 1;
	uint8_t res = 64;
	uint8_t my_excess_val=0;
	for(uint8_t i=0;i<64;++i, w>>=1){
		if(w&1)
			++my_excess_val;
		else
			--my_excess_val;
//if(i%8==0) std::cerr<<std::endl;
//std::cerr<<"("<<(int)i<<","<<(int)my_excess_val<<")";		
		if(excess_val==my_excess_val)
			res = i;
	}
//std::cerr<<std::endl;	
	return res;
}


inline uint16_t bit_magic::max_excess3(uint64_t x, uint16_t &b1Cnt){
	uint64_t sum = x-( (x>>1) & 0x5555555555555555ull);
	sum = (sum & 0x3333333333333333ull) + ((sum >> 2) & 0x3333333333333333ull);
	sum = (sum + (sum >> 4)) & 0x0f0f0f0f0f0f0f0full;
	sum *= 0x0101010101010101ull;
	b1Cnt = sum>>56;
	sum = (( (sum<<1) | 0x8080808080808080ULL) - 0x4038302820181008ULL); 
/*	uint8_t maxi = max_excess_8bit[(uint8_t)x] | 0x80;
	x>>=8;
	for(uint8_t i=1,m;i<8;++i, sum>>=8, x>>=8){
		if( (m = (((uint8_t)sum)+max_excess_8bit[(uint8_t)x])) > maxi ) maxi = m;
//		std::cerr<<"m="<<(int)(m&0x7F)<<" "<<(int)max_excess_8bit[(uint8_t)x]<<" sum="<<(int)((uint8_t)sum)<<std::endl;
	}
	return maxi&0x7F;
*/	uint8_t maxi[8] = {	static_cast<uint8_t>(max_excess_8bit[(uint8_t)x]|0x80), 
						static_cast<uint8_t>((sum)+		max_excess_8bit[(uint8_t)(x>>8)]),
						static_cast<uint8_t>((sum>>8)+	max_excess_8bit[(uint8_t)(x>>16)]),
						static_cast<uint8_t>((sum>>16)+	max_excess_8bit[(uint8_t)(x>>24)]),
						static_cast<uint8_t>((sum>>24)+	max_excess_8bit[(uint8_t)(x>>32)]),
						static_cast<uint8_t>((sum>>32)+	max_excess_8bit[(uint8_t)(x>>40)]),
						static_cast<uint8_t>((sum>>40)+	max_excess_8bit[(uint8_t)(x>>48)]),
						static_cast<uint8_t>((sum>>48)+	max_excess_8bit[(uint8_t)(x>>56)])};
	if(maxi[0]<maxi[1])maxi[0]=maxi[1];
	if(maxi[2]<maxi[3])maxi[2]=maxi[3];
	if(maxi[4]<maxi[5])maxi[4]=maxi[5];
	if(maxi[6]<maxi[7])maxi[6]=maxi[7];
	if(maxi[0]<maxi[2])maxi[0]=maxi[2];
	if(maxi[4]<maxi[6])maxi[4]=maxi[6];
	return (maxi[0]<maxi[4])?maxi[4]&0x7F:maxi[0]&0x7F;
}

// see page 11, Knuth TAOCP Vol 4 F1A
inline uint64_t bit_magic::b1Cnt(uint64_t x){
#ifdef __SSE4_2__	
	return __builtin_popcountll(x);
#else
	x = x-( (x>>1) & 0x5555555555555555ull);
	x = (x & 0x3333333333333333ull) + ((x >> 2) & 0x3333333333333333ull);
	x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0full;
	return  (0x0101010101010101ull*x >> 56);
#endif	
}

inline uint32_t bit_magic::b1Cnt32(uint32_t x){
	x = x-( (x>>1) & 0x55555555);
	x = (x & 0x33333333) + ((x>>2) & 0x33333333);
	return (0x10101010*x >>28)+(0x01010101*x >>28);
}

inline uint32_t bit_magic::b1CntNaive(uint64_t x){
	uint32_t res = 0;
	for(int i=0;i<64;++i){
		res += x&1;
		x>>=1;
	}
	return res;
}

/*
inline uint32_t bit_magic::b11Cnt(uint64_t x){
	// extract "10" 2bit blocks
	uint64_t ex10 = ((x^(x<<1))&0xAAAAAAAAAAAAAAAAULL)&x;
	// extract "10" 2bit blocks
	uint64_t ex01 = ((x^(x>>1))&0x5555555555555555ULL)&x;
	// extract "11" 2bit blocks
	uint64_t ex11 =  (x&(x<<1))&0xAAAAAAAAAAAAAAAAULL;
	ex11 |= (ex11>>1);
	x = (((ex11)+(ex10<<1))&(ex01))|(ex11&0x5555555555555555ULL);
	x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
	x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
	return  (0x0101010101010101ULL*x >> 56);
}*/

inline uint32_t bit_magic::b11Cnt(uint64_t x, uint64_t &c){
	// extract "11" 2bit blocks
	uint64_t ex11 =  (x&(x>>1))&0x5555555555555555ULL, t;
	// extract "10" 2bit blocks
	uint64_t ex10or01 = (ex11|(ex11<<1))^x;

	x = ex11 | ( (t=(ex11|(ex11<<1))+(((ex10or01<<1)&0x5555555555555555ULL)|c))&(ex10or01&0x5555555555555555ULL) );
	c = (ex10or01>>63) or (t < (ex11|(ex11<<1)));

	x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
	x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
	return  (0x0101010101010101ULL*x >> 56);
}

inline uint32_t bit_magic::b11Cnt(uint64_t x){
	// extract "11" 2bit blocks
	uint64_t ex11 =  (x&(x>>1))&0x5555555555555555ULL;
	// extract "10" 2bit blocks
	uint64_t ex10or01 = (ex11|(ex11<<1))^x;

	x = ex11 | ( ((ex11|(ex11<<1))+((ex10or01<<1)&0x5555555555555555ULL))&(ex10or01&0x5555555555555555ULL) );

	x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
	x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
	return  (0x0101010101010101ULL*x >> 56);
}

inline uint32_t bit_magic::b11CntS(uint64_t x, uint64_t &c){
	uint64_t ex11 =  (x&(x>>1))&0x5555555555555555ULL;
	uint64_t ex10or01 = (ex11|(ex11<<1))^x;
	ex11 += ( ((ex11|(ex11<<1))+(((ex10or01<<1)&0x5555555555555555ULL)|c)) & ((ex10or01&0x5555555555555555ULL)|ex11) );
	c = (x>>63) and (ex11>>62);
	x = ex11-( (ex11>>1) & 0x5555555555555555ULL);
	x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
	x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
	return  (0x0101010101010101ULL*x >> 56);
}

inline uint32_t bit_magic::b11CntS(uint64_t x){
	uint64_t ex11 =  (x&(x>>1))&0x5555555555555555ULL;
	uint64_t ex10or01 = (ex11|(ex11<<1))^x;
	ex11 += ( ((ex11|(ex11<<1))+((ex10or01<<1)&0x5555555555555555ULL)) & ((ex10or01&0x5555555555555555ULL)|ex11) );
	x = ex11-( (ex11>>1) & 0x5555555555555555ULL);
	x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
	x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
	return  (0x0101010101010101ULL*x >> 56);
}

inline uint32_t bit_magic::b11CntNaive(uint64_t x){
	uint32_t res = 0;
	for(int i=0;i<64;++i){
		if( (x&3) == 3 ){
			++res;
			++i;
			x>>=1;
		}
		x>>=1;
	}
	return res;
}

inline uint32_t bit_magic::eB11Cnt(uint64_t x){
	x = (x&(x>>1))&0x5555555555555555ULL;
	x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
	x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
	return  (0x0101010101010101ULL*x >> 56);
}

inline uint32_t bit_magic::b10Cnt(uint64_t x, uint64_t &c){
	uint32_t res = b1Cnt( (x ^ (( x<<1 ) | c) ) & (~x) );
	c = (x >> 63);
	return res;
}

inline uint64_t bit_magic::b10Map(uint64_t x, uint64_t c){
	return ( ( x ^ ( ( x << 1) | c ) ) & ( ~x ) );
}

inline uint32_t bit_magic::b01Cnt(uint64_t x, uint64_t &c){
	uint32_t res = b1Cnt( (x ^ (( x<<1 ) | c) ) & x );
	c = (x >> 63);
	return res;
}
inline uint64_t bit_magic::b01Map(uint64_t x, uint64_t c){
	return ( ( x ^ ( ( x << 1) | c ) ) &  x  );
}

inline uint32_t bit_magic::b10CntNaive(uint64_t x, uint64_t &c){
	uint32_t sum = 0, lastbit = c;
	for(uint32_t i=0; i<64; ++i){
		if( (x&1) == 0 and lastbit == 1){
			++sum;
		}
		lastbit = (x&1);
		x >>= 1;
	}
	c = lastbit;
	return sum;
}

inline uint32_t bit_magic::i1BP(uint64_t x, uint32_t i){
// TODO: Maybe case i<16 as a special case
//	assert(i>0 && i<=b1Cnt(x));
	/*register*/ uint64_t s = x, b;  // s = sum
	s = s-( (s>>1) & 0x5555555555555555ULL);
	s = (s & 0x3333333333333333ULL) + ((s >> 2) & 0x3333333333333333ULL);
	s = (s + (s >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
	s = 0x0101010101010101ULL*s;
	b = (s+PsOverflow[i]);//&0x8080808080808080ULL;// add something to the partial sums to cause overflow
//	b =DeBruijn64ToIndex[((b&-b)*DeBruijn64)>>58];
	i = (i-1)<<8;
	if(b&0x0000000080000000ULL)// byte <=3
		if(b&0x0000000000008000ULL)//byte <= 1
			if(b&0x0000000000000080ULL)
				return    Select256[ (x&0xFFULL) + i];
			else
				return 8 +Select256[(((x>>8)&0xFFULL)  + i - ((s&0xFFULL)<<8))&0x7FFULL];//byte 1;
		else//byte >1
			if(b&0x0000000000800000ULL)//byte <=2
				return 16+Select256[(((x>>16)&0xFFULL) + i - (s&0xFF00ULL))&0x7FFULL];//byte 2;
			else
				return 24+Select256[(((x>>24)&0xFFULL) + i - ((s>>8)&0xFF00ULL))&0x7FFULL];//byte 3;
	else//  byte > 3
		if(b&0x0000800000000000ULL)// byte <=5
			if(b&0x0000008000000000ULL)//byte <=4
				return 32+Select256[(((x>>32)&0xFFULL) + i - ((s>>16)&0xFF00ULL))&0x7FFULL];//byte 4;
			else
				return 40+Select256[(((x>>40)&0xFFULL) + i - ((s>>24)&0xFF00ULL))&0x7FFULL];//byte 5;
		else// byte >5
			if(b&0x0080000000000000ULL)//byte<=6
				return 48+Select256[(((x>>48)&0xFFULL) + i - ((s>>32)&0xFF00ULL))&0x7FFULL];//byte 6;
			else
				return 56+Select256[(((x>>56)&0xFFULL) + i - ((s>>40)&0xFF00ULL))&0x7FFULL];//byte 7;
	return 0;
}

inline uint32_t bit_magic::j1BP(uint64_t x, uint32_t j){
	register uint64_t s = x, b, l;  // s = sum   b = byte [1..8] multiplied by 8 in which the j 1 is located
	s = s-( (s>>1) & 0x5555555555555555ULL);
	s = (s & 0x3333333333333333ULL) + ((s >> 2) & 0x3333333333333333ULL);
	s = (s + (s >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
	s = 0x0101010101010101ULL*s;
	b = ((((((((j-1)*0x0101010101010101ULL)|0x8080808080808080ULL)-s)&0x8080808080808080ULL)>>7)*0x0101010101010101ULL)>>53)&0xFFFFFFFFFFFFFFF8ULL; // s <=_8 (j*L_8)
	// b>0 if there are at least j 1s in x
	l = j - (((s<<8)>>b)&0xFF);
	s = (((((((((x>>b)&0xFF)*0x0101010101010101ULL)&0x8040201008040201ULL)|0x8080808080808080ULL)-0x0101010101010101ULL)|(((x>>b)&0x80ULL)<<56))&0x8080808080808080ULL)>>7)*0x0101010101010101ULL;
//	std::cerr<<"j="<<j<<" in byte "<<b<<" l="<<l<<" s="<<s<<std::endl;
	return b + ((((((((l-1)*0x0101010101010101ULL)|0x8080808080808080ULL)-s)&0x8080808080808080ULL)>>7)*0x0101010101010101ULL)>>56);
}

inline uint32_t bit_magic::i1BPNaive(uint64_t x, uint32_t i){
	uint32_t pos = 0;
	while(x){
		i -= x&1;
		x>>=1;
		if(!i) break;
		++pos;
	}
	return pos;
}

// builtin version of sse version or
// 64 bit version of 32 bit proposial of
// http://www-graphics.stanford.edu/~seander/bithacks.html
inline const uint32_t bit_magic::l1BP(uint64_t x){
#ifdef __SSE4_2__
	if( x == 0 )
		return 0;
	return 63 - __builtin_clzll(x);
#else	
	register uint64_t t,tt; // temporaries
	if((tt = x >> 32)){ // l1BP >= 32
		if((t = tt >> 16)){ // l1BP >= 48
			return (tt = t >> 8) ? 56 + L1BP[tt] : 48 + L1BP[t];
		}
		else{ // l1BP < 48
			return (t = tt >> 8) ? 40 + L1BP[t] : 32 + L1BP[tt];
		}
	}
	else{ // l1BP < 32
		if((t = x >> 16)){ // l1BP >= 16
			return (tt = t >> 8) ? 24 + L1BP[tt] : 16 + L1BP[t];
		}
		else{ // l1BP < 16
			return (tt = x >> 8) ?  8 + L1BP[tt] : L1BP[x];
		}
	}
#endif	
}

inline const uint32_t bit_magic::l1BPNaive(uint64_t x){
	if(x&0x8000000000000000ULL)
		return 63;
	x<<=1;
	for(int i=62; i>=0; ++i)
		if(x&0x8000000000000000ULL)
			return i;
		else
			x<<=1;
	return 0;
}

// details see: http://citeseer.ist.psu.edu/leiserson98using.html
// or page 10, Knuth TAOCP Vol 4 F1A
inline const uint32_t bit_magic::r1BP(uint64_t x){
#ifdef __SSE4_2__
	if(x==0)
		return 0;
	return __builtin_ctzll(x);
#else	
	if(x&1) return 0;
	if(x&3) return 1;
	if(x&7) return 2;
	if(x&0x7F){ // in average every second random number x can be answered this way
//		return 0;
		return lookupr1BP[(x&0x7F)>>3]+3;
	}
	// x&-x equals x with only the lsb set
		//default:
	return DeBruijn64ToIndex[((x&-x)*DeBruijn64)>>58];
#endif
}

inline const uint32_t bit_magic::r1BPNaive(uint64_t x){
	if(x&1)
		return 0;
	x>>=1;
	for(int i=1;i<64;++i)
		if(x&1)
			return i;
		else
			x>>=1;
	return 63;
}

inline const uint32_t bit_magic::l11BP(uint64_t x){
	// extract "11" 2bit blocks
	uint64_t ex11 =  (x&(x>>1))&0x5555555555555555ULL;
//	std::cerr<<"ex11="<<ex11<<std::endl;
	// extract "10" 2bit blocks
	uint64_t ex10or01 = (ex11|(ex11<<1))^x;
//	uint64_t ex10or01 = (x^(x<<1))&0xAAAAAAAAAAAAAAAAULL;
//	ex10or01 = (ex10or01|(ex10or01>>1))&x;
	// extract "10" 2bit blocks
	ex11 += ( ((ex11|(ex11<<1))+((ex10or01<<1)&0x5555555555555555ULL)) & ((ex10or01&0x5555555555555555ULL)|ex11) );
    //
//	std::cerr<<"ex11="<<ex11<<std::endl;
//	uint32_t p = l1BP(ex11);// p % 2 == 0  and (x&(1<<p))==1
//	std::cerr<<"p="<<p<<std::endl;
	return l1BP(ex11);/*
	if( p == 0 and (x&1) == 0 )
		return 0;
	return p + (x>>(p+1)&1);
	*/
}
/*
inline const uint64_t bit_magic::all11BPs(uint64_t x, uint32_t c){
	uint64_t ex11 =  (x&(x>>1))&0x5555555555555555ULL;
	uint64_t ex10or01 = (ex11|(ex11<<1))^x;
	ex11 += ( ((ex11|(ex11<<1))+(((ex10or01<<1)&0x5555555555555555ULL)|c)) & ((ex10or01&0x5555555555555555ULL)|ex11) );
	return ex11;
}
*/

inline const uint64_t bit_magic::all11BPs(uint64_t x, bool &c){
	uint64_t ex11 =  (x&(x>>1))&0x5555555555555555ULL;
	uint64_t ex10or01 = (ex11|(ex11<<1))^x;
	ex11 += ( ((ex11|(ex11<<1))+(((ex10or01<<1)&0x5555555555555555ULL)|c)) & ((ex10or01&0x5555555555555555ULL)|ex11) );
	c = (x&(~ex11))>>63;
	return ex11;
}

inline const uint32_t bit_magic::i11BP(uint64_t x, uint32_t i, uint32_t c){
	uint64_t ex11 =  (x&(x>>1))&0x5555555555555555ULL;
	uint64_t ex10or01 = (ex11|(ex11<<1))^x;
	ex11 += ( ((ex11|(ex11<<1))+(((ex10or01<<1)&0x5555555555555555ULL)|c)) & ((ex10or01&0x5555555555555555ULL)|ex11) );
	return i1BP(ex11,i);
}

inline const uint32_t bit_magic::eI11BP(uint64_t x, uint32_t i){
	return i1BP((x&(x<<1))&0xAAAAAAAAAAAAAAAAULL, i);
}

inline void bit_magic::writeInt(uint64_t *word, uint64_t x, uint8_t offset, const uint8_t len){
	x &= bit_magic::Li1Mask[len];
	if( offset + len < 64 ){
		*word &=
			((bit_magic::All1Mask << (offset+len)) | bit_magic::Li1Mask[offset]); // mask 1..10..01..1
		*word |= (x << offset);
//		*word ^= ((*word ^ x) & (bit_magic::Li1Mask[len] << offset) );
//      surprisingly the above line is slower than the lines above		
	}else{
		*word &=
			((bit_magic::Li1Mask[offset]));  // mask 0....01..1
		*word |= (x << offset);
		if( (offset = (offset+len)&0x3F ) ){// offset+len > 64
			*(word+1) &= (~bit_magic::Li1Mask[offset]); // mask 1...10..0
//			*(word+1) &= bit_magic::Li0Mask[offset]; // mask 1...10..0
//          surprisingly the above line is slower than the line above 
			*(word+1) |= (x >> (len-offset));
		}
	}
}

inline void bit_magic::writeInt2(uint64_t *word, uint64_t x, uint8_t offset, const uint8_t len){
//	*word ^= ((x)& 111ULL << offset);
//	*word ^= ((*word ^ x)& bit_magic::Li1Mask[len] << offset);
//	if( offset + len > 64 ){
//		offset = offset + len - 64;
//		*(word+1) ^= (( *(word+1) ^ (x >> (len-offset)) ) & bit_magic::Li1Mask[offset]);
//	}
}

inline void bit_magic::writeIntAndMove(uint64_t* &word, uint64_t x, uint8_t &offset, const uint8_t len){
	x &= bit_magic::Li1Mask[len];
	if( offset + len < 64 ){
		*word &=
			((bit_magic::All1Mask << (offset+len)) | bit_magic::Li1Mask[offset]); // mask 1..10..01..1
		*word |= (x << offset);
		offset += len;
	}else{
		*word &=
			((bit_magic::Li1Mask[offset]));  // mask 0....01..1
		*word |= (x << offset);
		if( (offset= (offset+len))>64 ){// offset+len >= 64
			offset &= 0x3F;
			*(++word) &= (~bit_magic::Li1Mask[offset]); // mask 1...10..0
			*word |= (x >> (len-offset));
		}else{
			offset = 0;
			++word;
		}
	}
}

inline uint64_t bit_magic::readInt(const uint64_t *word, uint8_t offset, const uint8_t len){
	uint64_t w1 = (*word)>>offset;
	if( (offset+len) > 64 ){ // if offset+len > 64
		return w1 |  // w1 or w2 adepted:
			 ( (*(word+1) & bit_magic::Li1Mask[(offset+len)&0x3F])  // set higher bits zero
			  << (64-offset) ); // move bits to the left
	}else{
		return w1 & bit_magic::Li1Mask[len];
	}
}

inline uint64_t bit_magic::readIntAndMove(const uint64_t* &word, uint8_t &offset, const uint8_t len){
	uint64_t w1 = (*word)>>offset;
	if( (offset = (offset+len))>=64 ){  // if offset+len > 64
		if(offset==64){
			offset &= 0x3F;
			++word;
			return w1;
		}else{
			offset &= 0x3F;
			return w1 |
				(((*(++word)) & bit_magic::Li1Mask[offset]) << (len-offset));
		}
	}else{
		return w1 & bit_magic::Li1Mask[len];
	}
}

inline uint64_t bit_magic::readUnaryInt(const uint64_t* word, uint8_t offset){
	register uint64_t w = *word >> offset;
	if(w){
		return bit_magic::r1BP(w)/*+1*/;
	}else{
		if( 0!=(w=*(++word)) )
			return bit_magic::r1BP(w)+64-offset/*+1*/;
		uint64_t cnt=2;
		while( 0==(w=*(++word)) )
			++cnt;
		return bit_magic::r1BP(w)+(cnt<<6)/*cnt*64*/-offset/*+1*/;//+(64-offset)
	}
	return 0;
}

inline uint64_t bit_magic::readUnaryIntAndMove(const uint64_t* &word, uint8_t &offset){
	register uint64_t w = (*word) >> offset; // temporary variable is good for the performance
	if(w){
		uint8_t r = bit_magic::r1BP(w);
		offset = (offset + r+1)&0x3F;
		// we know that offset + r +1 <= 64, so if the new offset equals 0 increase word
		word += (offset==0);
		return r;
	}else{
		uint8_t rr=0;
		if( 0!=(w=*(++word)) ){
			rr = bit_magic::r1BP(w)+64-offset;
			offset = (offset+rr+1)&0x3F;
			word += (offset==0);
			return rr;
		}
		else{
			uint64_t cnt_1=1;
			while( 0==(w=*(++word)) )
				++cnt_1;
			rr = bit_magic::r1BP(w)+64-offset;
			offset = (offset+rr+1)&0x3F;
			word += (offset==0);
			return ((cnt_1)<<6) + rr;
		}
	}
	return 0;
}

inline void bit_magic::move_right(const uint64_t* &word, uint8_t &offset, const uint8_t len){
	if( (offset+=len)&0xC0 ){ // if offset >= 65
		offset&=0x3F;
		++word;
	}
}

inline void bit_magic::move_left(const uint64_t* &word, uint8_t &offset, const uint8_t len){
	if( (offset-=len)&0xC0  ){ // if offset-len<0
		offset&=0x3F;
		--word;
	}
}

inline uint64_t bit_magic::next(const uint64_t* word, uint64_t idx){
	word += (idx>>6);
	if( *word & ~Li1Mask[idx&0x3F] ){
		return (idx & ~((size_t)0x3F)) + r1BP(*word & ~Li1Mask[idx&0x3F]);
	}
	idx = (idx & ~((size_t)0x3F)) + 64;
	++word;
	while( *word==0 ){
		idx += 64;
		++word;
	}
	return idx + r1BP(*word);
}

inline uint64_t bit_magic::prev(const uint64_t* word, uint64_t idx){
	word += (idx>>6);
	if( *word & Li1Mask[(idx&0x3F)+1] ){
		return (idx & ~((size_t)0x3F)) + l1BP(*word & Li1Mask[(idx&0x3F)+1]);
	}
	idx = (idx & ~((size_t)0x3F)) - 64;
	--word;
	while( *word==0 ){
		idx -= 64;
		--word;
	}
	return idx + l1BP(*word);
}

} // end namespace sdsl

#endif
