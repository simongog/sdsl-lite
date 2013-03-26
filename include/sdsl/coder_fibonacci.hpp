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
/*! \file coder_fibonacci.hpp
    \brief coder_fibonacci.hpp contains the class sdsl::coder::fibonacci
	\author Simon Gog
 */
#ifndef SDSL_CODER_FIBONACCI_INCLUDED
#define SDSL_CODER_FIBONACCI_INCLUDED

#include "int_vector.hpp"

namespace sdsl
{

namespace coder
{

//! A class to encode and decode between Fibonacci and binary code.
class fibonacci
{
    public:
        //! Array contains precomputed values for the decoding of a number in the fibonacci system.
        static const uint64_t Fib2bin_0_95[(1<<12)*8];

        //! Array contains precomputed values for the decoding of a number in the fibonacci system
        /*! Entry x equals 0 if no 11 substring is in the binary representation of x otherwise
            Fib2binShift[x] contains the position (1..13) of the left 1 in the rightmost 11 substring in x.
        	E.g. Fib2binShift[3] = 2 and Fib2binShift[6] = 3.
         */
        static const uint8_t Fib2binShift[(1<<13)];

        //! Array contains precomputed values for the deconding of a prefix sum of fibonacci encoded integers
        /*! The 5 most significant bits contain information about how far to shift to get to the next encoded integer.
            If this 5 bits equal zero, there is no whole fibonacci number encoded in the 16 bits...
         */
        static const uint16_t Fib2bin_0_16_greedy[1<<16];

        typedef uint64_t size_type;

        static const uint8_t min_codeword_length = 2; // 11 represents 1 and is the code word with minimum length
        //! Get the number of bits that are necessary to encode the value w in fibonacci code.
        /*! \param w 64bit integer to get the length of its fibonacci encoding. Inclusive the terminating 1 of the code.
         */
        static uint8_t encoding_length(uint64_t w);
        //! Decode n Fibonacci encoded bits beginning at start_idx in the bitstring "data"
        /* \param data Bitstring
           \param start_idx Starting index of the decoding.
           \param n Number of values to decode from the bitstring.
           \param it Iterator
         */
        template<bool sumup, bool increment, class Iterator>
        static uint64_t decode(const uint64_t* data, const size_type start_idx, size_type n, Iterator it=(Iterator)NULL);

        template<bool sumup, bool increment, class Iterator>
        static uint64_t decode1(const uint64_t* data, const size_type start_idx, size_type n, Iterator it=(Iterator)NULL);



        //! Decode n Fibonacci encoded integers beginning at start_idx in the bitstring "data"  and return the sum of these values.
        /*! \param data Pointer to the beginning of the Fibonacci encoded bitstring.
            \param start_idx Index of the first bit to endcode the values from.
        	\param n Number of values to decode from the bitstring. Attention: There have to be at least n encoded values in the bitstring.
         */
        static uint64_t decode_prefix_sum(const uint64_t* data, const size_type start_idx, size_type n);

        //! Deocde n Fibonacci encoded integers beginning at start_idx and ending at end_idx (exclusive) in the bitstring "data" and return the sum of these values.
        /*! \sa decode_prefix_sum
          */
        static uint64_t decode_prefix_sum(const uint64_t* data, const size_type start_idx, const size_type end_idx, size_type n);

        template<class int_vector1, class int_vector2>
        static bool encode(const int_vector1& v, int_vector2& z);

        template<class int_vector>
        static uint64_t* raw_data(int_vector& v) {
            return v.m_data;
        };

        //! Encode one positive integer x to an int_vector at bit position start_idx.
        /*! \param x Positive integer to encode.
            \param z Raw data of vector to write the encoded form of x.
        	\param offset Start offset to write the encoded form of x in z. \f$0\leq offset< 64\f$.
         */
        static void encode(uint64_t x, uint64_t*& z, uint8_t& offset);

        template<class int_vector1, class int_vector2>
        static bool decode(const int_vector1& z, int_vector2& v);
};

inline uint8_t fibonacci::encoding_length(uint64_t w)
{
    if (w == 0) {
        return 93;
    }
    // This limit for the leftmost 1bit in the resulting fib code could be improved using a table
    uint8_t len_1 = bits::hi(w); // len-1 of the fib code
    while (++len_1 < (uint8_t)(sizeof(bits::Fib)/sizeof(bits::Fib[0])) && w >= bits::Fib[len_1]);
    return len_1+1;
}

template<class int_vector1, class int_vector2>
inline bool fibonacci::encode(const int_vector1& v, int_vector2& z)
{
    uint64_t z_bit_size = 0;
    register uint64_t w;
    const uint64_t zero_val = v.width() < 64 ? (1ULL)<<v.width() : 0;
    for (typename int_vector1::const_iterator it=v.begin(), end = v.end(); it != end; ++it) {
        if ((w=*it) == 0) {
            if (v.width() < 64) {
                w = zero_val;
            }
        }
        z_bit_size += encoding_length(w);
    }
    z.bit_resize(z_bit_size);
    if (z_bit_size & 0x3F) { // if z_bit_size % 64 != 0
        *(z.m_data + (z_bit_size>>6)) = 0; // initialize last word
    }
    uint64_t* z_data 	= z.m_data;
    uint8_t offset 		= 0;
    register uint64_t fibword_high = 0x0000000000000001ULL, fibword_low;
    register uint64_t t;
    for (typename int_vector1::const_iterator it=v.begin(), end = v.end(); it != end; ++it) {
        w = *it;
        if (w == 0) {
            w = zero_val;
        }
        int8_t len_1 = encoding_length(w)-1,j;
        fibword_low = 0x0000000000000001ULL;

        if (len_1 >= 64) { // length > 65
            fibword_high = 0x0000000000000001ULL;
            j = len_1-1;
            if (w == 0) { // handle special case
                fibword_high <<= 1;
                fibword_high |= 1;
                fibword_high <<= 1;
                w -= bits::Fib[len_1-1];
                j -= 2;
            }
            for (; j>63; --j) {
                fibword_high <<= 1;
                if (w >= (t=bits::Fib[j])) {
                    w -= t;
                    fibword_high |= 1;
                    if (w and j>64) {
                        fibword_high <<= 1;
                        --j;
                    } else {
                        fibword_high <<= (64-j);
                        break;
                    }
                }
            }
            j		= 64;
        } else {
            j = len_1-1;
        }

        for (; j >= 0; --j) {
            fibword_low <<= 1;
            if (w >= (t=bits::Fib[j])) {
                w -= t;
                fibword_low |= 1;
                if (w) {
                    fibword_low <<= 1;
                    --j;
                } else {
                    fibword_low <<= (j);
                    break;
                }
            }
        }
        if (len_1 >=64) {
            bits::write_int_and_move(z_data, fibword_low, offset, 64);
            bits::write_int_and_move(z_data, fibword_high, offset, len_1 - 63);
        } else {
            bits::write_int_and_move(z_data, fibword_low, offset, (len_1&0x3F) +1);
        }
    }
    z.width(v.width());
    return true;
}

inline void fibonacci::encode(uint64_t x, uint64_t*& z, uint8_t& offset)
{
    register uint64_t fibword_high = 0x0000000000000001ULL, fibword_low;
    register uint64_t t;
    int8_t len_1 = encoding_length(x)-1,j;
    fibword_low = 0x0000000000000001ULL;

    if (len_1 >= 64) { // length > 65
        fibword_high = 0x0000000000000001ULL;
        j = len_1-1;
        if (x == 0) { // handle special case
            fibword_high <<= 1;
            fibword_high |= 1;
            fibword_high <<= 1;
            x -= bits::Fib[len_1-1];
            j -= 2;
        }
        for (; j>63; --j) {
            fibword_high <<= 1;
            if (x >= (t=bits::Fib[j])) {
                x -= t;
                fibword_high |= 1;
                if (x and j>64) {
                    fibword_high <<= 1;
                    --j;
                } else {
                    fibword_high <<= (64-j);
                    break;
                }
            }
        }
        j	= 64;
    } else {
        j = len_1-1;
    }
    for (; j >= 0; --j) {
        fibword_low <<= 1;
        if (x >= (t=bits::Fib[j])) {
            x -= t;
            fibword_low |= 1;
            if (x) {
                fibword_low <<= 1;
                --j;
            } else {
                fibword_low <<= (j);
                break;
            }
        }
    }
    if (len_1 >=64) {
        bits::write_int_and_move(z, fibword_low, offset, 64);
        bits::write_int_and_move(z, fibword_high, offset, len_1 - 63);
    } else {
        bits::write_int_and_move(z, fibword_low, offset, (len_1&0x3F) +1);
    }
}

template<class int_vector1, class int_vector2>
bool fibonacci::decode(const int_vector1& z, int_vector2& v)
{
    uint64_t n = 0, carry = 0; // n = number of values to be decoded
    const uint64_t* data = z.data();
    // Determine size of v
    if (z.empty()) {// if z is empty we are ready with decoding
        v.width(z.width());
        v.resize(0);
        return true;
    }
    for (typename int_vector1::size_type i=0; i < (z.capacity()>>6)-1; ++i, ++data) {
        n += bits::b11Cnt(*data, carry);
    }
    if (z.capacity() != z.bit_size()) {
        n += bits::b11Cnt((*data) & bits::Li1Mask[z.bit_size()&0x3F], carry);
    } else {
        n += bits::b11Cnt(*data, carry);
    }
    std::cout<<"n="<<n<<std::endl;
    v.width(z.width()); v.resize(n);
    return decode<false, true>(z.data(), 0, n, v.begin());
}

template<bool sumup, bool increment, class Iterator>
inline uint64_t fibonacci::decode(const uint64_t* data, const size_type start_idx, size_type n, Iterator it)
{
    data += (start_idx >> 6);
    uint64_t w = 0, value = 0;
    int8_t buffered = 0; // bits buffered in w, in 0..64
    int8_t read = start_idx & 0x3F; // read bits in current *data 0..63
    int8_t shift = 0;
    uint32_t fibtable = 0;
    while (n) {// while not all values are decoded
        while (buffered < 13 and bits::b11Cnt(w) < n) {
            w |= (((*data)>>read)<<buffered);
            if (read >= buffered) {
                ++data;
                buffered += 64-read;
                read = 0;
            } else { // read < buffered
                read += 64-buffered;
                buffered = 64;
            }
        }
        value += Fib2bin_0_95[(fibtable<<12) | (w&0xFFF)];
        shift  = Fib2binShift[w&0x1FFF];
        if (shift > 0) {// if end of decoding
            w >>= shift;
            buffered -= shift;
            if (increment) *(it++) = value;
            if (!sumup and n!=1) value = 0;
            fibtable = 0;
            --n;
        } else { // not end of decoding
            w >>= 12;
            buffered -= 12;
            ++fibtable;
        }
    }
    return value;
}

template<bool sumup, bool increment, class Iterator>
inline uint64_t fibonacci::decode1(const uint64_t* data, const size_type start_idx, size_type n, Iterator it)
{
    data += (start_idx >> 6);
    uint64_t w = 0, value = 0;
    int8_t buffered = 0; // bits buffered in w, in 0..64
    int8_t read = start_idx & 0x3F; // read bits in current *data 0..63
    int8_t shift = 0;
    uint32_t fibtable = 0;
    uint8_t blocknr = (start_idx>>6)%9;
    while (n) {// while not all values are decoded
        while (buffered < 13 and bits::b11Cnt(w) < n) {
            w |= (((*data)>>read)<<buffered);
            if (read >= buffered) {
                ++blocknr;
                ++data;
                if (blocknr==8) {
                    ++data;
                    blocknr=0;
                }
                buffered += 64-read;
                read = 0;
            } else { // read < buffered
                read += 64-buffered;
                buffered = 64;
            }
        }
        value += Fib2bin_0_95[(fibtable<<12) | (w&0xFFF)];
        shift  = Fib2binShift[w&0x1FFF];
        if (shift > 0) {// if end of decoding
            w >>= shift;
            buffered -= shift;
            if (increment) *(it++) = value;
            if (!sumup)
                value = 0;
            fibtable = 0;
            --n;
        } else { // not end of decoding
            w >>= 12;
            buffered -= 12;
            ++fibtable;
        }
    }
    return value;
}

} // end namespace coder
} // end namespace sdsl

#endif
