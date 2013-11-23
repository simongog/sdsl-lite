/* sdsl - succinct data structures library
    Copyright (C) 2013 Simon Gog

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
/*! \file coder_elias_gamma.hpp
    \brief coder_elias_gamma.hpp contains the class sdsl::coder::elias_gamma
	\author Simon Gog
 */
#ifndef SDSL_CODER_ELIAS_GAMMA
#define SDSL_CODER_ELIAS_GAMMA

#include "int_vector.hpp"

namespace sdsl
{

namespace coder
{

//! A class to encode and decode between Elias-\f$\delta\f$ and binary code.
class elias_gamma
{
    public:
        typedef uint64_t size_type;

        static struct impl {
            //! Array contains precomputed values for the decoding of the prefix sum of Elias gamma encoded numbers.
            /*! The 8 most significant bits contain the length of decoded bits.
             *  The following 8 bits contain the number of decoded values.
             *  The last 16 bits contain the sum of the decoded values.
             */
            uint32_t prefixsum[1<<16];

            uint16_t prefixsum_8bit[(1<<8)*8];

            impl() {
                // initialize prefixsum
                for (uint64_t x=0; x < (1<<16); ++x) {
                    const uint64_t* w = &x; // copy of x
                    uint64_t value = 0;
                    uint16_t numbers = 0, offset=0, offset2=0;
                    while ((x >> offset) !=0) {
                        uint64_t len_1 = bits::read_unary(w, offset);
                        if (len_1 == 0) {
                            offset += 1;
                            value += 1;
                            ++numbers;
                        } else {
                            offset2 = offset + len_1 + 1;
                            if (offset2 + len_1 <= 16) {
                                value += bits::read_int(w, offset2, len_1) + (1ULL << len_1);
                                offset = offset2 + len_1;
                                ++numbers;
                            } else break;
                        }
                    }
                    uint32_t result=0;// the highest 8 bit equal the shift/offset, the second highest 8 bit equal the number of decoded values,
                    // and the last 16 bit equals value of decoded prefix sum
                    result = (offset << 24) | (numbers<<16) | value;
                    if (value>0)
                        assert(offset > 0  and numbers > 0 and offset<=16 and numbers <= 16);
                    prefixsum[x] = result;
                }
                // initialize prefixsum_8bit

                for (uint32_t maxi=1, idx=0; maxi<=8; ++maxi) {
                    for (uint64_t x=0; x < (1<<8); ++x) {
                        const uint64_t* w = &x; // copy of x
                        uint64_t value = 0;
                        uint32_t numbers = 0, offset=0, offset2=0;
                        while ((x >> offset) !=0 and numbers < maxi) {
                            uint64_t len_1 = bits::read_unary(w, offset);
                            if (len_1 == 0) {
                                offset += 1;
                                value += 1;
                                ++numbers;
                            } else {
                                offset2 = offset + len_1 + 1;
                                if (offset2 + len_1 <= 8) {
                                    value += bits::read_int(w, offset2, len_1) + (1ULL << len_1);
                                    offset = offset2 + len_1;
                                    ++numbers;
                                } else break;
                            }
                        }
                        uint16_t result=0;// the highest 8 bit equal the shift/offset, the second highest 8 bit equal the number of decoded values,
                        // and the last 16 bit equals value of decoded prefix sum
                        result = (offset << 12) | (numbers<<8) | value;
                        prefixsum_8bit[idx++] = result;
                    }
                }
            }
        } data;

        static const uint8_t min_codeword_length = 1; // 1 represents 1 and is the code word with minimum length
        static uint8_t encoding_length(uint64_t);
        //! Decode n Elias-delta encoded bits beginning at start_idx in the bitstring "data"
        /* \param data Bitstring
           \param start_idx Starting index of the decoding.
           \param n Number of values to decode from the bitstring.
           \param it Iterator to decode the values.
         */
        template<bool t_sumup, bool t_inc,class t_iter>
        static uint64_t decode(const uint64_t* data, const size_type start_idx, size_type n, t_iter it=(t_iter)nullptr);

        //! Decode n Elias gamma encoded integers beginning at start_idx in the bitstring "data"  and return the sum of these values.
        /*! \param data Pointer to the beginning of the Elias gamma encoded bitstring.
            \param start_idx Index of the first bit to endcode the values from.
        	\param n Number of values to decode from the bitstring. Attention: There have to be at least n encoded values in the bitstring.
         */
        static uint64_t decode_prefix_sum(const uint64_t* data, const size_type start_idx, size_type n);
        static uint64_t decode_prefix_sum(const uint64_t* data, const size_type start_idx, const size_type end_idx, size_type n);

        template<class int_vector>
        static bool encode(const int_vector& v, int_vector& z);
        template<class int_vector>
        static bool decode(const int_vector& z, int_vector& v);

        //! Encode one positive integer x to an int_vector at bit position start_idx.
        /* \param x Positive integer to encode.
           \param z Raw data of vector to write the encoded form of x.
           \param start_idx Beginning bit index to write the encoded form ox x in z.
        */
        static void encode(uint64_t x, uint64_t*& z, uint8_t& offset);

        template<class int_vector>
        static uint64_t* raw_data(int_vector& v) {
            return v.m_data;
        }
};

// \sa coder::elias_gamma::encoding_length
inline uint8_t elias_gamma::encoding_length(uint64_t w)
{
    uint8_t len_1 = w ? bits::hi(w) : 64;
    return 2*len_1 + 1;
}

template<class int_vector>
bool elias_gamma::encode(const int_vector& v, int_vector& z)
{
    typedef typename int_vector::size_type size_type;
    z.width(v.width());
    size_type z_bit_size = 0;
    uint64_t w;
    const uint64_t zero_val = v.width() < 64 ? (1ULL)<<v.width() : 0;
    for (typename int_vector::const_iterator it = v.begin(), end = v.end(); it != end; ++it) {
        if ((w=*it) == 0) {
            w = zero_val;
        }
        z_bit_size += encoding_length(w);
    }
    z.bit_resize(z_bit_size);   // Initial size of z
    if (z_bit_size & 0x3F) { // if z_bit_size % 64 != 0
        *(z.m_data + (z_bit_size>>6)) = 0; // initialize last word
    }
    z_bit_size = 0;
    uint64_t* z_data = z.m_data;
    uint8_t offset=0;
    size_type len_1; // TODO: change to uint8_t and test it
    for (typename int_vector::const_iterator it = v.begin(), end=v.end(); it != end; ++it) {
        w = *it;
        if (w == 0) {
            w = zero_val;
        }
        // (number of bits to represent w)-1
        if (!w) {
            len_1 = 64;
            bits::write_int_and_move(z_data, 0ULL, offset, 64);
            bits::write_int_and_move(z_data, 1ULL, offset, 1);
        } else {
            len_1 = bits::hi(w);
            bits::write_int_and_move(z_data, 1ULL << len_1, offset, len_1+1);
        }
        if (len_1) {
            bits::write_int_and_move(z_data, w, offset, len_1);
        }
    }
    return true;
}

inline void elias_gamma::encode(uint64_t x, uint64_t*& z, uint8_t& offset)
{
    uint8_t len_1 = 0;
    if (!x) {
        len_1 = 64;
        bits::write_int_and_move(z, 0, offset, 64);
        bits::write_int_and_move(z, 1, offset, 1);
    } else {
        len_1 = bits::hi(x);
        bits::write_int_and_move(z, 1ULL << len_1, offset, len_1+1);
    }
    if (len_1) {
        bits::write_int_and_move(z, x, offset, len_1);
    }
}

template<class int_vector>
bool elias_gamma::decode(const int_vector& z, int_vector& v)
{
    typename int_vector::size_type len_1, n = 0;
    const uint64_t* z_data	= z.data();
    const uint64_t* z_end	= z.data()+(z.bit_size()>>6);
    uint8_t offset 		= 0;
    while ((z_data < z_end) or (z_data==z_end and offset < (z.bit_size()&0x3F))) {
        len_1 = bits::read_unary_and_move(z_data, offset);
        if (len_1) {
            bits::move_right(z_data, offset, len_1);
        }
        ++n;
    }
    v.width(z.width());
    v.resize(n);
    return decode<false, true>(z.data(), 0, n, v.begin());
}

template<bool t_sumup, bool t_inc, class t_iter>
inline uint64_t elias_gamma::decode(const uint64_t* data, const size_type start_idx, size_type n, t_iter it)
{
    data += (start_idx >> 6);
    uint64_t value = 0;
    size_type i = 0;
    size_type len_1;
    uint8_t offset = start_idx & 0x3F;
    while (i++ < n) {// while not all values are decoded
        if (!t_sumup) value = 0;
        len_1 = bits::read_unary_and_move(data, offset); // read length of x-1
        if (!len_1) {
            value += 1;
        } else {
            value +=  bits::read_int_and_move(data, offset, len_1) + (len_1<64) * (1ULL << len_1);
        }
        if (t_inc) *(it++) = value;
    }
    return value;
}

} // end namespace coder
} // end namespace sdsl
#endif
