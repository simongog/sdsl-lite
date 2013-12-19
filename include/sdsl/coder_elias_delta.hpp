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
/*! \file coder_elias_delta.hpp
    \brief coder_elias_delta.hpp contains the class sdsl::coder::elias_delta
	\author Simon Gog
 */
#ifndef SDSL_CODER_ELIAS_DELTA
#define SDSL_CODER_ELIAS_DELTA

#include "int_vector.hpp"

namespace sdsl
{

namespace coder
{

//! A class to encode and decode between Elias-\f$\delta\f$ and binary code.
class elias_delta
{
    public:
        typedef uint64_t size_type;

        static struct impl {
            //! Array contains precomputed values for the decoding of the prefix sum of Elias delta encoded numbers.
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
                        uint64_t len_1_len = bits::read_unary(w, offset), len = 0;
                        if (len_1_len == 0) {
                            offset += 1;
                            value += 1;
                            ++numbers;
                        } else {
                            offset2 = offset + len_1_len +1;
                            len = bits::read_int(w, offset2, len_1_len) + (1ULL << len_1_len);
                            offset2 += len_1_len;
                            if (offset2 + len-1 <= 16) {
                                value += bits::read_int(w, offset2, len-1) + (1ULL << (len-1));
                                offset = offset2 + len - 1;
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
                            uint64_t len_1_len = bits::read_unary(w, offset), len = 0;
                            if (len_1_len == 0) {
                                offset += 1;
                                value += 1;
                                ++numbers;
                            } else {
                                offset2 = offset + len_1_len +1;
                                len = bits::read_int(w, offset2, len_1_len) + (1ULL << len_1_len);
                                offset2 += len_1_len;
                                if (offset2 + len-1 <= 8) {
                                    value += bits::read_int(w, offset2, len-1) + (1ULL << (len-1));
                                    offset = offset2 + len - 1;
                                    ++numbers;
                                } else break;
                            }
                        }
                        uint16_t result=0;// the highest 8 bit equal the shift/offset, the second highest 8 bit equal the number of decoded values,
                        // and the last 16 bit equals value of decoded prefix sum
                        result = (offset << 8) | (numbers<<4) | value;
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
        template<bool t_sumup, bool t_inc, class t_iter, class t_const_uint64_iter>
        static uint64_t decode(t_const_uint64_iter data, const size_type start_idx, size_type n, t_iter it=(t_iter)nullptr);

        //! Decode n Elias delta encoded integers beginning at start_idx in the bitstring "data"  and return the sum of these values.
        /*! \param data Pointer to the beginning of the Elias delta encoded bitstring.
            \param start_idx Index of the first bit to endcode the values from.
        	\param n Number of values to decode from the bitstring. Attention: There have to be at least n encoded values in the bitstring.
         */
        template<class t_const_uint64_iter>
        static uint64_t decode_prefix_sum(t_const_uint64_iter data, const size_type start_idx, size_type n);

        template<class t_const_uint64_iter>
        static uint64_t decode_prefix_sum(t_const_uint64_iter data, const size_type start_idx, const size_type end_idx, size_type n);

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

// \sa coder::elias_delta::encoding_length
inline uint8_t elias_delta::encoding_length(uint64_t w)
{
    uint8_t len_1 = w ? bits::hi(w) : 64;
    return len_1 + (bits::hi(len_1+1)<<1) + 1;
}

template<class int_vector>
bool elias_delta::encode(const int_vector& v, int_vector& z)
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
    size_type len, len_1_len; // TODO: change to uint8_t and test it
    for (typename int_vector::const_iterator it = v.begin(), end=v.end(); it != end; ++it) {
        w = *it;
        if (w == 0) {
            w = zero_val;
        }
        // (number of bits to represent w)
        len 		= w ? bits::hi(w)+1 : 65;
        // (number of bits to represent the length of w) -1
        len_1_len	= bits::hi(len);
        // Write unary representation for the length of the length of w
        bits::write_int_and_move(z_data, 1ULL << len_1_len, offset, len_1_len+1);
        if (len_1_len) {
            bits::write_int_and_move(z_data, len, offset, len_1_len);
            bits::write_int_and_move(z_data, w, offset, len-1);
        }
    }
    return true;
}

inline void elias_delta::encode(uint64_t x, uint64_t*& z, uint8_t& offset)
{
    uint8_t len, len_1_len;
    // (number of bits to represent w)
    len = x ? bits::hi(x)+1 : 65;
    // (number of bits to represent the length of w) - 1
    len_1_len	= bits::hi(len);
    // Write unary representation for the length of the length of w
    bits::write_int_and_move(z, 1ULL << len_1_len, offset, len_1_len+1);
    if (len_1_len) {
        bits::write_int_and_move(z, len, offset, len_1_len);
        bits::write_int_and_move(z, x, offset, len-1);
    }
}

template<class int_vector>
bool elias_delta::decode(const int_vector& z, int_vector& v)
{
    typename int_vector::size_type len_1_len, len, n = 0;
    const uint64_t* z_data	= z.data();
    const uint64_t* z_end	= z.data()+(z.bit_size()>>6);
    uint8_t offset 		= 0;
    while ((z_data < z_end) or (z_data==z_end and offset < (z.bit_size()&0x3F))) {
        len_1_len = bits::read_unary_and_move(z_data, offset);
        if (len_1_len) {
            len 	= bits::read_int_and_move(z_data, offset, len_1_len) + (1 << len_1_len);
            bits::move_right(z_data, offset, len-1);
        }
        ++n;
    }
    v.width(z.width());
    v.resize(n);
    return decode<false, true>(z.data(), 0, n, v.begin());
}

template<bool t_sumup, bool t_inc, class t_iter, class t_const_uint64_iter>
inline uint64_t elias_delta::decode(t_const_uint64_iter data, const size_type start_idx, size_type n, t_iter it)
{
    data += (start_idx >> 6);
    uint64_t value = 0;
    size_type i = 0;
    size_type len_1_len, len;
    uint8_t offset = start_idx & 0x3F;
    while (i++ < n) {// while not all values are decoded
        if (!t_sumup) value = 0;
        len_1_len = bits::read_unary_and_move(data, offset); // read length of length of x
        if (!len_1_len) {
            value += 1;
        } else {
            len 	=  bits::read_int_and_move(data, offset, len_1_len) + (1ULL << len_1_len);
            value	+= bits::read_int_and_move(data, offset, len-1) + (len-1<64) * (1ULL << (len-1));
        }
        if (t_inc) *(it++) = value;
    }
    return value;
}

template<class t_const_uint64_iter>
uint64_t elias_delta::decode_prefix_sum(t_const_uint64_iter data, const size_type start_idx, const size_type end_idx, size_type n)
{
    if (n==0)
        return 0;
    auto lastdata = data + ((end_idx+63) >> 6);
    data += (start_idx >> 6);
    uint64_t w = 0, value = 0;
    int16_t buffered = 0, read = start_idx & 0x3F;
    size_type i = 0;
    if (n + read <= 64) {
        if (((*data >> read)&bits::lo_set[n]) == bits::lo_set[n])
            return n;
    } else { // n+read > 64
        if ((*data >> read) == bits::lo_set[64-read]) {// all bits are set to 1
            value = 64-read;
            ++data;
            n -= (64-read);
            read = 0;
            while (n >= 64) {
                if (*data == 0xFFFFFFFFFFFFFFFFULL) {
                    value += 64;
                    ++data;
                    n -= 64;
                } else
                    goto start_decoding;
            }
            // 0 <= n <= 64
            if ((*data&bits::lo_set[n]) == bits::lo_set[n])
                return value + n;
        }
    }

start_decoding:
    while (i<n) {
        while (buffered < 64 and data < lastdata) {
fill_buffer:
            w |= (((*data)>>read)<<buffered);
            if (read >= buffered) {
                ++data;
                buffered += 64-read;
                read = 0;
            } else { // read buffered
                read += 64-buffered;
                buffered = 64;
            }
        }
        /*		if(w==0xFFFFFFFFFFFFFFFFULL){
        			i += 64;
        			value += 64;
        			if(i >= n)
        				return value - (i-n);
        			buffered = 0;
        			w = 0;
        		//	continue;
        			goto fill_buffer;
        		}
        */
//		uint32_t rbp = (w == 0xFFFFFFFFFFFFFFFFULL)?64:bits::lo(~w);
        uint32_t rbp = bits::lo(~w);
        if (rbp > 0) {
            i += rbp;
            value += rbp;
            if (i >= n) { // decoded too much
                return value - (i-n); // return corrected value
            }
            assert(buffered >= rbp);
            buffered -= rbp;
            w >>= rbp;
            if (buffered<16)
//				continue;
                goto fill_buffer;
        }
//		assert(w!=0xFFFFFFFFFFFFFFFFULL);
//		else
        {
            // i < n
begin_decode:
            uint32_t psum = elias_delta::data.prefixsum[w&0x0000FFFF];
            if (!psum or i+((psum>>16)&0x00FF) > n) {
                if (w==0) {// buffer is not full
                    w |= (((*data)>>read)<<buffered);
                    if (read >= buffered) {
                        ++data;
                        buffered += 64-read;
                        read = 0;
                    } else {
                        read += 64-buffered;
                        buffered = 64;
                    };
                    if (!w) {
                        w |= (((*data)>>read)<<buffered);
                        if (read >= buffered) {
                            ++data;
                            buffered += 64-read;
                            read = 0;
                        } else {
                            read += 64-buffered;
                            buffered = 64;
                        };
                    }
                }
//				assert(w>0);
                uint16_t len_1_len = bits::lo(w); // read length of length
                buffered -= (len_1_len+1);
                w >>= (len_1_len+1);
                if (len_1_len > buffered) {// buffer is not full
                    w |= (((*data)>>read)<<buffered);
                    if (read >= buffered) {
                        ++data;
                        buffered += 64-read;
                        read = 0;
                    } else {
                        read += 64-buffered;
                        buffered = 64;
                    };
                    if (len_1_len > buffered) {
                        w |= (((*data)>>read)<<buffered);
                        if (read >= buffered) {
                            ++data;
                            buffered += 64-read;
                            read = 0;
                        } else {
                            read += 64-buffered;
                            buffered = 64;
                        };
                    }
                }
//				assert(len_1_len <= buffered);
                uint16_t len_1 = (w&bits::lo_set[len_1_len]) + (1ULL << len_1_len) - 1;
                buffered -= len_1_len;
                w >>= len_1_len;
                if (len_1 > buffered) {// buffer is not full
                    w |= (((*data)>>read)<<buffered);
                    if (read >= buffered) {
                        ++data;
                        buffered += 64-read;
                        read = 0;
                    } else {
                        read += 64-buffered;
                        buffered = 64;
                    };
                    if (len_1 > buffered) {
                        w |= (((*data)>>read)<<buffered);
                        if (read >= buffered) {
                            ++data;
                            buffered += 64-read;
                            read = 0;
                        } else {
                            read += 64-buffered;
                            buffered = 64;
                        };
                    }
                }
//				if( len_1 > buffered ){
//					std::cerr<<"len_1="<<len_1<<" buffered = "<<buffered<<std::endl;
//				}
//				assert(len_1 <= buffered);
                value	+= (w&bits::lo_set[len_1]) + (len_1<64) * (1ULL << (len_1));
                buffered -= len_1;
                if (len_1 < 64) {
                    w >>= len_1;
                } else {
                    w = 0;
                }
                ++i;
                if (i==n)
                    return value;
                if (buffered>=16)
                    goto begin_decode;
            } else {
                value += (psum&0x0000FFFF);
                i += ((psum>>16)&0x00FF);
                if (i==n)
                    return value;
                buffered -= (psum>>24);
                w >>= (psum>>24);
                if (buffered>=16)
                    goto begin_decode;
            }
        }
//		std::cerr<<i<<" / "<<n<<std::endl;
    };
//	std::cerr<<value<<std::endl;
    return value;
}

template<class t_const_uint64_iter>
uint64_t elias_delta::decode_prefix_sum(t_const_uint64_iter data, const size_type start_idx, size_type n)
{
    if (n==0)
        return 0;
    data += (start_idx >> 6);
    uint64_t value = 0;
    size_type i = 0;
    uint8_t offset = start_idx & 0x3F;

    if (n < 24) {
        if (n + offset <= 64) {
            if (((*data >> offset)&bits::lo_set[n]) == bits::lo_set[n])
                return n;
        } else { // n+offset > 64
            if ((*data >> offset) == bits::lo_set[64-offset]) {// all bits are set to 1
                value = 64-offset;
                ++data;
                n -= (64-offset);
                offset = 0;
                while (n >= 64) {
                    if (*data == 0xFFFFFFFFFFFFFFFFULL) {
                        value += 64;
                        ++data;
                        n -= 64;
                    } else {
                        uint8_t temp = bits::lo(~(*data));
                        value += temp;
                        n -= temp;
                        offset = temp;
                        goto start_decoding;
                    }
                }
                // 0 <= n <= 64
                if ((*data&bits::lo_set[n]) == bits::lo_set[n])
                    return value + n;
            }
        }
    }

start_decoding:

    while (i < n) {// while not all values are decoded
        // n-i values to decode
        if (((*data>>offset)&0xF)==0xF) {
            uint8_t maxdecode = n-i > 63 ? 63 : n-i;
            uint8_t rbp = bits::lo(~bits::read_int(data, offset,maxdecode));
            i += rbp;
            value += rbp;
            if (rbp+offset>=64) {
                ++data;
                offset = (rbp+offset)&0x3F;
            } else {
                offset += rbp;
            }
            if (rbp == maxdecode)
                continue;
        }
        while (i < n) {
            uint32_t psum = elias_delta::data.prefixsum[bits::read_int(data, offset, 16)];
//			if( psum == 0 or i+((psum>>16)&0x00FF) > n ){ // value does not fit in 16 bits
            if (psum == 0) { // value does not fit in 16 bits
                goto decode_single;
            } else if (i+((psum>>16)&0x00FF) > n) { // decoded too much
                if (n-i <= 8) {
                    psum = elias_delta::data.prefixsum_8bit[bits::read_int(data, offset, 8) | ((n-i-1)<<8)];
                    if (psum > 0) {
                        value += (psum&0xF);
                        i += ((psum>>4)&0xF);
                        offset += (psum>>8);
                        if (offset>=64) {
                            offset&=0x3F;
                            ++data;
                        }
                    }
                }
                break;
            } else {
                value += (psum&0x0000FFFF);
                i += ((psum>>16)&0x00FF);
                offset += (psum>>24);
                if (offset>=64) {
                    offset&=0x3F;
                    ++data;
                }
            }
        }
        if (i<n) {
decode_single:
            i++;
            uint16_t len_1_len = bits::read_unary_and_move(data, offset); // read length of length of x
            uint16_t len_1 	=  bits::read_int_and_move(data, offset, len_1_len) + (1ULL << len_1_len) - 1;
            value	+= bits::read_int_and_move(data, offset, len_1) + (len_1<64) * (1ULL << (len_1));
//			std::cout<<"decode single ("<<len_1_len<<","<<len_1<<","<<value<<")"<<std::endl;
        }
    }
    return value;
}

} // end namespace coder
} // end namespace sdsl
#endif
