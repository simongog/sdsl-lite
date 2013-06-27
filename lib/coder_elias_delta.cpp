#include "sdsl/coder_elias_delta.hpp"

namespace sdsl
{

namespace coder
{

elias_delta::impl elias_delta::data;

uint64_t elias_delta::decode_prefix_sum(const uint64_t* data, const size_type start_idx, const size_type end_idx, size_type n)
{
    if (n==0)
        return 0;
    const uint64_t* lastdata = data + ((end_idx+63) >> 6);
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


uint64_t elias_delta::decode_prefix_sum(const uint64_t* data, const size_type start_idx, size_type n)
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
            }//else{

//		}
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

/*
 // naive decoding method
uint64_t elias_delta::decode_prefix_sum(const uint64_t *data, const size_type start_idx, size_type n){
	data += (start_idx >> 6);
	uint64_t value = 0;
	size_type i = 0;
	uint8_t offset = start_idx & 0x3F;//, maxdecode;
	while( i++ < n ){// while not all values are decoded
			uint16_t len_1_len = bits::read_unary_and_move(data, offset); // read length of length of x
			uint16_t len_1 	=  bits::read_int_and_move(data, offset, len_1_len) + (1ULL << len_1_len) - 1;
			value	+= bits::read_int_and_move(data, offset, len_1) + (1ULL << (len_1));
	}
	return value;
}
*/

} // end namespace sdsl
} // end namespace coder
