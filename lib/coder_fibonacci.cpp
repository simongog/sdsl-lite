#include "sdsl/coder_fibonacci.hpp"

namespace sdsl
{

namespace coder
{

fibonacci::impl fibonacci::data;

uint64_t fibonacci::decode_prefix_sum(const uint64_t* d, const size_type start_idx, size_type n)
{
    if (n==0)
        return 0;
//	return decode<true,false,int*>(data, start_idx, n);
    d += (start_idx >> 6);
    size_type i = 0;
    int32_t	bits_to_decode = 0;
    uint64_t w = 0, value = 0;
    int16_t buffered = 0, read = start_idx & 0x3F, shift = 0;
    uint16_t temp=0;
    uint64_t carry=0;
    i = bits::cnt11(*d & ~bits::lo_set[read], carry);
    if (i<n) {
        uint64_t oldcarry;
        w = 0;
        do {
            oldcarry = carry;
            i += (temp = bits::cnt11(*(d+(++w)), carry));
        } while (i<n);
        bits_to_decode += ((w-1)<<6) + bits::sel11(*(d+w), n-(i-temp), oldcarry) + 65 - read;
        w = 0;
    } else { // i>=n
        bits_to_decode = bits::sel11(*d >> read, n)+1;
    }
    if (((size_type)bits_to_decode) == n<<1)
        return n;
    if (((size_type)bits_to_decode) == (n<<1)+1)
        return n+1;
    i = 0;
//	while( bits_to_decode > 0 or buffered > 0){// while not all values are decoded
    do {
        while (buffered < 64 and bits_to_decode > 0) {
            w |= (((*d)>>read)<<buffered);
            if (read >= buffered) {
                ++d;
                buffered += 64-read;
                bits_to_decode -= (64-read);
                read = 0;
            } else { // read buffered
                read += 64-buffered;
                bits_to_decode -= (64-buffered);
                buffered = 64;
            }
            if (bits_to_decode < 0) {
                buffered += bits_to_decode;
                w &= bits::lo_set[buffered];
                bits_to_decode = 0;
            }
        }
        if (!i) { // try do decode multiple values
            if ((w&0xFFFFFF)==0xFFFFFF) {
                value += 12;
                w >>= 24;
                buffered -= 24;
                if ((w&0xFFFFFF)==0xFFFFFF) {
                    value += 12;
                    w >>= 24;
                    buffered -= 24;
                }
            }
            do {
                temp = fibonacci::data.fib2bin_16_greedy[w&0xFFFF];
                if ((shift=(temp>>11)) > 0) {
                    value += (temp & 0x7FFULL);
                    w >>= shift;
                    buffered -= shift;
                } else {
                    value += fibonacci::data.fib2bin_0_95[w&0xFFF];
                    w >>= 12;
                    buffered -= 12;
                    i = 1;
                    break;
                }
            } while (buffered>15);
        } else { // i > 0
            value += fibonacci::data.fib2bin_0_95[(i<<12) | (w&0xFFF)];
            shift  = fibonacci::data.fib2bin_shift[w&0x1FFF];
            if (shift > 0) { // if end of decoding
                w >>= shift;
                buffered -= shift;
                i = 0;
            } else { // not end of decoding
                w >>= 12;
                buffered -= 12;
                ++i;
            }
        }
    } while (bits_to_decode > 0 or buffered > 0);
    return value;
}

uint64_t fibonacci::decode_prefix_sum(const uint64_t* d, const size_type start_idx, SDSL_UNUSED const size_type end_idx, size_type n)
{
    return decode_prefix_sum(d, start_idx, n);
}

} // end namespace coder
} // end namespace sdsl
