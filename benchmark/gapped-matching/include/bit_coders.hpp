#pragma once

#include "sdsl/int_vector.hpp"
#include "bit_streams.hpp"

namespace coder
{

struct elias_gamma {
    template<typename T>
    static inline uint64_t encoded_length(const T& x)
    {
        uint8_t len_1 = sdsl::bits::hi(x);
        return len_1*2+1;
    }
    template<typename T>
    static void encode_check_size(bit_ostream& os,const T& x)
    {
        static_assert(std::numeric_limits<T>::is_signed == false,"can only encode unsigned integers");
        uint8_t len_1 = sdsl::bits::hi(x);
        os.expand_if_needed(len_1*2+1);
        os.put_unary_no_size_check(len_1);
        if (len_1) {
            os.put_int_no_size_check(x, len_1);
        }
    }
    template<typename T>
    static void encode(bit_ostream& os,const T& x)
    {
        static_assert(std::numeric_limits<T>::is_signed == false,"can only encode unsigned integers");
        uint8_t len_1 = sdsl::bits::hi(x);
        os.put_unary_no_size_check(len_1);
        if (len_1) {
            os.put_int_no_size_check(x, len_1);
        }
    }
    template<typename t_itr>
    static void encode(bit_ostream& os,t_itr begin,t_itr end)
    {
        uint64_t bits_required = 0;
        auto tmp = begin;
        while (tmp != end) {
            bits_required += encoded_length(*tmp);
            ++tmp;
        }
        os.expand_if_needed(bits_required);
        tmp = begin;
        while (tmp != end) {
            encode(os,*tmp);
            ++tmp;
        }
    }
    static uint64_t decode(const bit_istream& is)
    {
        uint64_t x;
        auto len = is.get_unary();
        x = 1ULL << len;
        if (len) x |= is.get_int(len);
        return x;
    }
    template<typename t_itr>
    static void decode(const bit_istream& is,t_itr it,size_t n)
    {
        for (size_t i=0; i<n; i++) {
            *it = decode(is);
            ++it;
        }
    }
};

struct elias_delta {
    template<typename T>
    static inline uint64_t encoded_length(const T& x)
    {
        uint8_t len_1 = sdsl::bits::hi(x);
        return len_1 + (sdsl::bits::hi(len_1+1)<<1) + 1;
    }
    template<typename T>
    static void encode_check_size(bit_ostream& os,const T& x)
    {
        static_assert(std::numeric_limits<T>::is_signed == false,"can only encode unsigned integers");
        // (number of sdsl::bits to represent x)
        uint8_t len         = sdsl::bits::hi(x)+1;
        // (number of sdsl::bits to represent the length of x) -1
        uint8_t len_1_len   = sdsl::bits::hi(len);
        uint8_t bits = len  + (len_1_len<<1) + 1;
        os.expand_if_needed(bits);

        // Write unary representation for the length of the length of x
        os.put_unary_no_size_check(len_1_len);
        if (len_1_len) {
            os.put_int_no_size_check(len, len_1_len);
            os.put_int_no_size_check(x, len-1);
        }
    }
    template<typename T>
    static void encode(bit_ostream& os,const T& x)
    {
        static_assert(std::numeric_limits<T>::is_signed == false,"can only encode unsigned integers");
        // (number of sdsl::bits to represent x)
        uint8_t len         = sdsl::bits::hi(x)+1;
        // (number of sdsl::bits to represent the length of x) -1
        uint8_t len_1_len   = sdsl::bits::hi(len);
        // Write unary representation for the length of the length of x
        os.put_unary_no_size_check(len_1_len);
        if (len_1_len) {
            os.put_int_no_size_check(len, len_1_len);
            os.put_int_no_size_check(x, len-1);
        }
    }
    template<typename t_itr>
    static void encode(bit_ostream& os,t_itr begin,t_itr end)
    {
        auto tmp = begin;
        uint64_t len = 0;
        while (tmp != end) {
            len += encoded_length(*tmp);
            ++tmp;
        }
        os.expand_if_needed(len);
        tmp = begin;
        while (tmp != end) {
            encode(os,*tmp);
            ++tmp;
        }
    }
    static uint64_t decode(const bit_istream& is)
    {
        uint64_t x = 1;
        auto len_1_len = is.get_unary();
        if (len_1_len) {
            auto len_1 = is.get_int(len_1_len);
            auto len = len_1 + (1ULL << len_1_len);
            x = is.get_int(len-1);
            x = x + (len-1<64) * (1ULL << (len-1));
        }
        return x;
    }
    template<typename t_itr>
    static void decode(const bit_istream& is,t_itr it,size_t n)
    {
        for (size_t i=0; i<n; i++) {
            *it = decode(is);
            ++it;
        }
    }
};

struct vbyte {
    template<typename T>
    static inline uint64_t encoded_length(const T& x)
    {
        const uint8_t vbyte_len[64] = {1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4,4,
                                       5,5,5,5,5,5,5,6,6,6,6,6,6,6,7,7,7,7,7,7,7,8,8,8,8,8,8,8
                                      };
        return 8*vbyte_len[sdsl::bits::hi(x)+1];
    }
    template<typename T>
    static void encode_check_size(bit_ostream& os,T y)
    {
        static_assert(std::numeric_limits<T>::is_signed == false,"can only encode unsigned integers");
        os.expand_if_needed(encoded_length(y));
        uint64_t x = y;
        uint8_t w = x & 0x7F;
        x >>= 7;
        while (x > 0) {
            w |= 0x80; // mark overflow bit
            os.put_int_no_size_check(w,8);
            w = x & 0x7F;
            x >>= 7;
        }
        os.put_int_no_size_check(w,8);
    }
    template<typename T>
    static void encode(bit_ostream& os,T y)
    {
        static_assert(std::numeric_limits<T>::is_signed == false,"can only encode unsigned integers");
        uint64_t x = y;
        uint8_t w = x & 0x7F;
        x >>= 7;
        while (x > 0) {
            w |= 0x80; // mark overflow bit
            os.put_int_no_size_check(w,8);
            w = x & 0x7F;
            x >>= 7;
        }
        os.put_int_no_size_check(w,8);
    }
    template<typename t_itr>
    static void encode(bit_ostream& os,t_itr begin,t_itr end)
    {
        uint64_t bits_required = 0;
        auto tmp = begin;
        while (tmp != end) {
            bits_required += encoded_length(*tmp);
            ++tmp;
        }
        os.expand_if_needed(bits_required);
        tmp = begin;
        while (tmp != end) {
            encode(os,*tmp);
            ++tmp;
        }
    }
    static uint64_t decode(const bit_istream& is)
    {
        uint64_t ww=0;
        uint8_t w=0;
        uint64_t shift=0;
        do {
            w = is.get_int(8);
            ww |= (((uint64_t)(w&0x7F))<<shift);
            shift += 7;
        } while ((w&0x80) > 0);
        return ww;
    }
    template<typename t_itr>
    static void decode(const bit_istream& is,t_itr it,size_t n)
    {
        for (size_t i=0; i<n; i++) {
            *it = decode(is);
            ++it;
        }
    }
};

template<class t_coder>
struct delta {
    template<typename t_itr>
    static inline uint64_t encoded_length(t_itr begin,t_itr end)
    {
        uint64_t bits_required = t_coder::encoded_length(*begin);
        auto prev = begin;
        auto cur = prev+1;
        while (cur != end) {
            bits_required += t_coder::encoded_length(*cur-*prev);
            prev = cur;
            ++cur;
        }
        return bits_required;
    }
    template<typename t_itr>
    static void encode(bit_ostream& os,t_itr begin,t_itr end)
    {
        uint64_t bits_required = encoded_length(begin,end);
        os.expand_if_needed(bits_required);
        auto prev = begin;
        auto cur = prev+1;
        t_coder::encode(os,*begin);
        while (cur != end) {
            t_coder::encode(os,*cur-*prev);
            prev = cur;
            ++cur;
        }
    }
    template<typename t_itr>
    static void decode(const bit_istream& is,t_itr it,size_t n)
    {
        auto prev = t_coder::decode(is);
        *it++ = prev;
        for (size_t i=1; i<n; i++) {
            *it = prev + t_coder::decode(is);
            prev = *it;
            ++it;
        }
    }
};

}

