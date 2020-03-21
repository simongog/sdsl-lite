#pragma once

#include "sdsl/bits.hpp"
#include "sdsl/int_vector.hpp"

namespace bit_magic
{

void print_word(const uint64_t w,const char* prefix)
{
    sdsl::bit_vector bv(64);
    bv.data()[0] = w;
    std::cout << prefix << bv << std::endl;
}

inline uint64_t next0(const uint64_t* word,uint64_t idx)
{
    word += (idx>>6);
    auto masked_inverse_word = ~(*word | sdsl::bits::lo_set[(idx&0x3F)+1]);
    if (masked_inverse_word) {
        return (idx & ~((size_t)0x3F)) + sdsl::bits::lo(masked_inverse_word);
    }
    idx = (idx & ~((size_t)0x3F)) + 64;
    ++word;
    while (*word==0xFFFFFFFFFFFFFFFFULL) {
        idx += 64;
        ++word;
    }
    return idx + sdsl::bits::lo(~(*word));
}

inline uint64_t next_Xth_zero(const uint64_t* word,uint64_t idx,uint64_t x)
{
    if (x==1) return next0(word,idx);
    word += (idx>>6);
    auto masked_inverse_word = ~(*word | sdsl::bits::lo_set[(idx&0x3F)+1]);
    auto zero_cnt = sdsl::bits::cnt(masked_inverse_word);
    if (zero_cnt >= x) {
        return (idx & ~((size_t)0x3F)) +  sdsl::bits::sel(masked_inverse_word,x);
    }
    idx = (idx & ~((size_t)0x3F)) + 64;
    x -= zero_cnt;
    ++word;
    zero_cnt = sdsl::bits::cnt(~*word);
    while (x > zero_cnt) {
        ++word;
        idx += 64;
        x -= zero_cnt;
        zero_cnt = sdsl::bits::cnt(~*word);
    }
    return idx + sdsl::bits::sel(~*word,x);
}

inline uint64_t next_Xth_one(const uint64_t* word,uint64_t idx,uint64_t x)
{
    word += (idx>>6);
    auto masked_word = *word & ~sdsl::bits::lo_set[(idx&0x3F)+1];
    auto one_cnt = sdsl::bits::cnt(masked_word);
    if (one_cnt >= x) {
        return (idx & ~((size_t)0x3F)) +  sdsl::bits::sel(masked_word,x);
    }
    idx = (idx & ~((size_t)0x3F)) + 64;
    x -= one_cnt;
    ++word;
    one_cnt = sdsl::bits::cnt(*word);
    while (x > one_cnt) {
        ++word;
        idx += 64;
        x -= one_cnt;
        one_cnt = sdsl::bits::cnt(*word);
    }
    return idx + sdsl::bits::sel(*word,x);;
}

}
