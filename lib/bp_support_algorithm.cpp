#include "sdsl/bp_support_algorithm.hpp"

namespace sdsl
{
excess::impl excess::data;

bit_vector
calculate_pioneers_bitmap(const bit_vector& bp, uint64_t block_size)
{
    bit_vector pioneer_bitmap(bp.size(), 0);

    std::stack<uint64_t> opening_parenthesis;
    uint64_t blocks = (bp.size()+block_size-1)/block_size;
    // calculate positions of findclose and findopen pioneers
    for (uint64_t block_nr = 0; block_nr < blocks; ++block_nr) {
        std::map<uint64_t, uint64_t> block_and_position; // for find_open and find_close
        std::map<uint64_t, uint64_t> matching_position;  // for find_open and find_close
        for (uint64_t i=0, j=block_nr*block_size; i < block_size and j < bp.size(); ++i, ++j) {
            if (bp[j]) {//opening parenthesis
                opening_parenthesis.push(j);
            } else { // closing parenthesis
                uint64_t position = opening_parenthesis.top();
                uint64_t blockpos = position/block_size;
                opening_parenthesis.pop();
                block_and_position[blockpos] = position;
                matching_position[blockpos]  = j; // greatest j is pioneer
            }
        }
        for (std::map<uint64_t, uint64_t>::const_iterator it = block_and_position.begin(),
             end = block_and_position.end(),
             mit = matching_position.begin(); it != end and it->first != block_nr; ++it, ++mit) {
            // opening and closing pioneers are symmetric
            pioneer_bitmap[it->second] = 1;
            pioneer_bitmap[mit->second] = 1;
        }
    }
    // assert that the sequence is balanced
    assert(opening_parenthesis.empty());
    return pioneer_bitmap;
}

bit_vector
calculate_pioneers_bitmap_succinct(const bit_vector& bp, uint64_t block_size)
{
    bit_vector pioneer_bitmap(bp.size(), 0);

    sorted_stack_support opening_parenthesis(bp.size());
    uint64_t cur_pioneer_block = 0, last_start = 0, last_j = 0, cur_block=0, first_index_in_block=0;
    // calculate positions of findclose and findopen pioneers
    for (uint64_t j=0, new_block=block_size; j < bp.size(); ++j, --new_block) {
        if (!(new_block)) {
            cur_pioneer_block = j/block_size;
            ++cur_block;
            first_index_in_block = j;
            new_block = block_size;
        }

        if (bp[j]) { // opening parenthesis
            if (/*j < bp.size() is not necessary as the last parenthesis is always a closing one*/
                new_block>1 and !bp[j+1]) {
                ++j; --new_block;
                continue;
            }
            opening_parenthesis.push(j);
        } else {
            assert(!opening_parenthesis.empty());
            uint64_t start = opening_parenthesis.top();
            opening_parenthesis.pop();
            if (start < first_index_in_block) {
                if ((start/block_size)==cur_pioneer_block) {
                    pioneer_bitmap[last_start] = pioneer_bitmap[last_j] = 0; // override false pioneer
                }
                pioneer_bitmap[start] = pioneer_bitmap[j] = 1;
                cur_pioneer_block	  = start/block_size;
                last_start			  = start;
                last_j			 	  = j;
            }
        }
    }
    // assert that the sequence is balanced
    assert(opening_parenthesis.empty());
    return pioneer_bitmap;
}

uint64_t
near_find_close(const bit_vector& bp, const uint64_t i,
                const uint64_t block_size)
{
    typedef bit_vector::difference_type difference_type;
    difference_type excess=1;

    const uint64_t end = ((i+1)/block_size+1)*block_size;
    const uint64_t l = (((i+1)+7)/8)*8;
    const uint64_t r = (end/8)*8;
    for (uint64_t j=i+1; j < std::min(end,l); ++j)	{
        if (bp[j])
            ++excess;
        else {
            --excess;
            if (excess == 0) {
                return j;
            }
        }
    }
    const uint64_t* b = bp.data();
    for (uint64_t j=l; j<r; j+=8) {
        if (excess <= 8) {
            assert(excess>0);
            uint32_t x = excess::data.min_match_pos_packed[((*(b+(j>>6)))>>(j&0x3F))&0xFF ];
            uint8_t p = (x >> ((excess-1)<<2))&0xF;
            if (p < 9) {
                return j+p;
            }
        }
        excess += excess::data.word_sum[((*(b+(j>>6)))>>(j&0x3F))&0xFF ];
    }
    for (uint64_t j=std::max(l,r); j < end; ++j)	{
        if (bp[j])
            ++excess;
        else {
            --excess;
            if (excess == 0) {
                return j;
            }
        }
    }
    return i;
}

uint64_t
near_find_closing(const bit_vector& bp, uint64_t i,
                  uint64_t closings,
                  const uint64_t block_size)
{
    typedef bit_vector::difference_type difference_type;
    difference_type excess=0;
    difference_type succ_excess=-closings;

    const uint64_t end = (i/block_size+1)*block_size;
    const uint64_t l = (((i)+7)/8)*8;
    const uint64_t r = (end/8)*8;
    for (uint64_t j=i; j < std::min(end,l); ++j)	{
        if (bp[j])
            ++excess;
        else {
            --excess;
            if (excess == succ_excess) {
                return j;
            }
        }
    }
    const uint64_t* b = bp.data();
    for (uint64_t j=l; j<r; j+=8) {
        if (excess-succ_excess <= 8) {
            uint32_t x = excess::data.min_match_pos_packed[((*(b+(j>>6)))>>(j&0x3F))&0xFF ];
            uint8_t p = (x >> (((excess-succ_excess)-1)<<2))&0xF;
            if (p < 9) {
                return j+p;
            }
        }
        excess += excess::data.word_sum[((*(b+(j>>6)))>>(j&0x3F))&0xFF ];
    }
    for (uint64_t j=std::max(l,r); j < end; ++j)	{
        if (bp[j])
            ++excess;
        else {
            --excess;
            if (excess == succ_excess) {
                return j;
            }
        }
    }
    return i-1;
}

uint64_t
near_fwd_excess(const bit_vector& bp, uint64_t i, bit_vector::difference_type rel, const uint64_t block_size)
{
    typedef bit_vector::difference_type difference_type;
    difference_type excess = rel;

    const uint64_t end = (i/block_size+1)*block_size;
    const uint64_t l = (((i)+7)/8)*8;
    const uint64_t r = (end/8)*8;
    for (uint64_t j=i; j < std::min(end,l); ++j) {
        excess += 1-2*bp[j];
        if (!excess) {
            return j;
        }
    }
    excess += 8;
    const uint64_t* b = bp.data();
    for (uint64_t j=l; j < r; j+=8) {
        if (excess >= 0 and  excess <= 16) {
            uint32_t x = excess::data.near_fwd_pos[(excess<<8) + (((*(b+(j>>6)))>>(j&0x3F))&0xFF) ];
            if (x < 8) {
                return j+x;
            }
        }
        excess -= excess::data.word_sum[((*(b+(j>>6)))>>(j&0x3F))&0xFF ];
    }
    excess -= 8;
    for (uint64_t j=std::max(l,r); j < end; ++j) {
        excess += 1-2*bp[j];
        if (!excess) {
            return j;
        }
    }
    return i-1;
}

uint64_t
near_rmq(const bit_vector& bp, uint64_t l, uint64_t r, bit_vector::difference_type& min_rel_ex)
{
    typedef bit_vector::difference_type difference_type;
    const uint64_t l8 = (((l+1)+7)/8)*8;
    const uint64_t r8 = (r/8)*8;
    difference_type excess = 0;
    difference_type min_pos=l;
    min_rel_ex = 0;
    for (uint64_t j=l+1; j < std::min(l8,r+1); ++j) {
        if (bp[j])
            ++excess;
        else {
            --excess;
            if (excess <= min_rel_ex) {
                min_rel_ex  = excess;
                min_pos 	= j;
            }
        }
    }

    const uint64_t* b = bp.data();
    for (uint64_t j=l8; j < r8; j+=8) {
        int8_t x = excess::data.min[(((*(b+(j>>6)))>>(j&0x3F))&0xFF)];
        if ((excess+x) <= min_rel_ex) {
            min_rel_ex = excess+x;
            min_pos    = j + excess::data.min_pos_max[(((*(b+(j>>6)))>>(j&0x3F))&0xFF)];
        }
        excess += excess::data.word_sum[((*(b+(j>>6)))>>(j&0x3F))&0xFF];
    }
    for (uint64_t j=std::max(l8,r8); j<r+1; ++j) {
        if (bp[j])
            ++excess;
        else {
            --excess;
            if (excess <= min_rel_ex) {
                min_rel_ex  = excess;
                min_pos 	= j;
            }
        }
    }
    return min_pos;
}

uint64_t
near_bwd_excess(const bit_vector& bp, uint64_t i, bit_vector::difference_type rel, const uint64_t block_size)
{
    typedef bit_vector::difference_type difference_type;
    difference_type excess = rel;
    const difference_type begin = ((difference_type)(i)/block_size)*block_size;
    const difference_type r = ((difference_type)(i)/8)*8;
    const difference_type l = ((difference_type)((begin+7)/8))*8;
    for (difference_type j=i+1; j >= /*begin*/std::max(r,begin); --j) {
        if (bp[j])
            ++excess;
        else
            --excess;
        if (!excess) return j-1;
    }

    excess += 8;
    const uint64_t* b = bp.data();
    for (difference_type j=r-8; j >= l; j-=8) {
        if (excess >= 0 and  excess <= 16) {
            uint32_t x = excess::data.near_bwd_pos[(excess<<8) + (((*(b+(j>>6)))>>(j&0x3F))&0xFF)];
            if (x < 8) {
                return j+x-1;
            }
        }
        excess += excess::data.word_sum[((*(b+(j>>6)))>>(j&0x3F))&0xFF];
    }
    excess -= 8;
    for (difference_type j=std::min(l, r); j > begin; --j) {
        if (bp[j])
            ++excess;
        else
            --excess;
        if (!excess) return j-1;
    }
    if (0==begin and -1==rel) {
        return -1;
    }
    return i+1;
}

uint64_t
near_find_open(const bit_vector& bp, uint64_t i, const uint64_t block_size)
{
    typedef bit_vector::difference_type difference_type;
    difference_type excess = -1;
    const difference_type begin = ((difference_type)(i-1)/block_size)*block_size;
    const difference_type r = ((difference_type)(i-1)/8)*8;
    const difference_type l = ((difference_type)((begin+7)/8))*8;
    for (difference_type j=i-1; j >= std::max(r,begin); --j) {
        if (bp[j]) {
            if (++excess == 0) {
                return j;
            }
        } else
            --excess;
    }
    const uint64_t* b = bp.data();
    for (difference_type j=r-8; j >= l; j-=8) {
        if (excess >= -8) {
            assert(excess<0);
            uint32_t x = excess::data.max_match_pos_packed[((*(b+(j>>6)))>>(j&0x3F))&0xFF];
            uint8_t p = (x >> ((-excess-1)<<2))&0xF;
            if (p < 9) {
                return j+p;
            }
        }
        excess += excess::data.word_sum[((*(b+(j>>6)))>>(j&0x3F))&0xFF];
    }
    for (difference_type j=std::min(l, r)-1; j >= begin; --j) {
        if (bp[j]) {
            if (++excess == 0) {
                return j;
            }
        } else
            --excess;
    }
    return i;
}


uint64_t
near_find_opening(const bit_vector& bp, uint64_t i, const uint64_t openings,const uint64_t block_size)
{
    typedef bit_vector::difference_type difference_type;
    difference_type excess = 0;
    difference_type succ_excess = openings;

    const difference_type begin = ((difference_type)(i)/block_size)*block_size;
    const difference_type r = ((difference_type)(i)/8)*8;
    const difference_type l = ((difference_type)((begin+7)/8))*8;
    for (difference_type j=i; j >= std::max(r,begin); --j) {
        if (bp[j]) {
            if (++excess == succ_excess) {
                return j;
            }
        } else
            --excess;
    }
    const uint64_t* b = bp.data();
    for (difference_type j=r-8; j >= l; j-=8) {
        if (succ_excess-excess <= 8) {
            assert(succ_excess-excess>0);
            uint32_t x = excess::data.max_match_pos_packed[((*(b+(j>>6)))>>(j&0x3F))&0xFF];
            uint8_t p = (x >> ((succ_excess-excess-1)<<2))&0xF;
            if (p < 9) {
                return j+p;
            }
        }
        excess += excess::data.word_sum[((*(b+(j>>6)))>>(j&0x3F))&0xFF];
    }
    for (difference_type j=std::min(l, r)-1; j >= begin; --j) {
        if (bp[j]) {
            if (++excess == succ_excess) {
                return j;
            }
        } else
            --excess;
    }
    return i+1;
}

uint64_t
near_enclose(const bit_vector& bp, uint64_t i, const uint64_t block_size)
{
    uint64_t opening_parentheses = 1;
    for (uint64_t j=i; j+block_size-1 > i and j>0; --j) {
        if (bp[j-1]) {
            ++opening_parentheses;
            if (opening_parentheses == 2) {
                return j-1;
            }
        } else
            --opening_parentheses;
    }
    return i;
}

uint64_t
near_rmq_open(const bit_vector& bp, const uint64_t begin, const uint64_t end)
{
    typedef bit_vector::difference_type difference_type;
    difference_type min_excess = end-begin+1, ex = 0;
    uint64_t result = end;

    const uint64_t l = ((begin+7)/8)*8;
    const uint64_t r = (end/8)*8;

    for (uint64_t k=begin; k < std::min(end,l); ++k) {
        if (bp[k]) {
            ++ex;
            if (ex <= min_excess) {
                result = k;
                min_excess = ex;
            }
        } else {
            --ex;
        }
    }
    const uint64_t* b = bp.data();
    for (uint64_t k = l; k < r; k+=8) {
        uint16_t x = excess::data.min_open_excess_info[((*(b+(k>>6)))>>(k&0x3F))&0xFF];
        int8_t ones = (x>>12);
        if (ones) {
            int8_t min_ex = (x&0xFF)-8;
            if (ex+min_ex <= min_excess) {
                result = k + ((x>>8)&0xF);
                min_excess = ex+min_ex;
            }
        }
        ex += ((ones<<1)-8);
    }
    for (uint64_t k=std::max(r,l); k < end; ++k) {
        if (bp[k]) {
            ++ex;
            if (ex <= min_excess) {
                result = k;
                min_excess = ex;
            }
        } else {
            --ex;
        }
    }
    if (min_excess <= ex)
        return result;
    return end;
}
}
