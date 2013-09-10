/* sdsl - succinct data structures library
    Copyright (C) 2009-2013 Simon Gog

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
/*! \file bp_support_algorithm.hpp
    \brief bp_support_algorithm.hpp contains algorithms for balanced parentheses sequences.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_BP_SUPPORT_ALGORITHM
#define INCLUDED_SDSL_BP_SUPPORT_ALGORITHM

#include "int_vector.hpp" // for bit_vector
#include <stack> // for calculate_pioneers_bitmap method
#include <map>   // for calculate_pioneers_bitmap method
#include <iomanip>
#include "sorted_stack_support.hpp"


namespace sdsl
{

// This structure contains lookup tables
struct excess {
    static struct impl {
        // Given an excess value x in [-8,8] and a 8-bit
        // word w interpreted as parentheses sequence.
        // near_fwd_pos[(x+8)<<8 | w] contains the minimal position
        // p in [0..7] where the excess value x is reached, or 8
        // if x is not reached in w.
        uint8_t near_fwd_pos[(8-(-8))*256];

        // Given an excess value of x in [-8,8] and a 8-bit
        // word w interpreted as parentheses sequence.
        // near_bwd_pos[(x+8)<<8 | w] contains the maximal position
        // p in [0..7] where the excess value x is reached, or 8
        // if x is not reached in w.
        uint8_t near_bwd_pos[(8-(-8))*256];

        // Given a 8-bit word w. word_sum[w] contains the
        // excess value of w.
        int8_t word_sum[256];

        // Given a 8-bit word w. min[w] contains the
        // minimal excess value in w.
        int8_t min[256];

        // Given a 8-bit word w. min_pos_max[w] contains
        // the maximal position p in w, where min[w] is
        // reached
        int8_t min_pos_max[256];

        // Given an excess value x in [1,8] and a 8-bit
        // word w interpreted as parentheses sequence.
        // min_match_pos_packed[w]:[(x-1)*4,x*4] contains
        // the minimal position, where excess value
        // -x is reached and 9, if there is no such position.
        uint32_t min_match_pos_packed[256];

        // Given an excess value x in [1,8] and a 8-bit
        // word w interpreted as parentheses sequence.
        // max_match_pos_packed[w]:[(x-1)*4,x*4] contains
        // the maximal position, where excess value
        // -x is reached and 9, if there is no such position.
        uint32_t max_match_pos_packed[256];

        // Given a 8-bit word w. x=min_and_info[w] contains
        // the following information.
        // * [0..7] the minimum excess value in w + 8 of an opening parenthesis
        // * [8..11] the maximal position of the minimal excess value
        // * [12..15] the number of ones in the word
        // if w != 0, and 17 for w=0.
        uint16_t min_open_excess_info[256];

        impl() {
            for (int32_t x = -8; x < 8; ++x) {
                for (uint16_t w=0; w < 256; ++w) {
                    uint16_t i = (x+8)<<8|w;
                    near_fwd_pos[i] = 8;
                    int8_t p=0;
                    int8_t excess = 0;
                    do {
                        excess += 1-2*((w&(1<<p))==0);
                        if (excess == x) {
                            near_fwd_pos[i] = p;
                            break;
                        }
                        ++p;
                    } while (p < 8);

                    near_bwd_pos[i] = 8;
                    p = 7;
                    excess = 0;
                    do {
                        excess += 1-2*((w&(1<<p))>0);
                        if (excess == x) {
                            near_bwd_pos[i] = p;
                            break;
                        }
                        --p;
                    } while (p > -1);
                }
            }
            int_vector<> packed_mins(1, 0, 32);
            int_vector<> packed_maxs(1, 0, 32);
            for (uint16_t w=0; w < 256; ++w) {
                int8_t excess = 0;
                int8_t rev_excess = 0;
                int32_t min_excess_of_open = 17;
                int32_t min_excess_of_open_pos = 0;
                uint32_t ones = 0;
                min[w] = 8;
                packed_mins[0] = 0x99999999U;
                packed_maxs[0] = 0x99999999U;
                packed_mins.width(4);
                packed_maxs.width(4);
                for (uint16_t p=0; p<8; ++p) {
                    ones += (w&(1<<p))!=0;
                    excess += 1-2*((w&(1<<p))==0);
                    if (excess <= min[w]) {
                        min[w] = excess;
                        min_pos_max[w] = p;
                    }
                    if (excess < 0 and packed_mins[-excess-1] == 9) {
                        packed_mins[-excess-1] = p;
                    }
                    if (w&(1<<p) and excess+8 <= min_excess_of_open) {
                        min_excess_of_open     = excess+8;
                        min_excess_of_open_pos = p;
                    }
                    rev_excess += 1-2*((w&(1<<(7-p)))>0);
                    if (rev_excess < 0 and packed_maxs[-rev_excess-1] == 9) {
                        packed_maxs[-rev_excess-1] = 7-p;
                    }
                }
                word_sum[w] = excess;
                packed_mins.width(32);
                min_match_pos_packed[w] = packed_mins[0];
                packed_maxs.width(32);
                max_match_pos_packed[w] = packed_maxs[0];
                min_open_excess_info[w] = (min_excess_of_open) |
                                          (min_excess_of_open_pos << 8) |
                                          (ones << 12);
            }
        }
    } data;
};

//! Calculate pioneers as defined in the paper of Geary et al. (CPM 2004)
/*! \param bp             The input balanced parentheses sequence.
 *  \param block_size     Block size.
 *  \param pioneer_bitmap Reference to the  pioneer bitmap, which is calculated.
 *  \par Time complexity
 *       \f$ \Order{n \log n} \f$, where \f$ n=\f$bp.size()
 *  \par Space complexity
 *       \f$ \Order{2n + min(block\_size, \frac{n}{block\_size} )\cdot \log n } \f$
 */
inline void calculate_pioneers_bitmap(const bit_vector& bp, bit_vector::size_type block_size, bit_vector& pioneer_bitmap)
{
    typedef bit_vector::size_type size_type;
    pioneer_bitmap.resize(bp.size());      // resize pioneers bitmap
    util::set_to_value(pioneer_bitmap, 0);  // initialize bitmap with zeros

    std::stack<size_type> opening_parenthesis;
    size_type blocks = (bp.size()+block_size-1)/block_size;
    // calculate positions of findclose and findopen pioneers
    for (size_type block_nr = 0; block_nr < blocks; ++block_nr) {
        std::map<size_type, size_type> block_and_position; // for find_open and find_close
        std::map<size_type, size_type> matching_position;  // for find_open and find_close
        for (size_type i=0, j=block_nr*block_size; i < block_size and j < bp.size(); ++i, ++j) {
            if (bp[j]) {//opening parenthesis
                opening_parenthesis.push(j);
            } else { // closing parenthesis
                size_type position = opening_parenthesis.top();
                size_type blockpos = position/block_size;
                opening_parenthesis.pop();
                block_and_position[blockpos] = position;
                matching_position[blockpos]  = j; // greatest j is pioneer
            }
        }
        for (std::map<size_type, size_type>::const_iterator it = block_and_position.begin(),
             end = block_and_position.end(),
             mit = matching_position.begin(); it != end and it->first != block_nr; ++it, ++mit) {
            // opening and closing pioneers are symmetric
            pioneer_bitmap[it->second] = 1;
            pioneer_bitmap[mit->second] = 1;
        }
    }
    // assert that the sequence is balanced
    assert(opening_parenthesis.empty());
}

//! Calculate pioneers as defined in the paper of Geary et al. (CPM 2004) with few extra space
/*! \param bp The balanced parentheses sequence for that the pioneers should be calculated.
 *  \param block_size Size of the blocks for which the pioneers should be calculated.
 *  \param pioneer_bitmap Reference to the resulting bit_vector.
 *  \par Time complexity
 *       \f$ \Order{n} \f$, where \f$ n=\f$bp.size()
 *  \par Space complexity
 *       \f$ \Order{2n + n} \f$ bits: \f$n\f$ bits for input, \f$n\f$ bits for output, and \f$n\f$ bits for a succinct stack.
 *  \pre The parentheses sequence represented by bp has to be balanced.
 */
template<class size_type>
void calculate_pioneers_bitmap_succinct(const bit_vector& bp, size_type block_size, bit_vector& pioneer_bitmap)
{
    pioneer_bitmap = bit_vector(bp.size(), 0);

    sorted_stack_support opening_parenthesis(bp.size());
    bit_vector::size_type cur_pioneer_block = 0, last_start = 0, last_j = 0, cur_block=0, first_index_in_block=0;
    // calculate positions of findclose and findopen pioneers
    for (bit_vector::size_type j=0, new_block=block_size; j < bp.size(); ++j, --new_block) {
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
            size_type start = opening_parenthesis.top();
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
}

//! Calculates matches (i.e. find_open for a closing parenthesis and find_close for an opening parenthesis)
/*! \param bp A bit_vector representing a balanced parentheses sequence.
 *  \param matches Reference to the result.
 *  \pre bp represents a balanced parentheses sequence.
 *  \par Time complexity
 *       \f$ \Order{n} \f$, where \f$ n=\f$bp.size()
 *  \par Space complexity
 *       \f$ \Order{n + 2n\log n } \f$
 */
template<class int_vector>
void calculate_matches(const bit_vector& bp, int_vector& matches)
{
    typedef bit_vector::size_type size_type;
    matches = int_vector(bp.size(), 0, bits::hi(bp.size())+1);
    std::stack<size_type> opening_parenthesis;
    for (size_type i=0; i < bp.size(); ++i) {
        if (bp[i]) {// opening parenthesis
            opening_parenthesis.push(i);
        } else { // closing parenthesis
            assert(!opening_parenthesis.empty());
            size_type position = opening_parenthesis.top();
            opening_parenthesis.pop();
            matches[i] = position;
            assert(matches[i]==position);
            matches[position] = i;
            assert(matches[position]==i);
        }
    }
    // assert that the sequence is balanced
    assert(opening_parenthesis.empty());
}

//! Calculates enclose answers for a balanced parentheses sequence.
/*! \param bp A bit_vector representing a balanced parentheses sequence.
 *  \param enclose Reference to the result.
 *  \pre bp represents a balanced parentheses sequence.
 *  \par Time complexity
 *       \f$ \Order{n} \f$, where \f$ n=\f$bp.size()
 *  \par Space complexity
 *       \f$ \Order{n + 2n\log n } \f$
 */
template<class int_vector>
void calculate_enclose(const bit_vector& bp, int_vector& enclose)
{
    typedef bit_vector::size_type size_type;
    enclose = int_vector(bp.size(), 0, bits::hi(bp.size())+1);
    std::stack<size_type> opening_parenthesis;
    for (size_type i=0; i < bp.size(); ++i) {
        if (bp[i]) {// opening parenthesis
            if (!opening_parenthesis.empty()) {
                size_type position = opening_parenthesis.top();
                enclose[i] = position;
                assert(enclose[i]==position);
            } else
                enclose[i] = bp.size();
            opening_parenthesis.push(i);
        } else { // closing parenthesis
            size_type position = opening_parenthesis.top();
            enclose[i] = position; // find open answer if i is a closing parenthesis
            opening_parenthesis.pop();
        }
    }
    // assert that the sequence is balanced
    assert(opening_parenthesis.empty());
}

inline bit_vector::size_type near_find_close(const bit_vector& bp, const bit_vector::size_type i, const bit_vector::size_type block_size)
{
    typedef bit_vector::size_type size_type;
    typedef bit_vector::difference_type difference_type;
    difference_type excess=1;

    const size_type end = ((i+1)/block_size+1)*block_size;
    const size_type l = (((i+1)+7)/8)*8;
    const size_type r = (end/8)*8;
    for (size_type j=i+1; j < std::min(end,l); ++j)	{
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
    for (size_type j=l; j<r; j+=8) {
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
    for (size_type j=std::max(l,r); j < end; ++j)	{
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

inline bit_vector::size_type near_find_closing(const bit_vector& bp, bit_vector::size_type i, bit_vector::size_type closings, const bit_vector::size_type block_size)
{
    typedef bit_vector::size_type size_type;
    typedef bit_vector::difference_type difference_type;
    difference_type excess=0;
    difference_type succ_excess=-closings;

    const size_type end = (i/block_size+1)*block_size;
    const size_type l = (((i)+7)/8)*8;
    const size_type r = (end/8)*8;
    for (size_type j=i; j < std::min(end,l); ++j)	{
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
    for (size_type j=l; j<r; j+=8) {
        if (excess-succ_excess <= 8) {
//			assert(excess>0);
            uint32_t x = excess::data.min_match_pos_packed[((*(b+(j>>6)))>>(j&0x3F))&0xFF ];
            uint8_t p = (x >> (((excess-succ_excess)-1)<<2))&0xF;
            if (p < 9) {
                return j+p;
            }
        }
        excess += excess::data.word_sum[((*(b+(j>>6)))>>(j&0x3F))&0xFF ];
    }
    for (size_type j=std::max(l,r); j < end; ++j)	{
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

//! This method searches the minimal parenthesis j, with j>=i, such that excess(j) = excess(i-1)+rel
/*! i>0
 */
inline bit_vector::size_type near_fwd_excess(const bit_vector& bp, bit_vector::size_type i, bit_vector::difference_type rel, const bit_vector::size_type block_size)
{
    typedef bit_vector::size_type size_type;
    typedef bit_vector::difference_type difference_type;
    difference_type excess = rel;

    const size_type end = (i/block_size+1)*block_size;
    const size_type l = (((i)+7)/8)*8;
    const size_type r = (end/8)*8;
    for (size_type j=i; j < std::min(end,l); ++j) {
        if (bp[j])
            --excess;
        else {
            ++excess;
        }
        if (!excess) {
            return j;
        }
    }
    excess += 8;
    const uint64_t* b = bp.data();
    for (size_type j=l; j < r; j+=8) {
        if (excess >= 0 and  excess <= 16) {
            uint32_t x = excess::data.near_fwd_pos[(excess<<8) + (((*(b+(j>>6)))>>(j&0x3F))&0xFF) ];
            if (x < 8) {
                return j+x;
            }
        }
        excess -= excess::data.word_sum[((*(b+(j>>6)))>>(j&0x3F))&0xFF ];
    }
    excess -= 8;
    for (size_type j=std::max(l,r); j < end; ++j) {
        if (bp[j])
            --excess;
        else {
            ++excess;
        }
        if (!excess) {
            return j;
        }
    }
    return i-1;
}

//! Calculate the position with minimal excess value in the interval [l..r].
/*! \param bp The bit_vector which represents the parentheses sequence
 *  \param l  The left border of the interval.
 *	\param r  The right border of the interval.
 *  \param min_rel_ex Reference to the relative minimal excess value with regards to excess(bp[l])
 */
inline bit_vector::size_type near_rmq(const bit_vector& bp, bit_vector::size_type l, bit_vector::size_type r, bit_vector::difference_type& min_rel_ex)
{
    typedef bit_vector::size_type size_type;
    typedef bit_vector::difference_type difference_type;
    const size_type l8 = (((l+1)+7)/8)*8;
    const size_type r8 = (r/8)*8;
    difference_type excess = 0;
    difference_type min_pos=l;
    min_rel_ex = 0;
    for (size_type j=l+1; j < std::min(l8,r+1); ++j) {
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
    for (size_type j=l8; j < r8; j+=8) {
        int8_t x = excess::data.min[(((*(b+(j>>6)))>>(j&0x3F))&0xFF)];
        if ((excess+x) <= min_rel_ex) {
            min_rel_ex = excess+x;
            min_pos    = j + excess::data.min_pos_max[(((*(b+(j>>6)))>>(j&0x3F))&0xFF)];
        }
        excess += excess::data.word_sum[((*(b+(j>>6)))>>(j&0x3F))&0xFF];
    }
    for (size_type j=std::max(l8,r8); j<r+1; ++j) {
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

//! This method searches the maximal parenthesis j, with \f$ j\leq i \f$, such that \f$ excess(j) = excess(i+1)+rel \f$ and
/*! i < bp.size()-1
 */
inline bit_vector::size_type near_bwd_excess(const bit_vector& bp, bit_vector::size_type i, bit_vector::difference_type rel, const bit_vector::size_type block_size)
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




//! Find the near opening parenthesis if it exists.
/*!
 * \param bp bit_vector containing the representation of the balanced parentheses sequence.
 * \param i  Position of the closing parenthesis for which we search the corresponding opening parenthesis.
 * \param block_size Number of entries to search for the corresponding opening parenthesis.
 * \return i if there is no near opening parenthesis, otherwise the position of the near opening parenthesis.
 * \pre We assert that \f$ bp[i]=0 \f$ holds, i.e. there is an closing parenthesis at position i.
 */
inline bit_vector::size_type near_find_open_naive(const bit_vector& bp, bit_vector::size_type i, const bit_vector::size_type block_size)
{
    typedef bit_vector::size_type size_type;
    size_type closing_parentheses = 1;
    for (size_type j=i; j+block_size-1 > i and j>0; --j)	{
        if (bp[j-1]) {
            --closing_parentheses;
            if (closing_parentheses == 0) {
                return j-1;
            }
        } else
            ++closing_parentheses;
    }
    return i;
}


inline bit_vector::size_type near_find_open(const bit_vector& bp, bit_vector::size_type i, const bit_vector::size_type block_size)
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


inline bit_vector::size_type near_find_opening(const bit_vector& bp, bit_vector::size_type i, const bit_vector::size_type openings,const bit_vector::size_type block_size)
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




//! Find the opening parenthesis of the enclosing pair if this parenthesis is near.
/*!
 * \param bp bit_vector containing the representation of the balanced parentheses sequence.
 * \param i Position of the opening parenthesis for which we search the position of the opening parenthesis of the enclosing parentheses pair.
 * \param block_size Number of entries to search for the corresponding opening parenthesis of the enclosing parentheses pair.
 * \return If no near enclose exists return i, otherwise the position of the opening parenthesis of the enclosing pair.
 * \pre We assert that \f$ bp[i]=1 \f$
 */
// TODO: implement a fast version using lookup-tables of size 8
inline bit_vector::size_type near_enclose(const bit_vector& bp, bit_vector::size_type i, const bit_vector::size_type block_size)
{
    typedef bit_vector::size_type size_type;
    size_type opening_parentheses = 1;
    for (size_type j=i; j+block_size-1 > i and j>0; --j) {
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

inline bit_vector::size_type near_rmq_open(const bit_vector& bp, const bit_vector::size_type begin, const bit_vector::size_type end)
{
    typedef bit_vector::size_type size_type;
    typedef bit_vector::difference_type difference_type;
    difference_type min_excess = end-begin+1, ex = 0;
    size_type result = end;

    const size_type l = ((begin+7)/8)*8;
    const size_type r = (end/8)*8;

    for (size_type k=begin; k < std::min(end,l); ++k) {
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
    const uint64_t* b = bp.data();// + (l>>6);
    for (size_type k = l; k < r; k+=8) {
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
    for (size_type k=std::max(r,l); k < end; ++k) {
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

inline bit_vector::size_type near_rmq_open_naive(const bit_vector& bp, const bit_vector::size_type begin, const bit_vector::size_type end)
{
    typedef bit_vector::size_type size_type;
    typedef bit_vector::difference_type difference_type;
    difference_type min_excess = end-begin+1, ex = 0;
    size_type result = end;
    for (size_type k=begin; k<end; ++k) {
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

}// end namespace sdsl

#endif


