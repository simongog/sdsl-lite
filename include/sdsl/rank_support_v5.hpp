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
/*! \file rank_support_v5.hpp
    \brief rank_support_v5.hpp contains rank_support_v5.5
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_RANK_SUPPORT_VFIVE
#define INCLUDED_SDSL_RANK_SUPPORT_VFIVE

#include "rank_support.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{

template<uint8_t, uint8_t>
struct rank_support_trait;

//! A class supporting rank queries in constant time.
/*! \par Space complexity
 *  \f$ 0.0625n\f$ bits for a bit vector of length n bits.
 *
 * The superblock size is 2048. Each superblock is subdivided into
 * 2048/(6*64) = 5 blocks (with some bit remaining). So absolute counts for
 * the superblock add 64/2048 bits on top of each supported bit. Since the
 * first of the 6 relative count values is 0, we can fit the remaining 5
 * (each of width log(2048)=11) in a 64 bit word. The relative counts add
 * another 64/2048 bits bits on top of each supported bit. In total this
 * results in 128/2048= 6.25% overhead.
 *
 * \tparam t_b       Bit pattern `0`,`1`,`10`,`01` which should be ranked.
 * \tparam t_pat_len Length of the bit pattern.
 *
 * @ingroup rank_support_group
 */
template<uint8_t t_b=1, uint8_t t_pat_len=1>
class rank_support_v5 : public rank_support
{
    private:
        static_assert(t_b == 1u or t_b == 0u or t_b == 10u or t_b == 11u, "rank_support_v5: bit pattern must be `0`,`1`,`10` or `01` or `11`");
        static_assert(t_pat_len == 1u or t_pat_len == 2u , "rank_support_v5: bit pattern length must be 1 or 2");
    public:
        typedef bit_vector bit_vector_type;
        typedef rank_support_trait<t_b, t_pat_len>  trait_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = t_pat_len };
    private:
//      basic block for interleaved storage of superblockrank and blockrank
        int_vector<64> m_basic_block;
    public:
        explicit rank_support_v5(const bit_vector* v = nullptr) {
            set_vector(v);
            if (v == nullptr) {
                return;
            } else if (v->empty()) {
                m_basic_block = int_vector<64>(2,0);
                return;
            }
            size_type basic_block_size = ((v->capacity() >> 11)+1)<<1;
            m_basic_block.resize(basic_block_size);   // resize structure for basic_blocks
            if (m_basic_block.empty())
                return;
            const uint64_t* data = m_v->data();
            size_type i, j=0;
            m_basic_block[0] = m_basic_block[1] = 0;

            uint64_t carry = trait_type::init_carry();
            uint64_t sum   = trait_type::args_in_the_word(*data, carry);
            uint64_t second_level_cnt = 0;
            uint64_t cnt_words=1;
            for (i = 1; i < (m_v->capacity()>>6) ; ++i, ++cnt_words) {
                if (cnt_words == 32) {
                    j += 2;
                    m_basic_block[j-1] = second_level_cnt;
                    m_basic_block[j]     = m_basic_block[j-2] + sum;
                    second_level_cnt = sum = cnt_words = 0;
                } else if ((cnt_words%6)==0) {
                    // pack the prefix sum for each 6x64bit block into the second_level_cnt
                    second_level_cnt |= sum<<(60-12*(cnt_words/6));//  48, 36, 24, 12, 0
                }
                sum += trait_type::args_in_the_word(*(++data), carry);
            }

            if ((cnt_words%6)==0) {
                second_level_cnt |= sum<<(60-12*(cnt_words/6));
            }
            if (cnt_words == 32) {
                j += 2;
                m_basic_block[j-1] = second_level_cnt;
                m_basic_block[j]   = m_basic_block[j-2] + sum;
                m_basic_block[j+1] = 0;
            } else {
                m_basic_block[j+1] = second_level_cnt;
            }
        }

        rank_support_v5(const rank_support_v5&) = default;
        rank_support_v5(rank_support_v5&&) = default;
        rank_support_v5& operator=(const rank_support_v5&) = default;
        rank_support_v5& operator=(rank_support_v5&&) = default;

        size_type rank(size_type idx) const {
            assert(m_v != nullptr);
            assert(idx <= m_v->size());
            const uint64_t* p = m_basic_block.data()
                                + ((idx>>10)&0xFFFFFFFFFFFFFFFEULL);// (idx/2048)*2
//                     ( prefix sum of the 6x64bit blocks | (idx%2048)/(64*6) )
            size_type result = *p
                               + ((*(p+1)>>(60-12*((idx&0x7FF)/(64*6))))&0x7FFULL)
                               + trait_type::word_rank(m_v->data(), idx);
            idx -= (idx&0x3F);
            uint8_t to_do = ((idx>>6)&0x1FULL)%6;
            --idx;
            while (to_do) {
                result += trait_type::full_word_rank(m_v->data(), idx);
                --to_do;
                idx-=64;
            }
            return result;
        }

        inline size_type operator()(size_type idx)const {
            return rank(idx);
        }
        size_type size()const {
            return m_v->size();
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            size_type written_bytes = 0;
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            written_bytes += m_basic_block.serialize(out, child, "cumulative_counts");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in, const bit_vector* v=nullptr) {
            set_vector(v);
            assert(m_v != nullptr); // supported bit vector should be known
            m_basic_block.load(in);
        }

        void set_vector(const bit_vector* v=nullptr) {
            m_v = v;
        }
        //! swap Operator
        void swap(rank_support_v5& rs) {
            if (this != &rs) { // if rs and _this_ are not the same object
                m_basic_block.swap(rs.m_basic_block);
            }
        }
};

}// end namespace sds

#endif // end file
