/* sdsl - succinct data structures library
    Copyright (C) 2008-2013 Simon Gog

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
/*! \file rank_support_v.hpp
    \brief rank_support_v.hpp contains rank_support_v.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_RANK_SUPPORT_V
#define INCLUDED_SDSL_RANK_SUPPORT_V

#include "rank_support.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{


//! A rank structure proposed by Sebastiano Vigna
/*! \par Space complexity
 *  \f$ 0.25n\f$ for a bit vector of length n bits.
 *
 * The superblock size is 512. Each superblock is subdivided into 512/64 = 8
 * blocks. So absolute counts for the superblock add 64/512 bits on top of each
 * supported bit. Since the first of the 8 relative count values is 0, we can
 * fit the remaining 7 (each of width log(512)=9) in a 64bit word.  The relative
 * counts add another 64/512 bits on top of each supported bit.
 * In total this results in 128/512=25% overhead.
 *
 * \tparam t_b       Bit pattern `0`,`1`,`10`,`01` which should be ranked.
 * \tparam t_pat_len Length of the bit pattern.
 *
 * \par Reference
 *    Sebastiano Vigna:
 *    Broadword Implementation of Rank/Select Queries.
 *    WEA 2008: 154-168
 *
 * @ingroup rank_support_group
 */
template<uint8_t t_b=1, uint8_t t_pat_len=1>
class rank_support_v : public rank_support
{
    private:
        static_assert(t_b == 1u or t_b == 0u or t_b == 10u or t_b == 11, "rank_support_v: bit pattern must be `0`,`1`,`10` or `01`");
        static_assert(t_pat_len == 1u or t_pat_len == 2u , "rank_support_v: bit pattern length must be 1 or 2");
    public:
        typedef bit_vector                          bit_vector_type;
        typedef rank_support_trait<t_b, t_pat_len>  trait_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = t_pat_len };
    private:
        // basic block for interleaved storage of superblockrank and blockrank
        int_vector<64> m_basic_block;
    public:
        explicit rank_support_v(const bit_vector* v = nullptr)
        {
            set_vector(v);
            if (v == nullptr) {
                return;
            } else if (v->empty()) {
                m_basic_block = int_vector<64>(2,0);   // resize structure for basic_blocks
                return;
            }
            size_type basic_block_size = ((v->capacity() >> 9)+1)<<1;
            m_basic_block.resize(basic_block_size);   // resize structure for basic_blocks
            if (m_basic_block.empty())
                return;
            const uint64_t* data = m_v->data();
            size_type i, j=0;
            m_basic_block[0] = m_basic_block[1] = 0;

            uint64_t carry = trait_type::init_carry();
            uint64_t sum   = trait_type::args_in_the_word(*data, carry);
            uint64_t second_level_cnt = 0;
            for (i = 1; i < (m_v->capacity()>>6) ; ++i) {
                if (!(i&0x7)) {// if i%8==0
                    j += 2;
                    m_basic_block[j-1] = second_level_cnt;
                    m_basic_block[j] 	= m_basic_block[j-2] + sum;
                    second_level_cnt = sum = 0;
                } else {
                    second_level_cnt |= sum<<(63-9*(i&0x7));//  54, 45, 36, 27, 18, 9, 0
                }
                sum += trait_type::args_in_the_word(*(++data), carry);
            }
            if (i&0x7) { // if i%8 != 0
                second_level_cnt |= sum << (63-9*(i&0x7));
                m_basic_block[j+1] = second_level_cnt;
            } else { // if i%8 == 0
                j += 2;
                m_basic_block[j-1] = second_level_cnt;
                m_basic_block[j]   = m_basic_block[j-2] + sum;
                m_basic_block[j+1] = 0;
            }
        }

        rank_support_v(const rank_support_v&)  = default;
        rank_support_v(rank_support_v&&) = default;
        rank_support_v& operator=(const rank_support_v&) = default;
        rank_support_v& operator=(rank_support_v&&) = default;


        size_type rank(size_type idx) const
        {
            assert(m_v != nullptr);
            assert(idx <= m_v->size());
            const uint64_t* p = m_basic_block.data()
                                + ((idx>>8)&0xFFFFFFFFFFFFFFFEULL); // (idx/512)*2
            if (idx&0x3F)  // if (idx%64)!=0
                return  *p + ((*(p+1)>>(63 - 9*((idx&0x1FF)>>6)))&0x1FF) +
                        trait_type::word_rank(m_v->data(), idx);
            else
                return  *p + ((*(p+1)>>(63 - 9*((idx&0x1FF)>>6)))&0x1FF);
        }

        inline size_type operator()(size_type idx)const
        {
            return rank(idx);
        }

        size_type size()const
        {
            return m_v->size();
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr,
                            std::string name="")const
        {
            size_type written_bytes = 0;
            structure_tree_node* child = structure_tree::add_child(v, name,
                                         util::class_name(*this));
            written_bytes += m_basic_block.serialize(out, child,
                             "cumulative_counts");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in, const int_vector<1>* v=nullptr)
        {
            set_vector(v);
            m_basic_block.load(in);
        }

        void set_vector(const bit_vector* v=nullptr)
        {
            m_v = v;
        }

        void swap(rank_support_v& rs)
        {
            if (this != &rs) { // if rs and _this_ are not the same object
                m_basic_block.swap(rs.m_basic_block);
            }
        }
};

}// end namespace sds

#endif // end file
