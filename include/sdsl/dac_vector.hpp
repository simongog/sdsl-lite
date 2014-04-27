/* sdsl - succinct data structures library
    Copyright (C) 2014 Simon Gog

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
/*! \file dac_vector.hpp
   \brief dac_vector.hpp contains a vector which stores the values with variable length codes.
   \author Simon Gog
*/
#ifndef SDSL_DAC_VECTOR
#define SDSL_DAC_VECTOR

#include "int_vector.hpp"
#include "rank_support_v5.hpp"
#include "iterators.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A generic immutable space-saving vector class for unsigned integers.
/*! The values of a dac_vector are immutable after the constructor call.
 *  The ,,escaping'' technique is used to encode values.
 *  This is defined as follows (see [1]):
 *  A k-bit integer is split into \f$K=\lceil k/(b-1)\rceil\f$ bits each and
 *  encoded into \f$K\f$ blocks of \f$ b \f$ bits each. All but the last block
 *  are marked with by a 1 in the most significant bit. Escaping with b=8 is
 *  also known as vbyte-coding (see [2]). A experimental study of using escaping
 *  for the LCP array is given in [3].
 *  \par Time complexity
 *        - \f$\Order{\log n/b}\f$ worst case, where b is the number of bits
          in a block
 *  \par References
 *       [1] F. Transier and P. Sanders: ,,Engineering Basic Search Algorithms
 *           of an In-Memory Text Search Engine'', ACM Transactions on
 *           Information Systems, Vol. 29, No.1, Article 2, 2010
 *       [2] H.E. Williams and J. Zobel: ,,Compressing integers for fast file
 *           access'', Computing Journal Vol 43, No.3, 1999
 *       [3] N. Brisboa, S. Ladra, G. Navarro: ,,Directly addressable variable-
 *           length codes'', Proceedings of SPIRE 2009.
 *
 * \tparam t_b    Split block size.
 * \tparam t_rank Rank structure to navigate between the different levels.

 */
template<uint8_t t_b = 4,
         typename t_rank = rank_support_v5<>
         >
class dac_vector
{
    private:
        static_assert(t_b > 0 , "dac_vector: t_b has to be larger than 0");
        static_assert(t_b < 64, "dac_vector: t_b has to be smaller than 64");
    public:
        typedef typename int_vector<>::value_type        value_type;
        typedef random_access_const_iterator<dac_vector> const_iterator;
        typedef const_iterator                           iterator;
        typedef const value_type                         const_reference;
        typedef const_reference                          reference;
        typedef const_reference*                         pointer;
        typedef const pointer                            const_pointer;
        typedef int_vector<>::size_type                  size_type;
        typedef ptrdiff_t                                difference_type;
        typedef t_rank                                   rank_support_type;
        typedef iv_tag                                   index_category;
    private:
        int_vector<t_b>   m_data;           // block data for every level
        bit_vector        m_overflow;       // mark non-end bytes
        rank_support_type m_overflow_rank;  // rank for m_overflow
        int_vector<64>    m_level_pointer_and_rank;
        uint8_t           m_max_level;      // maximum level < (log n)/b+1

        void copy(const dac_vector& v) {
            m_data                   = v.m_data;
            m_overflow               = v.m_overflow;
            m_overflow_rank          = v.m_overflow_rank;
            m_overflow_rank.set_vector(&m_overflow);
            m_level_pointer_and_rank = v.m_level_pointer_and_rank;
            m_max_level              = v.m_max_level;
        }

    public:
        dac_vector() {
            m_level_pointer_and_rank = int_vector<64>(4,0);
        }

        dac_vector(const dac_vector& v) {
            copy(v);
        }

        dac_vector(dac_vector&& v) {
            *this = std::move(v);
        }
        dac_vector& operator=(const dac_vector& v) {
            if (this != &v) {
                copy(v);
            }
            return *this;
        }

        dac_vector& operator=(dac_vector&& v) {
            if (this != &v) {
                m_data                   = std::move(v.m_data);
                m_overflow               = std::move(v.m_overflow);
                m_overflow_rank          = std::move(v.m_overflow_rank);
                m_overflow_rank.set_vector(&m_overflow);
                m_level_pointer_and_rank = std::move(v.m_level_pointer_and_rank);
                m_max_level              = std::move(v.m_max_level);
            }
            return *this;
        }

        //! Constructor for a Container of unsigned integers.
        /*! \param c A container of unsigned integers.
            \pre No two adjacent values should be equal.
          */
        template<class Container>
        dac_vector(const Container& c);

        //! Constructor for an int_vector_buffer of unsigned integers.
        template<uint8_t int_width>
        dac_vector(int_vector_buffer<int_width>& v_buf);

        //! The number of elements in the dac_vector.
        size_type size()const {
            return m_level_pointer_and_rank[2];
        }
        //! Return the largest size that this container can ever have.
        static size_type max_size() {
            return int_vector<>::max_size()/2;
        }

        //!    Returns if the dac_vector is empty.
        bool empty() const {
            return 0 == m_level_pointer_and_rank[2];
        }

        //! Swap method for dac_vector
        void swap(dac_vector& v) {
            m_data.swap(v.m_data);
            m_overflow.swap(v.m_overflow);
            util::swap_support(m_overflow_rank, v.m_overflow_rank,
                               &m_overflow, &(v.m_overflow));

            m_level_pointer_and_rank.swap(v.m_level_pointer_and_rank);
            std::swap(m_max_level, v.m_max_level);
        }

        //! Iterator that points to the first element of the dac_vector.
        const const_iterator begin()const {
            return const_iterator(this, 0);
        }


        //! Iterator that points to the position after the last element of the dac_vector.
        const const_iterator end()const {
            return const_iterator(this, size());
        }

        //! []-operator
        value_type operator[](size_type i)const {
            uint8_t level = 1;
            uint8_t offset = t_b;
            size_type result = m_data[i];
            const uint64_t* p = m_level_pointer_and_rank.data();
            uint64_t ppi = (*p)+i;
            while (level < m_max_level and m_overflow[ppi]) {
                p += 2;
                ppi = *p + (m_overflow_rank(ppi) - *(p-1));
                result |= (m_data[ppi] << (offset));
                ++level;
                offset += t_b;
            }
            return result;
        }

        //! Serializes the dac_vector to a stream.
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const;

        //! Load from a stream.
        void load(std::istream& in) {
            m_data.load(in);
            m_overflow.load(in);
            m_overflow_rank.load(in, &m_overflow);
            m_level_pointer_and_rank.load(in);
            read_member(m_max_level, in);
        }
};

template<uint8_t t_b, typename t_rank>
template<class Container>
dac_vector<t_b, t_rank>::dac_vector(const Container& c)
{
//  (1) Count for each level, how many blocks are needed for the representation
//      Running time: \f$ O(n \times \frac{\log n}{b}  \f$
//      Result is sorted in m_level_pointer_and_rank
    size_type n = c.size(), val=0;
    if (n == 0)
        return;
// initialize counter
    auto _size =  std::max(4*bits::hi(2), 2*(((bits::hi(n)+1)+t_b-1) / t_b));
    m_level_pointer_and_rank.resize(_size);
    for (size_type i=0; i < m_level_pointer_and_rank.size(); ++i)
        m_level_pointer_and_rank[i] = 0;
    m_level_pointer_and_rank[0] = n; // level 0 has n entries

    uint8_t level_x_2 = 0;
    for (size_type i=0; i < n; ++i) {
        val=c[i];
        val >>= t_b; // shift value b bits to the right
        level_x_2 = 2;
        while (val) {
            // increase counter for current level by 1
            ++m_level_pointer_and_rank[level_x_2];
            val >>= t_b; // shift value b bits to the right
            level_x_2 += 2; // increase level by 1
        }
    }
//  (2)    Determine maximum level and prefix sums of level counters
    m_max_level = 0;
    size_type sum_blocks = 0, last_block_size=0;
    for (size_type i=0, t=0; i < m_level_pointer_and_rank.size(); i+=2) {
        t = sum_blocks;
        sum_blocks += m_level_pointer_and_rank[i];
        m_level_pointer_and_rank[i] = t;
        if (sum_blocks > t) {
            ++m_max_level;
            last_block_size = sum_blocks - t;
        }
    }
    m_overflow = bit_vector(sum_blocks - last_block_size, 0);
    m_data.resize(sum_blocks);

    assert(last_block_size > 0);

//  (3)    Enter block and overflow data
    int_vector<64> cnt = m_level_pointer_and_rank;
    const uint64_t mask = bits::lo_set[t_b];

    for (size_type i=0, j=0; i < n; ++i) {
        val=c[i];
        j = cnt[0]++;
        m_data[ j ] =  val & mask;
        val >>= t_b; // shift value b bits to the right
        level_x_2 = 2;
        while (val) {
            m_overflow[j] = 1;
            // increase counter for current level by 1
            j = cnt[level_x_2]++;
            m_data[ j ] = val & mask;
            val >>= t_b; // shift value b bits to the right
            level_x_2 += 2; // increase level by 1
        }
    }

//  (4) Initialize rank data structure for m_overflow and precalc rank for
//      pointers
    util::init_support(m_overflow_rank, &m_overflow);
    for (size_type i=0; 2*i < m_level_pointer_and_rank.size() and
         m_level_pointer_and_rank[2*i] < m_overflow.size(); ++i) {
        m_level_pointer_and_rank[2*i+1] = m_overflow_rank(
                                              m_level_pointer_and_rank[2*i]);
    }

}

template<uint8_t t_b, typename t_rank>
template<uint8_t int_width>
dac_vector<t_b, t_rank>::dac_vector(int_vector_buffer<int_width>& v_buf)
{
//  (1) Count for each level, how many blocks are needed for the representation
//      Running time: \f$ O(n \times \frac{\log n}{b}  \f$
//      Result is sorted in m_level_pointer_and_rank
    size_type n = v_buf.size(), val=0;
    if (n == 0)
        return;
// initialize counter
    auto _size =  std::max(4*bits::hi(2), 2*(((bits::hi(n)+1)+t_b-1) / t_b));
    m_level_pointer_and_rank.resize(_size);
    for (size_type i=0; i < m_level_pointer_and_rank.size(); ++i)
        m_level_pointer_and_rank[i] = 0;
    m_level_pointer_and_rank[0] = n; // level 0 has n entries

    uint8_t level_x_2 = 0;
    for (size_type i=0; i < n; ++i) {
        val=v_buf[i];
        val >>= t_b; // shift value b bits to the right
        level_x_2 = 2;
        while (val) {
            // increase counter for current level by 1
            ++m_level_pointer_and_rank[level_x_2];
            val >>= t_b; // shift value b bits to the right
            level_x_2 += 2; // increase level by 1
        }
    }
//  (2)    Determine maximum level and prefix sums of level counters
    m_max_level = 0;
    size_type sum_blocks = 0, last_block_size=0;
    for (size_type i=0, t=0; i < m_level_pointer_and_rank.size(); i+=2) {
        t = sum_blocks;
        sum_blocks += m_level_pointer_and_rank[i];
        m_level_pointer_and_rank[i] = t;
        if (sum_blocks > t) {
            ++m_max_level;
            last_block_size = sum_blocks - t;
        }
    }
    m_overflow = bit_vector(sum_blocks - last_block_size, 0);
    m_data.resize(sum_blocks);

    assert(last_block_size > 0);

//  (3)    Enter block and overflow data
    int_vector<64> cnt = m_level_pointer_and_rank;
    const uint64_t mask = bits::lo_set[t_b];

    for (size_type i=0, j=0; i < n; ++i) {
        val=v_buf[i];
        j = cnt[0]++;
        m_data[ j ] =  val & mask;
        val >>= t_b; // shift value b bits to the right
        level_x_2 = 2;
        while (val) {
            m_overflow[j] = 1;
            // increase counter for current level by 1
            j = cnt[level_x_2]++;
            m_data[ j ] = val & mask;
            val >>= t_b; // shift value b bits to the right
            level_x_2 += 2; // increase level by 1
        }
    }

//  (4) Initialize rank data structure for m_overflow and precalc rank for
//      pointers
    util::init_support(m_overflow_rank, &m_overflow);
    for (size_type i=0; 2*i < m_level_pointer_and_rank.size() and
         m_level_pointer_and_rank[2*i] < m_overflow.size(); ++i) {
        m_level_pointer_and_rank[2*i+1] = m_overflow_rank(
                                              m_level_pointer_and_rank[2*i]);
    }
}

template<uint8_t t_b, typename t_rank>
dac_vector<>::size_type dac_vector<t_b, t_rank>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(
                                     v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_data.serialize(out, child, "data");
    written_bytes += m_overflow.serialize(out, child, "overflow");
    written_bytes += m_overflow_rank.serialize(out, child, "overflow_rank");
    written_bytes += m_level_pointer_and_rank.serialize(out,
                     child, "level_pointer_and_rank");
    written_bytes += write_member(m_max_level, out, child, "max_level");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

} // end namespace sdsl
#endif
