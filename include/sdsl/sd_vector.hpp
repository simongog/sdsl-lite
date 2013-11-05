/* sdsl - succinct data structures library
    Copyright (C) 2012-2013 Simon Gog

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
/*!\file sd_vector.hpp
   \brief sd_vector.hpp contains the sdsl::sd_vector class, and
          classes which support rank and select for sd_vector.
   \author Simon Gog
*/
#ifndef INCLUDED_SDSL_SD_VECTOR
#define INCLUDED_SDSL_SD_VECTOR

#include "int_vector.hpp"
#include "select_support_mcl.hpp"
#include "util.hpp"
#include "iterators.hpp"

//! Namespace for the succinct data structure library
namespace sdsl
{

// forward declaration needed for friend declaration
template<class t_hi_bit_vector= bit_vector,
         class t_select_1     = typename t_hi_bit_vector::select_1_type,
         class t_select_0     = typename t_hi_bit_vector::select_0_type>
class rank_support_sd;  // in sd_vector

// forward declaration needed for friend declaration
template<class t_hi_bit_vector= bit_vector,
         class t_select_1     = typename t_hi_bit_vector::select_1_type,
         class t_select_0     = typename t_hi_bit_vector::select_0_type>
class select_support_sd;  // in sd_vector

//! A bit vector which compresses very sparse populated bit vectors by
// representing the positions of 1 by the Elias-Fano representation for non-decreasing sequences
/*!
 * \par Other implementations of this data structure:
 *  - the sdarray of Okanohara and Sadakane
 *  - Sebastiano Vigna implemented a elias_fano class in this sux library.
 *
 * \par References
 *  - P. Elias: ,,Efficient storage and retrieval by content and address of static files'',
 *              Journal of the ACM, 1974
 *  - R. Fano: ,,On the number of bits required to implement an associative memory''.
 *             Memorandum 61. Computer Structures Group, Project MAC, MIT, 1971
 *  - D. Okanohara, K. Sadakane: ,,Practical Entropy-Compressed Rank/Select Dictionary'',
 *             Proceedings of ALENEX 2007.
 *
 *  \tparam t_hi_bit_vector Type of the bitvector used for the unary decoded differences of
 *                          the high part of the positions of the 1s.
 *  \tparam t_select_1      Type of the select structure which is used to select ones in HI.
 *  \tparam t_select_0      Type of the select structure which is used to select zeros in HI.
 */
template<class t_hi_bit_vector = bit_vector,
         class t_select_1     = typename t_hi_bit_vector::select_1_type,
         class t_select_0     = typename t_hi_bit_vector::select_0_type>
class sd_vector
{
    public:
        typedef bit_vector::size_type                   size_type;
        typedef size_type                               value_type;
        typedef bit_vector::difference_type             difference_type;
        typedef random_access_const_iterator<sd_vector> iterator;
        typedef bv_tag                                  index_category;
        typedef t_select_0                              select_0_support_type;
        typedef t_select_1                              select_1_support_type;

        friend class rank_support_sd<t_hi_bit_vector, select_1_support_type, select_0_support_type>;
        friend class select_support_sd<t_hi_bit_vector, select_1_support_type, select_0_support_type>;

        typedef rank_support_sd<t_hi_bit_vector, select_1_support_type, select_0_support_type> rank_1_type;
        typedef select_support_sd<t_hi_bit_vector, select_1_support_type, select_0_support_type> select_1_type;

        typedef t_hi_bit_vector hi_bit_vector_type;
    private:
        // we need this variables to represent the m ones of the original bit vector of size n
        size_type m_size = 0;  // length of the original bit vector
        uint8_t   m_wl   = 0;  // log n - log m, where n is the length of the original bit vector
        // and m is the number of ones in the bit vector, wl is the abbreviation
        // for ,,width (of) low (part)''

        int_vector<>          m_low;           // vector for the least significant bits of the positions of the m ones
        hi_bit_vector_type    m_high;          // bit vector that represents the most significant bit in permuted order
        select_1_support_type m_high_1_select; // select support for the ones in m_high
        select_0_support_type m_high_0_select; // select support for the zeros in m_high

        void copy(const sd_vector& v) {
            m_size = v.m_size;
            m_wl   = v.m_wl;
            m_low  = v.m_low;
            m_high = v.m_high;
            m_high_1_select = v.m_high_1_select;
            m_high_1_select.set_vector(&m_high);
            m_high_0_select = v.m_high_0_select;
            m_high_0_select.set_vector(&m_high);
        }

    public:
        const hi_bit_vector_type&    high          = m_high;
        const int_vector<>&          low           = m_low;
        const select_1_support_type& high_1_select = m_high_1_select;
        const select_0_support_type& high_0_select = m_high_0_select;

        sd_vector() { }

        sd_vector(const bit_vector& bv) {
            m_size = bv.size();
            size_type m = util::cnt_one_bits(bv);
            uint8_t logm = bits::hi(m)+1;
            uint8_t logn = bits::hi(m_size)+1;
            if (logm == logn) {
                --logm;    // to ensure logn-logm > 0
            }
            m_wl    = logn - logm;
            m_low = int_vector<>(m, 0, m_wl);
            bit_vector high = bit_vector(m + (1ULL<<logm), 0); //
            const uint64_t* bvp = bv.data();
            for (size_type i=0, mm=0,last_high=0,highpos=0; i < (bv.size()+63)/64; ++i, ++bvp) {
                size_type position = 64*i;
                uint64_t  w = *bvp;
                while (w) {  // process bit_vector word by word
                    uint8_t offset = bits::lo(w);
                    w >>= offset;   // note:  w >>= (offset+1) can not be applied for offset=63!
                    position += offset;
                    if (position >= bv.size()) // check that we have not reached the end of the bitvector
                        break;
                    // (1) handle high part
                    size_type cur_high = position >> m_wl;
                    highpos += (cur_high - last_high);   // write cur_high-last_high 0s
                    last_high = cur_high;
                    // (2) handle low part
                    m_low[mm++] = position; // int_vector truncates the most significant logm bits
                    high[highpos++] = 1;     // write 1 for the entry
                    position += 1;
                    w >>= 1;
                }
            }
            util::assign(m_high, high);
            util::init_support(m_high_1_select, &m_high);
            util::init_support(m_high_0_select, &m_high);
        }

        //! Accessing the i-th element of the original bit_vector
        /*! \param i An index i with \f$ 0 \leq i < size()  \f$.
        *   \return The i-th bit of the original bit_vector
        *   \par Time complexity
        *           \f$ \Order{t_{select0} + n/m} \f$, where m equals the number of zeros
        *    \par Remark
         *         The time complexity can be easily improved to
        *            \f$\Order{t_{select0}+\log(n/m)}\f$
        *        by using binary search in the second step.
        */
        value_type operator[](size_type i)const {
            size_type high_val = (i >> (m_wl));
            size_type sel_high = m_high_0_select(high_val + 1);
            size_type rank_low = sel_high - high_val;
            if (0 == rank_low)
                return 0;
            size_type val_low = i & bits::lo_set[ m_wl ]; // extract the low m_wl = log n -log m bits
            --sel_high; --rank_low;
            while (m_high[sel_high] and m_low[rank_low] > val_low) {
                if (sel_high > 0) {
                    --sel_high; --rank_low;
                } else
                    return 0;
            }
            return m_high[sel_high] and m_low[rank_low] == val_low;
        }

        //! Swap method
        void swap(sd_vector& v) {
            if (this != &v) {
                std::swap(m_size, v.m_size);
                std::swap(m_wl, v.m_wl);
                m_low.swap(v.m_low);
                m_high.swap(v.m_high);
                util::swap_support(m_high_1_select, v.m_high_1_select, &m_high, &v.m_high);
                util::swap_support(m_high_0_select, v.m_high_0_select, &m_high, &v.m_high);
            }
        }

        //! Returns the size of the original bit vector.
        size_type size()const {
            return m_size;
        }

        sd_vector& operator=(const sd_vector& v) {
            if (this != &v) {
                copy(v);
            }
            return *this;
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += write_member(m_wl, out, child, "wl");
            written_bytes += m_low.serialize(out, child, "low");
            written_bytes += m_high.serialize(out, child, "high");
            written_bytes += m_high_1_select.serialize(out, child, "high_1_select");
            written_bytes += m_high_0_select.serialize(out, child, "high_0_select");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            read_member(m_size, in);
            read_member(m_wl, in);
            m_low.load(in);
            m_high.load(in);
            m_high_1_select.load(in, &m_high);
            m_high_0_select.load(in, &m_high);
        }

        iterator begin() const {
            return iterator(this, 0);
        }

        iterator end() const {
            return iterator(this, size());
        }
};

//! Rank data structure for sd_vector
/*! \tparam t_hi_bit_vector Type of the bitvector used for the unary decoded differences of
 *                          the high part of the positions of the 1s.
 *  \tparam t_select_1      Type of the select structure which is used to select ones in HI.
 *  \tparam t_select_0      Type of the select structure which is used to select zeros in HI.
 */
template<class t_hi_bit_vector, class t_select_1, class t_select_0>
class rank_support_sd
{
    public:
        typedef bit_vector::size_type size_type;
        typedef sd_vector<t_hi_bit_vector, t_select_1, t_select_0> bit_vector_type;
    private:
        const bit_vector_type* m_v;

    public:

        explicit rank_support_sd(const bit_vector_type* v=nullptr) {
            set_vector(v);
        }

        size_type rank(size_type i)const {
            assert(m_v != nullptr);
            assert(i <= m_v->size());
            // split problem in two parts:
            // (1) find  >=
            size_type high_val = (i >> (m_v->m_wl));
            size_type sel_high = m_v->m_high_0_select(high_val + 1);
            size_type rank_low = sel_high - high_val; //
            if (0 == rank_low)
                return 0;
            size_type val_low = i & bits::lo_set[ m_v->m_wl ];
            // now since rank_low > 0 => sel_high > 0
            do {
                if (!sel_high)
                    return 0;
                --sel_high; --rank_low;
            } while (m_v->m_high[sel_high] and m_v->m_low[rank_low] >= val_low);
            return rank_low+1;
        }

        const size_type operator()(size_type i)const {
            return rank(i);
        }

        const size_type size()const {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr) {
            m_v = v;
        }

        rank_support_sd& operator=(const rank_support_sd& rs) {
            if (this != &rs) {
                set_vector(rs.m_v);
            }
            return *this;
        }

        void swap(rank_support_sd&) { }

        void load(std::istream&, const bit_vector_type* v=nullptr) {
            set_vector(v);
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            return serialize_empty_object(out, v, name, this);
        }
};

//! Select data structure for sd_vector
/*! \tparam t_hi_bit_vector Type of the bitvector used for the unary decoded differences of
 *                          the high part of the positions of the 1s.
 *  \tparam t_select_1      Type of the select structure which is used to select ones in HI.
 *  \tparam t_select_0      Type of the select structure which is used to select zeros in HI.
 */
template<class t_hi_bit_vector, class t_select_1, class t_select_0>
class select_support_sd
{
    public:
        typedef bit_vector::size_type size_type;
        typedef sd_vector<t_hi_bit_vector, t_select_1, t_select_0> bit_vector_type;
    private:
        const bit_vector_type* m_v;

    public:

        explicit select_support_sd(const bit_vector_type* v=nullptr) {
            set_vector(v);
        }

        //! Returns the position of the i-th occurrence in the bit vector.
        size_type select(size_type i)const {
            return m_v->m_low[i-1] +  // lower part of the number
                   ((m_v->m_high_1_select(i) + 1 - i)  << (m_v->m_wl));  // upper part
            //^-number of 0 before the i-th 1-^      ^-shift by m_wl
        }

        const size_type operator()(size_type i)const {
            return select(i);
        }

        const size_type size()const {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr) {
            m_v = v;
        }

        select_support_sd& operator=(const select_support_sd& rs) {
            if (this != &rs) {
                set_vector(rs.m_v);
            }
            return *this;
        }

        void swap(select_support_sd&) { }

        void load(std::istream&, const bit_vector_type* v=nullptr) {
            set_vector(v);
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            return serialize_empty_object(out, v, name, this);
        }
};

} // end namespace
#endif
