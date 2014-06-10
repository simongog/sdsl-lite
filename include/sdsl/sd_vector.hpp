/* sdsl - succinct data structures library
    Copyright (C) 2012-2014 Simon Gog

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
template<uint8_t t_b          = 1,
         class t_hi_bit_vector= bit_vector,
         class t_select_1     = typename t_hi_bit_vector::select_1_type,
         class t_select_0     = typename t_hi_bit_vector::select_0_type>
class rank_support_sd;  // in sd_vector

// forward declaration needed for friend declaration
template<uint8_t t_b          = 1,
         class t_hi_bit_vector= bit_vector,
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

        typedef rank_support_sd<0, t_hi_bit_vector, select_1_support_type, select_0_support_type> rank_0_type;
        typedef rank_support_sd<1, t_hi_bit_vector, select_1_support_type, select_0_support_type> rank_1_type;
        typedef select_support_sd<0, t_hi_bit_vector, select_1_support_type, select_0_support_type> select_0_type;
        typedef select_support_sd<1, t_hi_bit_vector, select_1_support_type, select_0_support_type> select_1_type;

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
        const uint8_t&               wl            = m_wl;
        const hi_bit_vector_type&    high          = m_high;
        const int_vector<>&          low           = m_low;
        const select_1_support_type& high_1_select = m_high_1_select;
        const select_0_support_type& high_0_select = m_high_0_select;

        sd_vector() { }

        sd_vector(const sd_vector& sd) {
            copy(sd);
        }

        sd_vector(sd_vector&& sd) {
            *this = std::move(sd);
        }

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

        //! Get the integer value of the binary string of length len starting at position idx.
        /*! \param idx Starting index of the binary representation of the integer.
         *  \param len Length of the binary representation of the integer. Default value is 64.
         *  \returns The integer value of the binary string of length len starting at position idx.
         *
         *  \pre idx+len-1 in [0..size()-1]
         *  \pre len in [1..64]
         */
        uint64_t get_int(size_type idx, const uint8_t len=64) const {
            uint64_t i = idx+len-1;
            uint64_t high_val = (i >> (m_wl));
            uint64_t sel_high = m_high_0_select(high_val + 1);
            uint64_t rank_low = sel_high - high_val;
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
            uint64_t res = 0;
            while (true) {
                while (!m_high[sel_high]) {
                    if (sel_high > 0 and(high_val << m_wl) >=idx) {
                        --sel_high; --high_val;
                    } else {
                        return res;
                    }
                }
                while (m_high[sel_high]) {
                    uint64_t val = (high_val << m_wl) + m_low[rank_low];
                    if (val >= idx) {
                        res |= 1ULL<<(val-idx);
                    } else {
                        return res;
                    }
                    if (sel_high > 0) {
                        --sel_high; --rank_low;
                    } else {
                        return res;
                    }
                }
            }
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

        sd_vector& operator=(sd_vector&& v) {
            if (this != &v) {
                m_size = v.m_size;
                m_wl   = v.m_wl;
                m_low  = std::move(v.m_low);
                m_high = std::move(v.m_high);
                m_high_1_select = std::move(v.m_high_1_select);
                m_high_1_select.set_vector(&m_high);
                m_high_0_select = std::move(v.m_high_0_select);
                m_high_0_select.set_vector(&m_high);
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

template<uint8_t t_b>
struct rank_support_sd_trait {
    typedef bit_vector::size_type size_type;
    static size_type adjust_rank(size_type r,size_type) {
        return r;
    }
};

template<>
struct rank_support_sd_trait<0> {
    typedef bit_vector::size_type size_type;
    static size_type adjust_rank(size_type r, size_type n) {
        return n - r;
    }
};

//! Rank data structure for sd_vector
/*! \tparam t_b             Bit pattern.
 *  \tparam t_hi_bit_vector Type of the bitvector used for the unary decoded differences of
 *                          the high part of the positions of the 1s.
 *  \tparam t_select_1      Type of the select structure which is used to select ones in HI.
 *  \tparam t_select_0      Type of the select structure which is used to select zeros in HI.
 */
template<uint8_t t_b, class t_hi_bit_vector, class t_select_1, class t_select_0>
class rank_support_sd
{
        static_assert(t_b == 1u or t_b == 0u , "rank_support_sd: bit pattern must be `0` or `1`");
    public:
        typedef bit_vector::size_type size_type;
        typedef sd_vector<t_hi_bit_vector, t_select_1, t_select_0> bit_vector_type;
        enum { bit_pat = t_b };
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
            size_type high_val = (i >> (m_v->wl));
            size_type sel_high = m_v->high_0_select(high_val + 1);
            size_type rank_low = sel_high - high_val; //
            if (0 == rank_low)
                return rank_support_sd_trait<t_b>::adjust_rank(0, i);
            size_type val_low = i & bits::lo_set[ m_v->wl ];
            // now since rank_low > 0 => sel_high > 0
            do {
                if (!sel_high)
                    return rank_support_sd_trait<t_b>::adjust_rank(0, i);
                --sel_high; --rank_low;
            } while (m_v->high[sel_high] and m_v->low[rank_low] >= val_low);
            return rank_support_sd_trait<t_b>::adjust_rank(rank_low+1, i);
        }

        size_type operator()(size_type i)const {
            return rank(i);
        }

        size_type size()const {
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

template<uint8_t t_b, class t_sd_vec>
struct select_support_sd_trait {
    typedef bit_vector::size_type size_type;
    static size_type select(size_type i, const t_sd_vec* v) {
        return v->low[i-1] +  // lower part of the number
               ((v->high_1_select(i) + 1 - i)  << (v->wl));  // upper part
        //^-number of 0 before the i-th 1-^    ^-shift by wl
    }
};

template<class t_sd_vec>
struct select_support_sd_trait<0, t_sd_vec> {
    typedef bit_vector::size_type size_type;
    static size_type select(size_type i, const t_sd_vec* v) {
        auto ones  = v->low.size();
        assert(0 < i and i <= v->size() - ones);
        size_type lb = 1, rb = ones+1;
        size_type r0 = 0;
        size_type pos = (size_type)-1;
        // rb exclusive
        // invariant: rank0(select_1(rb)) >= i
        while (lb < rb) {
            auto mid = lb + (rb-lb)/2;
            auto x = select_support_sd_trait<1, t_sd_vec>::select(mid, v);
            auto rank0 = x + 1 - mid;
            if (rank0 >= i) {
                rb = mid;
            } else {
                r0 = rank0;
                pos = x;
                lb = mid + 1;
            }
        }
        return pos + i - r0;
    }
};

//! Select data structure for sd_vector
/*! \tparam t_b             Bit pattern.
 *  \tparam t_hi_bit_vector Type of the bitvector used for the unary decoded differences of
 *                          the high part of the positions of the 1s.
 *  \tparam t_select_1      Type of the select structure which is used to select ones in HI.
 *  \tparam t_select_0      Type of the select structure which is used to select zeros in HI.
 */
template<uint8_t t_b, class t_hi_bit_vector, class t_select_1, class t_select_0>
class select_support_sd
{
    public:
        typedef bit_vector::size_type size_type;
        typedef sd_vector<t_hi_bit_vector, t_select_1, t_select_0> bit_vector_type;
        enum { bit_pat = t_b };

    private:
        const bit_vector_type* m_v;
    public:

        explicit select_support_sd(const bit_vector_type* v=nullptr) {
            set_vector(v);
        }

        //! Returns the position of the i-th occurrence in the bit vector.
        size_type select(size_type i)const {
            return select_support_sd_trait<t_b, bit_vector_type>::select(i, m_v);
        }

        size_type operator()(size_type i)const {
            return select(i);
        }

        size_type size()const {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr) {
            m_v = v;
        }

        select_support_sd& operator=(const select_support_sd& ss) {
            if (this != &ss) {
                set_vector(ss.m_v);
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


//! Select_0 data structure for sd_vector
/*! \tparam t_sd_vector sd_vector type
 *  \tparam t_rank_1    Rank support for high part of sd_vector
 */
template<typename t_sd_vector=sd_vector<>>
class select_0_support_sd
{
    public:
        typedef bit_vector::size_type size_type;
        typedef t_sd_vector           bit_vector_type;
        using rank_1 = typename t_sd_vector::rank_1_type;
        using sel0_type = typename t_sd_vector::select_0_type;
        typedef bit_vector           y_high_type;
        enum { bit_pat = 0 };

    private:
        const bit_vector_type* m_v;
        int_vector<>           m_pointer;
        int_vector<>           m_rank1;
    public:

        explicit select_0_support_sd(const bit_vector_type* v=nullptr) {
            set_vector(v);
            if (nullptr != m_v) {
                size_type rank_0 = 0; // rank0 in H
                const size_type bs = 1ULL << (m_v->wl);
                size_type z = 0;
                size_type rank1 = 0;// rank1 in H
                size_type zeros = m_v->size() - rank_1(m_v)(m_v->size()); // zeros in B
                m_pointer = int_vector<>(zeros/(64*bs)+1, 0, bits::hi(m_v->high.size()/64)+1);
                m_rank1   = int_vector<>(m_pointer.size(), 0, bits::hi(m_v->high.size())+1);
                uint64_t w=0;
                for (size_type i=0, sel0=1; i < m_v->high.size(); i+=64) {
                    size_type old_rank1 = rank1;
                    w = m_v->high.get_int(i, 64);
                    rank1 += bits::cnt(w);
                    rank_0 = (i+64)-rank1;
                    if (rank1 > 0 and (w>>63)&1) {
                        uint64_t pos = rank_0*bs + m_v->low[rank1-1]; // pos of last one (of previous block in B
                        z = pos + 1 - rank1;
                    } else {
                        z = rank_0*bs  - rank1;
                    }
                    while (sel0 <= z and sel0 <= zeros) {
                        m_pointer[(sel0-1)/(64*bs)] = i/64;
                        m_rank1[(sel0-1)/(64*bs)]   = old_rank1;
                        sel0 += 64*bs;
                    }
                }
            }
        }

        //! Returns the position of the i-th occurrence in the bit vector.
        size_type select(size_type i)const {
            const size_type bs = 1ULL << (m_v->wl);
            size_type j = m_pointer[(i-1)/(64*bs)]*64;// index into m_high
            size_type rank1 = m_rank1[(i-1)/(64*bs)]; // rank_1(j*bs*64) in B
            size_type pos = 0;
            size_type rank0 = 0;

            if (rank1 > 0 and (m_v->high[j-1])&1) {
                pos  = (j-rank1)*bs + m_v->low[rank1-1]; // starting position of current block
                rank0 = pos+1-rank1;
            } else {
                pos  = (j-rank1)*bs;// starting position of current block
                rank0 = pos-rank1;
            }
            uint64_t w = m_v->high.get_int(j, 64);
            do {
                uint64_t _rank1 = rank1 + bits::cnt(w);
                uint64_t _rank0 = 0;
                if (_rank1 > 0 and (w>>63)&1) {
                    pos = (j+64-_rank1)*bs + m_v->low[_rank1-1];
                    _rank0 = pos+1-_rank1;
                } else {
                    pos = (j+64-_rank1)*bs;
                    _rank0 = pos-_rank1;
                }
                if (_rank0 < i) {
                    j+=64;
                    w = m_v->high.get_int(j, 64);
                    rank1 = _rank1;
                } else {
                    break;
                }
            } while (true);
            // invariant i >zeros
            do {
                uint64_t _rank1 = rank1 + bits::lt_cnt[w&0xFFULL];
                uint64_t _rank0 = 0;
                if (_rank1 > 0 and (w>>7)&1) {
                    pos = (j+8-_rank1)*bs + m_v->low[_rank1-1];
                    _rank0 = pos+1-_rank1;
                } else {
                    pos = (j+8-_rank1)*bs;
                    _rank0 = pos-_rank1;
                }
                if (_rank0 < i) {
                    j+=8;
                    w >>= 8;
                    rank1 = _rank1;
                } else {
                    break;
                }
            } while (true);

            do {
                bool b = w&1ULL;
                w >>= 1; // zeros are shifted in
                ++j;
                if (0 == b) {
                    pos = (j-rank1)*bs;
                    size_type zeros = pos-rank1;
                    if (zeros >= i) {
                        pos = pos - (zeros-i) - 1;
                        break;
                    }
                } else {
                    pos = (j-1-rank1)*bs;
                    size_type one_pos = pos + m_v->low[rank1];
                    ++rank1;
                    size_type zeros = one_pos + 1 - rank1;
                    if (zeros >= i) {
                        pos = one_pos - (zeros-i) - 1;
                        break;
                    }
                }
                if (j%64==0) {
                    w = m_v->high.get_int(j,64);
                }
            } while (true);
            return pos;
        }

        size_type operator()(size_type i)const {
            return select(i);
        }

        size_type size()const {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr) {
            m_v = v;
        }

        select_0_support_sd& operator=(const select_0_support_sd& ss) {
            if (this != &ss) {
                m_pointer = ss.m_pointer;
                m_rank1   = ss.m_rank1;
                set_vector(ss.m_v);
            }
            return *this;
        }

        void swap(select_0_support_sd& ss) {
            m_pointer.swap(ss.m_pointer);
            m_rank1.swap(ss.m_rank1);
        }

        void load(std::istream& in, const bit_vector_type* v=nullptr) {
            m_pointer.load(in);
            m_rank1.load(in);
            set_vector(v);
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_pointer.serialize(out, child, "pointer");
            written_bytes += m_rank1.serialize(out, child, "rank1");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

};



} // end namespace
#endif
