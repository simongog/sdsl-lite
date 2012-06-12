/* sdsl - succinct data structures library
    Copyright (C) 2011 Simon Gog

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
/*! \file rrr_vector.hpp
   \brief rrr_vector.hpp contains the sdsl::rrr_vector class, and
          classes which support rank and select for rrr_vector.
   \author Simon Gog, Matthias Petri
*/
#ifndef SDSL_RRR_VECTOR
#define SDSL_RRR_VECTOR

#include "int_vector.hpp"
#include "bitmagic.hpp"
#include "util.hpp"
#include "rrr_helper.hpp" // for binomial helper class
#include <vector>
#include <algorithm> // for next_permutation
#include <iostream>

//! Namespace for the succinct data structure library
namespace sdsl
{


template<uint8_t b=1, uint16_t block_size=15, class wt_type=int_vector<> >  // forward declaration needed for friend declaration
class rrr_rank_support;                // in rrr_vector

template<uint8_t b=1, uint16_t block_size=15, class wt_type=int_vector<> >  // forward declaration needed for friend declaration
class rrr_select_support;                // in rrr_vector

//! A bit vector which compresses the input with the method from Raman, Raman, and Rao
/*!
    Recently, I discovered that Rasmus Pagh was the first who presented
	the design of this bitvector representation. For detail see the Technical Report
	,,Low redundancy in dictionaries with O(1) worst case lookup time''
    ftp://ftp.cs.au.dk/BRICS/Reports/RS/98/28/BRICS-RS-98-28.pdf, Section 2.

    This compact representation was presented by
	Rajeev Raman, V. Raman and S. Srinivasa Rao at SODA 2002 in the article:
	Succinct Indexable Dictionaries with Applications to representations
	of k-ary trees and multi-sets.

	In this version the block size can be adjust by the template parameter block_size!
	\sa sdsl::rrr_vector for a specialized version for block_size=15
*/

template<uint16_t block_size=15, class wt_type=int_vector<> >
class rrr_vector
{
    public:
        typedef bit_vector::size_type size_type;
        typedef bit_vector::value_type value_type;

        typedef rrr_rank_support<1, block_size, wt_type> rank_1_type; // typedef for default types for rank and select
        typedef rrr_rank_support<0, block_size, wt_type> rank_0_type;
        typedef rrr_select_support<1, block_size, wt_type> select_1_type;
        typedef rrr_select_support<0, block_size, wt_type> select_0_type;

//        template<uint16_t b, uint16_t block_size, class wt_type> Hm, the second and third argument is
//        friend class rrr_rank_support;  // fixed
        friend class rrr_rank_support<0, block_size, wt_type>;
        friend class rrr_rank_support<1, block_size, wt_type>;
        friend class rrr_select_support<0, block_size, wt_type>;
        friend class rrr_select_support<1, block_size, wt_type>;

        typedef rrr_helper<block_size> rrr_helper_type;
        typedef typename rrr_helper_type::number_type number_type;

        enum { rrr_block_size = block_size };
    private:
        size_type      m_size; // length of the original bit_vector
        uint16_t       m_sample_rate;
        wt_type        m_bt; // data structure, which stores the block types (bt). The block type equals the number
        // of ones in a block. Another option for this data structure is wt_huff
        bit_vector     m_btnr; // data structure, which stores the block type numbers of the blocks
        int_vector<>   m_btnrp; // sample pointers into btnr
        int_vector<>   m_rank;  // sample rank values
        bit_vector     m_invert; // specifies if a superblock (i.e. sample_rate blocks) have to be considered as inverted
        // i.e. 1 and 0 are swapped

        void copy(const rrr_vector& rrr) {
            m_size = rrr.m_size;
            m_sample_rate = rrr.m_sample_rate;
            m_bt = rrr.m_bt;
            m_btnr = rrr.m_btnr;
            m_btnrp = rrr.m_btnrp;
            m_rank = rrr.m_rank;
            m_invert = rrr.m_invert;
        }

    public:
        const wt_type& bt;
        const bit_vector& btnr;

        //! Default constructor
        rrr_vector(uint16_t sample_rate=32):m_size(0), m_sample_rate(sample_rate), bt(m_bt), btnr(m_btnr) {};

        //! Copy constructor
        rrr_vector(const rrr_vector& rrr):bt(m_bt), btnr(m_btnr) {
            copy(rrr);
        }

        //! Constructor
        /*!
        *  \param block_size Number of bits in one block. \f$ block\_size \in \{1,...,23\} \f$
        */
        rrr_vector(const bit_vector& bv, uint16_t sample_rate=32): m_sample_rate(sample_rate), bt(m_bt), btnr(m_btnr) {
            m_size = bv.size();
            int_vector<> bt_array;
            bt_array.set_int_width(bit_magic::l1BP(block_size)+1);
            bt_array.resize((m_size+block_size)/((size_type)block_size));

            // (1) calculate the block types and store them in m_bt
            size_type pos = 0, i = 0, x;
            size_type btnr_pos = 0;
            size_type sum_rank = 0;
            while (pos + block_size <= m_size) { // handle all blocks full blocks
                bt_array[ i++ ] = x = rrr_helper_type::get_bt(bv, pos, block_size);
                sum_rank += x;
                btnr_pos += rrr_helper_type::space_for_bt(x);
                pos += block_size;
            }
            if (pos < m_size) { // handle last not full block
                bt_array[ i++ ] = x = rrr_helper_type::get_bt(bv, pos, m_size - pos);
                sum_rank += x;
                btnr_pos += rrr_helper_type::space_for_bt(x);
            }
            util::assign(m_btnr, bit_vector(std::max(btnr_pos, (size_type)64), 0));      // max necessary for case: block_size == 1
            util::assign(m_btnrp, int_vector<>((bt_array.size()+m_sample_rate-1)/m_sample_rate, 0,  bit_magic::l1BP(btnr_pos)+1));
            util::assign(m_rank, int_vector<>((bt_array.size()+m_sample_rate-1)/m_sample_rate + 1, 0, bit_magic::l1BP(sum_rank)+1));
            util::assign(m_invert, bit_vector((bt_array.size()+m_sample_rate-1)/m_sample_rate, 0));

            // (2) calculate block type numbers and pointers into btnr and rank samples
            pos = 0; i = 0;
            btnr_pos= 0, sum_rank = 0;
            bool invert = false;
            while (pos + block_size <= m_size) {  // handle all full blocks
                if ((i % m_sample_rate) == 0) {
                    m_btnrp[ i/m_sample_rate ] = btnr_pos;
                    m_rank[ i/m_sample_rate ] = sum_rank;
                    // calculate invert bit for that superblock
                    if (i+m_sample_rate <= bt_array.size()) {
                        size_type gt_half_block_size = 0; // counter for blocks greater than half of the blocksize
                        for (size_type j=i; j < i+m_sample_rate; ++j) {
                            if (bt_array[j] > block_size/2)
                                ++gt_half_block_size;
                        }
                        if (gt_half_block_size > (m_sample_rate/2)) {
                            m_invert[ i/m_sample_rate ] = 1;
                            for (size_type j=i; j < i+m_sample_rate; ++j) {
                                bt_array[j] = block_size - bt_array[j];
                            }
                            invert = true;
                        } else {
                            invert = false;
                        }
                    } else {
                        invert = false;
                    }
                }
                uint16_t space_for_bt = rrr_helper_type::space_for_bt(x=bt_array[i++]);
                sum_rank += (invert ? (block_size - x) : x);
                if (space_for_bt) {
                    number_type bin = rrr_helper_type::decode_btnr(bv, pos, block_size);
                    number_type nr = rrr_helper_type::bin_to_nr(bin);
                    rrr_helper_type::set_bt(m_btnr, btnr_pos, nr, space_for_bt);
                }
                btnr_pos += space_for_bt;
                pos += block_size;
            }
            if (pos < m_size) { // handle last not full block
                if ((i % m_sample_rate) == 0) {
                    m_btnrp[ i/m_sample_rate ] = btnr_pos;
                    m_rank[ i/m_sample_rate ] = sum_rank;
                    m_invert[ i/m_sample_rate ] = 0; // default: set last block to not inverted
                    invert = false;
                }
                uint16_t space_for_bt = rrr_helper_type::space_for_bt(x=bt_array[i++]);
                sum_rank += invert ? (block_size - x) : x;
                if (space_for_bt) {
                    number_type bin = rrr_helper_type::decode_btnr(bv, pos, m_size-pos);
                    number_type nr = rrr_helper_type::bin_to_nr(bin);
                    rrr_helper_type::set_bt(m_btnr, btnr_pos, nr, space_for_bt);
                }
                btnr_pos += space_for_bt;
            }
            // for technical reasons add an additional element to m_rank
            m_rank[ m_rank.size()-1 ] = sum_rank; // sum_rank contains the total number of set bits in bv
            util::assign(m_bt, bt_array);
        }

        //! Swap method
        void swap(rrr_vector& rrr) {
            if (this != &rrr) {
                std::swap(m_size, rrr.m_size);
                std::swap(m_sample_rate, rrr.m_sample_rate);
                m_bt.swap(rrr.m_bt);
                m_btnr.swap(rrr.m_btnr);
                m_btnrp.swap(rrr.m_btnrp);
                m_rank.swap(rrr.m_rank);
                m_invert.swap(rrr.m_invert);
            }
        }

        //! Accessing the i-th element of the original bit_vector
        /*! \param i An index i with \f$ 0 \leq i < size()  \f$.
           \return The i-th bit of the original bit_vector
        */
        value_type operator[](size_type i)const {
            size_type bt_idx = i/block_size;
            uint16_t bt = m_bt[ bt_idx ];
            size_type sample_pos = bt_idx/m_sample_rate;
            if (m_invert[sample_pos])
                bt = block_size - bt;
            if (bt == 0 or bt == block_size) { // very effective optimization
                return bt>0;
            }
            uint16_t off = i % block_size; //i - bt_idx*block_size;
            size_type btnrp = m_btnrp[ sample_pos ];
            for (size_type j = sample_pos*m_sample_rate; j < bt_idx; ++j) {
                btnrp += rrr_helper_type::space_for_bt(m_bt[j]);
            }
            uint16_t btnrlen 	= rrr_helper_type::space_for_bt(bt);
            number_type	btnr	= rrr_helper_type::decode_btnr(m_btnr, btnrp, btnrlen);
            return rrr_helper_type::decode_bit(bt, btnr, off);
        }

        //! Assignment operator
        rrr_vector& operator=(const rrr_vector& rrr) {
            if (this != &rrr) {
                copy(rrr);
            }
            return *this;
        }

        //! Returns the size of the original bit vector.
        size_type size()const {
            return m_size;
        }

        //! Answers select queries
        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += util::write_member(m_size, out, child, "size");
            written_bytes += util::write_member(m_sample_rate, out, child, "sample_rate");
            written_bytes += m_bt.serialize(out, child, "bt");
            written_bytes += m_btnr.serialize(out, child, "btnr");
            written_bytes += m_btnrp.serialize(out, child, "btnrp");
            written_bytes += m_rank.serialize(out, child, "rank_samples");
            written_bytes += m_invert.serialize(out, child, "invert");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            util::read_member(m_size, in);
            util::read_member(m_sample_rate, in);
            m_bt.load(in);
            m_btnr.load(in);
            m_btnrp.load(in);
            m_rank.load(in);
            m_invert.load(in);
        }

        // Print information about the object to stdout.
        void print_info()const {
            size_type orig_bv_size = m_size; // size of the original bit vector in bits
            size_type rrr_size 	= util::get_size_in_bytes(*this)*8;
            size_type bt_size 	= util::get_size_in_bytes(m_bt)*8;
            size_type btnr_size 	= util::get_size_in_bytes(m_btnr)*8;
            size_type btnrp_and_rank_size 	= util::get_size_in_bytes(m_btnrp)*8 + util::get_size_in_bytes(m_rank)*8 + util::get_size_in_bytes(m_invert)*8;
            std::cout << "#block_size\tsample_rate\torig_bv_size\trrr_size\tbt_size\tbtnr_size\tbtnrp_and_rank_size" << std::endl;
            std::cout << (int)block_size << "\t" << m_sample_rate << "\t";
            std::cout << orig_bv_size << "\t" << rrr_size << "\t" << bt_size << "\t" << btnr_size << "\t"
                      << btnrp_and_rank_size << std::endl;
        }
};

template<uint8_t bit_pattern>
struct rrr_rank_support_trait {
    typedef bit_vector::size_type size_type;

    static size_type adjust_rank(size_type r, size_type n) {
        return r;
    }
};

template<>
struct rrr_rank_support_trait<0> {
    typedef bit_vector::size_type size_type;

    static size_type adjust_rank(size_type r, size_type n) {
        return n - r;
    }
};

//! rank_support for the rrr_vector class
/*! The first template parameter is the bit pattern of size one.
 *  The second one the block size and the third the array type
 *  that is used to store the block types.
 *  TODO: Test if the binary search can be speed up by
 *        saving the (n/2)-th rank value in T[0], the (n/4)-th in T[1],
 *        the (3n/4)-th in T[2],... for small number of rank values
 *    is this called hinted binary search???
 *    or is this called
 */
template< uint8_t b, uint16_t block_size, class wt_type>
class rrr_rank_support
{
    public:
        typedef rrr_vector<block_size, wt_type> bit_vector_type;
        typedef typename bit_vector_type::size_type size_type;
        typedef typename bit_vector_type::rrr_helper_type rrr_helper_type;
        typedef typename rrr_helper_type::number_type number_type;

    private:
        const bit_vector_type* m_v; //!< Pointer to the rank supported rrr_vector
        uint16_t m_sample_rate;  //!<    "     "   "      "

    public:
        //! Standard constructor
        /*! \param v Pointer to the rrr_vector, which should be supported
         */
        explicit rrr_rank_support(const bit_vector_type* v=NULL) {
            init(v);
        }

        //! Initialize the data structure with a rrr_vector, which should be supported
        void init(const bit_vector_type* v=NULL) {
            set_vector(v);
        }

        //! Answers rank queries
        /*! \param i Argument for the length of the prefix v[0..i-1], with \f$0\leq i \leq size()\f$.
           \returns Number of 1-bits in the prefix [0..i-1] of the original bit_vector.
           \par Time complexity
        		\f$ \Order{ sample\_rate of the rrr\_vector} \f$
        */
        const size_type rank(size_type i)const {
            size_type bt_idx = i/block_size;
            size_type sample_pos = bt_idx/m_sample_rate;
            size_type btnrp = m_v->m_btnrp[ sample_pos ];
            size_type rank  = m_v->m_rank[ sample_pos ];
            size_type diff_rank  = m_v->m_rank[ sample_pos+1 ] - rank;
            if (diff_rank == 0) {
                return  rrr_rank_support_trait<b>::adjust_rank(rank, i);
            } else if (diff_rank == block_size*m_sample_rate) {
                return  rrr_rank_support_trait<b>::adjust_rank(
                            rank + i - sample_pos*m_sample_rate*block_size, i);
            }
            const bool inv = m_v->m_invert[ sample_pos ];
            for (size_type j = sample_pos*m_sample_rate; j < bt_idx; ++j) {
                uint16_t r = m_v->m_bt[j];
                rank  += (inv ? block_size - r: r);
                btnrp += rrr_helper_type::space_for_bt(r);
            }
            uint16_t off = i % block_size;
            if (!off) {   // needed for special case: if i=size() is a multiple of block_size
                // the access to m_bt would cause a invalid memory access
                return rrr_rank_support_trait<b>::adjust_rank(rank, i);
            }
            uint16_t bt = inv ? block_size - m_v->m_bt[ bt_idx ] : m_v->m_bt[ bt_idx ];

            uint16_t btnrlen 	= rrr_helper_type::space_for_bt(bt);
            number_type	btnr	= rrr_helper_type::decode_btnr(m_v->m_btnr, btnrp, btnrlen);
            uint16_t popcnt 	= rrr_helper_type::decode_popcount(bt, btnr, off);
            return rrr_rank_support_trait<b>::adjust_rank(rank + popcnt, i);
        }

        //! Short hand for rank(i)
        const size_type operator()(size_type i)const {
            return rank(i);
        }

        //! Returns the size of the original vector
        const size_type size()const {
            return m_v->size();
        }

        //! Set the supported vector.
        void set_vector(const bit_vector_type* v=NULL) {
            m_v = v;
            if (v != NULL) {
                m_sample_rate = m_v->m_sample_rate;
            } else {
                m_sample_rate = 0;
            }
        }

        rrr_rank_support& operator=(const rrr_rank_support& rs) {
            if (this != &rs) {
                set_vector(rs.m_v);
                m_sample_rate = rs.m_sample_rate;
            }
            return *this;
        }

        void swap(rrr_rank_support& rs) {
            if (this != &rs) {
                std::swap(m_sample_rate, rs.m_sample_rate);
            }
        }

        bool operator==(const rrr_rank_support& rs)const {
            if (this == &rs)
                return true;
            return m_sample_rate == rs.m_sample_rate;
        }

        bool operator!=(const rrr_rank_support& rs)const {
            return !(*this == rs);
        }

        //! Load the data structure from a stream and set the supported vector.
        void load(std::istream& in, const bit_vector_type* v=NULL) {
            util::read_member(m_sample_rate, in);
            set_vector(v);
        }

        //! Serializes the data structure into a stream.
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += util::write_member(m_sample_rate, out, child, "sample_rate");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }
};


//! Select support for the rrr_vector class.
template< uint8_t b, uint16_t block_size, class wt_type>
class rrr_select_support
{
    public:
        typedef rrr_vector<block_size, wt_type> bit_vector_type;
        typedef typename bit_vector_type::size_type size_type;
        typedef typename bit_vector_type::rrr_helper_type rrr_helper_type;
        typedef typename rrr_helper_type::number_type number_type;

    private:
        const bit_vector_type* m_v; //!< Pointer to the rank supported rrr_vector
        uint16_t m_sample_rate;  //!<    "     "   "      "

        size_type select1(size_type i)const {
            if (m_v->m_rank[m_v->m_rank.size()-1] < i)
                return size();
            //  (1) binary search for the answer in the rank_samples
            size_type begin=0, end=m_v->m_rank.size()-1; // min included, max excluded
            size_type idx, rank;
            // invariant:  m_rank[end] >= i
            //             m_rank[begin] < i
            while (end-begin > 1) {
                idx  = (begin+end) >> 1; // idx in [0..m_rank.size()-1]
                rank = m_v->m_rank[idx];
                if (rank >= i)
                    end = idx;
                else { // rank < i
                    begin = idx;
                }
            }
            //   (2) linear search between the samples
            rank = m_v->m_rank[begin]; // now i>rank
            idx = begin * m_sample_rate; // initialize idx for select result
            size_type diff_rank  = m_v->m_rank[end] - rank;
            if (diff_rank == block_size*m_sample_rate) {// optimisiation for select<1>
                return idx*block_size + i-rank -1;
            }
            const bool inv = m_v->m_invert[ begin ];
            size_type btnrp = m_v->m_btnrp[ begin ];
            uint16_t bt = 0, btnrlen = 0; // temp variables for block_type and space for block type
            while (i > rank) {
                bt = m_v->m_bt[idx++]; bt = inv ? block_size-bt : bt;
                rank += bt;
                btnrp += (btnrlen=rrr_helper_type::space_for_bt(bt));
            }
            rank -= bt;
            number_type	btnr = rrr_helper_type::decode_btnr(m_v->m_btnr, btnrp-btnrlen, btnrlen);
            return (idx-1) * block_size + rrr_helper_type::decode_select(bt, btnr, i-rank);
        }

        size_type  select0(size_type i)const {
            if ((size() - m_v->m_rank[m_v->m_rank.size()-1]) < i) {
//				std::cout<<"i="<<i<<" size()-m_v->m_rank[m_v->m_rank.size()-1] = " << m_v->m_rank[m_v->m_rank.size()-1]<< std::endl;
                return size();
            }
            //  (1) binary search for the answer in the rank_samples
            size_type begin=0, end=m_v->m_rank.size()-1; // min included, max excluded
            size_type idx, rank;
            // invariant:  m_rank[end] >= i
            //             m_rank[begin] < i
            while (end-begin > 1) {
                idx  = (begin+end) >> 1; // idx in [0..m_rank.size()-1]
                rank = idx*block_size*m_sample_rate - m_v->m_rank[idx];
                if (rank >= i)
                    end = idx;
                else { // rank < i
                    begin = idx;
                }
            }
            //   (2) linear search between the samples
            rank = begin*block_size*m_sample_rate - m_v->m_rank[begin]; // now i>rank
            idx = begin * m_sample_rate; // initialize idx for select result
            if (m_v->m_rank[end] == m_v->m_rank[begin]) {      // only for select<0>
                return idx*block_size +  i-rank -1;
            }
            const bool inv = m_v->m_invert[ begin ];
            size_type btnrp = m_v->m_btnrp[ begin ];
            uint16_t bt = 0, btnrlen = 0; // temp variables for block_type and space for block type
            while (i > rank) {
                bt = m_v->m_bt[idx++]; bt = inv ? block_size-bt : bt;
                rank += (block_size-bt);
                btnrp += (btnrlen=rrr_helper_type::space_for_bt(bt));
            }
            rank -= (block_size-bt);
            number_type	btnr = rrr_helper_type::decode_btnr(m_v->m_btnr, btnrp-btnrlen, btnrlen);
            return (idx-1) * block_size + rrr_helper_type::template decode_select_bitpattern<0, 1>(bt, btnr, i-rank);
        }



    public:
        explicit rrr_select_support(const bit_vector_type* v=NULL) {
            init(v);
        }

        void init(const bit_vector_type* v=NULL) {
            set_vector(v);
        }

        //! Answers select queries
        size_type select(size_type i)const {
            return  b ? select1(i) : select0(i);
        }


        const size_type operator()(size_type i)const {
            return select(i);
        }

        const size_type size()const {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=NULL) {
            m_v = v;
            if (v != NULL) {
                m_sample_rate = m_v->m_sample_rate;
            } else {
                m_sample_rate = 0;
            }
        }

        rrr_select_support& operator=(const rrr_select_support& rs) {
            if (this != &rs) {
                set_vector(rs.m_v);
                m_sample_rate = rs.m_sample_rate;
            }
            return *this;
        }

        void swap(rrr_select_support& rs) {
            if (this != &rs) {
                std::swap(m_sample_rate, rs.m_sample_rate);
            }
        }

        bool operator==(const rrr_select_support& rs)const {
            if (this == &rs)
                return true;
            return m_sample_rate == rs.m_sample_rate;
        }

        bool operator!=(const rrr_select_support& rs)const {
            return !(*this == rs);
        }


        void load(std::istream& in, const bit_vector_type* v=NULL) {
            util::read_member(m_sample_rate, in);
            set_vector(v);
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += util::write_member(m_sample_rate, out, child, "sample_rate");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }
};

}// end namespace sdsl
#include "rrr_vector_15.hpp" // include specialization 

#endif
