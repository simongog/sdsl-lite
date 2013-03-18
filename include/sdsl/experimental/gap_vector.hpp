/* sdsl - succinct data structures library
    Copyright (C) 2012 Simon Gog

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
/*! \file gap_vector.hpp
   \brief gap_vector.hpp contains the sdsl::gap_vector class, and
          classes which support rank and select for gap_vector.
   \author Simon Gog
*/
#ifndef INCLUDED_SDSL_GAP_VECTOR
#define INCLUDED_SDSL_GAP_VECTOR

#include "int_vector.hpp"
#include "util.hpp"

//! Namespace for the succinct data structure library
namespace sdsl
{

template<bool t_b=true>// forward declaration needed for friend declaration
class gap_rank_support;  // in gap_vector

template<bool t_b=true>// forward declaration needed for friend declaration
class gap_select_support;  // in gap_vector

//! A bit vector which compresses very sparse populated bit vectors by representing either 1's or 0's by gap encoding
template<bool t_b=true>
class gap_vector
{
    public:
        typedef bit_vector::size_type size_type;
        typedef size_type value_type;

        friend class gap_rank_support<true>;
        friend class gap_select_support<true>;

        typedef gap_rank_support<true> rank_1_type;
        typedef gap_select_support<true> select_1_type;
    private:
        size_type m_size;
        int_vector<> m_position;

    public:
        gap_vector():m_size(0) {}

        gap_vector(const bit_vector& bv) {
            m_size = bv.size();
            if (m_size == 0)
                return;
            size_type ones = util::get_one_bits(bv);
            m_position = int_vector<>(ones, 0, bit_magic::l1BP(m_size)+1);
            const uint64_t* bvp = bv.data();
            for (size_type i=0, one_cnt=0; i < (bv.size()+63)/64; ++i, ++bvp) {
                if (*bvp) { // if there is a one in the word
                    for (size_type j=0; j<64 and 64*i+j < bv.size(); ++j) // check each bit of the word
                        if (bv[64*i+j]) {
                            m_position[one_cnt++] = 64*i+j;
                        }
                }
            }
        }

        //! Swap method
        void swap(gap_vector& v) {
            if (this != &v) {
                std::swap(m_size, v.m_size);
                m_position.swap(v.m_position);
            }
        }

        //! Accessing the i-th element of the original bit_vector
        /*! \param i An index i with \f$ 0 \leq i < size()  \f$.
           \return The i-th bit of the original bit_vector
           \par Time complexity
           		\f$ \Order{\log m} \f$, where m equals the number of zeros
        */
        value_type operator[](size_type i)const {
            // binary search the entries in m_position
            size_type lb=0, rb=m_position.size(), mid, pos;  // start interval [lb,rb)
            while (rb > lb) {
                mid = (lb+rb)/2; // then mid>=lb mid<rb
                pos = m_position[mid];
                if (i > pos) {
                    lb = mid+1;
                } else if (i < pos) {
                    rb = mid;
                } else { // i == pos
                    return 1;
                }
            }
            return 0;
        }

        //! Returns the size of the original bit vector.
        size_type size()const {
            return m_size;
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            size_type written_bytes = 0;
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            written_bytes += util::write_member(m_size, out, child, "size");
            written_bytes += m_position.serialize(out, child, "positions");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            util::read_member(m_size, in);
            m_position.load(in);
        }
};

template<bool t_b>
class gap_rank_support
{
    public:
        typedef bit_vector::size_type size_type;
        typedef gap_vector<t_b> bit_vector_type;
    private:
        const bit_vector_type* m_v;

    public:

        gap_rank_support(const bit_vector_type* v=NULL) {
            set_vector(v);
        }

        size_type rank(size_type i)const {
            // binary search the entries in m_position
            size_type lb=0, rb=m_v->m_position.size(), mid, pos=0;  // start interval [lb,rb)
            while (rb > lb) {
                mid = (lb+rb)/2; // then mid>=lb mid<rb
                pos = m_v->m_position[mid];
                if (i <= pos) {
                    rb = mid;
                } else {
                    lb = mid+1;
                }
            } // m_position[rb] >= i
            return rb;
        }

        const size_type operator()(size_type i)const {
            return rank(i);
        }

        const size_type size()const {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=NULL) {
            m_v = v;
        }

        gap_rank_support& operator=(const gap_rank_support& rs) {
            if (this != &rs) {
                set_vector(rs.m_v);
            }
            return *this;
        }

        void swap(gap_rank_support& rs) { }

        void load(std::istream& in, const bit_vector_type* v=NULL) {
            set_vector(v);
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }
};




template<bool t_b>
class gap_select_support
{
    public:
        typedef bit_vector::size_type size_type;
        typedef gap_vector<t_b> bit_vector_type;
    private:
        const bit_vector_type* m_v;

    public:
        gap_select_support(const bit_vector_type* v=NULL) {
            set_vector(v);
        }

        //! Returns the position of the i-th occurrence in the bit vector.
        size_type select(size_type i)const {
            return m_v->m_position[i-1];
        }

        const size_type operator()(size_type i)const {
            return select(i);
        }

        const size_type size()const {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=NULL) {
            m_v = v;
        }

        gap_select_support& operator=(const gap_select_support& rs) {
            if (this != &rs) {
                set_vector(rs.m_v);
            }
            return *this;
        }

        void swap(gap_select_support& rs) { }

        void load(std::istream& in, const bit_vector_type* v=NULL) {
            set_vector(v);
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }
};

}

#endif
