/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

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
/*! \file rank_support.hpp
    \brief rank_support.hpp contains classes that support a sdsl::bit_vector with constant time rank information.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_RANK_SUPPORT
#define INCLUDED_SDSL_RANK_SUPPORT

/** \defgroup rank_support_group Rank Support (RS)
 * This group contains data structures which support an sdsl::bit_vector with the rank method.
 */

#include "int_vector.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! The base class of classes supporting rank_queries for a sdsl::bit_vector in constant time.
/*!
*/
class rank_support
{
    protected:
        const bit_vector* m_v; //!< Pointer to the rank supported bit_vector
    public:
        typedef bit_vector::size_type size_type;

        //! Constructor
        /*! \param v The supported bit_vector.
         */
        rank_support(const bit_vector* v = nullptr);
        //! Copy constructor
        rank_support(const rank_support&) = default;
        rank_support(rank_support&&) = default;
        rank_support& operator=(const rank_support&) = default;
        rank_support& operator=(rank_support&&) = default;
        //! Destructor
        virtual ~rank_support() {}

        //! Answers rank queries for the supported bit_vector.
        /*!	\param i Argument for the length of the prefix v[0..i-1].
        	\returns Number of 1-bits in the prefix [0..i-1] of the supported bit_vector.
        	\note Method init has to be called before the first call of rank.
        	\sa init
         */
        virtual size_type rank(size_type i) const = 0;
        //! Alias for rank(i)
        virtual size_type operator()(size_type idx) const = 0;
        //! Serializes rank_support.
        /*! \param out Out-Stream to serialize the data to.
        */
        virtual size_type serialize(std::ostream& out, structure_tree_node* v, std::string name)const = 0;
        //! Loads the rank_support.
        /*! \param in In-Stream to load the rank_support data from.
            \param v The supported bit_vector.
         */
        virtual void load(std::istream& in, const bit_vector* v=nullptr) = 0;
        //! Sets the supported bit_vector to the given pointer.
        /*! \param v The new bit_vector to support.
         *  \note Method init has to be called before the next call of rank.
         *  \sa init, rank
         */
        virtual void set_vector(const bit_vector* v=nullptr) = 0;
};

inline rank_support::rank_support(const bit_vector* v)
{
    m_v = v;
}

//----------------------------------------------------------------------

template<uint8_t bit_pattern, uint8_t pattern_len>
struct rank_support_trait {
    typedef rank_support::size_type	size_type;

    static size_type args_in_the_word(uint64_t, uint64_t&) {
        return 0;
    }

    static uint32_t word_rank(const uint64_t*, size_type) {
        return 0;
    }

    static uint32_t full_word_rank(const uint64_t*, size_type) {
        return 0;
    }

    static uint64_t init_carry() {
        return 0;
    }
};

template<>
struct rank_support_trait<0,1> {
    typedef rank_support::size_type	size_type;

    static size_type args_in_the_word(uint64_t w, uint64_t&) {
        return bits::cnt(~w);
    }

    static uint32_t word_rank(const uint64_t* data, size_type idx) {
        return	bits::cnt((~*(data+(idx>>6))) & bits::lo_set[idx&0x3F]);
    }

    static uint32_t full_word_rank(const uint64_t* data, size_type idx) {
        return	bits::cnt((~*(data+(idx>>6))));
    }

    static uint64_t init_carry() {
        return 0;
    }
};

template<>
struct rank_support_trait<1,1> {
    typedef rank_support::size_type	size_type;

    static size_type args_in_the_word(uint64_t w, uint64_t&) {
        return bits::cnt(w);
    }

    static uint32_t word_rank(const uint64_t* data, size_type idx) {
        return	bits::cnt(*(data+(idx>>6)) & bits::lo_set[idx&0x3F]);
    }

    static uint32_t full_word_rank(const uint64_t* data, size_type idx) {
        return	bits::cnt(*(data+(idx>>6)));
    }

    static uint64_t init_carry() {
        return 0;
    }
};

template<>
struct rank_support_trait<10,2> {
    typedef rank_support::size_type	size_type;

    static size_type args_in_the_word(uint64_t w, uint64_t& carry) {
        return bits::cnt10(w, carry);
    }

    static uint32_t word_rank(const uint64_t* data, size_type idx) {
        data = data+(idx>>6);
        uint64_t carry = (idx>63) ? *(data-1)>>63 : 0;
        return	bits::cnt(bits::map10(*data, carry) & bits::lo_set[idx&0x3F]);
    }

    static uint32_t full_word_rank(const uint64_t* data, size_type idx) {
        data = data+(idx>>6);
        uint64_t carry = (idx>63) ? *(data-1)>>63 : 0;
        return	bits::cnt(bits::map10(*data, carry));
    }

    static uint64_t init_carry() {
        return 0;
    }
};

template<>
struct rank_support_trait<01,2> {
    typedef rank_support::size_type	size_type;

    static size_type args_in_the_word(uint64_t w, uint64_t& carry) {
        return bits::cnt01(w, carry);
    }

    static uint32_t word_rank(const uint64_t* data, size_type idx) {
        data = data+(idx>>6);
        uint64_t carry = (idx>63) ? *(data-1)>>63 : 0;
        return	bits::cnt(bits::map01(*data, carry) & bits::lo_set[idx&0x3F]);
    }

    static uint32_t full_word_rank(const uint64_t* data, size_type idx) {
        data = data+(idx>>6);
        uint64_t carry = (idx>63) ? *(data-1)>>63 : 0;
        return	bits::cnt(bits::map01(*data, carry));
    }

    static uint64_t init_carry() {
        return 1;
    }
};

}// end namespace sdsl

#include "rank_support_v.hpp"
#include "rank_support_v5.hpp"
#include "rank_support_scan.hpp"

#endif // end file 
