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
/*! \file rank_support_jmc.hpp
    \brief rank_support_jmc.hpp contains classes that support a sdsl::bit_vector with constant time rank information.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_RANK_SUPPORT_JMC
#define INCLUDED_SDSL_RANK_SUPPORT_JMC

#include "rank_support.hpp"
#include "int_vector.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A class supporting rank queries in constant time. The implementation is a lightweight version of the data structure proposed by Jacobson (1989), Munro (1996), and Clark (1996).
/*!
 *  - TODO: Space complexity
 * @ingroup rank_support_group
 */
class rank_support_jmc : public rank_support
{
    public:
        typedef bit_vector bit_vector_type;
    private:
        size_type m_logn;
        int_vector<0> m_superblockrank;
        int_vector<0> m_blockrank;
    public:
        rank_support_jmc(const int_vector<1>* v=NULL);
        rank_support_jmc(const rank_support_jmc& rs);
        ~rank_support_jmc();
        void init(const int_vector<1>* v=NULL);
        inline const size_type rank(size_type idx) const;
        inline const size_type operator()(size_type idx)const;
        size_type serialize(std::ostream& out)const;
        void load(std::istream& in, const int_vector<1>* v=NULL);
        void set_vector(const int_vector<1>* v);
        //! Assign Operator
        /*! Required for the Assignable Concept of the STL.
         */
        rank_support_jmc& operator=(const rank_support_jmc& rs);
        //! swap Operator
        /*! Swap two rank_support_jmc in constant time.
        	All members (the pointer to the supported SDSBitVector excluded) are swapped.

            Required for the Container Concept of the STL.
         */
        void swap(rank_support_jmc& rs);
        //! Equality Operator
        /*! Two rank_support_jmcs are equal if all member variables are equal.
         *
         * Required for the Equality Comparable Concept of the STL.
         * \sa operator!=
         */
        bool operator==(const rank_support_jmc& rs)const;
        //! Unequality Operator
        /*! Two rank_support_jmcs are not equal if any member variable are not equal.
         *
         * Required for the Equality Comparable Concept of the STL.
         * \sa operator==
         */
        bool operator!=(const rank_support_jmc& rs)const;
};

inline rank_support_jmc::rank_support_jmc(const int_vector<1>* v)
{
    m_logn = 0;
    init(v);
}

inline rank_support_jmc::rank_support_jmc(const rank_support_jmc& rs) : rank_support()
{
    set_vector(rs.m_v);
    m_superblockrank = rs.m_superblockrank;
    m_blockrank		 = rs.m_blockrank;
}

inline rank_support_jmc& rank_support_jmc::operator=(const rank_support_jmc& rs)
{
    if (this != &rs) {
        set_vector(rs.m_v);
        m_superblockrank = rs.m_superblockrank;
        m_blockrank		= rs.m_blockrank;
    }
    return *this;
}

inline void rank_support_jmc::swap(rank_support_jmc& rs)
{
    if (this != &rs) { // if rs and _this_ are not the same object
        std::swap(m_logn, rs.m_logn);
        m_superblockrank.swap(rs.m_superblockrank);
        m_blockrank.swap(rs.m_blockrank);
    }
}

inline void rank_support_jmc::init(const int_vector<1>* v)
{
    set_vector(v);
    if (m_v == NULL) return;
    if (m_v->empty()) {
        m_blockrank.set_int_width(1); m_superblockrank.set_int_width(1);
        m_blockrank.resize(1);		m_superblockrank.resize(1);
        m_blockrank[0] = 0; m_superblockrank[0] = 0;
        return;
    }
    m_blockrank.set_int_width(12);
    m_blockrank.resize((m_v->capacity()>>6) + (0==(m_v->size()&0x3F)));     // n/64 + 2*loglog 64
    m_superblockrank.set_int_width(m_logn);
    m_superblockrank.resize((m_blockrank.size()+63)>>6);

    m_blockrank[0]=0;
    m_superblockrank[0]=0;
    size_type cnt = 0, blockcnt = 0, wcnt = 0;
    const uint64_t* data = m_v->data();
    size_type i;
    for (i = 1; i < (m_v->capacity()>>6) ; ++i) {
        wcnt = bit_magic::b1Cnt(*data);
        ++data;
        blockcnt += wcnt;
        cnt 	 += wcnt;
        if ((i & 0x3F) == 0) {
            m_superblockrank[i>>6] = cnt;
            blockcnt = 0;
        }
        m_blockrank[i] = blockcnt;
    }
    if (0 == (m_v->size()&0x3F)) {
        wcnt = bit_magic::b1Cnt(*data);
        blockcnt += wcnt;
        cnt		 += wcnt;
        if ((i & 0x3F) == 0) {
            m_superblockrank[i>>6] = cnt;
            blockcnt = 0;
        }
        m_blockrank[i] = blockcnt;
    }
}

inline const rank_support_jmc::size_type rank_support_jmc::rank(size_type idx)const
{
    if ((idx & 0x3F) ==0)
        return m_blockrank[idx>>6]
               + m_superblockrank[idx>>12];
    return
        bit_magic::b1Cnt((*(m_v->data()+(idx>>6))&bit_magic::Li1Mask[idx & 0x3F]))
        + m_blockrank[idx>>6]
        + m_superblockrank[idx>>12];
}

inline const rank_support_jmc::size_type rank_support_jmc::operator()(size_type idx)const
{
    return rank(idx);
}

inline rank_support_jmc::size_type rank_support_jmc::serialize(std::ostream& out)const
{
    return    m_blockrank.serialize(out)
              + m_superblockrank.serialize(out);
}

inline void rank_support_jmc::load(std::istream& in, const int_vector<1>* v)
{
    set_vector(v);
    assert(m_v != NULL); // supported bit vector should be known
    m_blockrank.load(in);
    m_superblockrank.load(in);
}

inline rank_support_jmc::~rank_support_jmc()
{
}

inline void rank_support_jmc::set_vector(const int_vector<1>* v)
{
    if (v != NULL) {
        m_v = v;
        m_logn = bit_magic::l1BP(m_v->capacity())+1;
    }
}

inline bool rank_support_jmc::operator==(const rank_support_jmc& rs)const
{
    if (this == &rs)
        return true;
    return m_logn == rs.m_logn
           and m_superblockrank == rs.m_superblockrank
           and m_blockrank == rs.m_blockrank
           and *(m_v) == *(rs.m_v);
}

inline bool rank_support_jmc::operator!=(const rank_support_jmc& rs)const
{
    return !(*this == rs);
}

}// end namespace sds

#endif // end file 
