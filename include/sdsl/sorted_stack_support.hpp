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
/*! \file sorted_stack_support.hpp
    \brief sorted_stack_support.hpp contains a data structure for a stack which contains indices of a random access container and the elements of the indices are sorted in the stack order. This data structure was proposed by Johannes Fischer in the paper Wee LCP.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_SORTED_STACK_SUPPORT
#define INCLUDED_SDSL_SORTED_STACK_SUPPORT

#include "int_vector.hpp"
#include "bitmagic.hpp"
#include <vector>

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A stack class which contains indices of elements from an random access container and the elements are in sorted order on the stack.
/*! \par Space complexity
 *  \f$n\f$ bits
 */
class sorted_stack_support
{
    public:
        typedef int_vector<64>::size_type size_type;
    private:
        size_type m_n;   // Size of the supported vector.
        size_type m_cnt; // Counter for the indices on the stack.
        size_type m_top; // Topmost index of the stack.
        int_vector<64> m_stack; // Memory for the stack.

        inline size_type block_nr(size_type x) {
            return x/63;
        }; // TODO: maybe we can speed this up with bit hacks
        inline size_type block_pos(size_type x) {
            return x%63;
        }; // TODO: maybe we can speed this up with bit hacks
    public:
        //! Constructor
        /*! \param n Miximum that can be pushed onto the stack
         */
        sorted_stack_support(size_type n);

        sorted_stack_support(const sorted_stack_support& sis);
        ~sorted_stack_support() {};

        /*! Returns if the stack is empty.
         */
        bool empty() const {
            return 0==m_cnt;
        };

        /*! Returns the topmost index on the stack.
         * \pre empty()==false
         */
        const size_type top() const;

        /*! Pop the topmost index of the stack.
         */
        void pop();

        /*! Push the index x of vector vec onto the stack.
         * \par x Index of the value in vec which should be pushed onto the stack.
         * \pre top() < x
         */
        void push(size_type x);

        /*! Returns the number of element is the stack.
         */
        size_type size()const {
            return m_cnt;
        };

        size_type serialize(std::ostream& out)const;
        void load(std::istream& in);

        //! Assign Operator
        /*! Required for the Assignable Concept of the STL.
         */
        sorted_stack_support& operator=(const sorted_stack_support& sis);

        //! Equality Operator
        /*! Two sorted_stack_supports are equal if all member variables are equal.
         *
         * Required for the Equality Comparable Concept of the STL.
         * \sa operator!=
         */
        bool operator==(const sorted_stack_support& sis)const;
        //! Unequality Operator
        /*! Two sorted_stack_supports are not equal if any member variable are not equal.
         *
         * Required for the Equality Comparable Concept of the STL.
         * \sa operator==
         */
        bool operator!=(const sorted_stack_support& sis)const;

};

inline sorted_stack_support::sorted_stack_support(size_type n):m_n(n), m_cnt(0), m_top(0)
{
    m_stack = int_vector<64>(block_nr(m_n+1)+1, 0);
    m_stack[0] = 1;
}

inline sorted_stack_support::sorted_stack_support(const sorted_stack_support& sis):m_n(sis.m_n), m_cnt(sis.m_cnt), m_top(sis.m_top)
{
    m_stack = sis.m_stack;
}

inline sorted_stack_support& sorted_stack_support::operator=(const sorted_stack_support& sis)
{
    if (this != &sis) {
        m_n	 		= sis.m_n;
        m_cnt 		= sis.m_cnt;
        m_top		= sis.m_top;
        m_stack 	= sis.m_stack;
    }
    return *this;
}

inline bool sorted_stack_support::operator==(const sorted_stack_support& sis)const
{
    return m_n == sis.m_n and m_cnt == sis.m_cnt and m_top == sis.m_top and m_stack == sis.m_stack;
}

inline bool sorted_stack_support::operator!=(const sorted_stack_support& sis)const
{
    return !(*this==sis);
}

inline const sorted_stack_support::size_type sorted_stack_support::top()const
{
    return m_top-1;
}

inline void sorted_stack_support::push(size_type x)
{
    x += 1;
    ++m_cnt;	//< increment counter
    size_type bn = block_nr(x);
    m_stack[bn] ^= (1ULL << block_pos(x));
    if (bn > 0 and m_stack[bn-1] == 0) {
        m_stack[bn-1] = 0x8000000000000000ULL | m_top;
    }
    m_top = x;
}

inline void sorted_stack_support::pop()
{
    if (!empty()) {
        --m_cnt; //< decrement counter
        size_type bn = block_nr(m_top);
        uint64_t w = m_stack[ bn ];
        assert((w>>63) == 0);    // highest bit is not set, as the block contains no pointer
        w ^= (1ULL << block_pos(m_top));
        m_stack[ bn ] = w;
        if (w>0) {
            m_top = bn*63 + bit_magic::l1BP(w);
        } else { // w==0 and cnt>0
            assert(bn > 0);
            w = m_stack[ bn-1 ];
            if ((w>>63) == 0) { // highest bit is not set => the block contains no pointer
                assert(w>0);
                m_top = (bn-1)*63 + bit_magic::l1BP(w);
            } else { // block contains pointers
                m_stack[bn-1] = 0;
                m_top = w&0x7FFFFFFFFFFFFFFFULL;
            }
        }
    }
}

inline sorted_stack_support::size_type sorted_stack_support::serialize(std::ostream& out)const
{
    size_type written_bytes = 0;
    written_bytes += util::write_member(m_n, out);
    written_bytes += util::write_member(m_top, out);
    written_bytes += util::write_member(m_cnt, out);
    written_bytes += m_stack.serialize(out);
    return written_bytes;
}

inline void sorted_stack_support::load(std::istream& in)
{
    util::read_member(m_n, in);
    util::read_member(m_top, in);
    util::read_member(m_cnt, in);
    m_stack.load(in);
}

}// end namespace sdsl

#endif // end file 
