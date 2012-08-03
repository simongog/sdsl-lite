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
/*! \file rmq_succinct_sada.hpp
    \brief rmq_succinct_sada.hpp contains the class rmq_succinct_sada which supports range minimum or range maximum queries on a random access container in constant time and \f$4 n+o(n) bits\f$ space.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_RMQ_SUCCINCT_SADA
#define INCLUDED_SDSL_RMQ_SUCCINCT_SADA

#include "rmq_support.hpp"
#include "int_vector.hpp"
#include "algorithms_for_compressed_suffix_trees.hpp"
#include "bp_support_sada.hpp"
#include "rank_support_v.hpp"
#include "select_support_mcl.hpp"
#include "util.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{


template<class RandomAccessContainer = int_vector<>, bool Minimum = true, class Bp_support = bp_support_sada<>, class Rank_support10 = rank_support_v<10,2>, class Select_support10 = select_support_mcl<10,2> >
class rmq_succinct_sada;

template<class RandomAccessContainer = int_vector<>, class Bp_support = bp_support_sada<>, class Rank_support10 = rank_support_v<10,2>, class Select_support10 = select_support_mcl<10,2> >
struct range_maximum_support_sada {
    typedef rmq_succinct_sada<RandomAccessContainer, false, Bp_support, Rank_support10, Select_support10> type;
};

//! A class to support range minimum or range maximum queries on a random access container.
/*!
 * This class takes five template parameters:
 *  - RandomAccessContainer is the type of random access container, for which the range minimum/maximum support should be build,
 *  - minumum 				specifies whether the data structure should answer range minimum queries (mimumum=true) of range maximum queries (maximum=false), and
 *  - Bp_support			is the support structure for the balanced parentheses sequence of the balanced parentheses sequence of the extended Cartesian tree used internal in the class.
 *  - Rank_support10		is the rank structure which supports rank queries for the bit pattern "10".
 *	- Select_support10		is the select structure which supports select queries for the bit pattern "10".
 * \par Time complexity
 *		\f$ \Order{1} \f$ for the range minimum/maximum queries if the balanced parentheses support structure supports constant time operations.
 * \par Space complexity:
 *		\f$ 4n+o(n) \f$ bits for the data structure ( \f$ n=size() \f$ ).
 *
 * TODO: implement test
 */
template<class RandomAccessContainer, bool Minimum, class Bp_support, class Rank_support10, class Select_support10>
class rmq_succinct_sada
{
        bit_vector					m_ect_bp; 			//!< A bit vector which contains the balanced parentheses sequence of the extended Cartesian tree of the input container.
        Bp_support					m_ect_bp_support; 	//!< Support structure for the balanced parentheses sequence of the extended Cartesian tree.
        Rank_support10				m_ect_bp_rank10;		//!< A rank support (for bit pattern "10") which supports the balanced parentheses sequence of the extended Cartesian tree.
        Select_support10			m_ect_bp_select10; 	//!< A select support (for bit pattern "10") which supports the balanced parentheses sequence of the extended Cartesian tree.

    public:
        typedef typename RandomAccessContainer::size_type size_type;
        typedef typename RandomAccessContainer::size_type value_type;

        typedef Bp_support 			bp_support_type;
        typedef	Rank_support10 		rank_support10_type;
        typedef Select_support10	select_support10_type;

        const bit_vector&		 ect_bp;
        const Bp_support&		 ect_bp_support;
        const Rank_support10&	ect_bp_rank10;
        const Select_support10&	ect_bp_select10;

    private:

        typedef rmq_support_sct<RandomAccessContainer, Minimum> rmq_construct_helper_type;

        void _construct_bp_of_extended_cartesian_tree(size_type l, size_type r, size_type& bp_cnt, const rmq_construct_helper_type& rmq_helper) {
            if (r==(size_type)-1 or l > r)
                return;
            m_ect_bp[bp_cnt++] = 1; // write beginning of inner node
            size_type m = rmq_helper(l, r);
            _construct_bp_of_extended_cartesian_tree(l, m-1, bp_cnt, rmq_helper);
            m_ect_bp[bp_cnt++] = 1; // write leaf
            m_ect_bp[bp_cnt++] = 0;
            _construct_bp_of_extended_cartesian_tree(m+1, r, bp_cnt, rmq_helper);
            m_ect_bp[bp_cnt++] = 0; // write end of inner node
            assert(bp_cnt <= m_ect_bp.size());
        }

        void construct(const RandomAccessContainer* v) {
            if (v == NULL) {
                m_ect_bp = bit_vector(0); m_ect_bp_support = Bp_support();
                m_ect_bp_rank10 = Rank_support10(); m_ect_bp_select10 = Select_support10();
            } else {
                rmq_construct_helper_type rmq_helper(v);
                m_ect_bp.resize(4*v->size());
                size_type bp_cnt=0;
                _construct_bp_of_extended_cartesian_tree((size_type)0, v->size()-1, bp_cnt, rmq_helper);
                assert(bp_cnt == 4*v->size());
                m_ect_bp_support = Bp_support(&m_ect_bp);
                util::init_support(m_ect_bp_rank10, &m_ect_bp);
                util::init_support(m_ect_bp_select10, &m_ect_bp);
            }
        }

        void copy(const rmq_succinct_sada& rm) {
            m_ect_bp = rm.m_ect_bp;
            m_ect_bp_support = rm.m_ect_bp_support;
            m_ect_bp_support.set_vector(&m_ect_bp);
            m_ect_bp_rank10 = rm.m_ect_bp_rank10;
            m_ect_bp_rank10.set_vector(&m_ect_bp);
            m_ect_bp_select10 = rm.m_ect_bp_select10;
            m_ect_bp_select10.set_vector(&m_ect_bp);
        }

    public:

        //! Constructor
        rmq_succinct_sada(const RandomAccessContainer* v=NULL):ect_bp(m_ect_bp), ect_bp_support(m_ect_bp_support), ect_bp_rank10(m_ect_bp_rank10), ect_bp_select10(m_ect_bp_select10) {
            construct(v);
        }

        //! Copy constructor
        rmq_succinct_sada(const rmq_succinct_sada& rm) {
            if (this != &rm) { // if v is not the same object
                copy(rm);
            }
        }

        //! Destructor
        ~rmq_succinct_sada() { }

        rmq_succinct_sada& operator=(const rmq_succinct_sada& rm) {
            if (this != &rm) {
                copy(rm);
            }
            return *this;
        }

		//! Swap operator
        void swap(const rmq_sct& rm) {
            m_ect_bp.swap(rm.m_ect_bp);
			util::swap_support(m_ect_bp_support, rm.m_ect_bp_support, 
					          &m_ect_bp, &(rm.m_ect_bp));
			util::swap_support(m_ect_bp_rank10, rm.m_ect_bp_rank10, 
					          &m_ect_bp, &(rm.m_ect_bp));
			util::swap_support(m_ect_bp_select10, rm.m_ect_bp_select10,
					          &m_ect_bp, &(rm.m_ect_bp));
        }

        //! Range minimum/maximum query for the supported random access container v.
        /*!
         * \param l Leftmost position of the interval \f$[\ell..r]\f$.
         * \param r Rightmost position of the interval \f$[\ell..r]\f$.
         * \return The minimal index i with \f$\ell \leq i \leq r\f$ for which \f$ v[i] \f$ is minimal/maximal.
         * \pre
         *   - r < size()
         *   - \f$ \ell \leq r \f$
         * \par Time complexity
         *      \f$ \Order{1} \f$
         */
        size_type operator()(const size_type l, const size_type r)const {
            assert(l <= r); assert(r < size());
            if (l==r)
                return l;
            size_type x		= m_ect_bp_select10(l+1);
            size_type y		= m_ect_bp_select10(r+1);
            size_type z		= m_ect_bp_support.rmq(x, y);
            size_type f 	= z + 1 - 2*(m_ect_bp[z]);
            return m_ect_bp_rank10(f-1);
        }

        size_type size()const {
            return m_ect_bp.size()/4;
        }

        size_type serialize(std::ostream& out)const {
            size_type written_bytes = 0;
            written_bytes += m_ect_bp.serialize(out);
            written_bytes += m_ect_bp_support.serialize(out);
            written_bytes -= m_ect_bp_support.bp_select.serialize(out); // rmq_succinct_sada does not use the select support of bp_support
            written_bytes += m_ect_bp_rank10.serialize(out);
            written_bytes += m_ect_bp_select10.serialize(out);
            return written_bytes;
        }

        void load(std::istream& in) {
            m_ect_bp.load(in);
            m_ect_bp_support.load(in, &m_ect_bp);
            m_ect_bp_rank10.load(in, &m_ect_bp);
            m_ect_bp_select10.load(in, &m_ect_bp);
        }
};


}


#endif
