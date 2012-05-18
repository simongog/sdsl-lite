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
/*! \file rmq_support_sparse_table.hpp
    \brief rmq_support_sparse_table.hpp contains the class rmq_support_sparse_table which supports range minimum or range maximum queries on a random access container in constant time and \f$\Order{n\log^2 n} bits\f$ space. See paper LCP revisted of Bender and Farach.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_RMQ_SUPPORT_SPARSE_TABLE
#define INCLUDED_SDSL_RMQ_SUPPORT_SPARSE_TABLE

#include "rmq_support.hpp"
#include "int_vector.hpp"
#include "bitmagic.hpp"
#include <ostream>

//! Namespace for the succinct data structure library.
namespace sdsl
{


template<class RandomAccessContainer = int_vector<>, bool Minimum=true>
class rmq_support_sparse_table;


// see http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2002/n1406.pdf for a proposal for a better solution
template<class RandomAccessContainer = int_vector<> >
struct range_maximum_support_sparse_table {
    typedef rmq_support_sparse_table<RandomAccessContainer, false> type;
};

//! A class to support range minimum or range maximum queries on a random access container.
/*!
 * This class takes two template parameters. The first on is the type of random access container, for which
 * the range minimum/maximum support should be build. The second one specifies whether the data structure
 * should answer range minimum queries (Minimum=true) or range maximum queries (Maximum=false).
 *
 * \par Time complexity
 *		\f$ \Order{1} \f$ for the range minimum/maximum queries.
 * \par Space complexity:
 *      \f$ \Order{n\log^2 n} \f$ bits for the data structure ( \f$ n=size() \f$ ). We used bit compression to get a good result in practice.
 */
template<class RandomAccessContainer, bool Minimum>
class rmq_support_sparse_table
{
        const RandomAccessContainer* m_v; // pointer to the supported random access container
        bit_vector::size_type 		m_k; // size of m_table
        int_vector<>*	   			m_table;
        typedef min_max_trait<RandomAccessContainer, Minimum> mm_trait;

        void construct() {
            if (m_v == NULL)
                return;
            const size_type n = m_v->size();
            if (n < 2)  // for n<2 the queries could be answerd without any table
                return;
            size_type k=0;
            while (2*(1ULL<<k) < n) ++k;  // calculate maximal
            if (!(m_table == NULL))
                delete [] m_table;
            m_table = new int_vector<>[k];
            m_k = k;
            for (size_type i=0; i<k; ++i) {
                m_table[i] = int_vector<>(n-(1<<(i+1))+1, 0, i+1);
            }
            for (size_type i=0; i<n-1; ++i) {
                if (!mm_trait::compare((*m_v)[i], (*m_v)[i+1]))
                    m_table[0][i] = 1;
//				std::cerr<<i+m_table[0][i]<<" ";
            }
//			std::cerr<<std::endl;
            for (size_type i=1; i<k; ++i) {
                for (size_type j=0; j<m_table[i].size(); ++j) {
                    m_table[i][j] = mm_trait::compare((*m_v)[j+m_table[i-1][j]], (*m_v)[j+(1<<i)+m_table[i-1][j+(1<<i)]]) ? m_table[i-1][j] : (1<<i)+m_table[i-1][j+(1<<i)];
//					std::cerr<<j+m_table[i][j]<<" ";
                }
//				std::cerr<<std::endl;
            }
        }

        void copy(const rmq_support_sparse_table& rm) {
            m_v = rm.m_v;
            m_k = rm.m_k;
            if (m_table != NULL) {
                delete [] m_table;
                m_table = NULL;
            }
            if (m_k > 0) {
                m_table = new int_vector<>[m_k];
                for (size_type i=0; i<m_k; ++i)
                    m_table[i] = rm.m_table[i];
            }
        }

    public:
        typedef typename RandomAccessContainer::size_type size_type;
        typedef typename RandomAccessContainer::size_type value_type;

        rmq_support_sparse_table(const RandomAccessContainer* v=NULL):m_v(v), m_k(0), m_table(NULL) {
            construct();
        }

        //! Copy constructor
        rmq_support_sparse_table(const rmq_support_sparse_table& rm) {
            if (this != &rm) { // if v is not the same object
                copy(rm);
            }
        }


        ~rmq_support_sparse_table() {
            if (m_table != NULL)
                delete [] m_table;
        }

        rmq_support_sparse_table& operator=(const rmq_support_sparse_table& rm) {
            if (this != &rm) {
                copy(rm);
            }
            return *this;
        }

        void set_vector(const RandomAccessContainer* v) {
            m_v = v;
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
            if (l+1 == r)
                return mm_trait::compare((*m_v)[l],(*m_v)[r]) ? l : r;
            size_type k = bit_magic::l1BP(r-l);
            const size_type rr = r-(1<<k)+1;
            return mm_trait::compare((*m_v)[l+m_table[k-1][l]], (*m_v)[rr+m_table[k-1][rr]]) ? l+m_table[k-1][l] : rr+m_table[k-1][rr];
        }

        size_type size()const {
            if (m_v == NULL)
                return 0;
            else
                return m_v->size();
        }

        size_type serialize(std::ostream& out)const {
            size_type written_bytes = 0;
            out.write((char*)&m_k, sizeof(m_k));
            written_bytes += sizeof(m_k);
            if (m_k > 0) {
                assert(m_table != NULL);
                for (size_type i=0; i < m_k; ++i)
                    written_bytes += m_table[i].serialize(out);
            }
            return written_bytes;
        }

        void load(std::istream& in, const RandomAccessContainer* v) {
            set_vector(v);
            in.read((char*)&m_k, sizeof(m_k));
            if (m_k >0) {
                if (m_table != NULL)
                    delete [] m_table;
                m_table = new int_vector<>[m_k];
                for (size_type i=0; i < m_k; ++i)
                    m_table[i].load(in);
            }
        }
};

}// end namespace sds

#endif // end file 
