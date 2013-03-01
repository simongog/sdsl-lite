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
/*! \file bp_support_gg.hpp
    \brief bp_support_gg.hpp contains an implementation of a balanced parentheses support data structure.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_BP_SUPPORT_GG
#define INCLUDED_SDSL_BP_SUPPORT_GG

#include "int_vector.hpp"
#include "nearest_neighbour_dictionary.hpp"
#include "rank_support.hpp"
#include "select_support.hpp"
#include "algorithms.hpp"
#include "util.hpp"
#include "testutils.hpp"
#include <stack>
#include <map>
#include <set>
#include <utility>
#include <iostream>
#include <sstream> // for get_info method
#include <fstream>
#include <stdexcept>
#include <sstream>

namespace sdsl
{

//! A class that provides support for bit_vectors that represent a balanced parentheses sequence. Implementation was proposed by Geary et al. (CPM 2004) and extended by Ohlebusch and Gog (SPIRE 2009).
/*! This data structure supports the following methods on a bit_vector b that represents a balanced parentheses sequence:
 *    - rank
 *    - select
 *    - excess
 *    - find_open
 *    - find_close
 *    - enclose
 *    - rr_enclose
 *
 *  An opening parenthesis in the balanced parentheses sequence is represented by a 1 in the bit_vector
 *  and a closing parenthesis by a 0.
 *
 *  This class could be parametrized by four parameters:
 *    - NearestNeighbourDictionary is a class which supports rank and select with little space on sparse populated bit_vectors.
 *    - RankSupport is a class which support the rank operation on bit_vectors.
 *    - SelectSupport is a class which support the select operation on bit_vectors.
 *
 *  @ingroup bps
 */
template<class NearestNeighbourDictionary = nearest_neighbour_dictionary<30>, // TODO: make it possible that also rrr_vector can be used
         class RankSupport = rank_support_v<>,
         class SelectSupport = select_support_mcl<> >
class bp_support_gg
{
    public:
        typedef bit_vector::size_type size_type;
        typedef NearestNeighbourDictionary nnd_type;
        typedef RankSupport				   rank_support_type;
        typedef SelectSupport			   select_support_type;
        typedef bp_support_gg<nnd_type, RankSupport, select_support_bs<RankSupport> > bp_support_type;
    private:
        const bit_vector*	    m_bp;			  // the supported balanced parentheses sequence as bit_vector
        rank_support_type 		m_rank_bp;  	  // rank support for the balanced parentheses sequence => see excess() and rank()
        select_support_type		m_select_bp;      // select support for the balanced parentheses sequence => see select()

        nnd_type 				m_nnd; 			  // nearest neighbour dictionary for pioneers bit_vector

        bit_vector 				m_pioneer_bp;     // first level of recursion: balanced parentheses sequence of the pioneers
        bp_support_type* 		m_pioneer_bp_support;

        uint32_t m_block_size;
        size_type m_size;
        size_type m_blocks; // number of blocks

        void copy(const bp_support_gg& bp_support) {
            m_bp = bp_support.m_bp;
            m_rank_bp = bp_support.m_rank_bp;
            m_rank_bp.set_vector(m_bp);
            m_select_bp = bp_support.m_select_bp;
            m_select_bp.set_vector(m_bp);

            m_nnd = bp_support.m_nnd;

            m_block_size = bp_support.m_block_size;
            m_size = bp_support.m_size;
            m_blocks = bp_support.m_blocks;


            m_pioneer_bp = bp_support.m_pioneer_bp;
            if (bp_support.m_pioneer_bp_support == NULL) {
                if (m_pioneer_bp_support != NULL)
                    delete m_pioneer_bp_support;
                m_pioneer_bp_support = NULL;
            } else {
                if (m_pioneer_bp_support != NULL)
                    delete m_pioneer_bp_support;
                m_pioneer_bp_support = new bp_support_gg<nnd_type, RankSupport, select_support_bs<RankSupport> >(*(bp_support.m_pioneer_bp_support));
                assert(m_pioneer_bp_support != NULL);
                m_pioneer_bp_support->set_vector(&m_pioneer_bp);
            }
        }

    public:

        const RankSupport& bp_rank;
        const SelectSupport& bp_select;

        bp_support_gg():m_bp(NULL), m_pioneer_bp_support(NULL), m_block_size(840), m_size(0), m_blocks(0),bp_rank(m_rank_bp),bp_select(m_select_bp) {}

        //! Constructor
        // TODO: einen Konstruktor ohne das const bei bit_vector, damit bei calculate_pioneers_bitmap_succinct bp ueberschrieben werden kann?
        explicit bp_support_gg(const bit_vector* bp, uint32_t used_block_size = 840):m_bp(bp), m_pioneer_bp_support(NULL), m_block_size(used_block_size), m_size(bp==NULL?0:bp->size()), m_blocks((m_size+used_block_size-1)/used_block_size),bp_rank(m_rank_bp),bp_select(m_select_bp) {
            if (m_block_size<=2) {
                throw std::logic_error(util::demangle(typeid(this).name())+": block_size should be greater than 2!");
            }
            if (bp == NULL) { // -> m_bp == NULL
                return;
            }
            util::init_support(m_rank_bp, bp);
            util::init_support(m_select_bp, bp);
            {
                bit_vector pioneer;
                algorithm::calculate_pioneers_bitmap_succinct(*m_bp, m_block_size, pioneer);
				util::assign(m_nnd, nnd_type(pioneer));
            }

            m_pioneer_bp.resize(m_nnd.ones());
            if (m_nnd.ones() > 0  and m_nnd.ones() == m_bp->size()) { // m_bp != NULL see above
                throw std::logic_error(util::demangle(typeid(this).name())+": recursion in the construction does not terminate!");
            }

            for (size_type i=1; i<= m_nnd.ones(); ++i) { // replace this by an iterator!!! see TODO for the nnd data structure
                m_pioneer_bp[i-1] = (*m_bp)[m_nnd.select(i)];
            }

            if (m_bp->size() > 0) { // m_bp != NULL see above
                m_pioneer_bp_support = new bp_support_gg<nnd_type, RankSupport, select_support_bs<RankSupport> >(&m_pioneer_bp, m_block_size);
            }
        }

        //! Copy constructor
        bp_support_gg(const bp_support_gg& bp_support):m_pioneer_bp_support(NULL),bp_rank(m_rank_bp),bp_select(m_select_bp) {
            copy(bp_support);
        }

        //! Destructor
        ~bp_support_gg() {
            if (m_pioneer_bp_support != NULL)
                delete m_pioneer_bp_support;
        }

        //! Swap operator
        void swap(bp_support_gg& bp_support) {
            m_rank_bp.swap(bp_support.m_rank_bp);
            m_select_bp.swap(bp_support.m_select_bp);
            m_nnd.swap(bp_support.m_nnd);

            std::swap(m_block_size, bp_support.m_block_size);
            std::swap(m_size, bp_support.m_size);
            std::swap(m_blocks, bp_support.m_blocks);

            m_pioneer_bp.swap(bp_support.m_pioneer_bp);

            std::swap(m_pioneer_bp_support, bp_support.m_pioneer_bp_support);
            if (m_pioneer_bp_support != NULL) {
                m_pioneer_bp_support->set_vector(&m_pioneer_bp);
            }
            if (bp_support.m_pioneer_bp_support != NULL) {
                bp_support.m_pioneer_bp_support->set_vector(&(bp_support.m_pioneer_bp));
            }
        }

        //! Assignment operator
        bp_support_gg& operator=(const bp_support_gg& bp_support) {
            if (this != &bp_support) {
                copy(bp_support);
            }
            return *this;
        }

        void set_vector(const bit_vector* bp) {
            m_bp = bp;
            m_rank_bp.set_vector(bp);
            m_select_bp.set_vector(bp);
        }

        /*! Calculates the excess value at index i.
         * \param i The index of which the excess value should be calculated.
         */
        inline size_type excess(size_type i)const {
            return (m_rank_bp(i+1)<<1)-i-1;
        }

        /*! Returns the number of opening parentheses up to and including index i.
         * \pre{ \f$ 0\leq i < size() \f$ }
         */
        size_type rank(size_type i)const {
            return m_rank_bp(i+1);
        }

        /*! Returns the index of the i-th opening parenthesis.
         * \param i Number of the parenthesis to select.
         * \pre{ \f$1\leq i < rank(size())\f$ }
         * \post{ \f$ 0\leq select(i) < size() \f$ }
         */
        size_type select(size_type i)const {
            return m_select_bp(i);
        }

        /*! Calculate the index of the matching closing parenthesis to the parenthesis at index i.
         * \param i Index of an parenthesis. 0 <= i < size().
         * \return * i, if the parenthesis at index i is closing,
         *         * the position j of the matching closing parenthesis, if a matching parenthesis exists,
         *         * size() if no matching closing parenthesis exists.
         */
        size_type find_close(size_type i)const {
#ifdef SDSL_DEBUG_BP
            if (i >= m_size) {
                throw std::out_of_range("OUT_OF_RANGE: bp_support_gg::find_close");
            }
#endif
            if (!(*m_bp)[i]) {// if there is a closing parenthesis at index i return i
                return i;
            }
            size_type mi = 0; // match for i
            if ((mi=algorithm::near_find_closing(*m_bp, i+1, 1, m_block_size))==i) {
                const size_type i_ = m_nnd.rank(i+1)-1; // lemma that this gives us an opening pioneer
                assert(m_pioneer_bp[i_]==1); // assert that i2 is an opening parenthesis
                size_type mi_ = m_pioneer_bp_support->find_close(i_);    assert(m_pioneer_bp[mi_]==0);
                mi = m_nnd.select(mi_+1);  /* matching pioneer position in bp */ assert((*m_bp)[mi]==0);
                mi = (mi/m_block_size)*m_block_size;
//				size_type epb = excess(mi); // excess of first parenthesis in the pioneer block
                size_type epb2 = excess(mi-1); // excess of first parenthesis in the pioneer block
                const size_type ei = excess(i);  // excess at position i
                /* invariant: epb >= ei-1 */ //assert( epb+1 >= ei );
                /*				while( epb+1 != ei ){		 assert( mi < m_size );
                					if( (*m_bp)[++mi] )
                						++epb;
                					else
                						--epb;
                				}
                */
                return algorithm::near_find_closing(*m_bp, mi, epb2-ei+1, m_block_size);

            }
            return mi;
        }

        //! Calculate the matching opening parenthesis to the closing parenthesis at position i
        /*! \param i Index of a closing parenthesis.
          * \return * i, if the parenthesis at index i is closing,
          *         * the position j of the matching opening parenthesis, if a matching parenthesis exists,
          *         * size() if no matching closing parenthesis exists.
          */
        size_type find_open(size_type i)const {
#ifdef SDSL_DEBUG_BP
            if (i >= m_size) {
                std::stringstream ss;
                ss<<"i="<<i<<" m_size="<<m_size<<std::endl;
                throw std::out_of_range("OUT_OF_RANGE: bp_support_gg::find_open, "+ss.str());
            }
#endif
            if ((*m_bp)[i]) {// if there is a opening parenthesis at index i return i
                return i;
            }
            size_type mi = 0; // match for i
//			if( (mi=algorithm::near_find_open( *m_bp, i, m_block_size)) == i ){
            if ((mi=algorithm::near_find_opening(*m_bp, i-1, 1, m_block_size)) == i) {
//std::cerr<<"no near match found"<<" i="<<i<<" bp[i]="<<(*m_bp)[i]<<std::endl;
                const size_type i_ = m_nnd.rank(i); // lemma that this gives us an closing pioneer
                assert(m_pioneer_bp[i_]==0); // assert that i' is an opening parenthesis
                const size_type mi_ = m_pioneer_bp_support->find_open(i_); 		assert(m_pioneer_bp[mi_]==1);
                mi = m_nnd.select(mi_+1);  /* matching pioneer position in bp */ assert((*m_bp)[mi]==1);
                mi = (mi/m_block_size)*m_block_size + m_block_size - 1; 	assert(mi < i);
//				size_type epb = excess(mi); // excess of last parenthesis in the pioneer block
                size_type epb2 = excess(mi+1); // excess of last parenthesis in the pioneer block
                const size_type ei = excess(i);  // excess at position i
                /*invariant: epb >= ei+1*/ 	 //assert( epb >= ei+1 );
                /*				while( epb != ei ){			 assert( mi < m_size );
                					if( (*m_bp)[mi--] )
                						--epb;
                					else
                						++epb;
                				}
                				++mi;
                */
                return algorithm::near_find_opening(*m_bp, mi, epb2-ei+1-2*((*m_bp)[mi+1]), m_block_size);
            }
            return mi;
        }

        //! Calculate the index of the opening parenthesis corresponding to the closest matching parenthesis pair enclosing i.
        /*! \param i Index of an opening parenthesis.
         *  \return The index of the opening parenthesis corresponding to the closest matching parenthesis pair enclosing i,
         *          or size() if no such pair exists.
         */
        size_type enclose(size_type i)const {
#ifdef SDSL_DEBUG_BP
            if (i >= m_size) {
                throw std::out_of_range("OUT_OF_RANGE: bp_support_gg::enclose.");
            }
#endif
            if (!(*m_bp)[i]) { // if there is closing parenthesis at position i
                return find_open(i);
            }
            const size_type exi = excess(i);
            if (exi == 1)  // if i is not enclosed by a parentheses pair..
                return size();
            size_type ei; // enclose  for i
//			if( (ei=algorithm::near_enclose( *m_bp, i, m_block_size )) == i ){
            if ((ei=algorithm::near_find_opening(*m_bp, i-1, 1,  m_block_size)) == i) {
                const size_type i_ = m_nnd.rank(i); // next parenthesis in the pioneer bitmap
                size_type ei_; // enclose for i'
                ei_ = m_pioneer_bp_support->enclose(i_);
                assert(m_pioneer_bp[ei_]==1);
                ei = m_nnd.select(ei_+1);                                  assert((*m_bp)[ei]==1);
                ei = (ei/m_block_size)*m_block_size + m_block_size - 1;    assert(ei < i);
//				size_type epb = excess(ei); // excess of the last parenthesis in the pioneer block
                size_type epb2 = excess(ei+1); // excess of last parenthesis in the pioneer block
                /* invariant epb+1 >= exi */ //assert( epb+1 >= exi );
                /*				while( epb+2 != exi ){
                					if( (*m_bp)[ei--] )
                					   --epb;
                					else
                						++epb;
                				}
                				++ei;
                				*/
                return algorithm::near_find_opening(*m_bp, ei, epb2-exi+1+2*((*m_bp)[ei+1]==0), m_block_size);
            }
            return ei;
        }

        //! The range restricted enclose operation for parentheses pairs \f$(i,\mu(i))\f$ and \f$(j,\mu(j))\f$.
        /*! \param i First opening parenthesis.
         *  \param j Second opening parenthesis \f$ i<j \wedge findclose(i) < j \f$.
         *  \return The smallest index, say k, of an opening parenthesis such that findclose(i) < k < j and
         *  findclose(j) < findclose(k). If such a k does not exists, restricted_enclose(i,j) returns size().
         * \par Time complexity
         *     \f$ \Order{block\_size} \f$
         */
        size_type rr_enclose(const size_type i, const size_type j)const {
            assert(j < m_size);
            assert((*m_bp)[i]==1 and (*m_bp)[j]==1);
            const size_type mip1 = find_close(i)+1;
            if (mip1 >= j)
                return size();
            return rmq_open(mip1, j);
        }

        /*! Search the interval [l,r-1] for an opening parenthesis, say i, such that find_close(i) >= r.
         * \param l The left end (inclusive) of the interval to search for the result.
         * \param r The right end (exclusive) of the interval to search for the result.
         * \return the minimal opening parenthesis i with \f$ \ell \leq i < r \f$ and \f$ find_close(i) \geq r \f$;
         *         if no such i exists size() is returned.
         * \par Time complexity
         *     \f$ \Order{block\_size} \f$
         */
        size_type rmq_open(const size_type l, const size_type r)const {
            if (l >= r)
                return size();
            size_type		min_ex_pos = r;

            if (l/m_block_size == r/m_block_size) {
                min_ex_pos = algorithm::near_rmq_open(*m_bp, l, r);
            } else { // parentheses pair does not start in the same block
// muss nicht sein:				assert( l>=1 ); // l is at greater or equal than 1
                // note: l and r are not in the same block
                size_type		k, ex;	// helper variables
                size_type 		min_ex = excess(r) + 2*((*m_bp)[r]==0);// minimal excess


                // 1.2
                size_type l_ = m_nnd.rank(l);   //  l_ inclusive
                size_type r_ 	= m_nnd.rank(r);   // r_ exclusive

                size_type min_ex_pos_ = m_pioneer_bp_support->rmq_open(l_, r_);
                if (min_ex_pos_ < r_) {
                    k = m_nnd.select(min_ex_pos_ + 1);
                    min_ex = excess(k); min_ex_pos = k;
                } else {
                    // 1.1
                    k = algorithm::near_rmq_open(*m_bp, (r/m_block_size)*m_block_size, r);
                    if (k < r) {
                        assert(excess(k) < min_ex);
                        min_ex		= excess(k); min_ex_pos 	= k;
                    }
                }
                // 1.3
                k = algorithm::near_rmq_open(*m_bp, l, (l/m_block_size+1)*m_block_size);
                if (k < (l/m_block_size+1)*m_block_size and (ex=excess(k)) < min_ex) {
                    min_ex = ex; min_ex_pos = k;
                }
            }
            // 1.4
            if (min_ex_pos < r)
                return min_ex_pos;
            else
                return size();
        }
        /*
        		size_type _rmq_open(const size_type l, const size_type r)const{
        			if(l >= r)
        				return size();
        			size_type		min_ex_pos = r;

        			if( l/m_block_size == r/m_block_size ){
        				min_ex_pos = algorithm::near_rmq_open(*m_bp, l, r);
        			}
        			else{ // parentheses pair does not start in the same block
        // muss nicht sein				assert( l>=1 ); // l is at greater or equal than 1
        				// note: l and r are not in the same block
        				size_type		k, ex;	// helper variables
        				size_type 		min_ex = excess(r) + 2*((*m_bp)[r]==0);// minimal excess

        				// 1.3
        				k = algorithm::near_rmq_open(*m_bp, l, (l/m_block_size+1)*m_block_size );
        				if( k < (l/m_block_size+1)*m_block_size ){
        					ex = find_close(k);
        					if( ex >= r )
        						return k;
        				}

        				// 1.2
        				size_type l_ = m_nnd.rank( l ); //  l_ inclusive
        				size_type r_ 	= m_nnd.rank( r ); // r_ exclusive

        				size_type min_ex_pos_ = m_pioneer_bp_support->rmq_open(l_, r_);
        				if( min_ex_pos_ < r_ ){
        					k = m_nnd.select(min_ex_pos_ + 1 );
        					min_ex = excess(k); min_ex_pos = k;
        				}else{
        					// 1.1
        					k = algorithm::near_rmq_open(*m_bp, (r/m_block_size)*m_block_size, r);
        					if( k < r ){
        						assert(excess(k) < min_ex);
        						min_ex		= excess(k); min_ex_pos 	= k;
        					}
        				}
        			}
        			// 1.4
        			if( min_ex_pos < r )
        				return min_ex_pos;
        			else
        				return size();
        		}
        */

        //! The range restricted enclose operation
        /*! \param i Index of an opening parenthesis.
         *  \param j Index of an opening parenthesis \f$ i<j \wedge findclose(i) < j \f$.
         *  \return The minimal opening parenthesis, say k, such that \f$ findclose(i) < k < j\f$ and
         *  findclose(j) < findclose(k). If such a k does not exists, restricted_enclose(i,j) returns size().
         *  \par Time complexity
         *        \f$ \Order{size()}\f$ in the worst case.
        */
        size_type rr_enclose_naive(size_type i, size_type j)const {
            assert(j > i and j < m_size);
            assert((*m_bp)[i]==1 and (*m_bp)[j]==1);
            size_type mi = find_close(i); // matching parenthesis to i
            assert(mi > i and mi < j);
            assert(find_close(j) > j);
            size_type k = enclose(j);
            if (k == m_size or k < i) // there exists no opening parenthesis at position mi<k<j.
                return m_size;
            size_type kk;
            do {
                kk = k;
                k = enclose(k);
            } while (k != m_size and k > mi);
            return kk;
        }

        //! The range minimum query (rmq) returns the index of the parenthesis with minimal excess in the range \f$[l..r]\f$
        /*!	\param l The left border of the interval \f$[l..r]\f$ (\f$l\leq r\f$).
         *	\param r The right border of the interval \f$[l..r]\f$ (\f$l \leq r\f$).
         // TODO: bisher liefert rmq nicht immer das am weitesen rechts liegende minimum
         */
        size_type rmq(size_type l, size_type r)const {
            assert(l<=r);
            size_type m = rmq_open(l, r+1);
            if (m==size())
                return r;
            else if (m==l)
                return l;
            else { // m>l and m<=r
                assert(0 == (*m_bp)[m-1]);
                return m-1;
            }
        }

        //! The double enclose operation
        /*! \param i Index of an opening parenthesis.
         *  \param j Index of an opening parenthesis \f$ i<j \wedge findclose(i) < j \f$.
         *  \return The maximal opening parenthesis, say k, such that \f$ k<j \wedge k>findclose(j) \f$.
         *          If such a k does not exists, double_enclose(i,j) returns size().
         */
        size_type double_enclose(size_type i, size_type j)const {
            assert(j > i);
            assert((*m_bp)[i]==1 and (*m_bp)[j]==1);
            size_type k = rr_enclose(i, j);
            if (k == size())
                return enclose(j);
            else
                return enclose(k);
        }

        //! Return the number of zeros which procede position i in the balanced parentheses sequence.
        /*! \param i Index of an parenthesis.
         */
        size_type preceding_closing_parentheses(size_type i)const {
            assert(i < m_size);
            if (!i) return 0;
            size_type ones = m_rank_bp(i);
            if (ones) { // ones > 0
                assert(m_select_bp(ones) < i);
                return i - m_select_bp(ones) - 1;
            } else {
                return i;
            }
        }

        /*! The size of the supported balanced parentheses sequence.
         * \return the size of the supported balanced parentheses sequence.
         */
        size_type size() const {
            return m_size;
        }

        //! Serializes the bp_support_gg to a stream.
        /*!
         * \param out The outstream to which the data structure is written.
         * \return The number of bytes written to out.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += util::write_member(m_block_size, out, child, "block_size");
            written_bytes += util::write_member(m_size, out, child, "size");
            written_bytes += util::write_member(m_blocks, out, child, "block_cnt");

            written_bytes += m_rank_bp.serialize(out, child, "bp_rank");
            written_bytes += m_select_bp.serialize(out, child, "bp_select");
            written_bytes += m_nnd.serialize(out, child, "nearest_neighbour_dictionary");

            written_bytes += m_pioneer_bp.serialize(out, child, "pioneer_bp");
            if (m_bp != NULL and m_bp->size() > 0)
                written_bytes += m_pioneer_bp_support->serialize(out, child, "pioneer_bp_support");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load the bp_support_gg for a bit_vector v.
        /*!
         * \param in The instream from which the data strucutre is read.
         * \param bp Bit vector representing a balanced parentheses sequence that is supported by this data structure.
         */
        void load(std::istream& in, const bit_vector* bp) {
            m_bp = bp;
            util::read_member(m_block_size, in);
            util::read_member(m_size, in);
            util::read_member(m_blocks, in);

            m_rank_bp.load(in, m_bp);
            m_select_bp.load(in, m_bp);
            m_nnd.load(in);

            m_pioneer_bp.load(in);
            if (m_pioneer_bp_support != NULL) {
                delete m_pioneer_bp_support;
                m_pioneer_bp_support = NULL;
            }
            if (m_bp != NULL and m_bp->size() > 0) {
                m_pioneer_bp_support = new bp_support_gg<nnd_type, RankSupport, select_support_bs<RankSupport> >();
                m_pioneer_bp_support->load(in, &m_pioneer_bp);
            }
        }

        std::string get_info()const {
            std::stringstream ss;
            ss<<"number of parentheses: "<<m_bp->size()<<std::endl;

            return ss.str();
        }
};

}// end namespace




#endif
