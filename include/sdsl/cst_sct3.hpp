/* sdsl - succinct data structures library
    Copyright (C) 2010 Simon Gog

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
/*! \file cst_sct3.hpp
    \brief cst_sct3.hpp contains an implementation of an compressed suffix tree (CST).
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_CST_SCT3
#define INCLUDED_SDSL_CST_SCT3

#include "int_vector.hpp"
#include "algorithms.hpp"
#include "iterators.hpp"
#include "lcp_support_tree.hpp"
#include "lcp_wt.hpp"
#include "bp_support.hpp"
#include "csa_wt.hpp" // for std initialization of cst_sct3 
#include "cst_iterators.hpp"
#include "cst_sct.hpp" // for lcp_interval
#include "rank_support.hpp"
#include "select_support.hpp"
#include "testutils.hpp"
#include "util.hpp"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring> // for strlen
#include <iomanip>
#include <iterator>
#include <stack> // for the calculation of the balanced parantheses sequence
#include <ostream>

namespace sdsl
{

struct cst_tag; // forward declaration

template<class Int = int_vector<>::size_type>
struct bp_interval {
    Int i; 	//!< The left border of the lcp-interval \f$\ell-[left..right]\f$.
    Int j;	//!< The right border of the lcp-interval \f$\ell-[left..right]\f$.
    Int ipos;  // position of the i+1th opening parenthesis in the balanced parentheses sequence
    Int cipos; // position of the matching closing parenthesis of the i+1th opening parenthesis in the balanced parentheses sequence
    Int jp1pos; // position of the j+2th opening parenthesis in the balanced parentheses sequence

    //! Constructor
    /*!
     */
    bp_interval(Int i=0, Int j=0, Int ipos=0, Int cipos=0, Int jp1pos=0):i(i),j(j),ipos(ipos),cipos(cipos),jp1pos(jp1pos) {};


    bool operator<(const bp_interval& interval)const {
        if (i!=interval.i)
            return i<interval.i;
        return j<interval.j;
    }

    //! Equality operator.
    /*! Two lcp-intervals are equal if and only if all their corresponding member variables have the same values.
     */
    bool operator==(const bp_interval& interval)const {
        return i==interval.i and
               j==interval.j;
    }

    //! Inequality operator.
    /*! Two lcp-intervals are not equal if and only if not all their corresponding member variables have the same values.
      */
    bool operator!=(const bp_interval& interval)const {
        return !(*this==interval);
    }

    //! Assignment operator.
    /*! \param interval The interval which should be assigned to the current object.
    */
    bp_interval& operator=(const bp_interval& interval) {
        i = interval.i;
        j = interval.j;
        ipos = interval.ipos;
        cipos = interval.cipos;
        jp1pos = interval.jp1pos;
        return *this;
    }
};


template<class Int>
inline std::ostream& operator<<(std::ostream& os, const bp_interval<Int>& interval)
{
    os<<"-["<<interval.i<<","<<interval.j<<"]("<<interval.ipos<<","<<interval.cipos<<","<<interval.jp1pos<<")";
    return os;
}


//! A class for the Compressed Suffix Tree (CST) proposed by Ohlebusch and Gog.
/*! The CST is parameterized by
 *    - a (compressed) suffix array (accessible via member \f$csa\f$, default clas is sdsl::csa_wt),
 *    - a (compressed) longest common prefix array data structure (accessible via member \f$lcp\f$, default class is sdsl::lcp_support_tree), and
 *    - a support data structure for balanced parentheses sequences (accessible via member \f$bp\_support\f$, default class is sdsl::bp_support_sada).
 *    - a rank support data structure for the bitvector which indicates the first child property
 *
 *  It also contains a sdsl::bit_vector which represents the balanced parentheses sequence of the
 *  Super-Cartesian tree of the lcp array. This bit_vector can be accessed via the member \f$bp\f$.
 *  Another sdsl::bit_vector stores information, if a node is a first child of another node. This
 *  bit_vector can be access via the member first_child_bv and takes \f$n\f$ bits.
 *
 *  A node \f$v\f$ of the csa_sct is represented by an sdsl::bp_interval . The size
 *  of the sdsl::cst_sct3 is smaller than the size of a sdsl::cst_sada since the
 *  tree topology needs only \f$2n+n=3n\f$ bits in contrast to the \f$4n\f$ bits in sdsl::cst_sada.
 *
 *  Applications of the CST: The compressed suffix tree could be used for string matching and
 *  many other application in sequence analysis. 17 applications are in the book
 *  "Algorithms on Strings, Trees, and Sequences" of Dan Gusfield.
 *  @ingroup cst
 */
template<class Csa = csa_wt<>,                        // CSA type
         class Lcp = lcp_support_tree<lcp_wt<> >,     // LCP type
         class Bp_support = bp_support_sada<>,        // for the balanced parentheses
         class Rank_support = rank_support_v5<>       // for the 'first child' bit_vector
         >
class cst_sct3
{
    public:
        typedef typename Csa::value_type						value_type;	// STL Container requirement/TODO: ist das nicht gleich node type???
        typedef cst_dfs_const_forward_iterator<cst_sct3>		const_iterator;// STL Container requirement
        typedef cst_bottom_up_const_forward_iterator<cst_sct3>	const_bottom_up_iterator;
        typedef const value_type								const_reference;
        typedef const_reference									reference;
        typedef const_reference*								pointer;
        typedef const pointer									const_pointer;
        typedef typename Csa::size_type							size_type;		// STL Container requirement
        typedef size_type										cst_size_type;
        typedef ptrdiff_t  										difference_type; // STL Container requirement
        typedef Csa												csa_type;
        typedef typename Lcp::template type<cst_sct3>::lcp_type	lcp_type;
        typedef Bp_support										bp_support_type;
        typedef typename Csa::pattern_type						pattern_type;
        typedef typename Csa::char_type							char_type;
        typedef bp_interval<size_type>							node_type; //!< Type for the nodes in the tree
        typedef Rank_support									fc_rank_support_type;

		typedef typename Csa::alphabet_type::comp_char_type		comp_char_type;
		typedef typename Csa::alphabet_type::sigma_type			sigma_type;	

        typedef typename Csa::alphabet_category					alphabet_category;
        typedef cst_tag											index_category;
    private:
        Csa 					m_csa;
        lcp_type				m_lcp;
        bit_vector				m_bp;
        bp_support_type			m_bp_support;
        bit_vector				m_first_child; // implementation note: no rank structure is needed for the first_child bit_vector, except for id()
        fc_rank_support_type	m_first_child_rank;
        sigma_type				m_sigma; 
        size_type				m_nodes;

        void copy(const cst_sct3& cst) {
            m_csa 			= cst.m_csa;
            copy_lcp(m_lcp, cst.m_lcp, *this);
            m_bp  			= cst.m_bp;
            m_bp_support 	= cst.m_bp_support;
            m_bp_support.set_vector(&m_bp);
            m_first_child 	= cst.m_first_child;
            m_first_child_rank = cst.m_first_child_rank;
            m_first_child_rank.set_vector(&m_first_child);
            m_sigma			= cst.m_sigma;
            m_nodes			= cst.m_nodes;
        }

        // Get the first l index of a [i,j] interval.
        /* I.e. given an interval [i,j], the function returns the position of the smallest entry lcp[k] with \f$ i<k\leq j \f$
         * \par Time complexity
         * 	 \f$ \Order{1} \f$
         */
        inline size_type get_first_l_index(const node_type& node, size_type& kpos, size_type& ckpos)const {
            if (node.cipos > node.jp1pos) { // corresponds to m_lcp[i] <= m_lcp[j+1]
                ckpos 	= node.jp1pos-1;
            } else { // corresponds to m_lcp[i] > m_lcp[j+1]
                ckpos	= node.cipos-1;
            }
            assert(m_bp[ckpos]==0);
            kpos	= m_bp_support.find_open(ckpos);
            return m_bp_support.rank(kpos)-1;
        }

        // Get the ith l index of a node
        // if there exists no ith l-index return node.j+1
        size_type get_ith_l_index(const node_type& node, size_type i, size_type& kpos, size_type& ckpos)const {
            uint64_t children = 0;

            assert(i > 0);
            if (node.cipos > node.jp1pos) { // corresponds to m_lcp[i] <= m_lcp[j+1]
                ckpos	= node.jp1pos-1;
            } else { // corresponds to m_lcp[i] > m_lcp[j+1]
                ckpos	= node.cipos-1;
            }
            assert(m_bp[ckpos] == 0);   // at least the first l-index should be present, i.e. node is not leaf
            if (1 == i) {
                kpos	= m_bp_support.find_open(ckpos);
                return m_bp_support.rank(kpos)-1;
            } else { // i > 1
                size_type r = ckpos - m_bp_support.rank(ckpos); // numbers of closing parentheses - 1 = index of first child in m_first_child
                if (r+1 >= i) { // if there exist more than i l-indices
                    // check if m_first_child[r-i+1..r-1] consists of zeros
                    const uint64_t* p = m_first_child.data() + (r>>6);
                    uint8_t offset = r&0x3F;

                    uint64_t w = (*p) & bit_magic::Li1Mask[offset];
                    if (w) {
                        children = offset - bit_magic::l1BP(w) + 1;
                    } else if (m_first_child.data() == p) { // w==0 and we are in the first word
                        children = offset + 2; // da bit_magic::l1BP(w)=-1 sein muesste
                    } else {
                        children = offset + 2;
                        while (p > m_first_child.data()) {
                            w = *(--p);
                            if (0==w)
                                children += 64;
                            else {
                                children += (63-bit_magic::l1BP(w));
                                break;
                            }
                        }
                        children += (w==0);
                    }
                    if (i < children) {  // there exists an ith l-index
                        ckpos -= (i-1);
                        assert(m_bp[ckpos] == 0);
                        kpos   = m_bp_support.find_open(ckpos);
                        return m_bp_support.rank(kpos)-1;
                    }
                }
                // if i >= children
                kpos = node.jp1pos;
                if (kpos < m_bp.size())
                    ckpos = m_bp_support.find_close(kpos);
                else
                    ckpos = m_bp.size();
                return node.j+1;
            }
        }

        // Get the postion of the first l index of a l-[i,j] interval in the balanced parentheses sequence.
        /* \par Time complexity
         * 		\f$ \Order{1} \f$
         */
        inline size_type get_pos_of_closing_parenthesis_of_the_first_l_index(const node_type& node)const {
            if (node.cipos > node.jp1pos) { // corresponds to m_lcp[i] <= m_lcp[j+1]
                return node.jp1pos-1;
            } else { // corresponds to m_lcp[i] > m_lcp[j+1]
                return node.cipos-1;
            }
        }

        // Get the next smaller value.
        /* \par Time complexity
         *      \f$ \Order{1} \f$
         */
        // Returns n if there is no next smaller value in [i+1..n-1]
        inline size_type nsv(size_type i, size_type ipos)const { // possible optimization: calculate also position of nsv, i.e. next ( following position cipos
            size_type cipos = m_bp_support.find_close(ipos);
            size_type result = m_bp_support.rank(cipos);
            return result;
        }

        // Get the previous smaller value.
        /*
         * \par Time complexity
         *    \f$ \Order{\frac{\sigma}{w}} \f$, where w=64 is the word size, can be implemented in \f$\Order{1}\f$ with rank and select
         */
        inline size_type psv(size_type i, size_type ipos, size_type cipos, size_type& psvpos, size_type& psvcpos)const {
            if (cipos + m_sigma >= m_bp.size()) {  // if lcp[i]==0 => psv is the 0th index by definition
                psvpos = 0;
                psvcpos = m_bp.size()-1;
                return 0;
            }
            if (m_bp[cipos+1]) {
                psvpos = m_bp_support.enclose(ipos);
                psvcpos = m_bp_support.find_close(psvpos);
                return m_bp_support.rank(psvpos)-1;
            }

            size_type r0 = cipos - m_bp_support.rank(cipos); // index of clothing parenthesis in first_child bp
            size_type next_first_child = 0;
            const uint64_t* p = m_first_child.data() + (r0>>6);
            uint64_t w = (*p) >> (r0&0x3F);
            if (w) { // if w!=0
                next_first_child = cipos + bit_magic::r1BP(w);
                if (cipos == next_first_child and m_bp[next_first_child+1]) {
                    psvpos = m_bp_support.enclose(ipos);
                    psvcpos = m_bp_support.find_close(psvpos);
                    return m_bp_support.rank(psvpos)-1;
                }
            } else {
                cipos += 64-(r0&0x3F);
                ++p;
                while (!(w=*p)) { // while w==0
                    ++p;
                    cipos += 64;
                }
                next_first_child = cipos + bit_magic::r1BP(w);
            }
            if (!m_bp[next_first_child+1]) { // if next parenthesis is a closing one
                psvcpos = next_first_child+1;
                psvpos = m_bp_support.find_open(psvcpos);
                return m_bp_support.rank(psvpos)-1;
            } else {
                psvpos = m_bp_support.enclose(m_bp_support.find_open(next_first_child));
                psvcpos = m_bp_support.find_close(psvpos);
                return m_bp_support.rank(psvpos)-1;
            }
        }

        // Range minimum query based on the rr_enclose method.
        /* \par Time complexity
         *     \f$ \Order{\rrenclose} \f$
         */
        inline size_type rmq(size_type l, size_type r)const {
            size_type i 	= m_bp_support.select(l+1);
            size_type j 	= m_bp_support.select(r+1);
            size_type fc_i 	= m_bp_support.find_close(i);
            if (j < fc_i) { // i < j < find_close(j) < find_close(i)
                return l;
            } else { // i < find_close(i) < j < find_close(j)
                size_type ec = m_bp_support.rr_enclose(i,j);
                if (ec == m_bp_support.size()) {// no restricted enclosing pair found
                    return r;
                } else { // found range restricted enclosing pair
                    return m_bp_support.rank(ec)-1; // subtract 1, as the index is 0 based
                }
            }
        }

    public:
        const Csa& csa;       				//!< The compressed suffix array the suffix tree is based on.
        const lcp_type& lcp;       			//!< The lcp array the suffix tree is based on.
        const bit_vector& bp; 				//!< The balanced parentheses sequence of the Super-Cartesian tree the suffix tree is based on.
        const bp_support_type& bp_support;	//!< The balanced parentheses sequence support for member bp.

        const bit_vector& first_child_bv;
        const fc_rank_support_type& first_child_rank;

        /*! \defgroup cst_sct3_constructors Constructors of cst_sct3 */
        /* @{ */

        //! Default Constructor
        cst_sct3(): csa(m_csa), lcp(m_lcp), bp(m_bp), bp_support(m_bp_support), first_child_bv(m_first_child),
            first_child_rank(m_first_child_rank) {}

		//! Construct CST from file_map
        cst_sct3(tMSS& file_map, const std::string& dir, const std::string& id, bool build_only_bps=false);

        //! Copy constructor
        /*!
         *  \param cst The cst_sct3 which should be copied.
         *  \par Time complexity
         *       \f$ \Order{n} \f$, where \f$n=\f$cst_sct3.size()
         */
        cst_sct3(const cst_sct3& cst):csa(m_csa),lcp(m_lcp),bp(m_bp),bp_support(m_bp_support),first_child_bv(m_first_child),
            first_child_rank(m_first_child_rank) {
            copy(cst);
        }

        /* @} */

        //! Number of leaves of the suffix tree.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         */
        size_type size()const {
            return m_bp.size()>>1;
        }

        //! Returns the largest size that cst_sct3 can ever have.
        /*! Required for the Container Concept of the STL.
         *  \sa size
         */
        static size_type max_size() {
            return Csa::max_size();
        }

        //! Returns if the data strucutre is empty.
        /*! Required for the Container Concept of the STL.A
         * \sa size
         */
        bool empty()const {
            return m_csa.empty();
        }

        //! Swap method for cst_sct3
        /*! The swap method can be defined in terms of assignment.
        	This requires three assignments, each of which, for a container type, is linear
        	in the container's size. In a sense, then, a.swap(b) is redundant.
        	This implementation guaranties a run-time complexity that is constant rather than linear.
        	\param cst cst_sct3 to swap.

        	Required for the Assignable Conecpt of the STL.
          */
        void swap(cst_sct3& cst) {
            if (this != &cst) {
                m_csa.swap(cst.m_csa);
                m_bp.swap(cst.m_bp);
                util::swap_support(m_bp_support, cst.m_bp_support, &m_bp, &(cst.m_bp));
                m_first_child.swap(cst.m_first_child);
                util::swap_support(m_first_child_rank, cst.m_first_child_rank, &m_first_child, &(cst.m_first_child));
                std::swap(m_sigma, cst.m_sigma);
                std::swap(m_nodes, cst.m_nodes);
                // anything else has to be swapped before swapping lcp
                swap_lcp(m_lcp, cst.m_lcp, *this, cst);
            }
        }

        //! Returns a const_iterator to the first element of a depth first traversal of the tree.
        /*! Required for the STL Container Concept.
         *  \sa end
         */
        const_iterator begin()const {
            if (0 == m_bp.size())  // special case: tree is uninitialized
                return end();
            else if (2 == m_bp.size()) { // special case: the root is a leaf
                return const_iterator(this, root(), true, true);
            }
            return const_iterator(this, root(), false, true);
        };

        //! Returns a const_iterator to the element after the last element of a depth first traversal of the tree.
        /*! Required for the STL Container Concept.
         *  \sa begin.
         */
        const_iterator end()const {
            return const_iterator(this, root(), true, false);
        }

        //! Returns an iterator to the first element of a bottom-up traversal of the tree.
        const_bottom_up_iterator begin_bottom_up()const {
            if (0 == m_bp.size())  // special case: tree is uninitialized
                return end_bottom_up();
            return const_bottom_up_iterator(this, leftmost_leaf_in_the_subtree(root()));
        }

        //! Returns an iterator to the element after the last element of a bottom-up traversal of the tree.
        const_bottom_up_iterator end_bottom_up()const {
            return const_bottom_up_iterator(this, root(), false);
        }

        //! Assignment Operator.
        /*!
         *	Required for the Assignable Concept of the STL.
         */
        cst_sct3& operator=(const cst_sct3& cst);

        //! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const;

        //! Load from a stream.
        /*! \param in Inputstream to load the data structure from.
         */
        void load(std::istream& in);

        /*! \defgroup cst_sct3_tree_methods Tree methods of cst_sct3 */
        /* @{ */

        //! Return the root of the suffix tree.
        /*!
         * \par Time complexity O(1)
         *      \f$ \Order{1} \f$
         */
        node_type root() const {
            return node_type(0, size()-1, 0, m_bp.size()-1, m_bp.size());
        }

        //! Decide if a node is a leaf in the suffix tree.
        /*!
         * \param v A valid node of a cst_sct3.
         * \returns A boolean value indicating if v is a leaf.
         * \par Time complexity
         *      \f$ \Order{1} \f$
         */
        bool is_leaf(const node_type& v)const {
            return v.i==v.j;
        }

        //! Return the i-th leaf (1-based from left to right) of the suffix tree.
        /*!
         * \param i 1-based position of the leaf. \f$1\leq i\leq size()\f$.
         * \return The i-th leave.
         * \par Time complexity
         *      \f$ \Order{1} \f$
         * \pre \f$ 1 \leq i \leq size() \f$
         */
        node_type ith_leaf(size_type i)const {
            assert(i > 0 and i <= size());
            size_type ipos = m_bp_support.select(i);
            size_type jp1pos;
            if (i == size())
                jp1pos = m_bp.size();
            else if (m_bp[ipos+1])
                jp1pos = ipos+1;
            else
                jp1pos = m_bp_support.select(i+1);
            return node_type(i-1, i-1, ipos, m_bp_support.find_close(ipos), jp1pos);
        }

        //! Calculate the number of leaves in the subtree rooted at node v.
        /*! \param v A valid node of the suffix tree.
         *  \return The number of leaves in the subtree rooted at node v.
         *  \par Time complexity
         *       \f$ \Order{1} \f$
         *
         *  This method is used e.g. in the sdsl::algorithm::count<Cst> method.
         */
        size_type leaves_in_the_subtree(const node_type& v)const {
            return v.j-v.i+1;
        }

        //! Calculates the leftmost leaf in the subtree rooted at node v.
        /*! \param v A valid node of the suffix tree.
         * 	\return The leftmost leaf in the subtree rooted at node v.
         *	\par Time complexity
         *		\f$ \Order{1} \f$
         */
        node_type leftmost_leaf_in_the_subtree(const node_type& v)const {
            return ith_leaf(v.i+1);
        }

        //! Calculates the rightmost leaf in the subtree rooted at node v.
        /*! \param v A valid node of the suffix tree.
         * 	\return The rightmost leaf in the subtree rooted at node v.
         *	\par Time complexity
         *		\f$ \Order{1} \f$
         */
        node_type rightmost_leaf_in_the_subtree(const node_type& v)const {
            return ith_leaf(v.j+1);
        }

        //! Calculates the index of the leftmost leaf in the corresponding suffix array.
        /*! \param v A valid node of the suffix tree.
         * 	\return The index of the leftmost leaf in the corresponding suffix array.
         *	\par Time complexity
         *		\f$ \Order{1} \f$
         *  \par Note
         *  lb is an abbreviation for ,,left bound''
         */
        size_type lb(const node_type& v)const {
            return v.i;
        }

        //! Calculates the index of the rightmost leaf in the corresponding suffix array.
        /*! \param v A valid node of the suffix tree.
         * 	\return The index of the rightmost leaf in the corresponding suffix array.
         *	\par Time complexity
         *		\f$ \Order{1} \f$
         *  \par Note
         *   rb is an abbreviation for ,,right bound''
         */
        size_type rb(const node_type& v)const {
            return v.j;
        }

        //! Calculate the parent node of a node v.
        /*! \param v A valid node of the suffix tree.
         *  \return The parent node of v or the root if v==root().
         *  \par Time complexity
         *       \f$ \Order{1}\f$
         */
        node_type parent(const node_type& v) const {
            if (v.cipos > v.jp1pos) { // LCP[i] <= LCP[j+1]
                size_type psv_pos, psv_cpos, psv_v, nsv_v, nsv_p1pos;
                psv_v = psv(v.j+1, v.jp1pos, m_bp_support.find_close(v.jp1pos), psv_pos, psv_cpos);
                nsv_v = nsv(v.j+1, v.jp1pos)-1;
                if (nsv_v == size()-1) {
                    nsv_p1pos = m_bp.size();
                } else { // nsv_v < size()-1
                    nsv_p1pos = m_bp_support.select(nsv_v+2);
                }
                return node_type(psv_v, nsv_v, psv_pos, psv_cpos, nsv_p1pos);
            } else { // LCP[i] > LCP[j+1]
                size_type psv_pos, psv_cpos, psv_v;
                psv_v = psv(v.i, v.ipos, v.cipos, psv_pos, psv_cpos);
                return node_type(psv_v, v.j, psv_pos, psv_cpos, v.jp1pos);
            }
        }

        //! Returns the next sibling of node v.
        /*!
         * \param v A valid node v of the suffix tree.
         * \return The next (right) sibling of node v or root() if v has no next (right) sibling.
         * \par Time complexity
         *      \f$ \Order{1} \f$
         */
        node_type sibling(const node_type& v)const {
//Procedure:(1) Determine, if v has a right sibling.   Entsp. parent hat gleicht rechte Grenze wie v. Speziallfall rechter Rand ist hier schon behandelt!
            if (v.cipos < v.jp1pos) { // LCP[i] > LCP[j+1] => v has the same right border as parent(v) => no right sibling
                return root();
            }
//          (2)	There exists a right sibling, LCP[j+1] >= LCP[i] and j>i
            // Now it holds:  v.cipos > v.jp1pos
            size_type cjp1posm1 = m_bp_support.find_close(v.jp1pos)-1; // v.cipos-2 ???
            // m_bp[cjp1posm1] equals 1 =>  v is the last child 
            bool last_child = m_bp[cjp1posm1];
            // otherwise if m_bp[cjp1posm1] equals 0 => we don't know if it is the last child 
            if (!last_child) {
                size_type first_child_idx = cjp1posm1 - m_bp_support.rank(cjp1posm1);
                last_child = m_first_child[first_child_idx]; // if first_child indicator is true => the new sibling is the rightmost sibling
            }
            if (last_child) {
                size_type nsv_v = nsv(v.j+1, v.jp1pos)-1, nsv_p1pos;
                if (nsv_v == size()-1) {
                    nsv_p1pos = m_bp.size();
                } else {
                    nsv_p1pos = m_bp_support.select(nsv_v+2);
                }
                return node_type(v.j+1, nsv_v, v.jp1pos, m_bp_support.find_close(v.jp1pos), nsv_p1pos);
            } else {
                size_type new_j = m_bp_support.rank(m_bp_support.find_open(cjp1posm1))-2;
                return node_type(v.j+1, new_j, v.jp1pos, m_bp_support.find_close(v.jp1pos), m_bp_support.select(new_j+2));
            }
        }

        //! Get the ith child of a node v.
        /*!
         * 	\param v A valid tree node of the cst.
         *  \param i 1-based index of the child which should be returned. \f$i \geq 1\f$.
         *  \return The i-th child node of v or root() if v has no i-th child.
         *  \par Time complexity
         *    	\f$ \Order{\frac{\sigma}{w}} \f$, where w=64 is the word size, can be implemented in \f$\Order{1}\f$ with rank and select
         *  \pre \f$ 1 \leq i \leq degree(v) \f$
         */

        node_type ith_child(const node_type& v, size_type i)const {
            assert(i > 0);
            if (is_leaf(v))  // if v is a leave, v has no child
                return root();
            if (1 == i) {
                size_type k = 0, kpos = 0, k_find_close = 0;
                // v is not a leave: v has at least two children
                k = get_first_l_index(v, kpos, k_find_close);// get first l-index k and the position of k
                return node_type(v.i, k-1, v.ipos, v.cipos, kpos);
            } else { // i > 1
                size_type k1, kpos1, k_find_close1;
                k1 = get_ith_l_index(v, i-1, kpos1, k_find_close1);
                if (k1 == v.j+1)
                    return root();
                size_type k2, kpos2, k_find_close2;
                k2 = get_ith_l_index(v, i, kpos2, k_find_close2);
                return node_type(k1, k2-1, kpos1, k_find_close1, kpos2);
            }
        }

        //! Get the number of children of a node v.
        /*!
         *  \param v A valid node v of a cst_sct3.
         *  \returns The number of children of node v.
         *  \par Time complexity
         *    	\f$ \Order{\frac{\sigma}{w}} \f$, where w=64 is the word size, can be implemented in \f$\Order{1}\f$ with rank and select
         */
        size_type degree(const node_type& v)const {
            if (is_leaf(v))  // if v is a leave, v has no child
                return 0;
            // v is not a leave: v has at least two children
            size_type r = get_pos_of_closing_parenthesis_of_the_first_l_index(v);
            /*			if( m_bp[r-1] ){// if there exists no next l-index
            				return 2;
            			}
            */
            size_type r0 = r - m_bp_support.rank(r);
            const uint64_t* p = m_first_child.data() + (r0>>6);
            uint8_t offset = r0&0x3F;

            uint64_t w = (*p) & bit_magic::Li1Mask[offset];
            if (w) {
                return offset-bit_magic::l1BP(w)+1;
            } else if (m_first_child.data() == p) { // w==0 and we are in the first word
                return offset+2; // da bit_magic::l1BP(w)=-1 sein muesste
            } else {
                size_type res = offset+2;
                while (p > m_first_child.data()) {
                    w = *(--p);
                    if (0==w)
                        res += 64;
                    else {
                        return res + (63-bit_magic::l1BP(w));
                    }
                }
                return res + (w==0);
            }
        }

        // Returns the next sibling of node v.
        // Only for tests.
        node_type sibling_naive(const node_type& v)const {
            if (v==root())
                return root();
            node_type parent = this->parent(v);
            assert(parent != v);
            size_type nr = degree(parent);
            for (size_type i=1; i <= nr; ++i)
                if (ith_child(parent, i) == v and i!=nr)
                    return ith_child(parent, i+1);
            return root();
        }

        //! Get the child w of node v which edge label (v,w) starts with character c.
        /*!
         * 	\param v A valid tree node of the cst.
         *  \param c First character of the edge label from v to the desired child.
         *  \param char_pos Reference which will hold the position (0-based) of the matching char c in the sorted text/suffix array.
         *  \return The child node w which edge label (v,w) starts with c or root() if it does not exist.
         *  \par Time complexity
         *       \f$ \Order{(\saaccess+\isaaccess) \cdot \log\sigma + \lcpaccess} \f$
         */
        node_type child(const node_type& v, const char_type c, size_type& char_pos)const {
            if (is_leaf(v))  // if v is a leaf = (), v has no child
                return root();
            // else v = ( (     ))
            comp_char_type cc = m_csa.char2comp[c];
            if (cc==0 and c!=0) // TODO: aendere char2comp so ab, dass man diesen sonderfall nicht braucht
                return root();
            size_type char_ex_max_pos = m_csa.C[((size_type)1)+cc], char_inc_min_pos = m_csa.C[cc];

            size_type d			= depth(v);

//			(1) check the first child
            char_pos = get_char_pos(v.i, d, m_csa);
            if (char_pos >= char_ex_max_pos) {// the first character of the first child interval is lex. greater than c
                // => all other first characters of the child intervals are also greater than c => no solution
                return root();
            } else if (char_pos >= char_inc_min_pos) { // i.e. char_pos < char_ex_max_pos and char_pos >= char_inc_min_pos
                return ith_child(v, 1);
            }

            size_type child_cnt 	= degree(v);

//			(2) check the last child
            char_pos = get_char_pos(v.j, d, m_csa);
            if (char_pos < char_inc_min_pos) {// the first character of the last child interval is lex. smaller than c
                // =>	all other first characters of the child intervals are also smaller than c => no solution
                return root();
            } else if (char_pos < char_ex_max_pos) { // i.e. char_pos < char_ex_max_pos and char_pos >= char_inc_min_pos
                return ith_child(v, child_cnt);
            }

// 			(3) binary search for c in the children [2..children)
            size_type l_bound = 2, r_bound = child_cnt, mid, kpos, ckpos, l_index;
            while (l_bound < r_bound) {
                mid = (l_bound + r_bound) >> 1;

                l_index = get_ith_l_index(v, mid-1, kpos, ckpos);
                char_pos = get_char_pos(l_index, d, m_csa);

                if (char_inc_min_pos > char_pos) {
                    l_bound = mid+1;
                } else if (char_ex_max_pos <= char_pos) {
                    r_bound = mid;
                } else { // char_inc_min_pos <= char_pos < char_ex_max_pos => found child
                    // we know that the child is not the last child, see (2)
                    // find next l_index: we know that a new l_index exists: i.e. assert( 0 == m_bp[ckpos-1]);
                    size_type lp1_index = m_bp_support.rank(m_bp_support.find_open(ckpos-1))-1;
                    size_type jp1pos = m_bp.size();
                    if (lp1_index-1 < size()-1) {
                        jp1pos = m_bp_support.select(lp1_index+1);
                    }
                    return node_type(l_index, lp1_index-1, kpos, ckpos, jp1pos);
                }
            }
            return root();
        }

        //! Get the child w of node v which edge label (v,w) starts with character c.
        // \sa child(node_type v, const char_type c, size_type &char_pos)
        node_type child(const node_type& v, const char_type c) {
            size_type char_pos;
            return child(v, c, char_pos);
        }

        //! Returns the d-th character (1-based indexing) of the edge-label pointing to v.
        /*! \param v The node at which the edge path ends.
         * \param d The position (1-based indexing) of the requested character on the edge path from the root to v. \f$ d > 0 \wedge d < depth(v) \f$
         * \return The character at position d on the edge path from the root to v.
         * \par Time complexity
         *       \f$ \Order{ \log\sigma + (\saaccess+\isaaccess) } \f$
         * \pre \f$ 1 \leq d \leq depth(v)  \f$
         */
        char_type edge(const node_type& v, size_type d)const {
#ifndef NDEBUG
            if (d < 1 or d > depth(v)) {
                throw std::out_of_range("OUT_OF_RANGE_ERROR: "+util::demangle(typeid(this).name())+"::edge(node_type v, size_type d). d == 0 or d > depth(v)!");
            }
#endif
            size_type 	order 	= get_char_pos(v.i, d-1, m_csa);
            sigma_type 	c_begin	= 1, c_end = m_sigma+1, mid;
            while (c_begin < c_end) {
                mid = (c_begin+c_end)>>1;
                if (m_csa.C[mid] <= order) {
                    c_begin = mid+1;
                } else {
                    c_end = mid;
                }
            }
            return m_csa.comp2char[c_begin-1];
        }

        //! Calculate the lowest common ancestor (lca) of two nodes v and w of the suffix tree.
        /*!
         * \param v The first node for which the lca with the second node should be computed.
         * \param w The second node for which the lca with the first node should be computed.
         * \return A node that is the lowest common ancestor of v and w in the suffix tree.
         * \par Time complexity
         *      \f$ \Order{\rrenclose}\   \f$
         */

        node_type lca(node_type v, node_type w)const {
            if (v.i > w.i or (v.i == w.i and v.j < w.j)) {
                std::swap(v, w);
            }
            if (v.j >= w.j) { // v encloses w or v==w
                return v;
            } else { // v.i < v.j < w.i < w.j
                size_type min_index = rmq(v.i+1, w.j);
                size_type min_index_pos 	= m_bp_support.select(min_index+1);
                size_type min_index_cpos 	= m_bp_support.find_close(min_index_pos);

                if (min_index_cpos >= m_bp.size() - m_sigma) { // if lcp[min_index]==0 => return root
                    return root();
                }
                size_type new_j = nsv(min_index, min_index_pos)-1;
                size_type new_ipos, new_icpos;
                size_type new_i = psv(min_index, min_index_pos, min_index_cpos, new_ipos, new_icpos);
                size_type jp1pos = m_bp.size();
                if (new_j < size()-1) {
                    jp1pos = m_bp_support.select(new_j+2);
                }
                return node_type(new_i, new_j, new_ipos, new_icpos, jp1pos);
            }
        }

        //! Returns the string depth of node v.
        /*!
         * \param v A valid node of a cst_sct3.
         * \return The string depth of node v.
         * \par Time complexity
         *     \f$ \Order{1} \f$ for non-leaves and \f$\Order{t_{SA}}\f$ for leaves
         */
        size_type depth(const node_type& v)const {
            if (v.i == v.j) {
                return size()-m_csa[v.i];
            } else if (v == root()) {
                return 0;
            } else {
                size_type kpos, ckpos;
                size_type l = get_first_l_index(v, kpos, ckpos);
                return m_lcp[l];
            }
        }

        //! Returns the node depth of node v
        /*!
         *  \param v A valid node of a cst_sct3.
         *  \return The node depth of node v.
         *	\par Time complexity
         *		\f$ \Order{z} \f$, where \f$z\f$ is the resulting node depth.
         */
        // TODO: can be implemented in O(1) with o(n) space. See Jansson, Sadakane, Sung, SODA 2007, "Ultra-succinct Representation of Ordered Trees"
        size_type node_depth(node_type v)const {
            size_type d = 0;
            while (v != root()) {
                ++d;
                v = parent(v);
            }
            return d;
        }

        //! Compute the suffix link of node v.
        /*!
         * \param v A valid node of a cst_sct3.
         * \return The suffix link of node v.
         * \par Time complexity
         *      \f$ \Order{ \rrenclose } \f$
         */
        node_type sl(const node_type& v)const {
            if (v == root())
                return root();
            // get interval with first char deleted
            size_type i	 = m_csa.psi[v.i];
            if (is_leaf(v)) {
                if (v.i==0 and v.j==0) // if( v.l==1 )
                    return root();
                else
                    return ith_leaf(i+1);
            }
            size_type j	 = m_csa.psi[v.j];
            assert(i < j);
            size_type min_index = rmq(i+1, j); // rmq
            size_type min_index_pos 	= m_bp_support.select(min_index+1);
            size_type min_index_cpos 	= m_bp_support.find_close(min_index_pos);
            if (min_index_cpos >= m_bp.size() - m_sigma) { // if lcp[min_index]==0 => return root
                return root();
            }
            size_type new_j = nsv(min_index, min_index_pos)-1;
            size_type new_ipos, new_icpos;
            size_type new_i = psv(min_index, min_index_pos, min_index_cpos, new_ipos, new_icpos);
            size_type jp1pos = m_bp.size();
            if (new_j < size()-1) {
                jp1pos = m_bp_support.select(new_j+2);
            }
            return node_type(new_i, new_j, new_ipos, new_icpos, jp1pos);
        }


        //! Compute the Weiner link of node v and character c.
        /*
         *  \param v A valid not of a cst_sct3.
         *  \param c The character which should be prepended to the string of the current node.
         *	\return root() if the Weiner link of (v, c) does not exist, otherwise the Weiner link is returned.
         *  \par Time complexity
         *		\f$ \Order{ t_{rank\_bwt} } \f$
         *
         */
        node_type wl(const node_type& v, const char_type c) const {
            size_type c_left	= m_csa.rank_bwt(v.i, c);
            size_type c_right	= m_csa.rank_bwt(v.j+1, c);
            if (c_left == c_right)  // there exists no Weiner link
                return root();
            if (c_left+1 == c_right)
                return ith_leaf(m_csa.C[m_csa.char2comp[c]] + c_left + 1);
            else {
                size_type left	= m_csa.C[m_csa.char2comp[c]] + c_left;
                size_type right	= m_csa.C[m_csa.char2comp[c]] + c_right - 1;
                assert(left < right);

                size_type ipos = m_bp_support.select(left+1);
                size_type jp1pos = m_bp.size();
                if (right < size()-1) {
                    jp1pos = m_bp_support.select(right+2);
                }
                return node_type(left, right, ipos,
                                 m_bp_support.find_close(ipos), jp1pos);
            }
        }

        //! Computes the suffix number of a leaf node v.
        /*! \param v A valid leaf node of a cst_sct3.
         *  \return The suffix array value corresponding to the leaf node v.
         *  \par Time complexity
         *		\f$ \Order{ \saaccess } \f$
         */
        size_type sn(const node_type& v)const {
            assert(is_leaf(v));
            return m_csa[v.i];
        }

        //! Computes a unique identification number for a node of the suffx tree in the range [0..nodes()-1]
        /*!
         *	\param v A valid node of a cst_sct3.
         *  \return A unique identification number for the node v in the range [0..nodes()-1]
         *  \par Time complexity
         *		\f$ \Order{1} \f$
         */
        size_type id(const node_type& v)const {
            if (is_leaf(v)) { // return id in the range from 0..csa.size()-1
                return v.i;
            }
            size_type ckpos; // closing parentheses of the l-index
            if (v.cipos > v.jp1pos) { // corresponds to m_lcp[i] <= m_lcp[j+1]
                ckpos 	= v.jp1pos-1;
            } else { // corresponds to m_lcp[i] > m_lcp[j+1]
                ckpos	= v.cipos-1;
            }
            assert(m_bp[ckpos]==0);
            size_type r0ckpos = ckpos-m_bp_support.rank(ckpos); // determine the rank of the closing parenthesis
            return size()+m_first_child_rank(r0ckpos);
        }

        //! Computes the node for such that id(v)=id.
        /*!
         *	\param id An id in the range [0..nodes()-1].
         *  \return A node v of the CST such that id(v)=id.
         *  \par Time complexity
         *		\f$ \Order{1} \f$ for leaves and \f$ \Order{\log size()} \f$ for inner nodes
         *  \sa id(node_type v)
         */
        node_type inv_id(size_type id) {
            if (id < size()) {  // the corresponding node is a leaf
                return ith_leaf(id+1);
            } else { // the corresponding node is a inner node
                // (1) get index of the closing parenthesis in m_first_child
                size_type r0ckpos = 0;
                {
                    //binary search for the position of the (id-size()+1)-th set bit in
                    id = id-size()+1;
                    size_type lb = 0, rb = m_bp.size(); // lb inclusive, rb exclusive
                    // invariant: arg(lb) < id, arg(rb) >= id
                    while (rb-lb > 1) {
                        size_type mid = lb + (rb-lb)/2;
                        size_type arg = m_first_child_rank(mid); // ones in the prefix [0..mid-1]
                        if (arg < id) {
                            lb = mid;
                        } else { // arg >= id
                            rb = mid;
                        }
                    }
                    r0ckpos = lb;
                }
                // (2) determine position clpos of the r0clpos-th closing parentheses in the parentheses sequence
                size_type ckpos = 0;
                {
                    // binary search for the position of the (r0ckpos+1)-th closing parenthesis
                    size_type lb = 0, rb = m_bp.size(); // lb inclusive, rb exclusive
                    // invariant: arg(lb) < r0ckpos+1,  arg(rb) >= r0ckpos+1
                    while (rb-lb > 1) {
                        size_type mid = lb + (rb-lb)/2;
                        size_type arg = mid - m_bp_support.rank(mid-1);  // zeros in the prefix [0..mid-1]
                        if (arg < r0ckpos+1) {
                            lb = mid;
                        } else { // arg >= x
                            rb = mid;
                        }
                    }
                    ckpos = lb;
                }
//				if ( m_bp[ckpos] ){
//					std::cerr<<"m_bp[ckpos] should be zero! id=" << id << std::endl;
//					std::cerr<<"r0ckpos="<<r0ckpos<<" rank_0(ckpos)="<< ckpos - m_bp_support.rank(ckpos-1)  << std::endl;
//				}
                if (ckpos == m_bp.size()-1) {
                    return root();
                }
                if (m_bp[ckpos+1]) {  // jp1pos < cipos
//					std::cout<<"case1"<<std::endl;
                    size_type jp1pos= ckpos+1;
                    size_type j 	= m_bp_support.rank(jp1pos-1)-1;
                    size_type kpos  = m_bp_support.find_open(ckpos);
                    size_type ipos	= m_bp_support.enclose(kpos);
                    size_type cipos = m_bp_support.find_close(ipos);
                    size_type i		= m_bp_support.rank(ipos-1);
                    return node_type(i, j, ipos, cipos, jp1pos);
                } else { //
//					std::cout<<"case2"<<std::endl;
                    size_type cipos = ckpos+1;
                    size_type ipos  = m_bp_support.find_open(cipos);
                    size_type i     = m_bp_support.rank(ipos-1);
                    size_type j     = nsv(i, ipos)-1;
                    size_type jp1pos= m_bp.size();
                    if (j != size()-1) {
                        jp1pos = m_bp_support.select(j+2);
                    }
                    return node_type(i, j, ipos, cipos, jp1pos);
                }
            }
        }

        //! Get the number of nodes of the suffix tree.
        size_type nodes()const {
            return m_nodes;
        }

        //! Get the node in the suffix tree which corresponds to the lcp-interval [lb..rb]
        /* \param lb Left bound of the lcp-interval [lb..rb] (inclusive).
         * \param rb Right bound of the lcp-interval [lb..rb] (inclusive).
         * \param l  Depth of the lcp-interval.
         *\ return The node in the suffix tree corresponding lcp-interval [lb..rb]
         */
        node_type node(size_type lb, size_type rb, size_type l=0) const {
            size_type ipos = m_bp_support.select(lb+1);
            size_type jp1pos;
            if (rb == size()-1) {
                jp1pos = m_bp.size();
            } else {
                jp1pos = m_bp_support.select(rb+2);
            }
            return node_type(lb, rb, ipos, m_bp_support.find_close(ipos), jp1pos);
        }

        //! Maps an index i to the position in TLCP where LCP[i] can be found
        /*!
         * \param i The index in the LCP array
         * \return The corresponding position in the TLCP array
         */
        size_type tlcp_idx(size_type i) const {
            size_type ipos 	= m_bp_support.select(i+1);
            size_type cipos = m_bp_support.find_close(ipos);
            return m_first_child_rank.rank(((ipos+cipos-1)>>1)-i);
        }
        /* @} */
};

// == template functions ==


template<class Csa, class Lcp, class Bp_support, class Rank_support>
cst_sct3<Csa, Lcp, Bp_support, Rank_support>::cst_sct3(tMSS& file_map, const std::string& dir, const std::string& id, bool build_only_bps):csa(m_csa), lcp(m_lcp), bp(m_bp), bp_support(m_bp_support), first_child_bv(m_first_child), first_child_rank(m_first_child_rank)
{
    write_R_output("cst", "construct BPS", "begin", 1, 0);
    int_vector_file_buffer<> lcp_buf(file_map["lcp"].c_str());
    m_nodes = algorithm::construct_supercartesian_tree_bp_succinct_and_first_child(lcp_buf, m_bp, m_first_child) + m_bp.size()/2;
    write_R_output("cst", "construct BPS", "end", 1, 0);
    write_R_output("cst", "construct BPSS", "begin", 1, 0);
    util::init_support(m_bp_support, &m_bp);
    util::init_support(m_first_child_rank, &m_first_child);
    write_R_output("cst", "construct BPSS", "end", 1, 0);

    if (!build_only_bps) {
        write_R_output("cst", "construct CLCP", "begin", 1, 0);
        construct_lcp(m_lcp, *this, file_map, dir, id);
        write_R_output("cst", "construct CLCP", "end", 1, 0);
    }
    if (!build_only_bps) {
        util::load_from_file(m_csa, file_map[util::class_to_hash(m_csa)].c_str());
    }
    m_sigma = degree(root());
}

template<class Csa, class Lcp, class Bp_support, class Rank_support>
typename cst_sct3<Csa, Lcp, Bp_support, Rank_support>::size_type cst_sct3<Csa, Lcp, Bp_support, Rank_support>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_csa.serialize(out, child, "csa");
    written_bytes += m_lcp.serialize(out, child, "lcp");
    written_bytes += m_bp.serialize(out, child, "bp");
    written_bytes += m_bp_support.serialize(out, child, "bp_support");
    written_bytes += m_first_child.serialize(out, child, "mark_child");
    written_bytes += m_first_child_rank.serialize(out, child, "mark_child_rank");
    written_bytes += util::write_member(m_sigma, out, child, "sigma");
    written_bytes += util::write_member(m_nodes, out, child, "node_cnt");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<class Csa, class Lcp, class Bp_support, class Rank_support>
void cst_sct3<Csa, Lcp, Bp_support, Rank_support>::load(std::istream& in)
{
    m_csa.load(in);
    load_lcp(m_lcp, in, *this);
    m_bp.load(in);
    m_bp_support.load(in, &m_bp);
    m_first_child.load(in);
    m_first_child_rank.load(in, &m_first_child);
    util::read_member(m_sigma, in);
    util::read_member(m_nodes, in);
#ifdef SDSL_DEBUG
    assert(algorithm::check_bp_support(m_bp, m_bp_support));
    std::cerr<<"checked bp_support"<<std::endl;
#endif
}

template<class Csa, class Lcp, class Bp_support, class Rank_support>
cst_sct3<Csa, Lcp, Bp_support, Rank_support>& cst_sct3<Csa, Lcp, Bp_support, Rank_support>::operator=(const cst_sct3& cst)
{
    if (this != &cst) {
        copy(cst);
    }
    return *this;
}





} // end namespace sdsl


#endif
