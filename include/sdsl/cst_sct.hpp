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
/*! \file cst_sct.hpp
    \brief cst_sct.hpp contains an implementation of an compressed suffix tree (CST).
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_CST_SCT
#define INCLUDED_SDSL_CST_SCT

#include "int_vector.hpp"
#include "algorithms.hpp"
#include "iterators.hpp"
#include "lcp_support_sada.hpp"
#include "bp_support.hpp"
#include "csa_wt.hpp" // for std initialization of cst_sct 
#include "csa_bitcompressed.hpp"
#include "csa_wt.hpp"
#include "cst_iterators.hpp"
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

//! A struct for the representation of an lcp-interval \f$\ell-[left..right]\f$.
/*! lcp intervals were introduced by Abouelhoda, Kurtz, and Ohlebusch 2004 in the article
 *  ,,Replacing Suffix Trees with Enhanced Suffix Arrays''.
 *
 *  Given a lcp array \f$LCP\f$ of length \f$n\f$. An interval \f$[left..right]\f$, where \f$1\leq i<j\leq n\f$
 *  is called a lcp-interval of lcp-value \f$\ell\f$ (denoted by \f$\ell-[left..right]\f$) if
 *    -# \f$LCP[left]<\ell\f$
 *    -# \f$LCP[k]\geq \ell\f$ for all \f$k\f$ with \f$left+1\leq k \leq right\f$
 *    -# \f$LCP[k] = \ell\f$ for all \f$k\f$ with \f$left+1\leq k \leq right\f$
 *    -# \f$LCP[right+1]<\ell\f$
 *
 *  In addition to that, we define the lcp-inteval if \f$left=right\f$. In this case \f$\ell=LCP[left]\f$.
 *
 *  Lcp-intervals are used to represent nodes of the compressed suffix tree sdsl::cst_sct.
 */
template<class Int = size_t>
struct lcp_interval {
    Int l; 		//!< The \f$\ell\f$ value of the lcp-interval \f$\ell-[left..right]\f$.
    Int left; 	//!< The left border of the lcp-interval \f$\ell-[left..right]\f$.
    Int right;	//!< The right border of the lcp-interval \f$\ell-[left..right]\f$.

    //! Constructor
    /*!
     *  \param l The \f$\ell\f$ value of the interval lcp-interval \f$\ell-[left..right]\f$.
     *  \param left The left border of the interval lcp-interval \f$\ell-[left..right]\f$.
     *  \param right The right border of the interval lcp-interval \f$\ell-[left..right]\f$.
     */
    lcp_interval(Int l=0, Int left=0, Int right=0):l(l),left(left),right(right) {};


    bool operator<(const lcp_interval& interval)const {
        if (l!=interval.l)
            return l<interval.l;
        if (left!=interval.left)
            return left<interval.left;
        return right<interval.right;

    }

    //! Equality operator.
    /*! Two lcp-intervals are equal if and only if all their corresponding member variables have the same values.
     */
    bool operator==(const lcp_interval& interval)const {
        return l==interval.l and
               left==interval.left and
               right==interval.right;
    }

    //! Inequality operator.
    /*! Two lcp-intervals are not equal if and only if not all their corresponding memeber variables have the same values.
      */
    bool operator!=(const lcp_interval& interval)const {
        return !(*this==interval);
    }

    //! Assignment operator.
    /*! \param interval The interval which should be assigned to the current object.
    */
    lcp_interval& operator=(const lcp_interval& interval) {
        l = interval.l;
        left = interval.left;
        right = interval.right;
        return *this;
    }
};

template<class Int>
inline std::ostream& operator<<(std::ostream& os, const lcp_interval<Int>& interval)
{
    os<<interval.l<<"-["<<interval.left<<","<<interval.right<<"]";
    return os;
}



//! A class for the Compressed Suffix Tree (CST) proposed by Ohlebusch and Gog.
/*! The CST is parameterized by
 *    - a (compressed) suffix array (accessible via member \f$csa\f$, default clas is sdsl::csa_wt),
 *    - a (compressed) longest common prefix array data structure (accessible via member \f$lcp\f$, default class is sdsl::lcp_support_sada), and
 *    - a support data structure for balanced parentheses sequences (accessible via member \f$bp\_support\f$, default class is sdsl::bp_support_sada).
 *
 *  It also contains a sdsl::bit_vector which represents the balanced parentheses sequence of the
 *  Super-Cartesian tree of the lcp array. This bit_vector can be accessed via the member \f$bp\f$.
 *
 *  A node \f$v\f$ of the csa_sct is represented by an sdsl::lcp_interval . The size
 *  of the sdsl::cst_sct is smaller than the size of a sdsl::cst_sada since the
 *  tree topology needs only \f$2n\f$ bits in contrast to the \f$4n\f$ bits in sdsl::cst_sada.
 *  See the diagram below for an example which gives us the sizes
 *  for the compressed suffix arrays for the file english.50MB from the pizza and chili website.
 *
 *	\image html ./example/size_of_csts.png "Size comparison between cst_sct and cst_sada"
 *	\image latex ./example/size_of_csts.eps  "Size comparison between cst_sct and cst_sada" width=0.7\textwidth
 *
 *  Applications of the CST: The compressed suffix tree could be used for string matching and
 *  many other application in sequence analysis. 17 applications are in the book
 *  "Algorithms on Strings, Trees, and Sequences" of Dan Gusfield.
 *  @ingroup cst
 */
template<class Csa = csa_wt<>, class Lcp = lcp_support_sada<Csa>, class Bp_support = bp_support_sada<> >
class cst_sct
{
    public:
        typedef typename Csa::value_type							 value_type;	// STL Container requirement/TODO: ist das nicht gleich node type???
        typedef cst_dfs_const_forward_iterator<cst_sct>	 const_iterator;// STL Container requirement
        typedef const_iterator 										 iterator;		// STL Container requirement
        typedef cst_bottom_up_const_forward_iterator<cst_sct>		 const_bottom_up_iterator;
        typedef  const_bottom_up_iterator							 bottom_up_iterator;
        typedef const value_type									 const_reference;
        typedef const_reference										 reference;
        typedef const_reference*									 pointer;
        typedef const pointer										 const_pointer;
        typedef typename Csa::size_type								 size_type;		// STL Container requirement
        typedef size_type											 cst_size_type;
        typedef ptrdiff_t  											 difference_type; // STL Container requirement
        typedef Csa													 csa_type;
        typedef Lcp													 lcp_type;
        typedef Bp_support											 bp_support_type;
        typedef typename Csa::pattern_type							 pattern_type;
        typedef typename Csa::char_type								 char_type;
        typedef lcp_interval<size_type>								 node_type; //!< Type for the nodes in the tree

        typedef cst_tag												 index_category;
    private:
        Csa 			m_csa;
        Lcp 			m_lcp;
        bit_vector		m_bp;
        Bp_support		m_bp_support;
        size_type		m_nodes;

        void copy(const cst_sct& cst) {
            m_csa = cst.m_csa;
//			lcp_trait<Csa, Lcp>::copy_lcp(m_lcp, cst.m_lcp, m_csa);
            copy_lcp(m_lcp, cst.m_lcp, *this);
            m_bp = cst.m_bp;
            m_bp_support = cst.m_bp_support;
            m_bp_support.set_vector(&m_bp);
            m_nodes = cst.m_nodes;
        }

        // Get the first l index of a l-[i,j] interval.
        /* I.e. given an interval [i,j], the function returns the position of the smallest entry lcp[k] with \f$ i<k\leq j \f$
         * Time complexity: O(1) + O(lcp access)
         */ // TODO: if interval is last child interval -> we don't need the m_lcp comparisson
        size_type get_first_l_index(size_type i, size_type j, size_type& kpos, size_type& find_close_k)const {
            if (j+1 != m_csa.size() and m_lcp[i] <= m_lcp[j+1]) {
                size_type jp1pos = m_bp_support.select(j+2);
                assert(m_bp[jp1pos] == 1);  // opening parenthesis
                assert(jp1pos > 0 and m_bp[jp1pos-1] == 0); // closing parenthesis
                find_close_k = jp1pos-1;
                kpos 	= m_bp_support.find_open(find_close_k);
                return m_bp_support.rank(kpos)-1;
            } else {
                size_type ipos = m_bp_support.select(i+1);
                assert(m_bp[ipos] == 1);  // opening parenthesis
                ipos = m_bp_support.find_close(ipos);
                assert(m_bp[ipos] == 0 and m_bp[ipos-1] == 0);   // closing parenthesis
                find_close_k = ipos-1;
                kpos	= m_bp_support.find_open(find_close_k);
                return m_bp_support.rank(kpos)-1;
            }
        }

        // Get the postion of the first l index of a l-[i,j] interval in the balanced parentheses sequence.
        /* \par Time complexity
         * 		\f$ \Order{\lcpaccess} \f$
         */
        size_type get_pos_of_closing_parenthesis_of_the_first_l_index(size_type i, size_type j)const {
            if (j+1 != m_csa.size() and m_lcp[i] <= m_lcp[j+1]) {
                size_type jp1pos = m_bp_support.select(j+2);
                assert(m_bp[jp1pos] == 1);  // opening parenthesis
                assert(jp1pos > 0 and m_bp[jp1pos-1] == 0); // closing parenthesis
                return jp1pos-1;
            } else {
                size_type ipos = m_bp_support.select(i+1);
                assert(m_bp[ipos] == 1);  // opening parenthesis
                ipos = m_bp_support.find_close(ipos);
                assert(m_bp[ipos] == 0 and m_bp[ipos-1] == 0);   // closing parenthesis
                return ipos-1;
            }
        }

        // Get the next smaller value.
        /* \par Time complexity
         *      \f$ \Order{1} \f$
         */
        size_type nsv(size_type i)const {
            size_type pos =  m_bp_support.select(i+1);
            size_type rank = m_bp_support.find_close(pos);
            size_type result = m_bp_support.rank(rank);
            return result;
        }

        // Get the previous smaller value.
        /* \pre lcp[i] > 0!
         * \par Time complexity
         *     \f$ \Order{\sigma \cdot \lcpaccess} \f$
         * \sa psv
         */
        size_type psv_linear(size_type i, size_type lcpi)const {
            assert(lcpi>0);
            size_type pos =  m_bp_support.select(i+1);
            size_type ec = m_bp_support.enclose(pos);
            size_type j = m_bp_support.rank(ec);
            while (lcpi == m_lcp[j-1]) {
                ec = m_bp_support.enclose(ec);
                j = m_bp_support.rank(ec);
            }
            return j-1;
        }

        // Get the previous smaller value.
        /* \pre lcp[i] > 0!
         * \par Time complexity
         *    \f$ \Order{\log\sigma \cdot \lcpaccess} \f$
         */
        size_type psv(size_type i, size_type lcpi)const {
            assert(lcpi!=0);// if(m_lcp[i]==0) return 0;
            size_type begin = m_bp_support.find_close(m_bp_support.select(i+1)); // calculate NSV
            size_type end = m_bp_support.rank(begin)+1;                          //     "      " // exclusive
            if (end > m_csa.size()) {
                assert(end == m_csa.size()+1);
                end = m_bp.size();
            } else {
                end = m_bp_support.select(end);
            }
            if (++begin < end) {
                size_type result = m_bp_support.rank(m_bp_support.find_open(end-1))-1;
                if (m_lcp[result] < lcpi) { // TODO: explain this condition
                    // binary search
                    --end;
                    while (begin!=end) {
                        size_type mid, temp;
                        mid = (begin+end)>>1; // begin <= mid < end
                        if (m_lcp[temp = m_bp_support.rank(m_bp_support.find_open(mid))-1] >= lcpi) {
                            begin = mid+1;
                        } else {
                            result = temp;
                            end = mid;
                        }
                    }
                    return result;
                }
            }

            return m_bp_support.rank(m_bp_support.enclose(end))-1;
        }

        // Range minimum query based on the rr_enclose method.
        /* \par Time complexity
         *     \f$ \Order{\rrenclose} \f$
         */
        size_type rmq(size_type l, size_type r)const {
            size_type i 	= m_bp_support.select(l+1);
            size_type j 	= m_bp_support.select(r+1);
            size_type fc_i 	= m_bp_support.find_close(i);
            if (j < fc_i) { // i < j < find_close(j) < find_close(i)
                return l;
            } else { // i < find_close(i) < j < find_close(j)
                size_type ec = m_bp_support.rr_enclose(i,j);
                if (ec == m_bp_support.size()) {// no restricted enclosing pair found
                    return r;
                } else { // found range restriced enclosing pair
                    return m_bp_support.rank(ec)-1; // subtract 1, as the index is 0 based
                }
            }
        }

    public:
        const csa_type& csa;       			//!< The compressed suffix array the suffix tree is based on.
        const lcp_type& lcp;       			//!< The lcp array the suffix tree is based on.
        const bit_vector& bp; 			//!< The balanced parentheses sequence of the Super-Cartesian tree the suffix tree is based on.
        const bp_support_type& bp_support;	//!< The balanced parentheses sequence support for member bp.

        /*! \defgroup cst_sct_constructors Constructors of cst_sct */
        /* @{ */

        //! Default Constructor
        cst_sct(): csa(m_csa), lcp(m_lcp), bp(m_bp), bp_support(m_bp_support) {}

        // Constructor for the cst_sct taking a string for that the compressed suffix tree should be calculated.
        /*
         * \param str Text for which the \f$ \CST \f$ should be constructed. The text should be terminated by a zero byte.
         * \pre The text has to be terminated by a zero byte.
         * \sa cst_sct::cst_sct(const Csa &csa, const unsigned char *str)
         */
//		cst_sct(const unsigned char *str);

        template<uint8_t int_width, class size_type_class, uint8_t int_width_1, class size_type_class_1, uint8_t int_width_2, class size_type_class_2>
        cst_sct(const std::string& cst_file_name,
                int_vector_file_buffer<int_width, size_type_class>& lcp_buf,
                int_vector_file_buffer<int_width_1, size_type_class_1>& sa_buf,
                int_vector_file_buffer<int_width_2, size_type_class_2>& isa_buf,
                std::string dir="./",
                bool build_only_bps=false);

        cst_sct(tMSS& file_map, const std::string& dir, const std::string& id, bool build_only_bps=false);

        //! Copy constructor
        /*!
         *  \param cst The cst_sct which should be copied.
         *  \par Time complexity
         *       \f$ \Order{n} \f$, where \f$n=\f$cst_sct.size()
         */
        cst_sct(const cst_sct& cst):csa(m_csa),lcp(m_lcp),bp(m_bp),bp_support(m_bp_support) {
            copy(cst);
        }

        /* @} */

        void construct(tMSS& file_map, const std::string& dir, const std::string& id, bool build_only_bps=false);

        //! Default Destructor
        ~cst_sct() {}

        //! Number of leafs of the suffix tree.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         */
        size_type size()const {
            return m_csa.size();
        }

        //! Returns the largest size that cst_sct can ever have.
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

        //! Swap method for cst_sct
        /*! The swap method can be defined in terms of assignment.
        	This requires three assignments, each of which, for a container type, is linear
        	in the container's size. In a sense, then, a.swap(b) is redundant.
        	This implementation guaranties a run-time complexity that is constant rather than linear.
        	\param cst cst_sct to swap.

        	Required for the Assignable Conecpt of the STL.
          */
        void swap(cst_sct<Csa, Lcp>& cst);

        //! Returns a const_iterator to the first element.
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

        //! Returns a const_iterator to the element after the last element.
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
        cst_sct& operator=(const cst_sct& cst);

        //! Equality Operator
        /*! Two Instances of cst_sct are equal if
         *  all their members are equal.
         *  \par Required for the Equality Comparable Concept of the STL.
         *  \sa operator!=
         */
        bool operator==(const cst_sct& cst)const;

        //! Unequality Operator
        /*! Two Instances of cst_sct are equal if
         *  not all their members are equal.
         *  \par Required for the Equality Comparable Concept of the STL.
         *  \sa operator==
         */
        bool operator!=(const cst_sct& cst)const;

        //! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out) const;

        //! Load from a stream.
        /*! \param in Inputstream to load the data structure from.
         */
        void load(std::istream& in);

        /*! \defgroup cst_sct_tree_methods Tree methods of cst_sct */
        /* @{ */

        //! Return the root of the suffix tree.
        /*!
         * \par Time complexity O(1)
         *      \f$ \Order{1} \f$
         */
        node_type root() const {
            return node_type(0, 0, m_csa.size()-1);
        }

        //! Decide if a node is a leaf in the suffix tree.
        /*!
         * \param v A valid node of a cst_sct.
         * \returns A boolean value indicating if v is a leaf.
         * \par Time complexity
         *      \f$ \Order{1} \f$
         */
        bool is_leaf(const node_type& v)const {
            return v.right==v.left;
        }

        //! Return the i-th leaf (1-based from left to right) of the suffix tree.
        /*!
         * \param i 1-based position of the leaf. \f$1\leq i\leq csa.size()\f$.
         * \return The i-th leave.
         * \par Time complexity
         *      \f$ \Order{\saaccess} \f$
         * \pre \f$ 1 \leq i \leq csa.size() \f$
         */
        node_type ith_leaf(size_type i)const {
            assert(i > 0 and i <= m_csa.size());
            return node_type(m_csa.size()-m_csa[i-1], i-1, i-1);
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
            return v.right-v.left+1;
        }

        //! Calculates the leftmost leaf in the subtree rooted at node v.
        /*! \param v A valid node of the suffix tree.
         * 	\return The leftmost leaf in the subtree rooted at node v.
         *	\par Time complexity
         *		\f$ \Order{\saaccess} \f$
         */
        node_type leftmost_leaf_in_the_subtree(const node_type& v)const {
            return node_type(m_csa.size()-m_csa[v.left], v.left, v.left);
        }

        //! Calculates the rightmost leaf in the subtree rooted at node v.
        /*! \param v A valid node of the suffix tree.
         * 	\return The rightmost leaf in the subtree rooted at node v.
         *	\par Time complexity
         *		\f$ \Order{\saaccess} \f$
         */
        node_type rightmost_leaf_in_the_subtree(const node_type& v)const {
            return node_type(m_csa.size()-m_csa[v.right], v.right, v.right);
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
            return v.left;
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
            return v.right;
        }

        //! Calculate the parent node of a node v.
        /*! \param v A valid node of the suffix tree.
         *  \return The parent node of v or the root if v==root().
         *  \par Time complexity
         *       \f$ \Order{\lcpaccess \cdot \log\sigma}\f$
         */
        node_type parent(const node_type& v) const {
#ifndef NDEBUG
//			std::cout<<"parent interval of "<<v.l<<"-["<<v.left<<" "<<v.right<<"]"<<std::endl;
#endif
            if (v.l <= 1) {// if v.l < 1 => root()
                return root();  // return root node
            } else { // l > 1
                if (v.right+1 == m_csa.size()) {
                    // search for the previous smaller value
                    size_type lcpi = m_lcp[v.left];
                    if (lcpi==0)
                        return root();
                    assert(lcpi!=0);
                    return node_type(lcpi, psv(v.left, lcpi), v.right);
                } else if (v.left == 0) {
                    // search for the next smaller value
                    return node_type(m_lcp[v.right+1],0, nsv(v.right+1)-1);
                } else {
                    // TODO: speed up, only two cases?
                    size_type lcpl = m_lcp[v.left], lcpr = m_lcp[v.right+1];
                    if (lcpl < lcpr) {
                        return node_type(lcpr, v.left, nsv(v.right+1)-1);
                    } else if (lcpl > lcpr) {
                        return node_type(lcpl, psv(v.left, lcpl), v.right);
                    } else {
                        if (lcpl==0)
                            return root();
                        return node_type(lcpl, psv(v.left, lcpl), nsv(v.right+1)-1);
                    }
                }
            }
        }


        //! Get the ith child of a node v.
        /*!
         * 	\param v A valid tree node of the cst.
         *  \param i 1-based index of the child which should be returned. \f$i \geq 1\f$.
         *  \return The i-th child node of v or root() if v has no i-th child.
         *  \par Time complexity
         *		  \f$ \Order{i \cdot \lcpaccess}\f$, \f$ i \leq \sigma \f$
         *  \pre \f$ 1 \leq i \leq degree(v) \f$
         *  \note Could be implemented in time \f$\Order{\log i \cdot \lcpaccess}\f$ with binary search.
         */
        node_type ith_child(const node_type& v, size_type i)const {
            if (is_leaf(v))  // if v is a leave, v has no child
                return root();
            size_type k = 0, kpos = 0, k_find_close = 0;
            // v is not a leave: v has at least two children
            k = get_first_l_index(v.left, v.right, kpos, k_find_close);// get first l-index k and the position of k
            assert(m_lcp[k]==v.l);
            size_type left = v.left, right = k-1;
            size_type j=1;
            while (j < i) {
                // get next l-index m and the position of m
                size_type m    = m_bp_support.rank(m_bp_support.find_open(k_find_close-j))-1;
                ++j;
                if (m==k or v.l != m_lcp[m]) { // no next l-index exists
                    if (i!=j)
                        return root();
                    left 	= k;
                    right 	= v.right;
                    break;
                }
                left  = k;
                right = m-1;
                k = m;
            }
            if (left==right) { // child interval is a leaf
                assert(m_csa.size()-m_csa[left] > v.l);
                return node_type(m_csa.size()-m_csa[left],left,right);
            } else {
                // find l of the first child interval
                size_type k1 = get_first_l_index(left, right, kpos, kpos);
                assert(m_lcp[k1] > v.l);
                return node_type(m_lcp[k1], left, right);
            }
            return root();
        }

        //! Get the number of children of a node v.
        /*!
         *  \param v A valid node v of a cst_sct.
         *  \returns The number of children of node v.
         *  \par Time complexity
         *       \f$ \Order{\lcpaccess \cdot \log\sigma} \f$
         */
        size_type degree(const node_type& v)const {
            if (is_leaf(v))  // if v is a leave, v has no child
                return 0;
            // v is not a leave: v has at least two children
            size_type r = get_pos_of_closing_parenthesis_of_the_first_l_index(v.left, v.right);
            assert(m_lcp[m_bp_support.rank(m_bp_support.find_open(r))-1]==v.l);
            size_type begin = m_bp_support.preceding_closing_parentheses(r);
            assert(m_csa.sigma>0);
            if (begin > (size_type)m_csa.sigma-1) begin = m_csa.sigma-1;
            begin = r-begin;
            size_type end = r, mid;
            while (end != begin) {
                mid = (begin+end)>>1;  // begin <= mid < end
                if (m_lcp[m_bp_support.rank(m_bp_support.find_open(mid))-1] > v.l)
                    begin = mid+1;
                else
                    end = mid;
            }
            return r-end+2;
        }

        //! Returns the next sibling of node v.
        /*!
         * \param v A valid node v of the suffix tree.
         * \return The next (right) sibling of node v or root() if v has no next (right) sibling.
         * \par Time complexity
         *      \f$ \Order{\lcpaccess} \f$
         */
        node_type sibling(const node_type& v)const {
            if (v == root())
                return root();
            if (v.right == m_csa.size()-1) // there exists no next sibling if the interval touches the right border
                return root();
            // now let [p..q] be the interval [v.left, v.right]
            size_type lcp_p   = m_lcp[v.left];
            size_type lcp_qp1 = m_lcp[v.right+1]; // is defined as v.right < m_csa.size()-1
            if (lcp_p > lcp_qp1)
                return root();
            else { // search next l-index
                size_type new_q;
                size_type qp1_closing_pos = m_bp_support.find_close(m_bp_support.select(v.right+2));
                if (m_bp[qp1_closing_pos-1]) {// if there exists no closing parenthesis left to qp1_closing_pos
                    new_q = nsv(v.right+1);   // i.e. we could not find a next l-index, only a smaller value: last sibling
                    assert(new_q > v.right+1);
                } else {
                    new_q = m_bp_support.rank(m_bp_support.find_open(qp1_closing_pos-1))-1;
                    assert(new_q > v.right+1);
                    if (m_lcp[new_q] > lcp_qp1) { // no next l-index exists: find next smaller value, last sibling
                        new_q = nsv(v.right+1);
                    } else {
                        assert(m_lcp[new_q] == lcp_qp1);
                    }
                }
                assert(new_q > v.right+1);
                if (new_q-1 == v.right+1) { // sibling is a singlton interval
                    return ith_leaf(new_q/*v.right+2*/);
                } else {
                    // get l-index m of the sibling
                    size_type m_pos, m_close_pos;
                    size_type m = get_first_l_index(v.right+1, new_q-1, m_pos, m_close_pos);
                    return node_type(m_lcp[m], v.right+1, new_q-1);
                }

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

        //! Get the child w of v which edge label (v,w) starts with character c.
        /*!
         * 	\param v A valid tree node of the cst.
         *  \param c First character of the edge label from v to the desired child.
         *  \param char_pos Reference which will hold the position (0-based) of the matching char c in the sorted text/suffix array.
         *  \return The child node w which edge label (v,w) starts with c or root() if it does not exist.
         *  \par Time complexity
         *       \f$ \Order{ \sigma \cdot (\saaccess+\isaaccess+\lcpaccess) } \f$
         */
        node_type child_linear(const node_type& v, unsigned char c, size_type& char_pos)const {
            if (is_leaf(v))  // if v is a leave, it has no child
                return root();
            size_type k, kpos, k_find_close;
            // v is not a leave: v has at least two children
            k = get_first_l_index(v.left, v.right, kpos, k_find_close);// get first l-index k
            assert(m_lcp[k]==v.l);
            size_type left = v.left, right = k-1;
            char_pos = m_csa(m_csa[left]+v.l); // char_pos = invSA[SA[left]+v.l+1-1]
            uint16_t cc = m_csa.char2comp[c];
            size_type char_ex_max_pos = m_csa.C[cc+1], char_inc_min_pos = m_csa.C[cc];
            while (char_pos < char_ex_max_pos) {
                if (char_pos >= char_inc_min_pos) {
                    goto found_child;
                }
                // get next l-index
                size_type mpos = m_bp_support.find_open((m_bp_support.find_close(kpos)-1));
                size_type m    = m_bp_support.rank(mpos)-1;
                if (m==k or v.l != m_lcp[m]) { // no next l-index exists
                    left 	= k;
                    right 	= v.right;
                    char_pos = m_csa(m_csa[left]+v.l);
                    if (char_pos < char_inc_min_pos or char_pos >= char_ex_max_pos)  // check if character is two small
                        return root();
                    goto found_child;
                }
                left  = k;
                right = m-1;
                k = m;
                kpos = mpos;
                assert(m_csa(m_csa[left]+v.l) > char_pos);
                char_pos = m_csa(m_csa[left]+v.l);
            }
            return root();
found_child:
            if (left==right) { // child interval is a leaf
                assert(m_csa.size()-m_csa[left] > v.l);
                return node_type(m_csa.size()-m_csa[left],left,right);
            } else {
                // find l of the first child interval
                size_type k1pos, k1 = get_first_l_index(left, right, k1pos, k_find_close);
                assert(m_lcp[k1] > v.l);
                return node_type(m_lcp[k1], left, right);
            }
        }



        //! Get the child w of node v which edge label (v,w) starts with character c.
        /*!
         * 	\param v A valid tree node of the cst.
         *  \param c First character of the edge label from v to the desired child.
         *  \param char_pos Reference which will hold the position (0-based) of the matching char c in the sorted text/suffix array.
         *  \return The child node w which edge label (v,w) starts with c or root() if it does not exist.
         *  \par Time complexity
         *       \f$ \Order{(\lcpaccess+\saaccess+\isaaccess) \cdot \log\sigma} \f$
         */
        // TODO const unsigned char c durch char_type ersetzen
        node_type child(const node_type& v, const unsigned char c, size_type& char_pos)const {
            if (is_leaf(v))  // if v is a leaf = (), v has no child
                return root();
            // else v = ( (     ))
            uint16_t cc = m_csa.char2comp[c];
            if (cc==0 and c!=0) // TODO: aendere char2comp so ab, dass man diesen sonderfall nicht braucht
                return root();
            size_type char_ex_max_pos = m_csa.C[cc+1], char_inc_min_pos = m_csa.C[cc];
            // v is not a leave: v has at least two children
            size_type k, kpos, k_find_close;
            k = get_first_l_index(v.left, v.right, kpos, k_find_close);
            assert(m_lcp[k]==v.l);
            size_type left = v.left, right = k-1;
            char_pos = m_csa(m_csa[left]+v.l);
            if (char_pos >= char_ex_max_pos) {// c is lex. greater than the first character of the first child interval
                return root();
            } else if (char_pos >= char_inc_min_pos) { // found interval [left, right]
                goto found_child;
            }
            {
                size_type k_char_pos = m_csa(m_csa[k]+v.l);
                {
                    size_type begin = m_bp_support.preceding_closing_parentheses(k_find_close);
                    assert(m_csa.sigma>0);
                    if (begin > (size_type)m_csa.sigma-1) begin = m_csa.sigma-1;
                    begin = k_find_close - begin;
                    size_type end = k_find_close, m_find_close, m;
                    while (end != begin) {
                        m_find_close = (begin+end)>>1;  // begin <= m_find_close < end
                        m = m_bp_support.rank(m_bp_support.find_open(m_find_close))-1;
                        if (m_lcp[m] > v.l or /*m_lcp[m]==v.l and */ (char_pos = m_csa(m_csa[m]+v.l)) >= char_ex_max_pos) {
                            begin = m_find_close+1;
                        } else { // m_lcp[m] == v.l and char_pos < char_ex_max_pos
                            end = m_find_close;
                            k = m;
                            k_char_pos  = char_pos;
                        }
                    }
                    assert(m_lcp[k] == v.l);
                    assert(m_bp_support.find_close(m_bp_support.select(k+1)) == end);
                    if (k_char_pos < char_ex_max_pos and k_char_pos >= char_inc_min_pos) {// child interval exists
                        left = k;
                        char_pos = k_char_pos;
                        assert(end > 0);
                        if (m_bp[end-1] == 1) { // last child interval
                            right = v.right;
                        } else {
                            size_type k1 = m_bp_support.rank(m_bp_support.find_open(end-1))-1;
                            if (m_lcp[k1] > v.l) {// last child interval
                                right = v.right;
                            } else {
                                right = k1-1;
                            }
                        }
                    } else {
                        return root();
                    }
                }
            }
found_child:
            if (left==right) { // child interval is a leaf
                assert(m_csa.size()-m_csa[left] > v.l);
                return node_type(m_csa.size()-m_csa[left],left,right);
            } else {
                // find l of the first child interval
                size_type k1pos, k1 = get_first_l_index(left, right, k1pos, k_find_close);
                assert(m_lcp[k1] > v.l);
                return node_type(m_lcp[k1], left, right);
            }
        }

        //! Get the child w of node v which edge label (v,w) starts with character c.
        // \sa child(node_type v, const unsigned char c, size_type &char_pos)
        node_type child(const node_type& v, const unsigned char c) {
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
        unsigned char edge(const node_type& v, size_type d)const {
#ifndef NDEBUG
            if (d < 1 or d > depth(v)) {
                throw std::out_of_range("OUT_OF_RANGE_ERROR: "+util::demangle(typeid(this).name())+"cst_sct<>::edge(node_type v, size_type d). d == 0 or d > depth(v)!");
            }
#endif
            size_type order = m_csa(m_csa[v.left]+d-1);
            uint16_t c_begin = 1, c_end = m_csa.sigma+1, mid;
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
            if (v.left > w.left or (v.left == w.left and v.right < w.right)) {
                std::swap(v, w);
            }
            assert(v.left < w.left or (v.left==w.left and v.right >= w.right));
            assert(!(v.left < w.left and v.right > w.left and v.right < w.right));   // assert that v and w do not overlapp
            if (v.right >= w.right) { // v encloses w or v==w
                return v;
            } else { // v.left < v.right < w.left < w.right
                size_type min_index = rmq(v.left+1, w.right);
                size_type l = m_lcp[min_index];
                if (l==0)
                    return root();
                else
                    return node_type(l, psv(min_index, l), nsv(min_index)-1);
            }
        }

        //! Returns the string depth of node v.
        /*!
         * \param v A valid node of a cst_sct.
         * \returns The string depth of node v.
         * \par Time complexity
         *     \f$ \Order{1} \f$
         */
        size_type depth(const node_type& v)const {
            return v.l;
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
         * \param v A valid node of a cst_sct.
         * \return The suffix link of node v.
         * \par Time complexity
         *      \f$ \Order{ \rrenclose + \lcpaccess \cdot \log\sigma} \f$
         */
        node_type sl(const node_type& v)const {
            if (v == root())
                return root();
            // get interval with first char deleted
            size_type left	 = m_csa.psi[v.left];
            if (is_leaf(v)) {
                if (v.l==1)
                    return root();
                else
                    return node_type(v.l-1, left, left);
            }
            size_type right	 = m_csa.psi[v.right];
            assert(left < right);
            // rmq
            size_type min_index = rmq(left+1, right);
            size_type l = m_lcp[min_index];
            if (l==0)
                return root();
            else
                return node_type(l, psv(min_index, l), nsv(min_index)-1);
        }

        //! Compute the Weiner link of node v and character c.
        /*
         *  \param v A valid not of a cst_sct.
         *  \param c The character which should be prepended to the string of the current node.
         *	\return root() if the Weiner link of (v, c) does not exist, otherwise the Weiner link is returned.
         *  \par Time complexity
         *		\f$ \Order{ t_{rank\_bwt}\cdot t_{rmq} } \f$
         */
        node_type wl(const node_type& v, const unsigned char c) const {
            size_type c_left	= m_csa.rank_bwt(v.left, c);
            size_type c_right	= m_csa.rank_bwt(v.right+1, c);
            if (c_left == c_right)  // there exists no Weiner link
                return root();
            if (c_left+1 == c_right)
                return ith_leaf(m_csa.C[m_csa.char2comp[c]] + c_left + 1);
            else {
                size_type left	= m_csa.C[m_csa.char2comp[c]] + c_left;
                size_type right	= m_csa.C[m_csa.char2comp[c]] + c_right - 1;
                assert(left < right);
                // rmq
                size_type min_index = rmq(left+1, right);
                size_type l = m_lcp[min_index];
                return node_type(l, left, right);
            }
        }

        //! Compute the suffix number of a leaf node v.
        /*! \param v A valid leaf node of a cst_sct.
         *  \return The suffix array value corresponding to the leaf node v.
         *  \par Time complexity
         *		\f$ \Order{ \saaccess } \f$
         */
        size_type sn(const node_type& v)const {
            assert(is_leaf(v));
            return m_csa[v.left];
        }

        //! Computes a unique identification number for a node of the suffx tree in the range [0..2*size()-1]
        /*!
         *	\param v A valid node of a cst_sct.
         *  \return A unique identification number for the node v in the range [0..2*size()-1]
         *  \par Time complexity
         *		\f$ \Order{\lcpaccess} \f$
         *  \par Changes in future
         *		This method will someday hopefully return a value in the range [0..nodes()-1]
         */
        size_type id(const node_type& v)const {
            if (is_leaf(v)) { // return i in the range from 0..csa.size()-1
                return v.left;
            }
            size_type kpos, find_close_k;
            return (m_bp.size()>>1) + get_first_l_index(v.left, v.right, kpos, find_close_k);
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
            return node_type(l, lb, rb);
        }

        //! Print some infos about the size of the compressed suffix tree
        void print_info()const {
            std::cout << "# size of cst in bytes per character" << std::endl;
            size_type cst_size = util::get_size_in_bytes(*this);
            std::cout << ((double)cst_size)/csa.size() << " # = "<< ((double)cst_size)/(1<<20) <<" MB"<<std::endl;
            std::cout<< ((double)util::get_size_in_bytes(csa))/cst_size << " # ratio of csa size" << std::endl;
            std::cout<< ((double)util::get_size_in_bytes(lcp))/cst_size << " # ratio of lcp size" << std::endl;
            std::cout<< ((double)util::get_size_in_bytes(bp))/cst_size << " # ratio of bp size" << std::endl;
            std::cout<< ((double)util::get_size_in_bytes(bp_support))/cst_size << " # ratio of bp_support size" << std::endl;
            std::cout<< 0 << " # ratio of bp_rank_10 size" << std::endl;
            std::cout<< 0 << " # ratio of bp_select_10 size" << std::endl;
        }

        /* @} */

};

// == template functions ==

template<class Csa, class Lcp, class Bp_support>
template<uint8_t int_width, class size_type_class, uint8_t int_width_1, class size_type_class_1, uint8_t int_width_2, class size_type_class_2>
cst_sct<Csa, Lcp, Bp_support>::cst_sct(const std::string& csa_file_name,
                                       int_vector_file_buffer<int_width, size_type_class>& lcp_buf,
                                       int_vector_file_buffer<int_width_1, size_type_class_1>& sa_buf,
                                       int_vector_file_buffer<int_width_2, size_type_class_2>& isa_buf,
                                       std::string dir,
                                       bool build_only_bps
                                      ):csa(m_csa), lcp(m_lcp), bp(m_bp), bp_support(m_bp_support)
{
    std::string id =  util::to_string(util::get_pid())+"_"+util::to_string(util::get_id()).c_str();
    write_R_output("cst", "construct BPS", "begin", 1, 0);
    {
        bit_vector m_first_child;
        m_nodes = algorithm::construct_supercartesian_tree_bp_succinct_and_first_child(lcp_buf, m_bp, m_first_child) + m_bp.size()/2;
    }
    write_R_output("cst", "construct BPS", "end", 1, 0);

    write_R_output("cst", "construct BPSS", "begin", 1, 0);
    m_bp_support = Bp_support(&m_bp);
    write_R_output("cst", "construct BPSS", "end", 1, 0);

    write_R_output("cst", "construct CLCP", "begin", 1, 0);
    construct_lcp(m_lcp, *this, lcp_buf, isa_buf);
    write_R_output("cst", "construct CLCP", "end", 1, 0);

    util::load_from_file(m_csa, csa_file_name.c_str());
}

template<class Csa, class Lcp, class Bp_support>
cst_sct<Csa, Lcp, Bp_support>::cst_sct(tMSS& file_map, const std::string& dir, const std::string& id, bool build_only_bps):csa(m_csa), lcp(m_lcp), bp(m_bp), bp_support(m_bp_support)
{
    construct(file_map, dir, id, build_only_bps);
}

template<class Csa, class Lcp, class Bp_support>
void cst_sct<Csa, Lcp, Bp_support>::construct(tMSS& file_map, const std::string& dir, const std::string& id, bool build_only_bps)
{
    write_R_output("cst", "construct BPS", "begin", 1, 0);
    {
        bit_vector m_first_child;
        int_vector_file_buffer<> lcp_buf(file_map["lcp"].c_str());
        m_nodes = algorithm::construct_supercartesian_tree_bp_succinct_and_first_child(lcp_buf, m_bp, m_first_child) + m_bp.size()/2;
    }
    write_R_output("cst", "construct BPS", "end", 1, 0);

    write_R_output("cst", "construct CLCP", "begin", 1, 0);
    construct_lcp(m_lcp, *this, file_map, dir, id);
    write_R_output("cst", "construct CLCP", "end", 1, 0);

    write_R_output("cst", "construct BPSS", "begin", 1, 0);
    m_bp_support = Bp_support(&m_bp);
    write_R_output("cst", "construct BPSS", "end", 1, 0);

    util::load_from_file(m_csa, file_map["csa"].c_str());
    std::cout<<"debug = "<<m_csa.size()<<std::endl;
    std::cout<<"csa file = "<<file_map["csa"]<<std::endl;
}



template<class Csa, class Lcp, class Bp_support>
typename cst_sct<Csa, Lcp, Bp_support>::size_type cst_sct<Csa, Lcp, Bp_support>::serialize(std::ostream& out) const
{
    size_type written_bytes = 0;
    written_bytes += m_csa.serialize(out);
    written_bytes += m_lcp.serialize(out);
    written_bytes += m_bp.serialize(out);
    written_bytes += m_bp_support.serialize(out);
    out.write((char*) &m_nodes, sizeof(m_nodes));
    written_bytes += sizeof(m_nodes);
    return written_bytes;
}

template<class Csa, class Lcp, class Bp_support>
void cst_sct<Csa, Lcp, Bp_support>::load(std::istream& in)
{
    m_csa.load(in);
//	lcp_trait<Csa, Lcp>::load_lcp(m_lcp, in, m_csa);
    load_lcp(m_lcp, in, *this);
    m_bp.load(in);
    m_bp_support.load(in, &m_bp);
    in.read((char*) &m_nodes, sizeof(m_nodes));
#ifdef SDSL_DEBUG
    assert(algorithm::check_bp_support(m_bp, m_bp_support));
    std::cerr<<"checked bp_support"<<std::endl;
#endif
}

template<class Csa, class Lcp, class Bp_support>
cst_sct<Csa, Lcp, Bp_support>& cst_sct<Csa, Lcp, Bp_support>::operator=(const cst_sct& cst)
{
    if (this != &cst) {
        copy(cst);
    }
    return *this;
}





} // end namespace sdsl


#endif
