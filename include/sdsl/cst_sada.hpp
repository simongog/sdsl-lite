/* sdsl - succinct data structures library
    Copyright (C) 2009-2013 Simon Gog

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
/*! \file cst_sada.hpp
    \brief cst_sada.hpp contains an implementation of Sadakane's CST.
    \author Simon Gog, Uwe Baier
*/
#ifndef INCLUDED_SDSL_CST_SADA
#define INCLUDED_SDSL_CST_SADA

#include "int_vector.hpp"
#include "iterators.hpp"
#include "lcp_support_sada.hpp"
#include "select_support_mcl.hpp"
#include "sorted_stack_support.hpp"
#include "bp_support.hpp"
#include "bp_support_sada.hpp"
#include "csa_sada.hpp" // for std initialization of cst_sada
#include "cst_iterators.hpp"
#include "cst_sct3.hpp"
#include "sdsl_concepts.hpp"
#include "construct.hpp"
#include "suffix_tree_helper.hpp"
#include "suffix_tree_algorithm.hpp"
#include "util.hpp"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring> // for strlen
#include <iomanip>
#include <iterator>

namespace sdsl
{

//! A class for the Compressed Suffix Tree (CST) proposed by Sadakane.
/*!
 * \tparam t_csa       Type of a CSA (member of this type is accessible via
 *                     member `csa`, default class is sdsl::csa_sada).
 * \tparam t_lcp       Type of a LCP structure (member is accessible via member
 *                     `lcp`, default class is sdsl::lcp_support_sada),
 * \tparam t_bp_support Type of a BPS structure (member accessible via member
 *                      `bp_support`, default class is sdsl::bp_support_sada),
 * \tparam t_rank_10    Type of a rank structure for the 2-bit pattern `10`
 *                      (accessible via member `bp_rank_10`, default class is
 *                      sdsl::rank_support_v5)
 * \tparam t_select_10  Type of a select structure for the 2-bit pattern `10`
 *                      (accessible via member \f$bp\_select\_10\f$, default
 *                      class is sdsl::select_support_mcl).
 *
 * It also contains a sdsl::bit_vector which represents the balanced
 * parentheses sequence of the suffix tree. This bit_vector can be accessed
 * via member `bp`.
 *
 * A node `v` of the `csa_sada` is represented by an integer `i` which
 * corresponds to the position of the opening parenthesis of the parentheses
 * pair \f$(i,\mu(i))\f$ that corresponds to `v` in `bp`.
 *
 * \par Reference
 *  Kunihiko Sadakane:
 *  Compressed Suffix Trees with Full Functionality.
 *  Theory Comput. Syst. 41(4): 589-607 (2007)
 *
 * @ingroup cst
 */
template<class t_csa = csa_sada<>,
         class t_lcp = lcp_support_sada<>,
         class t_bp_support = bp_support_sada<>,
         class t_rank_10 = rank_support_v5<10,2>,
         class t_select_10 = select_support_mcl<10,2>
         >
class cst_sada
{
        static_assert(std::is_same<typename index_tag<t_csa>::type, csa_tag>::value,
                      "First template argument has to be a compressed suffix array.");
    public:
        typedef cst_dfs_const_forward_iterator<cst_sada>          const_iterator;
        typedef cst_bottom_up_const_forward_iterator<cst_sada>    const_bottom_up_iterator;
        typedef typename t_csa::size_type                         size_type;
        typedef ptrdiff_t                                         difference_type;
        typedef t_csa                                             csa_type;
        typedef typename t_lcp::template type<cst_sada>           lcp_type;
        typedef typename t_csa::char_type                         char_type;
        typedef typename t_csa::string_type                       string_type;
        typedef size_type                                         node_type; //!< Type for the nodes  in the tree.
        typedef t_bp_support                                      bp_support_type;
        typedef t_rank_10                                         rank_10_type;
        typedef t_select_10                                       select_10_type;

        typedef typename t_csa::alphabet_type::comp_char_type     comp_char_type;
        typedef typename t_csa::alphabet_type::sigma_type         sigma_type;

        typedef typename t_csa::alphabet_category                 alphabet_category;
        typedef cst_tag                                           index_category;
    private:
        t_csa           m_csa; // suffix array
        lcp_type        m_lcp; // lcp information
        bit_vector      m_bp;  // balanced parentheses sequence for suffix tree
        bp_support_type m_bp_support; // support for the balanced parentheses sequence
        rank_10_type    m_bp_rank10;  // rank_support for leaves, i.e. "10" bit pattern
        select_10_type  m_bp_select10;// select_support for leaves, i.e. "10" bit pattern

        /* Get the number of leaves that are in the subtree rooted at the first child of v +
         * number of leafs in the subtrees rooted at the children of parent(v) which precede v in the tree.
         */
        size_type inorder(node_type v)const
        {
            return m_bp_rank10(m_bp_support.find_close(v+1)+1);
        }

        void copy(const cst_sada& cst)
        {
            m_csa           = cst.m_csa;
            copy_lcp(m_lcp, cst.m_lcp, *this);
            m_bp            = cst.m_bp;
            m_bp_support    = cst.m_bp_support;
            m_bp_support.set_vector(&m_bp);
            m_bp_rank10     = cst.m_bp_rank10;
            m_bp_rank10.set_vector(&m_bp);
            m_bp_select10   = cst.m_bp_select10;
            m_bp_select10.set_vector(&m_bp);
        }

    public:
        const t_csa&           csa          = m_csa;
        const lcp_type&        lcp          = m_lcp;
        const bit_vector&      bp           = m_bp;
        const bp_support_type& bp_support   = m_bp_support;
        const rank_10_type&    bp_rank_10   = m_bp_rank10;
        const select_10_type&  bp_select_10 = m_bp_select10;

//! Default constructor
        cst_sada() { }


//! Copy constructor
        cst_sada(const cst_sada& cst)
        {
            copy(cst);
        }

//! Move constructor
        cst_sada(cst_sada&& cst)
        {
            *this = std::move(cst);
        }

//! Construct CST from file_map
        cst_sada(cache_config& config)
        {
            {
                auto event = memory_monitor::event("bps-dfs");
                int_vector_buffer<> lcp(cache_file_name(conf::KEY_LCP, config));

                const bool o_par = true;
                const bool c_par = !o_par;

                //trim bps to maximal size of tree
                m_bp.resize( 4 * lcp.size() );

                if (lcp.size() > 0) {
                    //run from back to front of lcp, enumerate intervals and count
                    // opening parentheses per position i
                    sorted_stack_support stack( lcp.size()+1 );
                    stack.push( 0 ); //for lcp[n+1]
                    size_type p = m_bp.size() - 1;
                    for (size_type i = lcp.size() - 1; i > 0; --i) {
                        //compute number of opening parentheses at position i
                        size_type co = 1;  //for singleton interval
                        size_type x = lcp[i]+1; //to indicate start and end of lcp-array
                        while (stack.top() > x) {
                            stack.pop(); ++co;
                        }
                        if (stack.top() < x) {
                            stack.push(x);
                        }
                        //encode number of opening parenthesis at i as unary number
                        m_bp[p--] = o_par;
                        while (--co > 0)	m_bp[p--] = c_par;
                    }
                    //handle last value lcp[0] separate, since it virtually is a -1, but in real is a 0
                    m_bp[p--] = o_par; //code last number of opening parenthesis
                    while (stack.size() > 1) { //remove all elements except the zero from stack for next run
                        stack.pop();
                        m_bp[p--] = c_par;	//move k to first bit before unary number
                    }
                

                    //run from front to back of lcp, enumerate intervals,
                    //write opening parentheses and leave out closing parentheses
                    size_type q = 0;
                    for (size_type i = 1; i < lcp.size(); ++i) {
                        //compute number of opening parentheses at position i-1 using
                        //the unary coding from the last step
                        size_type co = 0;
                        do {
                            ++co;
                        } while (m_bp[++p] == c_par);

                        //compute number of closing parentheses at position i-1
                        size_type cc = 1; //for singleton interval
                        size_type x = lcp[i]+1;
                        while (stack.top() > x) {
                            stack.pop();	++cc;
                        }
                        if (stack.top() < x) {
                            stack.push(x);
                        }
                        //write sequence for position i-1
                        while (co-- > 0)	m_bp[q++] = o_par;
                        while (cc-- > 0)	m_bp[q++] = c_par;
                    }
                    //handle last value lcp[n+1] separate
                    m_bp[q++] = o_par;
                    while (!stack.empty()) {
                        m_bp[q++] = c_par;
                        stack.pop();
                    }

                    //trim bps to correct size and stop
                    m_bp.resize(q);
                }
            }
            {
                auto event = memory_monitor::event("bpss-dfs");
                util::assign(m_bp_support, bp_support_type(&m_bp));
                util::init_support(m_bp_rank10,   &m_bp);
                util::init_support(m_bp_select10, &m_bp);
            }
            {
                auto event = memory_monitor::event("clcp");
                cache_config tmp_config(false, config.dir, config.id, config.file_map);
                construct_lcp(m_lcp, *this, tmp_config);
                config.file_map = tmp_config.file_map;
            }
            {
                auto event = memory_monitor::event("load csa");
                load_from_cache(m_csa,std::string(conf::KEY_CSA)+"_"+util::class_to_hash(m_csa), config);
            }
        }

//! Number of leaves in the suffix tree.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         */
        size_type size()const
        {
            return m_csa.size();
        }

//! Returns the maximal length of text for that a suffix tree can be build.
        /*! Required for the Container Concept of the STL.
         *  \sa size
         */
        static size_type max_size()
        {
            return t_csa::max_size();
        }

//! Returns if the data structure is empty.
        /*! Required for the Container Concept of the STL.
         * \sa size
         */
        bool empty()const
        {
            return m_csa.empty();
        }

//! Swap method for cst_sada
        /*! The swap method can be defined in terms of assignment.
            This requires three assignments, each of which, for a container type, is linear
            in the container's size. In a sense, then, a.swap(b) is redundant.
            This implementation guaranties a run-time complexity that is constant rather than linear.
            \param cst cst_sada to swap.

            Required for the Assignable Concept of the STL.
          */
        void swap(cst_sada& cst)
        {
            if (this != &cst) {
                m_csa.swap(cst.m_csa);
                m_bp.swap(cst.m_bp);
                util::swap_support(m_bp_support, cst.m_bp_support, &m_bp, &(cst.m_bp));
                util::swap_support(m_bp_rank10, cst.m_bp_rank10, &m_bp, &(cst.m_bp));
                util::swap_support(m_bp_select10, cst.m_bp_select10, &m_bp, &(cst.m_bp));
                // anything else has to be swapped before swapping lcp
                swap_lcp(m_lcp, cst.m_lcp, *this, cst);
            }
        }

//! Returns a const_iterator to the first element.
        /*! Required for the STL Container Concept.
         *  \sa end
         */
        const_iterator begin()const
        {
            if (0 == m_bp.size())  // special case: tree is uninitialized
                return end();
            return const_iterator(this, root(), false, true);
        }

        //! Returns a const_iterator to the first element of a depth first traversal of the subtree rooted at node v.
        const_iterator begin(const node_type& v)const
        {
            if (0 == m_bp.size() and root()==v)
                return end();
            return const_iterator(this, v, false, true);
        }

//! Returns a const_iterator to the element after the last element.
        /*! Required for the STL Container Concept.
         *  \sa begin.
         */
        const_iterator end()const
        {
            return const_iterator(this, root(), true, false);
        }

        //! Returns a const_iterator to the element past the end of a depth first traversal of the subtree rooted at node v.
        const_iterator end(const node_type& v)const
        {
            if (root() == v)
                return end();
            return ++const_iterator(this, v, true, true);
        }

//! Returns an iterator to the first element of a bottom-up traversal of the tree.
        const_bottom_up_iterator begin_bottom_up()const
        {
            if (0 == m_bp.size())  // special case: tree is uninitialized
                return end_bottom_up();
            return const_bottom_up_iterator(this, leftmost_leaf(root()));
        }

//! Returns an iterator to the element after the last element of a bottom-up traversal of the tree.
        const_bottom_up_iterator end_bottom_up()const
        {
            return const_bottom_up_iterator(this, root(), false);
        }

//! Assignment Operator.
        /*!
         *    Required for the Assignable Concept of the STL.
         */
        cst_sada& operator=(const cst_sada& cst)
        {
            if (this != &cst) {
                copy(cst);
            }
            return *this;
        }

//! Assignment Move Operator.
        /*!
         *    Required for the Assignable Concept of the STL.
         */
        cst_sada& operator=(cst_sada&& cst)
        {
            if (this != &cst) {
                m_csa           = std::move(cst.m_csa);
                move_lcp(m_lcp, cst.m_lcp, *this);
                m_bp            = std::move(cst.m_bp);
                m_bp_support    = std::move(cst.m_bp_support);
                m_bp_support.set_vector(&m_bp);
                m_bp_rank10     = std::move(cst.m_bp_rank10);
                m_bp_rank10.set_vector(&m_bp);
                m_bp_select10   = std::move(cst.m_bp_select10);
                m_bp_select10.set_vector(&m_bp);
            }
            return *this;
        }

//! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_csa.serialize(out, child, "csa");
            written_bytes += m_lcp.serialize(out, child, "lcp");
            written_bytes += m_bp.serialize(out, child, "bp");
            written_bytes += m_bp_support.serialize(out, child, "bp_support");
            written_bytes += m_bp_rank10.serialize(out, child, "bp_rank_10");
            written_bytes += m_bp_select10.serialize(out, child, "bp_select_10");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

//! Load from a stream.
        /*! \param in Inputstream to load the data structure from.
         */
        void load(std::istream& in)
        {
            m_csa.load(in);
            load_lcp(m_lcp, in, *this);
            m_bp.load(in);
            m_bp_support.load(in, &m_bp);
            m_bp_rank10.load(in, &m_bp);
            m_bp_select10.load(in, &m_bp);
        }

        /*! \defgroup cst_sada_tree_methods Tree methods of cst_sada */
        /* @{ */

//! Return the root of the suffix tree.
        /*!
         * \par Time complexity
         *   \f$ \Order{1} \f$
         */
        node_type root() const
        {
            return 0;
        }

//! Decide if a node is a leaf in the suffix tree.
        /*!
        * \param v A valid node of a cst_sada.
        * \returns A boolean value indicating if v is a leaf.
        * \par Time complexity
        *      \f$ \Order{1} \f$
        */
        bool is_leaf(node_type v)const
        {
            assert(m_bp[v]==1);  // assert that v is a valid node of the suffix tree
            // if there is a closing parenthesis at position v+1, the node is a leaf
            return !m_bp[v+1];
        }

//! Return the i-th leaf (1-based from left to right) of the suffix tree.
        /*!
         * \param i 1-based position of the leaf. \f$1\leq i\leq csa.size()\f$.
         * \return The i-th leave.
         * \par Time complexity
         *     \f$ \Order{1} \f$
         * \pre \f$ 1 \leq i \leq csa.size() \f$
         */
        node_type select_leaf(size_type i)const
        {
            assert(i > 0 and i <= m_csa.size());
            // -1 as select(i) returns the position of the 0 of pattern 10
            return m_bp_select10.select(i)-1;
        }

//! Returns the depth of node v.
        /*!
         * \param v A valid node of the suffix tree.
         * \return The depth of the node.
         * \par Time complexity
         *    \f$ \Order{\lcpaccess \vee \saaccess} \f$
         */
        size_type depth(node_type v)const
        {
            if (v == root())  // if v is the root
                return 0;

            if (is_leaf(v)) { // if v is a leave
                size_type i = m_bp_rank10(v); // get the index in the suffix array
                return m_csa.size() - m_csa[i];
            }
            assert(inorder(v)>0);
            return m_lcp[inorder(v)];
        }

//! Returns the node depth of node v.
        /*!
         * \param v A valid node of a cst_sada.
         * \return The node depth of node v.
         * \par Time complexity
         *   \f$ \Order{1} \f$
         */
        size_type node_depth(node_type v)const
        {
            // -2 as the root() we assign depth=0 to the root
            return (m_bp_support.rank(v)<<1)-v-2;
        }

//! Calculate the number of leaves in the subtree rooted at node v.
        /*! \param v A valid node of the suffix tree.
         *  \return The number of leaves in the subtree rooted at node v.
         *  \par Time complexity
         *      \f$ \Order{1} \f$
         *
         *  This method is used e.g. in the count method.
         */
        size_type size(node_type v)const
        {
            size_type r = m_bp_support.find_close(v);
            return m_bp_rank10(r+1) - m_bp_rank10(v);
        }

//! Calculates the leftmost leaf in the subtree rooted at node v.
        /*! \param v A valid node of the suffix tree.
         *  \return The leftmost leaf in the subtree rooted at node v.
         *  \par Time complexity
         *    \f$ \Order{1} \f$
         */
        node_type leftmost_leaf(const node_type v)const
        {
            return m_bp_select10(m_bp_rank10(v)+1)-1;
        }

//! Calculates the rightmost leaf in the subtree rooted at node v.
        /*!\param v A valid node of the suffix tree.
         * \return The rightmost leaf in the subtree rooted at node v.
         * \par Time complexity
         *   \f$ \Order{1} \f$
         */
        node_type rightmost_leaf(const node_type v)const
        {
            size_type r = m_bp_support.find_close(v);
            return m_bp_select10(m_bp_rank10(r+1))-1;
        }

//!Calculates the index of the leftmost leaf in the corresponding suffix array.
        /*!\param v A valid node of the suffix tree.
         * \return The index of the leftmost leaf in the corresponding suffix array.
         * \par Time complexity
         *  \f$ \Order{1} \f$
         * \par Note
         * lb is an abbreviation for ,,left bound''
         */
        size_type lb(const node_type v)const
        {
            return m_bp_rank10(v);
        }

//! Calculates the index of the rightmost leaf in the corresponding suffix array.
        /*! \param v A valid node of the suffix tree.
         *  \return The index of the rightmost leaf in the corresponding suffix array.
         *  \par Time complexity
         *    \f$ \Order{1} \f$
         *  \par Note
         *   rb is an abbreviation for ,,right bound''
         */
        size_type rb(const node_type v)const
        {
            size_type r = m_bp_support.find_close(v);
            return m_bp_rank10(r+1)-1;
        }

//! Calculate the parent node of a node v.
        /*! \param v A valid node of the suffix tree.
         *  \return The parent node of v or root() if v equals root().
         *  \par Time complexity
         *       \f$ \Order{1} \f$
         */
        node_type parent(node_type v) const
        {
            assert(m_bp[v]==1); // assert a valid node
            if (v == root())
                return root();
            else {
                return m_bp_support.enclose(v);
            }
        }

        //! Return a proxy object which allows iterating over the children of a node
        /*! \param v A valid node of the suffix tree.
         *  \return The proxy object of v containing all children
         *  \par Time complexity
         *     \f$ \Order{1}\f$
         */
        cst_node_child_proxy<cst_sada> children(node_type v) const
        {
            return cst_node_child_proxy<cst_sada>(this,v);
        }


//! Returns the next sibling of node v.
        /*!
         * \param v A valid node v of the suffix tree.
         * \return The next (right) sibling of node v or root() if v has no next (right) sibling.
         * \par Time complexity
         *    \f$ \Order{1} \f$
         */
        node_type sibling(node_type v)const
        {
            if (v==root())
                return root();
            node_type sib = m_bp_support.find_close(v)+1;
            if (m_bp[sib])
                return sib;
            else
                return root();
        }

//! Get the child w of node v which edge label (v,w) starts with character c.
        /*
         * \param v A valid tree node of the cst.
         * \param c First character of the edge label from v to the desired child.
         * \param char_pos Reference which will hold the position (0-based) of the matching char c in the sorted text/suffix array.
         * \return The child node w which edge label (v,w) starts with c or root() if it does not exist.
         * \par Time complexity
         *   \f$ \Order( (\saaccess+\isaaccess) \cdot \sigma + \lcpaccess) \f$
         * \par Note
         *   With range median mimimum queries (RMMQ) one can code this operation in \f$\log \sigma \f$ time
         */
        node_type child(node_type v, const char_type c, size_type& char_pos)const
        {
            if (is_leaf(v))  // if v is a leaf = (), v has no child
                return root();
            // else v = ( (     ))
            comp_char_type cc = m_csa.char2comp[c];
            if (cc==0 and c!=0) // TODO: aendere char2comp so ab, dass man diesen sonderfall nicht braucht
                return root();
            size_type char_ex_max_pos = m_csa.C[cc+1], char_inc_min_pos = m_csa.C[cc];

            size_type d = depth(v);  // time complexity: \lcpaccess
            size_type res = v+1;
            while (true) {
                if (is_leaf(res)) {
                    char_pos = get_char_pos(m_bp_rank10(res), d, m_csa);
                } else {
                    char_pos = get_char_pos(inorder(res), d, m_csa);
                }
                if (char_pos >= char_ex_max_pos)  // if the current char is lex. greater than the searched char: exit
                    return root();
                if (char_pos >= char_inc_min_pos)  // if the current char is lex. equal with the
                    return res;
                res = m_bp_support.find_close(res)+1;
                if (!m_bp[res]) // closing parenthesis: there exists no next child
                    return root();
            }
        }

        //! Get the child w of node v which edge label (v,w) starts with character c.
        // \sa child(node_type v, const char_type c, size_type &char_pos)
        node_type child(node_type v, const char_type c) const
        {
            size_type char_pos;
            return child(v, c, char_pos);
        }

        //! Get the i-th child of a node v.
        /*!
         * \param v A valid tree node of the cst.
         * \param i 1-based Index of the child which should be returned. \f$i \geq 1\f$.
         * \return The i-th child node of v or root() if v has no i-th child.
         * \par Time complexity
         *   \f$ \Order{i} \f$ for \f$  i \leq \sigma \f$
         *  \pre \f$ 1 \leq i \leq degree(v) \f$
         */
        node_type select_child(node_type v, size_type i)const
        {
            if (is_leaf(v))  // if v is a leave, v has no child
                return root();
            size_type res = v+1;
            while (i > 1) {
                res = m_bp_support.find_close(res)+1;
                if (!m_bp[res]) {// closing parenthesis: there exists no next child
                    return root();
                }
                --i;
            }
            return res;
        }

//! Returns the d-th character (1-based indexing) of the edge-label pointing to v.
        /*!\param v The node at which the edge path ends.
         * \param d The position (1-based indexing) of the requested character on the edge path from the root to v. \f$ d > 0 \wedge d <= depth(v) \f$
         * \return The character at position d on the edge path from the root to v.
         * \par Time complexity
         *       \f$ \Order{ \log\sigma + (\saaccess+\isaaccess) } \f$
         * \pre \f$ 1 \leq d \leq depth(v)  \f$
         */
        char_type edge(node_type v, size_type d)const
        {
            assert(1 <= d);
            assert(d <= depth(v));
            size_type i = 0;// index of the first suffix in the subtree of v
            if (is_leaf(v)) { // if v is a leave
                i = m_bp_rank10(v); // get the index in the suffix array
            } else {
                i = inorder(v);
            }
            size_type      order    = get_char_pos(i, d-1, m_csa);
            size_type    c_begin    = 1, c_end = ((size_type)m_csa.sigma)+1, mid;
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
         *    \f$ \Order{\rrenclose}\   \f$
         */
        node_type lca(node_type v, node_type w)const
        {
            assert(m_bp[v] == 1 and m_bp[w] == 1);
            if (v > w) {
                std::swap(v,w);
            } else if (v==w) {
                return v;
            }
            if (v == root())
                return root();
            return m_bp_support.double_enclose(v, w);
        }

//! Compute the suffix link of node v.
        /*!
         * \param v A valid node of a cst_sada.
         * \return The suffix link of node v.
         * \par Time complexity
         *   \f$ \Order{ 1 } \f$
         */
        node_type sl(node_type v)const
        {
            if (v == root())
                return root();
            // get leftmost leaf in the tree rooted at v
            size_type left        = m_bp_rank10(v);
            if (is_leaf(v)) {
                return select_leaf(m_csa.psi[left]+1);
            }
            // get the rightmost leaf in the tree rooted at v
            size_type right         = m_bp_rank10(m_bp_support.find_close(v))-1;
            assert(left < right);
            node_type left_leaf = select_leaf(m_csa.psi[left]+1);
            node_type right_leaf= select_leaf(m_csa.psi[right]+1);
            return lca(left_leaf, right_leaf);
        }


        //! Compute the suffix link of node v applied a number of times consecutively.
        /*!
         * \param v A valid node of a cst_sada
         * \return The suffix link ^ i of node v.
         * \par Time complexity
         *   \f$ \Order{ i } \f$
         */
        node_type sl(node_type v, size_type i)const
        {
            if (v == root())
                return root();
            // get leftmost leaf in the tree rooted at v
            size_type left        = m_bp_rank10(v);
            if (is_leaf(v)) {
                return select_leaf(get_char_pos(left, i, m_csa)+1);
            }
            // get the rightmost leaf in the tree rooted at v
            size_type right         = m_bp_rank10(m_bp_support.find_close(v))-1;
            assert(left < right);
            node_type left_leaf = select_leaf(get_char_pos(left, i, m_csa)+1);
            node_type right_leaf= select_leaf(get_char_pos(right, i, m_csa)+1);
            return lca(left_leaf, right_leaf);
        }

//! Compute the Weiner link of node v and character c.
        /*
         * \param v A valid node of a cst_sada.
         * \param c The character which should be prepended to the string of the current node.
         *   \return root() if the Weiner link of (v, c) does not exist, otherwise the Weiner link is returned.
         * \par Time complexity
         *    \f$ \Order{ t_{rank\_bwt} + t_{lca}}\f$
         */
        node_type wl(node_type v, const char_type c) const
        {
            // get leftmost leaf in the tree rooted at v
            size_type left        = m_bp_rank10(v);
            // get the rightmost leaf in the tree rooted at v
            size_type right = is_leaf(v) ? left : m_bp_rank10(m_bp_support.find_close(v))-1;

            size_type c_left    = m_csa.bwt.rank(left, c);
            size_type c_right    = m_csa.bwt.rank(right+1, c);

            if (c_left == c_right)  // there exists no Weiner link
                return root();
            if (c_left+1 == c_right)
                return select_leaf(m_csa.C[m_csa.char2comp[c]] + c_left + 1);
            else {
                size_type left    = m_csa.C[m_csa.char2comp[c]] + c_left;
                size_type right    = m_csa.C[m_csa.char2comp[c]] + c_right - 1;
                assert(left < right);
                node_type left_leaf = select_leaf(left+1);
                node_type right_leaf= select_leaf(right+1);
                return lca(left_leaf, right_leaf);
            }
        }

//! Compute the suffix number of a leaf node v.
        /*!\param v A valid leaf node of a cst_sada.
         * \return The suffix array value corresponding to the leaf node v.
         * \par Time complexity
         *   \f$ \Order{ \saaccess } \f$
         */
        size_type sn(node_type v)const
        {
            assert(is_leaf(v));
            // count the leaves left to leaf v
            return m_csa[m_bp_rank10(v)];
        }

//! Computes a unique identification number for a node of the suffix tree in the range [0..nodes()-1]
        /*!
         *\param v A valid node of a cst_sada.
         * \return A unique identification number for the node v in the range [0..nodes()-1]
         * \par Time complexity
         *     \f$ \Order{1} \f$
         * \sa inv_id(size_type id)
         */
        size_type id(node_type v)const
        {
            // v+1 is < m_bp.size(), as v is the position of an open parenthesis
            if (m_bp[v+1]) {    // case (a) inner node
                return size() + (m_bp_support.rank(v) - 1) - m_bp_rank10(v);
            } else {            // case (b) leaf
                return m_bp_rank10(v);
            }
        }

//! Computes the node for such that id(v)=id.
        /*!
         * \param id An id in the range [0..nodes()-1].
         * \return A node v of the CST such that id(v)=id.
         * \par Time complexity
         *        \f$ \Order{1} \f$ for leaves and \f$ \log n \f$ for inner nodes
         * \sa id(node_type v)
         */
        size_type inv_id(size_type id)
        {
            if (id < size()) {  // the corresponding node is a leaf
                return select_leaf(id+1);
            } else { // the corresponding node is a inner node
                id = id + 1 - size();
                // solved by binary search; TODO: can be done in constant time by using a select structure on the bitpattern 11
                size_type lb = 0, rb = m_bp.size(); // lb inclusive, rb exclusive
                // invariant: arg(lb) < id, arg(rb)>= id
                while (rb-lb > 1) {
                    size_type mid = lb + (rb-lb)/2; // mid \in [0..m_bp.size()-1]
                    if (m_bp[mid] == 0 and m_bp[mid-1] == 1) {  // if we are ``half on a leaf''
                        ++mid; //we step one to the right to include it
                    }
                    // get the number of open inner nodes before position mid, i.e. arg(mid)
                    size_type mid_id = m_bp_support.rank(mid-1) - m_bp_rank10(mid);  // Note: mid-1 is valid of mid is of type ``size_type'' as us the parameter of rank
                    if (mid_id < id) {
                        lb = mid;
                    } else { // mid_id >= x
                        rb = mid;
                    }
                }
                return lb;
            }
        }

//! Get the number of nodes of the suffix tree.
        /*
         *  \return The number of nodes of the suffix tree.
         *  \par Time complexity
         *    \f$ \Order{1} \f$
         */
        size_type nodes()const
        {
            return m_bp.size()>>1;
        }

//! Get the node in the suffix tree which corresponds to the lcp-interval [lb..rb]
        /* \param lb Left bound of the lcp-interval [lb..rb] (inclusive).
         * \param rb Right bound of the lcp-interval [lb..rb] (inclusive).
         *\ return The node in the suffix tree corresponding lcp-interval [lb..rb]
         */
        node_type node(size_type lb, size_type rb) const
        {
            return lca(select_leaf(lb+1), select_leaf(rb+1));
        }

//! Get the number of children of a node v.
        /*!
         *  \param v A valid node v of a cst_sada.
         *  \returns The number of children of node v.
         *  \par Time complexity
         *       \f$ \Order{\sigma} \f$
         */
        size_type degree(node_type v)const
        {
            size_type res = 0;
            v = v+1;
            while (m_bp[v]) { // found open parentheses
                ++res;
                v = m_bp_support.find_close(v)+1;
            }
            return res;
        }

//! Maps an index i to the position in TLCP where LCP[i] can be found
        /*!
         * \param i The index in the LCP array
         * \return The corresponding position in the TLCP array
         */
        size_type tlcp_idx(size_type i) const
        {
            size_type ii = 0;
            if (i > 0) {
                size_type ipos   = m_bp_select10(i) - 1;  // -1 as select returns the position of the zero
                size_type ip1pos = m_bp_select10(i+1) - 1;// "  "    "        "    "     "     "   "   "
                ii    = m_bp_support.double_enclose(ipos, ip1pos);
            }
            ii = m_bp_support.find_close(ii);
            // all right, as bp[ii] = 0
            return ii - m_bp_support.rank(ii) - m_bp_rank10(ii);
        }
        /* @} */
};

} // end namespace sdsl
#endif
