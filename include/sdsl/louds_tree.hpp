/* sdsl - succinct data structures library
    Copyright (C) 2012 Simon Gog
*/
/*! \file louds_tree.hpp
    \brief louds_tree.hpp contains a classes for the succinct tree representation LOUDS (level order unary degree sequence).
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_LOUDS_TREE
#define INCLUDED_SDSL_LOUDS_TREE

#include "int_vector.hpp"
#include <ostream>

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A class for the node representation of louds_tree
class louds_node
{
    public:
        typedef bit_vector::size_type size_type;
    private:
        size_type m_nr; 	// node number
        size_type m_pos;  // position in the LOUDS
    public:
        const size_type& nr;
        const size_type& pos;

        louds_node(size_type f_nr=0, size_type f_pos=0):m_nr(f_nr), m_pos(f_pos),nr(m_nr),pos(m_pos) {}

        bool operator==(const louds_node& v)const {
            return m_nr == v.m_nr and m_pos ==v.m_pos;
        }

        bool operator!=(const louds_node& v)const {
            return !(v==*this);
        }
};

std::ostream& operator<<(std::ostream& os, const louds_node& v);

//! A tree class based on the level order unary degree sequence (LOUDS) representation.
/*!
 * \tparam BitVector The bit vector representation used for LOUDS.
 * \tparam SelectSupport1 A select_support on 1-bits; it is needed for the child(v,i) operation.
 * \tparam SelectSupport0 A select_support on 0-bits; it is needed for the parent operation.
 *
 * Disadvantages of louds: No efficient support for subtree size.
*/
template<class BitVector = bit_vector, class SelectSupport1 = typename BitVector::select_1_type, class SelectSupport0 = typename BitVector::select_0_type>
class louds_tree
{
    public:
        typedef bit_vector::size_type 	size_type;
        typedef	louds_node				node_type;
        typedef	BitVector				bit_vector_type;
        typedef SelectSupport1			select_1_type;
        typedef SelectSupport0			select_0_type;

    private:
        bit_vector_type m_bv;  			// bit vector for the LOUDS sequence
        select_1_type   m_bv_select1;	// select support for 1-bits on m_bv
        select_0_type   m_bv_select0;	// select support for 0-bits on m_bv
    public:
        const bit_vector_type& bv;      // const reference to the LOUDS sequence

        //! Constructor for a cst and a root node for the traversal
        template<class Cst, class CstBfsIterator>
        louds_tree(const Cst& cst, const CstBfsIterator begin, const CstBfsIterator end):bv(m_bv) {
            bit_vector tmp_bv(4*cst.leaves_in_the_subtree(*begin) , 0); // resize the bit_vector to the maximal
            // possible size 2*2*#leaves in the tree
            size_type pos = 0;
            for (CstBfsIterator it = begin; it != end;) {
                tmp_bv[pos++] = 1;
                size_type size = it.size();
                ++it;
                pos += it.size()+1-size;
            }
            tmp_bv.resize(pos);
            // TODO do the BFS traversal
            util::assign(m_bv, tmp_bv);
            m_bv_select1.init(&m_bv);
            m_bv_select0.init(&m_bv);
        }

        //! Returns the root node
        node_type root() const {
            return louds_node(0, 0);
        }

        //! Returns the number of nodes in the tree.
        size_type nodes()const {
            return m_bv.size()/2;
        }

        //! Indicates if a node is a leaf.
        /*! \param v A node.
         */
        bool is_leaf(const node_type& v) const {
            // node is the last leaf        or  has no children, so m_bv[v.pos]==1
            return (v.pos+1 == m_bv.size()) or m_bv[v.pos+1];
        }

        //! Returns the number of children of a node.
        /*!
         *  \param v A node.
         */
        size_type degree(const node_type& v) const {
            if (is_leaf(v)) {  // handles boundary cases
                return 0;
            }
            // position of the next node  - node position - 1
            return m_bv_select1(v.nr+2) - v.pos - 1;
        }

        //! Returns the i-child of a node.
        /*!
         * \param v	The parent node.
         * \param i Index of the child. Indexing starts at 1.
         * \pre \f$ i \in [1..degree(v)] \f$
         */
        node_type child(const node_type& v, size_type i)const {
            size_type pos   = v.pos+i; // go to the position of the child's zero
            // (#bits = pos+1) - (#1-bits = v.nr+1)
            size_type zeros = pos+1 - (v.nr+1);
            return louds_node(zeros, m_bv_select1(zeros+1));
        }

        //! Returns the parent of a node v or root() if v==root().
        node_type parent(const node_type& v)const {
            if (v == root()) {
                return root();
            }
            size_type zero_pos   = m_bv_select0(v.nr);
            size_type parent_nr  = (zero_pos+1) - v.nr - 1;
            return node_type(parent_nr, m_bv_select1(parent_nr+1));
        }

        //! Returns an unique id for each node in [0..size()-1]
        size_type id(const node_type& v)const {
            return v.nr;
        }
};

}// end namespace sdsl

#endif // end file 
