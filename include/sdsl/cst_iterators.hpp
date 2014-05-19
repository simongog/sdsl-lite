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
/*! \file cst_iterators.hpp
    \brief cst_iterators.hpp contains iterator classes for traversing (compressed) suffix arrays.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_CST_ITERATORS
#define INCLUDED_SDSL_CST_ITERATORS

namespace sdsl
{

//! An forward iterator for (compressed) suffix trees.
/*! The cst_dfs_const_forward_iterator iterates through the nodes of a (compressed) suffix tree in
    depth first search (dfs) order. Note, that
	  - each inner node is visited twice, and
	  - each leaf node is visited only once.

	If the current node is a inner node, the method visit() returns 1 for the first visit and 2 for
    the second one.

	\par Time complexity
		\f$\Order{1}\f$ for all methods
	\par Space complexity
		\f$\Order{1} \f$

	This iterator is the standard iterator for the classes
		- sdsl::cst_sada,
		- sdsl::cst_sct3
 */
// TODO: implement operator--
template<class Cst, uint32_t cache_size=128>
class cst_dfs_const_forward_iterator: public std::iterator<std::forward_iterator_tag, typename Cst::node_type>
{
    public:
        typedef typename Cst::node_type value_type;
        typedef const value_type const_reference;
        typedef typename Cst::size_type size_type;
        typedef cst_dfs_const_forward_iterator<Cst> iterator;
        typedef typename Cst::node_type node_type;
    private:
        const Cst* m_cst;
        node_type m_v;
        bool	 m_visited;
        bool 	 m_valid;
        node_type* m_stack_cache;
        uint32_t m_stack_size;


        inline node_type parent() {
            --m_stack_size; // decrease stack size
            if (m_stack_cache != nullptr and m_stack_size < cache_size) {
                return m_stack_cache[m_stack_size];
            } else
                return m_cst->parent(m_v);
        }

        inline node_type first_child() {
            if (m_stack_cache != nullptr and m_stack_size < cache_size)   // push node to the stack
                m_stack_cache[ m_stack_size ] = m_v;
            m_stack_size++;
            return m_cst->select_child(m_v, 1);
        }

        //! Default constructor
        cst_dfs_const_forward_iterator():m_cst(nullptr),m_visited(false),m_valid(false), m_stack_cache(nullptr)
        {}
    public:

        //! Constructor
        cst_dfs_const_forward_iterator(const Cst* cst, const value_type node, bool visited=false, bool valid=true):m_visited(visited), m_valid(valid), m_stack_cache(nullptr) {
            m_cst = cst;
            m_v = node;
            if (m_cst == nullptr) {
                m_valid = false;
            } else if (m_v == m_cst->root() and !m_visited and m_valid) { // if the iterator equal cst.begin()
                m_stack_cache = new node_type[cache_size];
                m_stack_size  = 0;
//			std::cerr<<"#creating stack "<<m_cst->lb(m_v)<<" "<<m_cst->rb(m_v)<<std::endl;
            }
        }

        //! Copy Constructor
//	cst_dfs_const_forward_iterator(const cst_dfs_const_forward_iterator &it):m_cst(it.cst),m_v(it.m_v),m_valid(m_valid), m_stack_cache(nullptr),m_stack_size(0){
//	}

        ~cst_dfs_const_forward_iterator() {
            if (m_stack_cache != nullptr) {
//			std::cerr<<"#deleting stack "<<m_cst->lb(m_v)<<" "<<m_cst->rb(m_v)<<std::endl;
                delete [] m_stack_cache;
            }
        }

        //! Returns how often the current node was already visited.
        uint8_t visit()const {
            return 1+(uint8_t)m_visited;
        }

        void skip_subtree() {
            if (m_valid) {
                if (!m_visited) {
                    m_visited = true;
                }
            }
        }

        //! Method for dereferencing the iterator.
        const_reference operator*()const {
            return m_v;
        }

        //! Prefix increment of the iterator.
        iterator& operator++() {
            if (!m_valid)
                return *this;
            if (m_v == m_cst->root() and m_visited) {
                m_valid = false;
                return *this;
            }
            value_type w;
            if (!m_visited) { // go down, if possible
                if (m_cst->is_leaf(m_v)) {
                    w = m_cst->sibling(m_v);  // determine sibling of leaf v
                    if (w == m_cst->root()) { // if there exists no right sibling of the leaf v
//					w = m_cst->parent(m_v);
                        w = parent();
                        m_visited = true; // go up
                    }
                } else { // v is not a leaf => go down the tree
                    w = first_child();
                }
            } else { //
                w = m_cst->sibling(m_v);
                if (w == m_cst->root()) {   // if there exists no right sibling
                    w = parent();
                } else {
                    m_visited = false;
                }
            }
            m_v = w;
            return *this;
        }

        //! Postfix increment of the iterator.
        void operator++(int) {
            ++(*this);
        }

        //! Equality operator.
        bool operator==(const iterator& it)const {
            return (it.m_visited == m_visited) // visited status is equal
                   and (it.m_valid == m_valid) // valid status is equal => for end() iterator
                   and (it.m_v == m_v)    // nodes are equal
                   and (it.m_cst == m_cst);  // iterator belongs to the same cst
        }

        //! Inequality operator.
        bool operator!=(const iterator& it)const {
            return !(*this==it);
        }

};

//! A forward iterator for a bottom up traversal of a suffix tree
template<class Cst>
class cst_bottom_up_const_forward_iterator: public std::iterator<std::forward_iterator_tag, typename Cst::node_type>
{
    public:
        typedef typename Cst::node_type value_type;
        typedef const value_type const_reference;
        typedef typename Cst::size_type size_type;
        typedef cst_bottom_up_const_forward_iterator<Cst> iterator;
    private:
        const Cst* m_cst;
        typename Cst::node_type m_v;
        bool 	 m_valid;
    public:

        //! Default constructor
        cst_bottom_up_const_forward_iterator():m_cst(nullptr),m_valid(false) {}

        //! Constructor
        cst_bottom_up_const_forward_iterator(const Cst* cst, const value_type node, bool valid=true):m_valid(valid) {
            m_cst = cst;
            m_v = node;
            if (m_cst == nullptr)
                m_valid = false;
        }

        //! Method for dereferencing the iterator.
        const_reference operator*()const {
            return m_v;
        }

        //! Prefix increment of the iterator.
        iterator& operator++() {
            if (!m_valid)
                return *this;
            if (m_v == m_cst->root()) {
                m_valid = false;
                return *this;
            }
            value_type w = m_cst->sibling(m_v);
            if (w == m_cst->root()) {   // if no next right sibling exist
                m_v = m_cst->parent(m_v);    // go to parent
            } else { // if next right sibling exist
                m_v = m_cst->leftmost_leaf(w);   // go to leaftmost leaf in the subtree of w
            }
            return *this;
        }

        //! Postfix increment of the iterator.
        iterator operator++(int) {
            iterator it = *this;
            ++(*this);
            return it;
        }

        //! Equality operator.
        bool operator==(const iterator& it)const {
            return (it.m_valid == m_valid) // valid status is equal => for end() iterator
                   and (it.m_v == m_v)    // nodes are equal
                   and (it.m_cst == m_cst);  // iterator belongs to the same cst
        }

        //! Inequality operator.
        bool operator!=(const iterator& it)const {
            return !(*this==it);
        }

};

//! A forward iterator for a breath first traversal of a tree
/*!
 *	\tparam Cst 	A class which fulfills the CST concept
 *  \tparam Queue	A queue for the traversal. Note that for large data,
 *                  you should use an external implementation of a queue.
 */
template<class Cst, class Queue = std::queue<typename Cst::node_type> >
class cst_bfs_iterator: public std::iterator<std::forward_iterator_tag, typename Cst::node_type>
{
    public:
        typedef typename Cst::node_type 		value_type;
        typedef const value_type 				const_reference;
        typedef typename Cst::size_type 		size_type;
        typedef cst_bfs_iterator<Cst, Queue> 	iterator;
        typedef Queue 							queue_type;
    private:
        const Cst* 	m_cst;   // Pointer to the cst.
        queue_type	m_queue; //
        bool		m_valid; // State of the iterator.

    public:

        //! Constructor
        /*!
         * \param cst	Pointer to the compressed suffix tree.
         * \param node  Root node of the traversal.
         * \param valid State of the iterator.
         * \param end   If valid=true and end=true, we get the end() iterator otherwise ``end'' has no effect.
         */
        cst_bfs_iterator(const Cst* cst, const value_type node, bool valid=true, bool end_it=false) {
            m_cst = cst;
            m_valid = valid;
            if (m_cst != nullptr and !end_it) {
                m_queue.push(node);
            }
        }

        //! Returns the current number of nodes in the queue.
        size_type size()const {
            return m_queue.size();
        }

        //! Method for dereferencing the iterator.
        const_reference operator*()const {
            return m_queue.front();
        }

        //! Prefix increment of the iterator.
        iterator& operator++() {
            if (!m_valid)
                return *this;
            if (m_queue.empty()) {
                m_valid = false;
                return *this;
            }
            value_type v = m_queue.front();
            m_queue.pop();
            value_type child = m_cst->select_child(v, 1);
            while (m_cst->root() != child) {
                m_queue.push(child);
                child = m_cst->sibling(child);
            }
            return *this;
        }

        //! Postfix increment of the iterator.
        iterator operator++(int) {
            iterator it = *this;
            ++(*this);
            return it;
        }

        //! Equality operator.
        bool operator==(const iterator& it)const {
            if (m_queue.size() != it.m_queue.size()) {   // if the queue size is different
                return false;                            // the state of the to iterator are different
            }
            if (m_queue.empty()) {  // if the queue is empty, we have to check if they are valid and
                return it.m_valid == m_valid and it.m_cst == m_cst; // belong to the same cst
            }
            return (it.m_valid == m_valid) // valid status is equal => for end() iterator
                   and (it.m_cst == m_cst) // iterator belongs to the same cst
                   and (it.m_queue.front() == m_queue.front())  // front element and
                   and (it.m_queue.back() == m_queue.back());  //  back element are the same.
        }

        //! Inequality operator.
        bool operator!=(const iterator& it)const {
            return !(*this==it);
        }

};


} // end namespace sdsl

#endif
