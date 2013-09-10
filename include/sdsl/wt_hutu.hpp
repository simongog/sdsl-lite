/* sdsl - succinct data structures library
    Copyright (C) 2013 Simon Gog, Timo Beller and Markus Brenner

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
/*! \file wt_hutu.hpp
    \brief wt_hutu.hpp contains a class for a Hu-Tucker shaped wavelet tree
           over byte sequences.
    \author Simon Gog, Markus Brenner
*/
#ifndef INCLUDED_SDSL_WT_HUTU
#define INCLUDED_SDSL_WT_HUTU

#include "wt_pc.hpp"
#include <vector>

//! Namespace for the succinct data structure library.
namespace sdsl
{

// forward declaration
struct hutu_shape;

//! A Hu-Tucker-shaped wavelet tree.
/*!
 *  \tparam t_bitvector   Underlying bitvector structure.
 *  \tparam t_rank        Rank support for pattern `1` on the bitvector.
 *  \tparam t_select      Select support for pattern `1` on the bitvector.
 *  \tparam t_select_zero Select support for pattern `0` on the bitvector.
 *  \tparam t_dfs_shape   Layout of the tree structure in memory. Set 0
 *                        for BFS layout and 1 fro DFS layout.
 *  \par Space complexity
 *     Almost \f$n H_0 + 2|\Sigma|\log n\f$ bits, where \f$n\f$ is the size of
 *     the vector the wavelet tree was build for.
 *
 *   @ingroup wt
 */
template<class t_bitvector   = bit_vector,
         class t_rank        = typename t_bitvector::rank_1_type,
         class t_select      = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type,
         class t_tree_strat  = byte_tree<> >
using wt_hutu = wt_pc<hutu_shape,
      t_bitvector,
      t_rank,
      t_select,
      t_select_zero,
      t_tree_strat>;

// Hu Tucker shape for wt_pc
template<class t_wt>
struct _hutu_shape {
        typedef typename t_wt::size_type size_type;
        enum { lex_ordered = 1 };

        //! Node class used by the leftist heap
        template <class t_element>
        struct heap_node {
            t_element* item;   // pointer to the represented item
            heap_node* left, *right, *parent;  // pointer to left/right child, parent
            int64_t rank; // rank of the heap node
            //! Constructor
            heap_node(t_element* it=nullptr) : item(it), left(nullptr),
                right(nullptr), parent(nullptr),
                rank(0) { }
            //! Less then operator
            bool operator< (const heap_node& other) {
                return *item < *(other.item);
            }
        };

        // Implementation of a leftist heap as needed in the first phase of
        // Hu-Tucker Code construction
        template <class t_element>
        class l_heap
        {
            private:
                heap_node<t_element>* m_root;  // pointer to the root

                // fixes node information after the deletion of elements
                void fix_node(heap_node<t_element>* item) {
                    if (item != nullptr) {
                        if (!item->left || !item->right) { // if node has only one child
                            // only go on fixing if the node information needs to be changed
                            if (item->rank != 0) {
                                item->rank = 0;
                                if (item->parent) fix_node(item->parent);
                            }
                        } else { // node information has to be adapted
                            int64_t nn = (item->left->rank > item->right->rank) ? item->right->rank : item->left->rank;
                            if (item->rank != nn && item->parent != 0) {
                                item->rank = nn;
                                fix_node(item->parent);
                            }
                        }
                    }
                }

                // helper function to remove the data structure from memory
                void free_node(heap_node<t_element>* item) {
                    if (item->left) {
                        free_node(item->left);
                        delete item->left;
                        item->left = nullptr;
                    }
                    if (item->right) {
                        free_node(item->right);
                        delete item->right;
                        item->right = nullptr;
                    }
                }

                // internal merge function
                heap_node<t_element>* merge(heap_node<t_element>* h1, heap_node<t_element>* h2) {
                    if (!h1) return h2;
                    if (!h2) return h1;
                    if (*(h1->item) < *(h2->item))  return merge1(h1, h2);
                    else                            return merge1(h2, h1);
                }
                // internal merge function
                heap_node<t_element>* merge1(heap_node<t_element>* h1, heap_node<t_element>* h2) {
                    if (!h1->left) { // if h1 has no children, the merge is simple
                        h1->left = h2;
                        h2->parent = h1; // adjust the parent pointer
                    } else {
                        h1->right = merge(h1->right, h2);
                        if (h1->right) h1->right->parent = h1;

                        if (h1->left->rank < h1->right->rank) {
                            heap_node<t_element>* tmp = h1->left;
                            h1->left = h1->right;
                            h1->right = tmp;
                        }
                        h1 -> rank = h1 -> right -> rank + 1;
                    }
                    return h1;
                }

            public:

                //! Default constructor
                l_heap() : m_root(nullptr) { }

                //! Indicates if the heap is empty
                bool empty() const {
                    return (m_root==nullptr);
                }

                //! Get the smallest element
                /*! \return The smallest element in the heap
                 *          or nullptr if it does not exist.
                 */
                heap_node<t_element>* find_min() const {
                    return m_root;
                }

                //! Get the second smallest element
                /*! \return The second smallest element in the heap
                 *         or nullptr if it does not exist.
                 */
                heap_node<t_element>* find_snd_min() const {
                    if (m_root == nullptr) return nullptr;
                    if (m_root->left == nullptr) return m_root->right;
                    if (m_root->right == nullptr) return m_root->left;

                    if (m_root->left->operator< (*m_root->right)) return m_root->left;
                    else return m_root->right;
                }

                //! Insert an element into the heap
                /*! \param  x Element that is inserted into the heap.
                 *  \return The new generated heap node.
                 */
                heap_node<t_element>* insert(t_element* x) {
                    heap_node<t_element>* n = new heap_node<t_element>(x);
                    l_heap<t_element> lh;
                    lh.m_root = n;
                    merge(&lh);
                    return n;
                }

                //! Delete the smallest element in the heap
                void delete_min() {
                    heap_node<t_element>* old_root = m_root;
                    m_root = merge(m_root->left, m_root->right);
                    if (m_root) m_root->parent = nullptr;
                    delete old_root;
                }

                // deletes an arbitrary element from the heap
                // this function assumes, that item is an element of the heap
                void delete_element(heap_node<t_element>* item) {
                    if (item != nullptr) {
                        if (m_root == item) { // deleting the root is trivial
                            delete_min();
                        } else {
                            // otherwise we have to adapt the parent node and
                            // the children of item
                            heap_node<t_element>* h1 = merge(item->left,item->right);
                            if (h1) h1->parent = item->parent;
                            if (item == item->parent->left) {
                                item->parent->left = h1;
                            } else if (item == item->parent->right) {
                                item->parent->right = h1;
                            }
                            // fix node information considering rank
                            fix_node(item->parent);
                            delete item; // remove the item from memory
                        }
                    }
                }

                // public merge function
                void merge(l_heap<t_element>* rhs) {
                    m_root = merge(m_root, rhs->m_root);
                    rhs->m_root = nullptr;
                }

                // removes the whole data structure from memory
                void free_memory() {
                    if (m_root != nullptr) {
                        free_node(m_root);
                        delete m_root;
                        m_root = nullptr;
                    }
                }
        };


        // forward declaration of node classes
        struct ht_node;

        // Master node as used in the first phase of the Hu-Tucker algorithm
        struct m_node {
            // min sum of the two min elements of the hpq this node points to
            size_type min_sum;
            int64_t i; // position of the left node in the working sequence
            int64_t j; // position of the right node in the working sequence
            // pointer to the corresponding heap element (used for deletion)
            heap_node<m_node>* qel;
            l_heap<ht_node>* myhpq;  // pointer to the hpq

            ht_node* lt;  // pointer to the left- and rightmost leafs of the hpq
            ht_node* rt;  // need for merge operations

            m_node() : qel(0), myhpq(0), lt(0), rt(0) { }

            bool operator<(const m_node other) {
                if (min_sum != other.min_sum) {
                    return min_sum < other.min_sum;
                }
                if (i != other.i) {
                    return i < other.i;
                }
                return j < other.j;
            }

            bool operator> (const m_node other) {
                return other < *this;
            }
        };

        // Hu-Tucker node as used in the first phase of the Hu-Tucker algorithm
        struct ht_node {
            int64_t   pos;   // position of the node
            uint64_t  c;     // the represented letter
            size_type w;     // frequency of the node
            bool      t;     // whether the node is a leaf
            int64_t   level; // level in the tree

            // pointer to the two master nodes
            // (as a node can belong to up to two hpqs)
            m_node*  mpql;
            m_node*  mpqr; // only mpql is used for inner nodes
            // pointer to the two heap nodes (as a node can belong to up to two hpqs)
            heap_node<ht_node>* ql;
            heap_node<ht_node>* qr; // only ql is used for inner nodes
            ht_node* left;  // left child
            ht_node* right; // right child

            ht_node() : mpql(0), mpqr(0), ql(0), qr(0),
                left(nullptr), right(nullptr) { }

            bool operator< (const ht_node& other) {
                if (w != other.w) {
                    return w < other.w;
                }
                return pos < other.pos;
            }

            bool operator> (const ht_node& other) {
                return other < *this;
            }
        };


        template<class t_rac>
        static void
        construct_tree(t_rac& C, std::vector<pc_node>& temp_nodes) {
            //create a leaf for every letter
            std::vector<ht_node> node_vector;
            for (size_t i = 0; i < C.size(); i++) {
                if (C[i]) {
                    ht_node n;
                    n.c = (uint64_t)i;
                    n.w = C[i];
                    n.t = true;
                    n.pos = node_vector.size();
                    node_vector.push_back(n);
                }
            }
            if (node_vector.size() == 1) {
                // special case of an alphabet of size 1:
                // just instantly create the tree and return it
                temp_nodes.emplace_back(pc_node(node_vector[0].w,
                                                (size_type)node_vector[0].c));
                return;
            }
            size_type       sigma = node_vector.size();
            std::vector<ht_node>  T(sigma); // physical leaves
            std::vector<ht_node*> A(sigma); // the current working sequence
            // Priority Queues, containing the Huffman Sequences
            std::vector<l_heap<ht_node>> HPQ(sigma);
            l_heap<m_node>  MPQ;      // Master Priority Queue

            // init T, A, HPQs and MPQ
            T[0] = node_vector[0];
            A[0] = &T[0];

            // initialization needed for every leaf
            for (size_type i = 1; i < sigma; i++) {
                T[i] = node_vector[i];
                A[i] = &T[i];

                T[i - 1].qr = HPQ[i - 1].insert(&T[i - 1]);
                T[i].ql = HPQ[i - 1].insert(&T[i]);

                m_node* m = new m_node();
                m->min_sum = T[i - 1].w + T[i].w;
                m->i = i - 1;
                m->j = i;
                m->lt = &T[i-1];
                m->rt = &T[i];
                m->myhpq = &HPQ[i - 1];

                m->qel = MPQ.insert(m);

                T[i-1].mpqr = m;
                T[i].mpql = m;
            }

            // main action loop
            for (size_type k = 1; k < sigma; k++) {
                m_node* m = MPQ.find_min()->item;
                ht_node* l = A[m->i];
                ht_node* r = A[m->j];
                int64_t lpos = m->i;
                int64_t rpos = m->j;

                l_heap<ht_node>* n_hpq = nullptr;
                ht_node* n_rt = nullptr;
                ht_node* n_lt = nullptr;

                // create a new master priority queue
                m_node* n_m = new m_node();
                // delete old nodes from all hpqs
                if (l->t) {
                    if (l->mpql) l->mpql->myhpq->delete_element(l->ql);
                    l->ql = nullptr;
                    if (l->mpqr) l->mpqr->myhpq->delete_element(l->qr);
                    l->qr = nullptr;
                } else {
                    m->myhpq->delete_element(l->ql);
                    l->ql = nullptr;
                }
                if (r->t) {
                    if (r->mpql) r->mpql->myhpq->delete_element(r->ql);
                    l->ql = nullptr;

                    if (r->mpqr) r->mpqr->myhpq->delete_element(r->qr);
                    r->qr = nullptr;
                } else {
                    m->myhpq->delete_element(r->ql);
                    r->ql = nullptr;
                }
                // handle the merge of hpqs
                if (l->t && r ->t) {
                    // both nodes are leaves
                    l_heap<ht_node>* h1 = nullptr;
                    l_heap<ht_node>* h2 = nullptr;
                    l_heap<ht_node>* h3 = nullptr;
                    if (l -> mpql) {
                        n_lt = l->mpql->lt;
                        if (n_lt == l) n_lt = nullptr;
                        if (n_lt) n_lt -> mpqr = n_m;

                        h1 = l->mpql->myhpq;
                        h2 = l->mpqr->myhpq;

                        h1 -> merge(h2);
                        MPQ.delete_element(l->mpql->qel);
                        MPQ.delete_element(l->mpqr->qel);
                        delete l->mpql;
                        delete l->mpqr;

                    } else {
                        h1 = l->mpqr->myhpq;
                        h2 = l->mpqr->myhpq;
                        n_lt = nullptr;

                        MPQ.delete_element(l->mpqr->qel);
                        delete l->mpqr;
                    }
                    if (r->mpqr) {
                        n_rt = r->mpqr->rt;
                        if (n_rt == r) n_rt = nullptr;
                        if (n_rt) n_rt -> mpql = n_m;

                        h3 = r->mpqr->myhpq;
                        h1->merge(h3);
                        MPQ.delete_element(r->mpqr->qel);
                        delete r->mpqr;

                        n_hpq = h1;
                        if (n_rt) n_rt -> mpql = n_m;
                    } else {
                        n_rt = nullptr;
                        n_hpq = h1;
                    }
                } else if (l->t) { // the left node is a leaf
                    if (l->mpql) {
                        n_lt = l->mpql->lt;
                        if (n_lt) n_lt->mpqr = n_m;
                        n_rt = l->mpqr->rt;
                        if (n_rt) n_rt ->mpql = n_m;

                        l -> mpql ->myhpq -> merge(l->mpqr->myhpq);
                        n_hpq=l->mpql->myhpq;
                        MPQ.delete_element(l->mpql->qel);
                        MPQ.delete_element(l->mpqr->qel);
                        delete l->mpql;
                        delete l->mpqr;
                    } else {
                        n_lt = nullptr;
                        n_rt = l->mpqr->rt;
                        if (n_rt) n_rt->mpql = n_m;

                        n_hpq = l->mpqr->myhpq;
                        MPQ.delete_element(l->mpqr->qel);
                        delete l->mpqr;
                    }
                } else if (r->t) { // right node is a leaf
                    if (r->mpqr) {
                        n_lt = r->mpql->lt;
                        if (n_lt) n_lt->mpqr = n_m;
                        n_rt = r->mpqr->rt;
                        if (n_rt) n_rt->mpql = n_m;

                        r -> mpql ->myhpq -> merge(r->mpqr->myhpq);
                        n_hpq=r->mpql->myhpq;
                        MPQ.delete_element(r->mpql->qel);
                        MPQ.delete_element(r->mpqr->qel);
                        delete r->mpql;
                        delete r->mpqr;
                    } else {
                        n_lt = r->mpql->lt;
                        if (n_lt) n_lt->mpqr = n_m;
                        n_rt = nullptr;

                        n_hpq = r->mpql->myhpq;
                        MPQ.delete_element(r->mpql->qel);
                        delete r->mpql;
                    }
                } else {
                    // merge of two inner nodes
                    // no need to merge hpqs
                    MPQ.delete_element(m->qel);

                    n_hpq = m->myhpq;
                    n_lt = m->lt;
                    n_rt = m->rt;

                    if (n_lt) n_lt->mpqr = n_m;
                    if (n_rt) n_rt->mpql = n_m;

                    delete m;
                }

                // create a new node with the information gained above
                ht_node* new_node = new ht_node();
                new_node -> c = ' ';
                new_node -> w = l->w + r->w;
                new_node -> t = false;
                new_node -> pos = lpos;
                new_node -> left = l;
                new_node -> right = r;
                // insert node to the correct hpq
                new_node -> ql = n_hpq->insert(new_node);

                // update working sequence
                A[lpos] = new_node;
                A[rpos] = nullptr;
                // update information in the new master node and reinsert it to mpq
                ht_node* tmp_min = n_hpq->find_min()->item;
                heap_node<ht_node>* tmpsnd = n_hpq->find_snd_min();
                if (tmpsnd) {
                    ht_node* tmp_snd = n_hpq->find_snd_min()->item;
                    n_m->min_sum = tmp_min->w + tmp_snd->w;

                    if (tmp_min -> pos < tmp_snd->pos) {
                        n_m->i = tmp_min -> pos;
                        n_m->j = tmp_snd -> pos;
                    } else {
                        n_m->i = tmp_snd -> pos;
                        n_m->j = tmp_min -> pos;
                    }
                    n_m->qel = MPQ.insert(n_m);
                    n_m->myhpq = n_hpq;
                    n_m->lt = n_lt;
                    n_m->rt = n_rt;
                } else {
                    // free the last remaining hpq
                    n_hpq->free_memory();
                    delete n_m;
                }
            }

            // level assignment and deletion of unneeded nodes
            assign_level(A[0], 0);

            // reconstruction phase using the stack algorithm
            ht_node* stack[sigma];

            for (size_type i = 0; i < sigma; i++) {
                stack[i] = nullptr;
                temp_nodes.emplace_back(pc_node(T[i].w, (size_type)T[i].c));
                T[i].pos = i;
            }

            int64_t spointer = -1;
            uint64_t qpointer = 0; // use the Array T as a stack
            int64_t max_nodes = sigma;
            while (qpointer < sigma or spointer >= 1LL) {
                if (spointer >= 1LL and
                    (stack[spointer]->level == stack[spointer-1]->level)) {
                    ht_node* n_node = new ht_node();
                    max_nodes++;
                    n_node->t = false;
                    n_node->left = stack[spointer-1];
                    n_node->right = stack[spointer];
                    n_node->level = stack[spointer]->level-1;
                    n_node->w = stack[spointer]->w + stack[spointer-1]->w;
                    n_node->c = '|';

                    n_node->pos = temp_nodes.size();
                    temp_nodes[stack[spointer-1]->pos].parent = temp_nodes.size();
                    temp_nodes[stack[spointer]->pos].parent = temp_nodes.size();
                    temp_nodes.emplace_back(pc_node(n_node->w, 0,
                                                    pc_node::undef,
                                                    stack[spointer-1]->pos,
                                                    stack[spointer]->pos));

                    if (!stack[spointer-1]->t) delete stack[spointer-1];
                    if (!stack[spointer]->t) delete stack[spointer];

                    stack[--spointer] = n_node;
                } else {
                    stack[++spointer] = &T[qpointer++];
                }
            }
            delete stack[0];
        }

        static void assign_level(ht_node* n, int64_t lvl) {
            if (n) {
                n->level = lvl;
                assign_level(n->left, lvl + 1);
                assign_level(n->right, lvl + 1);

                if (!n->t) {
                    delete n;
                }
            }
        }
};

struct hutu_shape {
    template<class t_wt>
    using type = _hutu_shape<t_wt>;
};



}// end namespace sdsl

#endif // end file
