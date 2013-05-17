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
/*! \file wt_hutu.hpp
    \brief wt_hutu.hpp contains a class for the wavelet tree of byte sequences which is in Hu Tucker shape.
    \author Simon Gog, Timo Beller and Markus Brenner
*/
#ifndef INCLUDED_SDSL_WT_HUTU
#define INCLUDED_SDSL_WT_HUTU

#include "int_vector.hpp"
#include "rank_support_v.hpp"
#include "rank_support_v5.hpp"
#include "select_support_mcl.hpp"
#include "rrr_vector.hpp"
#include "bits.hpp"
#include "util.hpp"
#include "wt_helper.hpp"
#include "wt.hpp"
#include <algorithm> // for std::swap
#include <stdexcept>
#include <vector>
#include <utility> // for pair
#include <deque>
#include <queue>
#include <iostream>


#ifdef SDSL_DEBUG
#define SDSL_DEBUG_WT_HUTU
#endif


//#define SDSL_DEBUG_WAVELET_TREE

#ifdef SDSL_DEBUG_WAVELET_TREE
#include "testutils.hpp"
#endif

//! Namespace for the succinct data structure library.
namespace sdsl
{

// predeclarations of node classes
template<class size_type> class HTNode;

template<class size_type> class MNode;

// node as used by the leftist heap
template <class Element>
class HeapNode
{
    public:
        Element* item;   // pointer to the represented item

        HeapNode* left;  // pointer to left child
        HeapNode* right; // pointer to right child
        HeapNode* parent; // pointer to parent node (needed to delete elements)
        int rank; // rank of the heap node

        // constructor
        HeapNode() {
            item = NULL;
            left = NULL;
            right = NULL;
            parent = NULL;
            rank = 0;
        }

        // comparison operators needed for the heap structure
        bool operator< (HeapNode<Element> other) {
            return item->operator<(* other.item);
        }
        bool operator> (HeapNode<Element> other) {
            return *item < other.*item;
        }
};

// Implementation of a leftist heap as needed in the first phase of HuTucker Code construction
template <class Element>
class LHeap
{
    private:
        HeapNode<Element> * root; // pointer to the root

        // fixes node information after the deletion of elements
        void fixNode(HeapNode<Element> * item) {
            if (item) {
                if (!item->left || !item->right) { // if node has only one child
                    if (item->rank != 0) { // only go on fixing if the node information needs to be changed
                        item->rank = 0;
                        if (item->parent) fixNode(item->parent);
                    }
                } else { // node information has to be adapted
                    int nn = (item->left->rank > item->right->rank) ? item->right->rank : item->left->rank;
                    if (item->rank != nn && item->parent != 0) {
                        item->rank = nn;
                        fixNode(item->parent);
                    }
                }
            }
        }

        // helper function to remove the data structure from memory
        void freeNode(HeapNode<Element> * item) {
            if (item->left) {
                freeNode(item->left);
                delete item->left;
                item->left = NULL;
            }
            if (item->right) {
                freeNode(item->right);
                delete item->right;
                item->right = NULL;
            }
        }

        // internal merge function
        HeapNode<Element> * merge(HeapNode<Element> *h1, HeapNode<Element> *h2) {
            if (!h1) return h2;
            if (!h2) return h1;

            if (*(h1->item) < *(h2->item))   return merge1(h1, h2);
            else                            return merge1(h2, h1);
        }
        // internal merge function
        HeapNode<Element> * merge1(HeapNode<Element> *h1, HeapNode<Element> *h2) {
            if (!h1->left) { // if h1 has no children, the merge is simple
                h1->left = h2;
                h2->parent = h1; // adjust the parent pointer
            } else {
                h1->right = merge(h1->right, h2);
                if (h1->right) h1->right->parent = h1;

                if (h1->left->rank < h1->right->rank) {
                    HeapNode<Element> *tmp = h1->left;
                    h1->left = h1->right;
                    h1->right = tmp;
                }

                h1 -> rank = h1 -> right -> rank + 1;
            }
            return h1;
        }

    public:
        // constructor
        LHeap() {
            root = NULL;
        }

        bool isEmpty() {
            return (root==NULL);
        }

        // the minimum is always at the root
        HeapNode<Element> * findMin() {
            return root;
        }

        // returns the second smallest node, which is one of
        // the two children of the root (if existent)
        HeapNode<Element> * findSndMin() {
            if (!root) return NULL;
            if (root->left == NULL) return root->right;
            if (root->right == NULL) return root->left;

            if (root->left->operator< (*root->right)) return root->left;
            else return root->right;
        }

        // inserts a node into the heap
        HeapNode<Element> * insert(Element* x) {
            HeapNode<Element> * n = new HeapNode<Element>();
            n->item = x;

            LHeap<Element> lh;
            lh.root = n;
            merge(&lh);
            return n;
        }

        // removes the minimum (the root) from the heap
        void deleteMin() {
            HeapNode<Element> *oldRoot = root;
            root = merge(root->left, root->right);
            if (root) root->parent = 0;
            delete oldRoot;
        }

        // deletes an arbitrary element from the heap
        // this function assumes, that item is an element of the heap
        void deleteElement(HeapNode<Element> * item) {
            if (item) {
                if (root == item) { // trivial case: we want to delete the root
                    deleteMin();
                } else {
                    // otherwise we have to adapt the parent node and
                    // the children of item
                    HeapNode<Element> * h1 = merge(item->left, item->right);
                    if (h1) h1->parent = item->parent;
                    if (item == item->parent->left) {
                        item->parent->left = h1;
                    } else if (item == item->parent->right) {
                        item->parent->right = h1;
                    }
                    // fix node information considering rank
                    fixNode(item->parent);
                    delete item; // remove the item from memory
                }
            }
        }

        // public merge function
        void merge(LHeap<Element> *rhs) {
            root = merge(root, rhs->root);
            rhs->root = NULL;
        }

        // removes the whole data structure from memory
        void freeMemory() {
            if (root) {
                freeNode(root);
                delete root;
                root = NULL;
            }
        }
};

// Master node as used in the first phase of the HuTucker algorithm
template <class size_type>
class MNode
{
    public:
        size_type MinSum; // minimal sum of the two minimal elements of the hpq this node points to
        int i; // position of the left node in the working sequence
        int j; // position of the right node in the working sequence
        HeapNode<MNode<size_type> > * qel; // pointer to the corresponding heap element (used for deletion)
        LHeap<HTNode<size_type> > * myhpq; // pointer to the hpq

        HTNode<size_type> * lt; // pointer to the left- and rightmost leafs of the hpq
        HTNode<size_type> * rt; // need for merge operations

        MNode() {
            qel = 0;
            myhpq = 0;
            lt = 0;
            rt = 0;
        }

        bool operator< (const MNode other) {
            if (MinSum < other.MinSum) return true;
            else if (other.MinSum < MinSum) return false;
            else if (i < other.i) return true;
            else if (other.i < i) return false;
            else return j < other.j;
        }
        bool operator> (const MNode other) {
            if (MinSum > other.MinSum) return true;
            else if (other.MinSum > MinSum) return false;
            else if (i > other.i) return true;
            else if (other.i > i) return false;
            else return j > other.j;
        }
        bool operator== (const MNode other) {
            return MinSum == other.MinSum;
        }
};

// HuTucker node as used in the first phase of the HuTucker algorithm
template <class size_type>
class HTNode
{
    public:
        int pos;     // position of the node
        char c;      // the represented letter
        size_type w; // frequency of the node
        bool t;      // wether the node is a leaf
        int level;   // level in the tree

        MNode<size_type> *mpql; // pointer to the two master nodes (as a node can belong to up to two hpqs)
        MNode<size_type> *mpqr; // only mpql is used for inner nodes

        HeapNode<HTNode> * ql; // pointer to the two heap nodes (as a node can belong to up to two hpqs)
        HeapNode<HTNode> * qr; // only ql is used for inner nodes

        HTNode<size_type> *left;  // left child
        HTNode<size_type> *right; // right child

        HTNode() {
            mpql = 0;
            mpqr = 0;
            ql = 0;
            qr = 0;
            left = NULL;
            right = NULL;
        }

        bool operator< (const HTNode other) {
            if (w < other.w) return true;
            if (other.w < w) return false;
            else {
                return pos < other.pos;
            }

            return w < other.w;
        }
        bool operator> (const HTNode other) {
            if (w > other.w) return true;
            if (other.w < w) return false;
            else {
                return pos > other.pos;
            }
        }
        bool operator== (const HTNode other) {
            return w == other.w;
        }
};


const int_vector<>::size_type ZOO[2] = {0, (int_vector<>::size_type)-1};

template<class size_type>
struct _node_ht {
    size_type 	tree_pos; 		// pointer into the bit_vector, which represents the wavelet tree
    size_type 	tree_pos_rank;	// pre-calculated rank for the prefix up to but not including tree_pos
    uint16_t	parent;			// pointer to the parent
    uint16_t	child[2];		// pointer to the children

    _node_ht(size_type tree_pos=0, size_type tree_pos_rank=0, uint16_t parent=_undef_node,
             uint16_t child_left=_undef_node, uint16_t child_right=_undef_node):
        tree_pos(tree_pos), tree_pos_rank(tree_pos_rank), parent(parent) {
        child[0] = child_left;
        child[1] = child_right;
    }

    _node_ht& operator=(const _node_ht& v) {
        if (this != &v) {
            tree_pos 		= v.tree_pos;
            tree_pos_rank 	= v.tree_pos_rank;
            parent			= v.parent;
            child[0] 		= v.child[0];
            child[1] 		= v.child[1];
        }
        return *this;
    }

    size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
        structure_tree_node* st_child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += write_member(tree_pos, out);
        written_bytes += write_member(tree_pos_rank, out);
        written_bytes += write_member(parent, out);
        out.write((char*)child, 2*sizeof(child[0]));
        written_bytes += 2*sizeof(child[0]);
        structure_tree::add_size(st_child, written_bytes);
        return written_bytes;
    }

    void load(std::istream& in) {
        read_member(tree_pos, in);
        read_member(tree_pos_rank, in);
        read_member(parent, in);
        in.read((char*) child, 2*sizeof(child[0]));
    }
};

//! A Wavelet Tree class for byte sequences.
/*!
 * A wavelet tree is build for a vector of characters over the alphabet \f$\Sigma\f$.
 * This class should be used only for small alphabets \f$\Sigma \ll n\f$ (see int_wavelet_tree for a wavelet tree for big alphabets).
 * The wavelet tree \f$wt\f$ consists of a tree of bitvectors and provides three efficient methods:
 *   - The "[]"-operator: \f$wt[i]\f$ returns the ith symbol of vector for which the wavelet tree was build for.
 *   - The rank method: \f$wt.rank(i,c)\f$ returns the number of occurences of symbol \f$c\f$ in the prefix [0..i-1] in the vector for which the wavelet tree was build for.
 *   - The select method: \f$wt.select(j,c)\f$ returns the index \f$i\in [0..size()-1]\f$ of the jth occurence of symbol \f$c\f$.
 *
 *	\par Space complexity
 *		\f$\Order{n H_0 + 2|\Sigma|\log n}\f$ bits, where \f$n\f$ is the size of the vector the wavelet tree was build for.
 *
 *   @ingroup wt
 */

template<class BitVector 		 = bit_vector,
              class RankSupport 		 = typename BitVector::rank_1_type,
              class SelectSupport	 = typename BitVector::select_1_type,
              class SelectSupportZero = typename BitVector::select_0_type,
              bool dfs_shape=0 >
class wt_hutu
{
    public:
        typedef int_vector<>::size_type size_type;
        typedef unsigned char           value_type;
        typedef BitVector               bit_vector_type;
        typedef RankSupport             rank_1_type;
        typedef SelectSupport           select_1_type;
        typedef SelectSupportZero       select_0_type;
        typedef wt_tag                  index_category;
        typedef byte_alphabet_tag       alphabet_category;
        enum { lex_ordered=1 };

    private:
#ifdef WT_HUTU_CACHE
        mutable value_type m_last_access_answer;
        mutable size_type  m_last_access_i;
        mutable size_type  m_last_access_rl;
#endif

        size_type 			m_size;
        size_type 			m_sigma; 		//<- \f$ |\Sigma| \f$
        bit_vector_type		m_tree;			// bit vector to store the wavelet tree
        RankSupport			m_tree_rank;	// rank support for the wavelet tree bit vector
        SelectSupport		m_tree_select1;	// select support for the wavelet tree bit vector
        SelectSupportZero	m_tree_select0;

        _node_ht<size_type> 	m_nodes[511]; // nodes for the HuTucker tree structure
        uint16_t			m_c_to_leaf[256]; // map symbol c to a leaf in the tree structure
        // if m_c_to_leaf[c] == _undef_node the char does not exists in the text
        uint64_t			m_path[256];     // path information for each char:
        // the bits at position 0..55 hold path information and
        // bits 56..63 the length of the path in binary representation

        void copy(const wt_hutu& wt) {
            m_size 			= wt.m_size;
            m_sigma 		= wt.m_sigma;
            m_tree			= wt.m_tree;
            m_tree_rank 	= wt.m_tree_rank;
            m_tree_rank.set_vector(&m_tree);
            m_tree_select1	= wt.m_tree_select1;
            m_tree_select1.set_vector(&m_tree);
            m_tree_select0	= wt.m_tree_select0;
            m_tree_select0.set_vector(&m_tree);
            for (size_type i=0; i < 511; ++i)
                m_nodes[i] = wt.m_nodes[i];
            for (size_type i=0; i<256; ++i)
                m_c_to_leaf[i] = wt.m_c_to_leaf[i];
            for (size_type i=0; i<256; ++i) {
                m_path[i] = wt.m_path[i];
            }
        }

        // insert a character into the wavelet tree, see constuct method
        void insert_char(uint8_t old_chr, size_type* tree_pos, size_type times, bit_vector& f_tree) {
            uint32_t path_len = (m_path[old_chr]>>56);
            uint64_t p = m_path[old_chr];
            for (uint32_t node=0, l=0; l<path_len; ++l, p >>= 1) {
                if (p&1) {
                    f_tree.set_int(tree_pos[node], 0xFFFFFFFFFFFFFFFFULL,times);
                }
                tree_pos[node] += times;
                node = m_nodes[node].child[p&1];
            }
        }


        // constructs the Hu-Tucker tree, writes the node to the given array and returns the number of nodes
        size_type construct_hutucker_tree(size_type* C, _node_ht<size_type>* tmp_nodes) {
            //create a leaf for every letter
            std::vector<HTNode<size_type> > node_vector;
            for (int i = 0; i < 256; i++) {
                if (C[i]) {
                    HTNode<size_type> n;
                    n.c = (char)i;
                    n.w = C[i];
                    n.t = true;
                    n.pos = node_vector.size();
                    node_vector.push_back(n);
                }
            }

            if (node_vector.size() == 1) {
                // special case of an alphabet of size 1:
                // just instantly create the tree and return it
                tmp_nodes[0] = _node_ht<size_type>(node_vector[0].w, (size_type)node_vector[0].c);
                return 1;
            }

            size_type sigma = node_vector.size();
            // physical Leafs
            HTNode<size_type> T[sigma];
            // the current working sequence
            HTNode<size_type> *A[sigma];
            // Priority Queues, containing the Huffman Sequences
            LHeap<HTNode<size_type> > HPQ[sigma];
            // Master Priority Queue
            LHeap<MNode<size_type> > MPQ;

            // init T, A, HPQs and MPQ
            T[0] = node_vector[0];
            A[0] = &T[0];

            // initialization needed for every leaf
            for (size_type i = 1; i < sigma; i++) {
                T[i] = node_vector[i];
                A[i] = &T[i];

                T[i - 1].qr = HPQ[i - 1].insert(&T[i - 1]);
                T[i].ql = HPQ[i - 1].insert(&T[i]);

                MNode<size_type> * m = new MNode<size_type>();
                m->MinSum = T[i - 1].w + T[i].w;
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
                MNode<size_type> * m = MPQ.findMin()->item;
                HTNode<size_type> * l = A[m->i];
                HTNode<size_type> * r = A[m->j];
                int lpos = m->i;
                int rpos = m->j;

                LHeap<HTNode<size_type> > * n_hpq = NULL;
                HTNode<size_type> * n_rt = NULL;
                HTNode<size_type> * n_lt = NULL;

                // create a new master priority queue
                MNode<size_type> * n_m = new MNode<size_type>();
                // delete old nodes from all hpqs
                if (l->t) {
                    if (l->mpql) l->mpql->myhpq->deleteElement(l->ql);
                    l->ql = NULL;
                    if (l->mpqr) l->mpqr->myhpq->deleteElement(l->qr);
                    l->qr = NULL;
                } else {
                    m->myhpq->deleteElement(l->ql);
                    l->ql = NULL;
                }
                if (r->t) {
                    if (r->mpql) r->mpql->myhpq->deleteElement(r->ql);
                    l->ql = NULL;

                    if (r->mpqr) r->mpqr->myhpq->deleteElement(r->qr);
                    r->qr = NULL;
                } else {
                    m->myhpq->deleteElement(r->ql);
                    r->ql = NULL;
                }
                // handle the merge of hpqs
                if (l->t && r ->t) {
                    // both nodes are leaves
                    LHeap<HTNode<size_type> > * h1 = NULL;
                    LHeap<HTNode<size_type> > * h2 = NULL;
                    LHeap<HTNode<size_type> > * h3 = NULL;
                    if (l -> mpql) {
                        n_lt = l->mpql->lt;
                        if (n_lt == l) n_lt = NULL;
                        if (n_lt) n_lt -> mpqr = n_m;

                        h1 = l->mpql->myhpq;
                        h2 = l->mpqr->myhpq;

                        h1 -> merge(h2);
                        MPQ.deleteElement(l->mpql->qel);
                        MPQ.deleteElement(l->mpqr->qel);
                        delete l->mpql;
                        delete l->mpqr;

                    } else {
                        h1 = l->mpqr->myhpq;
                        h2 = l->mpqr->myhpq;
                        n_lt = NULL;

                        MPQ.deleteElement(l->mpqr->qel);
                        delete l->mpqr;
                    }
                    if (r->mpqr) {
                        n_rt = r->mpqr->rt;
                        if (n_rt == r) n_rt = NULL;
                        if (n_rt) n_rt -> mpql = n_m;

                        h3 = r->mpqr->myhpq;
                        h1->merge(h3);
                        MPQ.deleteElement(r->mpqr->qel);
                        delete r->mpqr;

                        n_hpq = h1;
                        if (n_rt) n_rt -> mpql = n_m;
                    } else {
                        n_rt = NULL;
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
                        MPQ.deleteElement(l->mpql->qel);
                        MPQ.deleteElement(l->mpqr->qel);
                        delete l->mpql;
                        delete l->mpqr;
                    } else {
                        n_lt = NULL;
                        n_rt = l->mpqr->rt;
                        if (n_rt) n_rt->mpql = n_m;

                        n_hpq = l->mpqr->myhpq;
                        MPQ.deleteElement(l->mpqr->qel);
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
                        MPQ.deleteElement(r->mpql->qel);
                        MPQ.deleteElement(r->mpqr->qel);
                        delete r->mpql;
                        delete r->mpqr;
                    } else {
                        n_lt = r->mpql->lt;
                        if (n_lt) n_lt->mpqr = n_m;
                        n_rt = NULL;

                        n_hpq = r->mpql->myhpq;
                        MPQ.deleteElement(r->mpql->qel);
                        delete r->mpql;
                    }
                } else {
                    // merge of two inner nodes
                    // no need to merge hpqs
                    MPQ.deleteElement(m->qel);

                    n_hpq = m->myhpq;
                    n_lt = m->lt;
                    n_rt = m->rt;

                    if (n_lt) n_lt->mpqr = n_m;
                    if (n_rt) n_rt->mpql = n_m;

                    delete m;
                }

                // create a new node with the information gained above
                HTNode<size_type> * new_node = new HTNode<size_type>();
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
                A[rpos] = NULL;
                // update information in the new master node and reinsert it to the mpq
                HTNode<size_type> * tmp_min = n_hpq->findMin()->item;
                HeapNode<HTNode<size_type> > * tmpsnd = n_hpq->findSndMin();
                if (tmpsnd) {
                    HTNode<size_type> * tmp_snd = n_hpq->findSndMin()->item;
                    n_m->MinSum = tmp_min->w + tmp_snd->w;

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
                    n_hpq->freeMemory();
                    delete n_m;
                }
            }

            // level assignment and deletion of unneeded nodes
            assignLevel(A[0], 0);

            // reconstruction phase using the stack algorithm
            HTNode<size_type>* stack[sigma];

            size_type node_count=0;
            for (size_type i = 0; i < sigma; i++) {
                stack[i] = NULL;
                tmp_nodes[i] = _node_ht<size_type>(T[i].w, (size_type)T[i].c);
                T[i].pos = i;
                node_count++;
            }

            int spointer = -1;
            unsigned int qpointer = 0; // use the Array T as a stack
            int max_nodes = sigma;
            while (qpointer < sigma || spointer >= 1) {
                if (spointer >= 1 && (stack[spointer]->level == stack[spointer-1]->level)) {
                    HTNode<size_type>* n_node = new HTNode<size_type>();
                    max_nodes++;
                    n_node->t = false;
                    n_node->left = stack[spointer-1];
                    n_node->right = stack[spointer];
                    n_node->level = stack[spointer]->level-1;
                    n_node->w = stack[spointer]->w + stack[spointer-1]->w;
                    n_node->c = '|';

                    n_node->pos = node_count;
                    tmp_nodes[stack[spointer-1]->pos].parent = node_count;
                    tmp_nodes[stack[spointer]->pos].parent = node_count;
                    tmp_nodes[node_count++] = _node_ht<size_type>(n_node->w, 0, _undef_node, stack[spointer-1]->pos, stack[spointer]->pos);

                    if (!stack[spointer-1]->t) delete stack[spointer-1];
                    if (!stack[spointer]->t) delete stack[spointer];

                    stack[--spointer] = n_node;
                } else {
                    stack[++spointer] = &T[qpointer++];
                }
            }
            delete stack[0];
            return node_count;
        }

        void assignLevel(HTNode<size_type> *n, int lvl) {
            if (n) {
                n->level = lvl;
                assignLevel(n->left, lvl + 1);
                assignLevel(n->right, lvl + 1);

                if (!n->t) {
                    delete n;
                }
            }
        }

        // calculates the HuTucker tree and returns the size of the WT bit vector
        size_type construct_wavelet_tree(size_type* C) {
            _node_ht<size_type> temp_nodes[2*m_sigma-1];
            size_type node_cnt = construct_hutucker_tree(C, temp_nodes);
            // Convert HuTucker tree into breadth first search order in memory and
            // calculate tree_pos values
            m_nodes[0] = temp_nodes[node_cnt-1];  // insert root at index 0
            size_type tree_size = 0;
            node_cnt = 1;
            uint16_t last_parent = _undef_node;
            std::deque<size_type> q;
            q.push_back(0);
            while (!q.empty()) {
                size_type idx;
                if (!dfs_shape) {
                    idx = q.front(); q.pop_front();
                } else {
                    idx = q.back(); q.pop_back();
                }
                size_type frq = m_nodes[idx].tree_pos; // frq_sum was stored in tree_pos
                m_nodes[idx].tree_pos = tree_size;
                if (m_nodes[idx].child[0] != _undef_node)  // if node is not a leaf
                    tree_size += frq;					   // add frequency, as leaves have size 0
                if (idx > 0) { // node is not the root
                    if (last_parent != m_nodes[idx].parent)
                        m_nodes[m_nodes[idx].parent].child[0] = idx;
                    else
                        m_nodes[m_nodes[idx].parent].child[1] = idx;
                    last_parent = m_nodes[idx].parent;
                }
                if (m_nodes[idx].child[0] != _undef_node) { // if node is not a leaf
                    for (size_type k=0; k<2; ++k) {			// add children to tree
                        m_nodes[node_cnt] = temp_nodes[ m_nodes[idx].child[k] ];
                        m_nodes[node_cnt].parent = idx;
                        q.push_back(node_cnt);
                        m_nodes[idx].child[k] = node_cnt++;
                    }
                }
            }

            // initialize m_c_to_leaf
            for (size_type i=0; i<256; ++i)
                m_c_to_leaf[i] = _undef_node; // if c is not in the alphabet m_c_to_leaf[c] = _undef_node
            for (size_type i=0; i < 2*sigma-1; ++i) {
                if (m_nodes[i].child[0] == _undef_node) 				// if node is a leaf
                    m_c_to_leaf[(uint8_t)m_nodes[i].tree_pos_rank] = i; // calculate value
            }
            // initialize path information
            // Note: In the case of a bfs search order,
            // we can classify nodes as rigth child and left child with an easy criterion:
            //       node is a left child, if node%2==1
            //       node is a rigth child, if node%2==0
            for (size_type c=0; c<256; ++c) {
                if (m_c_to_leaf[c] != _undef_node) { // if char exists in the alphabet
                    size_type node = m_c_to_leaf[c];
                    uint64_t w = 0; // path
                    uint64_t l = 0; // path len
                    while (node != 0) { // while node is not the root
                        w <<= 1;
                        if (m_nodes[m_nodes[node].parent].child[1] == node)
                            w |= 1ULL;
                        ++l;
                        node = m_nodes[node].parent; // go up the tree
                    }
                    if (l > 56) {
                        std::cerr<<"HuTucker tree has max depth > 56!!! ERROR"<<std::endl;
                        throw std::logic_error("HuTucker tree size is greater than 56!!!");
                    }
                    m_path[c] = w | (l << 56);
                } else {
                    m_path[c] = 0; // i.e. len is also 0, good for special case in rank()
                }
            }
            return tree_size;
        }

        void construct_init_rank_select() {
            util::init_support(m_tree_rank, &m_tree);
            util::init_support(m_tree_select0, &m_tree);
            util::init_support(m_tree_select1, &m_tree);
        }

        void construct_precalc_node_ranks() {
            for (size_type i=0; i<2*m_sigma-1; ++i) {
                if (m_nodes[i].child[0] != _undef_node)  // if node is not a leaf
                    m_nodes[i].tree_pos_rank = m_tree_rank(m_nodes[i].tree_pos);
            }
        }

        // recursive internal version of the method interval_symbols
        void _interval_symbols(size_type i, size_type j, size_type& k,
                               std::vector<unsigned char>& cs,
                               std::vector<size_type>& rank_c_i,
                               std::vector<size_type>& rank_c_j, uint16_t node) const {
            // invariant: j>i
            size_type i_new = (m_tree_rank(m_nodes[node].tree_pos + i) - m_nodes[node].tree_pos_rank);
            size_type j_new = (m_tree_rank(m_nodes[node].tree_pos + j) - m_nodes[node].tree_pos_rank);
            i -= i_new; j -= j_new;
            // goto left child
            if (i != j) {
                uint16_t node_new = m_nodes[node].child[0];
                // if node is not a leaf
                if (m_nodes[node_new].child[0] != _undef_node) {
                    _interval_symbols(i, j, k, cs, rank_c_i, rank_c_j, node_new);
                } else {
                    rank_c_i[k] = i;
                    rank_c_j[k] = j;
                    cs[k++] = m_nodes[node_new].tree_pos_rank;
                }
            }
            // goto right child
            if (i_new!=j_new) {
                uint16_t node_new = m_nodes[node].child[1];
                // if node is not a leaf
                if (m_nodes[node_new].child[0] != _undef_node) {
                    _interval_symbols(i_new, j_new, k, cs, rank_c_i, rank_c_j, node_new);
                } else {
                    rank_c_i[k] = i_new;
                    rank_c_j[k] = j_new;
                    cs[k++] = m_nodes[node_new].tree_pos_rank;
                }
            }
        }


    public:

        const size_type& sigma;
        const bit_vector_type& tree;

        // Default constructor
        wt_hutu():m_size(0),m_sigma(0), sigma(m_sigma),tree(m_tree) {};



        //! Constructor
        /*!
         *	\param rac Reference to the vector (or unsigned char array) for which the wavelet tree should be build.
         *	\param size Size of the prefix of the vector (or unsigned char array) for which the wavelet tree should be build.
         *	\par Time complexity
         *		\f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         */
        template<typename RandomAccessContainer>
        wt_hutu(const RandomAccessContainer& rac, size_type size):m_size(size), m_sigma(0), sigma(m_sigma), tree(m_tree) {
            construct(rac, size);
        }

        template<uint8_t w>
        wt_hutu(const int_vector<w>& rac):m_size(rac.size()), m_sigma(0), sigma(m_sigma), tree(m_tree) {
            construct(rac, rac.size());
        }


        template<typename RandomAccessContainer>
        void construct(const RandomAccessContainer& rac, size_type size) {
            m_size = size;
            if (m_size == 0)
                return;
            // O(n + |\Sigma|\log|\Sigma|) algorithm for calculating node sizes
            size_type C[256] = {0};
            //  1. Count occurences of characters
            for (size_type i=0; i < size; ++i) {
                ++C[rac[i]];
            }
            // 2. Calculate effective alphabet size
            calculate_effective_alphabet_size(C, m_sigma);
            // 3. Generate HuTucker tree
            size_type tree_size = construct_wavelet_tree(C);
            // 4. Generate wavelet tree bit sequence m_tree
            bit_vector tmp_tree(tree_size, 0);  // initialize bit_vector for the tree
            //  Calculate starting position of wavelet tree nodes
            size_type tree_pos[511];
            for (size_type i=0; i < 2*sigma-1; ++i) {
                tree_pos[i] = m_nodes[i].tree_pos;
            }
            uint8_t old_chr = rac[0], times = 0;
            for (size_type i=0; i < m_size; ++i) {
                uint8_t chr = rac[i];
                if (chr	!= old_chr) {
                    insert_char(old_chr, tree_pos, times, tmp_tree);
                    times = 1;
                    old_chr = chr;
                } else { // chr == old_chr
                    ++times;
                    if (times == 64) {
                        insert_char(old_chr, tree_pos, times, tmp_tree);
                        times = 0;
                    }
                }
            }
            if (times > 0) {
                insert_char(old_chr, tree_pos, times, tmp_tree);
            }
            util::assign(m_tree, tmp_tree);
            // 5. Initialize rank and select data structures for m_tree
            construct_init_rank_select();
            // 6. Finish inner nodes by precalculating the tree_pos_rank values
            construct_precalc_node_ranks();
        }

        wt_hutu(int_vector_file_buffer<8>& rac, size_type size):m_size(size), m_sigma(0), sigma(m_sigma), tree(m_tree) {
            construct(rac, size);
        }

        //! Construct the wavelet tree from a random access container
        /*! \param rac A random access container
         *	\param size The length of the prefix of the random access container, for which the wavelet tree should be build
         */
        void construct(int_vector_file_buffer<8>& rac, size_type size) {
            m_size = size;
            if (m_size == 0)
                return;
            // O(n + |\Sigma|\log|\Sigma|) algorithm for calculating node sizes
            size_type C[256] = {0};
            // 1. Count occurences of characters
            calculate_character_occurences(rac, m_size, C);
            // 2. Calculate effective alphabet size
            calculate_effective_alphabet_size(C, m_sigma);
            // 3. Generate HuTucker tree
            size_type tree_size = construct_wavelet_tree(C);

            // 4. Generate wavelet tree bit sequence m_tree

            bit_vector tmp_tree(tree_size, 0);  // initialize bit_vector for the tree
            //  Calculate starting position of wavelet tree nodes
            size_type tree_pos[511];
            for (size_type i=0; i < 2*sigma-1; ++i) {
                tree_pos[i] = m_nodes[i].tree_pos;
            }
            rac.reset();
            if (rac.int_vector_size < size) {
                throw std::logic_error("wt_huff::construct: stream size is smaller than size!");
                return;
            }
            for (size_type i=0, r_sum=0, r = rac.load_next_block(); r_sum < m_size;) {
                if (r_sum + r > size) {  // read not more than size chars in the next loop
                    r = size-r_sum;
                }
                uint8_t old_chr = rac[i-r_sum], times = 0;
                for (; i < r_sum+r; ++i) {
                    uint8_t chr = rac[i-r_sum];
                    if (chr	!= old_chr) {
                        insert_char(old_chr, tree_pos, times, tmp_tree);
                        times = 1;
                        old_chr = chr;
                    } else { // chr == old_chr
                        ++times;
                        if (times == 64) {
                            insert_char(old_chr, tree_pos, times, tmp_tree);
                            times = 0;
                        }
                    }
                }
                if (times > 0) {
                    insert_char(old_chr, tree_pos, times, tmp_tree);
                }
                r_sum += r; r = rac.load_next_block();
            }
            util::assign(m_tree, tmp_tree);
            // 5. Initialize rank and select data structures for m_tree
            construct_init_rank_select();
            // 6. Finish inner nodes by precalculating the tree_pos_rank values
            construct_precalc_node_ranks();
        }


        //! Copy constructor
        wt_hutu(const wt_hutu& wt):sigma(m_sigma), tree(m_tree) {
            copy(wt);
        }

        //! Assignment operator
        wt_hutu& operator=(const wt_hutu& wt) {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Swap operator
        void swap(wt_hutu& wt) {
            if (this != &wt) {
                std::swap(m_size, wt.m_size);
                std::swap(m_sigma,  wt.m_sigma);
                m_tree.swap(wt.m_tree);
                util::swap_support(m_tree_rank, wt.m_tree_rank, &m_tree, &(wt.m_tree));

                util::swap_support(m_tree_select1, wt.m_tree_select1, &m_tree, &(wt.m_tree));
                util::swap_support(m_tree_select0, wt.m_tree_select0, &m_tree, &(wt.m_tree));

                for (size_type i=0; i < 511; ++i)
                    std::swap(m_nodes[i], wt.m_nodes[i]);
                for (size_type i=0; i<256; ++i)
                    std::swap(m_c_to_leaf[i], wt.m_c_to_leaf[i]);
                for (size_type i=0; i<256; ++i)
                    std::swap(m_path[i], wt.m_path[i]);
            }
        }

        //! Returns the size of the original vector.
        size_type size()const {
            return m_size;
        }

        //! Returns whether the wavelet tree contains no data.
        bool empty()const {
            return m_size == 0;
        }

        //! Recovers the ith symbol of the original vector.
        /*! \param i The index of the symbol in the original vector. \f$i \in [0..size()-1]\f$
         *  \return The ith symbol of the original vector.
         *  \par Time complexity
         *		\f$ \Order{H_0} \f$ on average, where \f$ H_0 \f$ is the zero order entropy
         *      of the sequence.
         */
        value_type operator[](size_type i)const { // TODO: Maybe it is good to integrate a cache here
            // which stores how many of the next symbols are equal
            // with the current char
            assert(i < size());
            size_type node = 0; // start at root node
            while (m_nodes[node].child[0] != _undef_node) { // while node is not a leaf
                if (m_tree[ m_nodes[node].tree_pos + i]) {  // goto the right child
                    i = m_tree_rank(m_nodes[node].tree_pos + i) - m_nodes[node].tree_pos_rank;
                    node = m_nodes[node].child[1];
                } else { // goto the left child
                    i -= (m_tree_rank(m_nodes[node].tree_pos + i) - m_nodes[node].tree_pos_rank);
                    node = m_nodes[node].child[0];
                }
            }
            return m_nodes[node].tree_pos_rank;
        };

        //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
         *  \param c The symbol to count the occurences in the prefix.
         *	\return The number of occurences of symbol c in the prefix [0..i-1] of the supported vector.
         *  \par Time complexity
         *		\f$ \Order{H_0} \f$
         */
        size_type rank(size_type i, value_type c)const {
            assert(i <= size());
            uint64_t p = m_path[c];
            uint32_t path_len = (m_path[c]>>56); // equals zero if char was not present in the original text or m_sigma=1
            if (!path_len and 1 == m_sigma) {    // if m_sigma == 1 return result immediately
                if (m_c_to_leaf[c] == _undef_node) { // if character does not exist return 0
                    return 0;
                }
                return std::min(i, m_size);
            }
            size_type result = i & ZOO[path_len>0]; // important: result has type size_type and ZOO has type size_type
            uint32_t node=0;
            for (uint32_t l=0; l<path_len and result; ++l, p >>= 1) {
                if (p&1) {
                    result 	= (m_tree_rank(m_nodes[node].tree_pos+result) -  m_nodes[node].tree_pos_rank);
                } else {
                    result -= (m_tree_rank(m_nodes[node].tree_pos+result) -  m_nodes[node].tree_pos_rank);
                }
                node = m_nodes[node].child[p&1]; // goto child
            }
            return result;
        };


        //! Calculates for symbol c, how many symbols smaller and bigger c occure in wt[i..j-1].
        /*!
         *  \param i The start index (inclusive) of the interval.
         *  \param j The end index (exclusive) of the interval.
         *  \param c The symbol to count the occurences in the interval.
         *  \param smaller Reference that will contain the number of symbols smaller than c in wt[i..j-1].
         *  \param bigger Reference that will contain the number of symbols bigger than c in wt[i..j-1].
         *  \return The number of occurences of symbol c in wt[0..i-1].
         *
         *  \par Precondition
         *       \f$ i \leq j and i < n and j \leq n \f$
         *       \f$ c must exist in wt \f$
         */
        size_type bounds(size_type i, size_type j, value_type c, size_type& smaller, size_type& bigger)const {
            assert(i<=j and i < size() and j <= size());
            smaller = 0;
            bigger = 0;
            if (1==m_sigma) {
                return i;
            }
            if (i==j) {
                return rank(i,c);
            }
            uint64_t p = m_path[c];
            uint32_t path_len = (m_path[c]>>56); // equals zero if char was not present in the original text
            assert(path_len>0);
            size_type res1 = i;
            size_type res2 = j;
            uint32_t node=0;
            for (uint32_t l=0; l<path_len; ++l, p >>= 1) {
                if (p&1) {
                    size_type r1_1 = (m_tree_rank(m_nodes[node].tree_pos+res1)-m_nodes[node].tree_pos_rank);
                    size_type r1_2 = (m_tree_rank(m_nodes[node].tree_pos+res2)-m_nodes[node].tree_pos_rank);

                    smaller += res2 - r1_2 - res1 + r1_1;

                    res1 = r1_1;
                    res2 = r1_2;
                } else {
                    size_type r1_1 = (m_tree_rank(m_nodes[node].tree_pos+res1)-m_nodes[node].tree_pos_rank);
                    size_type r1_2 = (m_tree_rank(m_nodes[node].tree_pos+res2)-m_nodes[node].tree_pos_rank);

                    bigger += r1_2 - r1_1;

                    res1 -= r1_1;
                    res2 -= r1_2;
                }
                node = m_nodes[node].child[p&1];
            }
            return res1;
        };


        //! Calculates how many occurrences of symbol wt[i] are in the prefix [0..i-1] of the original sequence.
        /*!
         *  \param i The index of the symbol.
         *  \param c Reference that will contain symbol wt[i].
         *  \return The number of occurrences of symbol wt[i] in the prefix [0..i-1]
         *  \par Time complexity
         *        \f$ \Order{H_0} \f$
         */
        size_type inverse_select(size_type i, value_type& c)const {
            assert(i < size());
            uint32_t node=0;
            while (m_nodes[node].child[0] != _undef_node) { // while node is not a leaf
                if (m_tree[m_nodes[node].tree_pos + i]) { // if bit is set at position goto right child
                    i 	= (m_tree_rank(m_nodes[node].tree_pos + i) -  m_nodes[node].tree_pos_rank);
                    node = m_nodes[node].child[1];
                } else { // goto left child
                    i -= (m_tree_rank(m_nodes[node].tree_pos + i) -  m_nodes[node].tree_pos_rank);
                    node = m_nodes[node].child[0];
                }
            }
            c = m_nodes[node].tree_pos_rank;
            return i;
        }

        //! Calculates the ith occurence of the symbol c in the supported vector.
        /*!
         *  \param i The ith occurence. \f$i\in [1..rank(size(),c)]\f$.
         *  \param c The symbol c.
         *  \par Time complexity
         *		\f$ \Order{H_0} \f$
         */
        size_type select(size_type i, value_type c)const {
            assert(i > 0);
            assert(i <= rank(size(), c));
            uint16_t node = m_c_to_leaf[c];
            if (node == _undef_node) { // if c was not present in the original text
                return m_size;		   // -> return a position right to the end
            }
            if (m_sigma == 1) {
                return std::min(i-1,m_size);
            }
            size_type result = i-1;		// otherwise
            uint64_t p = m_path[c];
            uint32_t path_len = (p>>56);
            p <<= (64-path_len); // Note: path_len > 0, since we have handled m_sigma = 1.

            for (uint32_t l=0; l<path_len; ++l, p <<= 1) {
                if ((p & 0x8000000000000000ULL)==0) { // node was a left child
                    node = m_nodes[node].parent;
                    result = m_tree_select0(m_nodes[node].tree_pos-m_nodes[node].tree_pos_rank + result + 1)
                             - m_nodes[node].tree_pos;
                } else { // node was a right child
                    node = m_nodes[node].parent;
                    result = m_tree_select1(m_nodes[node].tree_pos_rank + result + 1)
                             - m_nodes[node].tree_pos;
                }
            }
            return result;
        };


        //! Calculates for each symbol c in wt[i..j-1], how many times c occurs in wt[0..i-1] and wt[0..j-1].
        /*!
         *	\param i The start index (inclusive) of the interval.
         *	\param j The end index (exclusive) of the interval.
         *	\param k Reference that will contain the number of different symbols in wt[i..j-1].
         *  \param cs Reference to a vector that will contain in cs[0..k-1] all symbols that occur in wt[i..j-1] in ascending order.
         *  \param rank_c_i Reference to a vector which equals rank_c_i[p] = rank(i,cs[p]), for \f$ 0 \leq p < k \f$
         *  \param rank_c_j Reference to a vector which equals rank_c_j[p] = rank(j,cs[p]), for \f$ 0 \leq p < k \f$
         *	\par Time complexity
         *		\f$ \Order{\min{\sigma, k \log \sigma}} \f$
         *
         *  \par Precondition
         *       \f$ i \leq j and i < n and j \leq n \f$
         *       \f$ cs.size() \geq \sigma \f$
         *       \f$ rank_c_i.size() \geq \sigma \f$
         *       \f$ rank_c_j.size() \geq \sigma \f$
         */
        void interval_symbols(size_type i, size_type j, size_type& k,
                              std::vector<unsigned char>& cs,
                              std::vector<size_type>& rank_c_i,
                              std::vector<size_type>& rank_c_j) const {
            assert(i<=j and i < size() and j <= size());
            if (i==j) {
                k = 0;
            } else if (1==m_sigma) {
                k = 1;
                cs[0] = m_nodes[0].tree_pos_rank;
                rank_c_i[0] = std::min(i,m_size);
                rank_c_j[0] = std::min(j,m_size);
            } else if ((j-i)==1) {
                k = 1;
                rank_c_i[0] = inverse_select(i, cs[0]);
                rank_c_j[0] = rank_c_i[0]+1;
            } else if ((j-i)==2) {
                rank_c_i[0] = inverse_select(i, cs[0]);
                rank_c_i[1] = inverse_select(i+1, cs[1]);
                if (cs[0]==cs[1]) {
                    k = 1;
                    rank_c_j[0] = rank_c_i[0]+2;
                } else {
                    k = 2;
                    if (cs[0] > cs[1]) {
                        std::swap(cs[0],cs[1]);
                        std::swap(rank_c_i[0],rank_c_i[1]);
                    }
                    rank_c_j[0] = rank_c_i[0]+1;
                    rank_c_j[1] = rank_c_i[1]+1;
                }
            } else {
                k = 0;
                _interval_symbols(i, j, k, cs, rank_c_i, rank_c_j, 0);
            }
        }



        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += write_member(m_sigma, out, child, "sigma");
            written_bytes += m_tree.serialize(out, child, "tree");
            written_bytes += m_tree_rank.serialize(out, child, "tree_rank");
            written_bytes += m_tree_select1.serialize(out, child, "tree_select_1");
            written_bytes += m_tree_select0.serialize(out, child, "tree_select_0");
            for (size_type i=0; i < 511; ++i) {
                written_bytes += m_nodes[i].serialize(out);
            }
            out.write((char*) m_c_to_leaf, 256*sizeof(m_c_to_leaf[0]));
            written_bytes += 256*sizeof(m_c_to_leaf[0]); // add written bytes from previous loop
            out.write((char*) m_path, 256*sizeof(m_path[0]));
            written_bytes += 256*sizeof(m_path[0]); // add written bytes from previous loop
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            read_member(m_size, in);
            read_member(m_sigma, in);
            m_tree.load(in);
            m_tree_rank.load(in, &m_tree);
            m_tree_select1.load(in, &m_tree);
            m_tree_select0.load(in, &m_tree);
            for (size_type i=0; i < 511; ++i) {
                m_nodes[i].load(in);
            }
            in.read((char*) m_c_to_leaf, 256*sizeof(m_c_to_leaf[0]));
            in.read((char*) m_path, 256*sizeof(m_path[0]));
        }

#ifdef MEM_INFO
        //! Print some infos about the size of the compressed suffix tree
        void mem_info(std::string label="")const {
            if (label=="")
                label="wt_hutu";
            size_type bytes = util::get_size_in_bytes(*this);
            std::cout << "list(label = \""<<label<<"\", size = "<< bytes/(1024.0*1024.0) <<"\n,";
            m_tree.mem_info("data"); std::cout<<",";
            m_tree_rank.mem_info("rank"); std::cout<<",";
            m_tree_select1.mem_info("select 1"); std::cout<<",";
            m_tree_select0.mem_info("select 0"); std::cout << ")\n";
            // TODO: add m_nodes, m_c_to_leaf, m_path?
        }
#endif
};

typedef wt_hutu<rrr_vector<>,
        rrr_vector<>::rank_1_type,
        rrr_vector<>::select_1_type,
        rrr_vector<>::select_0_type, 0> wt_hutu_rrr;

}// end namespace sdsl

#endif // end file
