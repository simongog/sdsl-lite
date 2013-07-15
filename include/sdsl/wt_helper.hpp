#ifndef INCLUDED_SDSL_WT_HELPER
#define INCLUDED_SDSL_WT_HELPER

#include "int_vector.hpp"
#include <algorithm>
#include <limits>
#include <deque>
#include <queue>

namespace sdsl
{


//! Count for each character the number of occurrences in rac[0..size-1]
/*!
 * \param C An array of size 256, which contains for each character the number of occurrences in rac[0..size-1]
 */
template<class t_file_buffer,class t_rac>
void calculate_character_occurences(t_file_buffer& text, const int_vector_size_type size, t_rac& C)
{
    C = t_rac(256, 0);
    text.reset();
    if (text.int_vector_size < size) {
        throw std::logic_error("calculate_character_occurrences: stream size is smaller than size!");
        return;
    }
    for (int_vector_size_type i=0, r_sum=0, r = text.load_next_block(); r_sum < size;) {
        if (r_sum + r > size) {  // read not more than size chars in the next loop
            r = size-r_sum;
        }
        for (; i < r_sum+r; ++i) {
            uint64_t c = text[i-r_sum];
            if (c >= C.size()) { C.resize(c+1, 0); }
            ++C[c];
        }
        r_sum += r; r = text.load_next_block();
    }
}


template<class t_rac, class sigma_type>
void calculate_effective_alphabet_size(const t_rac& C, sigma_type& sigma)
{
    sigma = std::count_if(begin(C),end(C),[](decltype(*begin(C)) &x) {
        return x > 0;
    });
}

struct pc_node {
    uint64_t  freq;     // frequency of symbol sym
    uint64_t  sym;      // symbol
    uint64_t  parent;   // pointer to the parent
    uint64_t  child[2]; // pointer to the children

    static const uint64_t undef = 0xFFFFFFFFFFFFFFFF; // max uint64_t value

    pc_node(uint64_t freq=0, uint64_t sym=0, uint64_t parent=undef,
            uint64_t child_left=undef, uint64_t child_right=undef);

    pc_node& operator=(const pc_node& v);
};

// Strategy class for tree representation of a WT
struct byte_tree {
    using alphabet_category = byte_alphabet_tag;
    using value_type = uint8_t;
    using node_type = uint16_t; // node is represented by index in m_nodes
    static const node_type undef = 0xFFFF; // max uint16_t value

    struct _node {
        uint64_t bv_pos      = 0;      // pointer into the bit_vector, which represents the wavelet tree
        uint64_t bv_pos_rank = 0;      // pre-calculated rank for the prefix up to but not including bv_pos
        uint16_t parent      = undef;      // pointer to the parent
        uint16_t child[2]    = {undef,undef};  // pointer to the children

        _node(uint64_t bv_pos=0, uint64_t bv_pos_rank=0, uint16_t parent=undef,
              uint16_t child_left=undef, uint16_t child_right=undef):
            bv_pos(bv_pos), bv_pos_rank(bv_pos_rank), parent(parent) {
            child[0] = child_left;
            child[1] = child_right;
        }

        _node& operator=(const _node& v) {
            if (this != &v) {
                bv_pos      = v.bv_pos;
                bv_pos_rank = v.bv_pos_rank;
                parent        = v.parent;
                child[0]      = v.child[0];
                child[1]      = v.child[1];
            }
            return *this;
        }

        _node& operator=(const pc_node& v) {
            bv_pos      = v.freq;
            bv_pos_rank = v.sym;
            parent        = v.parent;
            child[0]      = v.child[0];
            child[1]      = v.child[1];
            return *this;
        }


        uint64_t serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* st_child = structure_tree::add_child(v, name, util::class_name(*this));
            uint64_t written_bytes = 0;
            written_bytes += write_member(bv_pos, out);
            written_bytes += write_member(bv_pos_rank, out);
            written_bytes += write_member(parent, out);
            out.write((char*)child, 2*sizeof(child[0]));
            written_bytes += 2*sizeof(child[0]);
            structure_tree::add_size(st_child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            read_member(bv_pos, in);
            read_member(bv_pos_rank, in);
            read_member(parent, in);
            in.read((char*) child, 2*sizeof(child[0]));
        }
    };

    std::vector<_node> m_nodes;          // nodes for the prefix code tree structure
    node_type          m_c_to_leaf[256]; // map symbol c to a leaf in the tree structure
    // if m_c_to_leaf[c] == undef the char does
    // not exists in the text
    uint64_t           m_path[256];      // path information for each char; the bits at position
    // 0..55 hold path information; bits 56..63 the length
    // of the path in binary representation

    void copy(const byte_tree& bt) {
        m_nodes = bt.m_nodes;
        for (int32_t i=0; i<256; ++i)
            m_c_to_leaf[i] = bt.m_c_to_leaf[i];
        for (int32_t i=0; i<256; ++i)
            m_path[i] = bt.m_path[i];
    }

    byte_tree() {}

    template<class size_type>
    byte_tree(const std::vector<pc_node>& temp_nodes, bool dfs_shape, size_type& tree_size) {
        m_nodes.resize(temp_nodes.size());
        m_nodes[0] = temp_nodes.back(); // insert root at index 0
        tree_size = 0;
        size_t node_cnt = 1;
        uint16_t last_parent = undef;
        std::deque<size_type> q;
        q.push_back(0);
        while (!q.empty()) {
            size_type idx;
            if (!dfs_shape) {
                idx = q.front(); q.pop_front();
            } else {
                idx = q.back(); q.pop_back();
            }
            // frq_sum is store in bv_pos value
            uint64_t frq = m_nodes[idx].bv_pos;
            m_nodes[idx].bv_pos = tree_size;
            if (m_nodes[idx].child[0] != undef) // if node is not a leaf
                tree_size += frq;               // add frequency
            if (idx > 0) { // node is not the root
                if (last_parent != m_nodes[idx].parent)
                    m_nodes[m_nodes[idx].parent].child[0] = idx;
                else
                    m_nodes[m_nodes[idx].parent].child[1] = idx;
                last_parent = m_nodes[idx].parent;
            }
            if (m_nodes[idx].child[0] != undef) { // if node is not a leaf
                for (size_type k=0; k<2; ++k) {       // add children to tree
                    m_nodes[node_cnt] = temp_nodes[ m_nodes[idx].child[k] ];
                    m_nodes[node_cnt].parent = idx;
                    q.push_back(node_cnt);
                    m_nodes[idx].child[k] = node_cnt++;
                }
            }
        }
        // initialize m_c_to_leaf
        for (size_type i=0; i<256; ++i)
            m_c_to_leaf[i] = undef; // if c is not in the alphabet m_c_to_leaf[c] = undef
        for (size_type i=0; i < m_nodes.size(); ++i) {
            if (m_nodes[i].child[0] == undef)               // if node is a leaf
                m_c_to_leaf[(uint8_t)m_nodes[i].bv_pos_rank] = i; // calculate value
        }
        // initialize path information
        // Note: In the case of a bfs search order,
        // we can classify nodes as right child and left child with an easy criterion:
        //   node is a left child, if node%2==1
        //   node is a right child, if node%2==0
        for (size_type c=0; c<256; ++c) {
            if (m_c_to_leaf[c] != undef) { // if char exists in the alphabet
                size_type node = m_c_to_leaf[c];
                uint64_t w = 0; // path
                uint64_t l = 0; // path len
                while (node != 0) { // while node is not the root
                    w <<= 1;
                    if (m_nodes[m_nodes[node].parent].child[1] == node) // if the node is a right child
                        w |= 1ULL;
                    ++l;
                    node = m_nodes[node].parent; // go up the tree
                }
                if (l > 56) {
                    throw std::logic_error("Code depth greater than 56!!!");
                }
                m_path[c] = w | (l << 56);
            } else {
                m_path[c] = 0;// i.e. len is  0, good for special case in rank
            }
        }
    }

    template<class t_rank_type>
    void init_node_ranks(const t_rank_type& rank) {
        for (uint64_t i=0; i<m_nodes.size(); ++i) {
            if (m_nodes[i].child[0] != undef)  // if node is not a leaf
                m_nodes[i].bv_pos_rank = rank.rank(m_nodes[i].bv_pos);
        }
    }

    byte_tree(const byte_tree& bt) {
        copy(bt);
    }

    void swap(byte_tree& bt) {
        std::swap(m_nodes, bt.m_nodes);
        for (uint32_t i=0; i<256; ++i)
            std::swap(m_c_to_leaf[i], bt.m_c_to_leaf[i]);
        for (uint32_t i=0; i<256; ++i)
            std::swap(m_path[i], bt.m_path[i]);
    }

    byte_tree& operator=(const byte_tree& bt) {
        if (this != &bt) {
            copy(bt);
        }
        return *this;
    }

    //! Serializes the data structure into the given ostream
    uint64_t serialize(std::ostream& out, structure_tree_node* v=nullptr,
                       std::string name="") const {
        structure_tree_node* child = structure_tree::add_child(
                                         v, name, util::class_name(*this));
        uint64_t written_bytes = 0;
        uint64_t m_nodes_size = m_nodes.size();
        write_member(m_nodes_size, out, child, "m_nodes.size()");
        serialize_vector(m_nodes, out, child, "m_nodes");
        out.write((char*) m_c_to_leaf, 256*sizeof(m_c_to_leaf[0]));
        written_bytes += 256*sizeof(m_c_to_leaf[0]);// bytes from previous loop
        out.write((char*) m_path, 256*sizeof(m_path[0]));
        written_bytes += 256*sizeof(m_path[0]);// bytes from previous loop
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    //! Loads the data structure from the given istream.
    void load(std::istream& in) {
        uint64_t m_nodes_size = 0;
        read_member(m_nodes_size, in);
        m_nodes = std::vector<_node>(m_nodes_size);
        load_vector(m_nodes, in);
        in.read((char*) m_c_to_leaf, 256*sizeof(m_c_to_leaf[0]));
        in.read((char*) m_path, 256*sizeof(m_path[0]));
    }

    //! Get corresponding leaf for symbol c.
    inline node_type c_to_leaf(value_type c)const {
        return m_c_to_leaf[c];
    }
    //! Return the root node of the tree.
    inline node_type root()const {
        return 0;
    }
    //! Return the parent node of v.
    inline node_type parent(node_type v)const {
        return m_nodes[v].parent;
    }
    //! Return left (i=0) or right (i=1) child node of v.
    inline node_type child(node_type v, uint8_t i)const {
        return m_nodes[v].child[i];
    }

    inline uint64_t bit_path(value_type c)const {
        return m_path[c];
    }

    inline uint64_t bv_pos(node_type v)const {
        return m_nodes[v].bv_pos;
    }

    inline uint64_t bv_pos_rank(node_type v)const {
        return m_nodes[v].bv_pos_rank;
    }
};

// Strategy class for tree representation of a WT
struct int_tree {
    using alphabet_category = int_alphabet_tag;
    using value_type = uint64_t;
};



} // end namespace sdsl
#endif
