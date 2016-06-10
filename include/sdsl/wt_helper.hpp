#ifndef INCLUDED_SDSL_WT_HELPER
#define INCLUDED_SDSL_WT_HELPER

#include "int_vector.hpp"
#include <algorithm>
#include <limits>
#include <deque>
#include <queue>
#include <vector>
#include <utility>
#include <array>

namespace sdsl
{

typedef std::array<int_vector<>::size_type, 2> range_type;
typedef std::vector<range_type>         range_vec_type;

//! Empty range check
/*! \param r Range to check
 *  \returns True if the range is empty, false otherwise.
 */
bool empty(const range_type& r);

//! Size of a range
/*! \param r Range to check
 *  \returns True if the range is empty, false otherwise.
 */
int_vector<>::size_type size(const range_type& r);

//! Count for each character the number of occurrences in rac[0..size-1]
/*!
 * \param C An array of size 256, which contains for each character the number of occurrences in rac[0..size-1]
 */
template<class t_file_buffer,class t_rac>
void calculate_character_occurences(t_file_buffer& text, const int_vector_size_type size, t_rac& C)
{
    C = t_rac();
    if (text.size() < size) {
        throw std::logic_error("calculate_character_occurrences: stream size is smaller than size!");
        return;
    }
    for (int_vector_size_type i=0; i < size; ++i) {
        uint64_t c = text[i];
        if (c >= C.size()) { C.resize(c+1, 0); }
        ++C[c];
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

    enum :uint64_t {undef = 0xFFFFFFFFFFFFFFFFULL}; // max uint64_t value

    pc_node(uint64_t freq=0, uint64_t sym=0, uint64_t parent=undef,
            uint64_t child_left=undef, uint64_t child_right=undef);

    pc_node& operator=(const pc_node& v);
};

template<class t_tree_strat_fat>
struct _node {
    using node_type = typename t_tree_strat_fat::node_type;
    typedef uint64_t size_type;
    uint64_t bv_pos      = 0;             // pointer into the bit_vector, which represents the wavelet tree
    uint64_t bv_pos_rank = 0;             // pre-calculated rank for the prefix up to but not including bv_pos
    node_type parent      = t_tree_strat_fat::undef; // pointer to the parent
    node_type child[2]    = {t_tree_strat_fat::undef,t_tree_strat_fat::undef}; // pointer to the children

    _node(uint64_t bv_pos=0, uint64_t bv_pos_rank=0, node_type parent=t_tree_strat_fat::undef,
          node_type child_left=t_tree_strat_fat::undef, node_type child_right=t_tree_strat_fat::undef):
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

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
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

// TODO: version of _byte_tree for lex_ordered tree shapes
//       m_c_to_leaf can be compressed and
//       m_path is only needed for sigma chars

// Strategy class for tree representation of a WT
template<bool t_dfs_shape, class t_wt>
struct _byte_tree {
    using alphabet_category = byte_alphabet_tag;
    using value_type = uint8_t;
    using node_type = uint16_t; // node is represented by index in m_nodes
    using data_node = _node<_byte_tree>;
    enum :uint16_t {undef       = 0xFFFF}; // max uint16_t value
    enum :uint32_t {fixed_sigma = 256};
    enum :uint8_t {int_width   = 8};      // width of the input integers



    std::vector<data_node> m_nodes;              // nodes for the prefix code tree structure
    node_type          m_c_to_leaf[fixed_sigma]; // map symbol c to a leaf in the tree structure
    // if m_c_to_leaf[c] == undef the char does
    // not exists in the text
    uint64_t           m_path[fixed_sigma];      // path information for each char; the bits at position
    // 0..55 hold path information; bits 56..63 the length
    // of the path in binary representation

    void copy(const _byte_tree& bt) {
        m_nodes = bt.m_nodes;
        for (uint32_t i=0; i<fixed_sigma; ++i)
            m_c_to_leaf[i] = bt.m_c_to_leaf[i];
        for (uint32_t i=0; i<fixed_sigma; ++i)
            m_path[i] = bt.m_path[i];
    }

    _byte_tree() {}

    _byte_tree(const std::vector<pc_node>& temp_nodes, uint64_t& bv_size, const t_wt*) {
        m_nodes.resize(temp_nodes.size());
        m_nodes[0] = temp_nodes.back(); // insert root at index 0
        bv_size = 0;
        size_t node_cnt = 1;
        node_type last_parent = undef;
        std::deque<node_type> q;
        q.push_back(0);
        while (!q.empty()) {
            node_type idx;
            if (!t_dfs_shape) {
                idx = q.front(); q.pop_front();
            } else {
                idx = q.back(); q.pop_back();
            }
            // frq_sum is store in bv_pos value
            uint64_t frq = m_nodes[idx].bv_pos;
            m_nodes[idx].bv_pos = bv_size;
            if (m_nodes[idx].child[0] != undef) // if node is not a leaf
                bv_size += frq;               // add frequency
            if (idx > 0) { // node is not the root
                if (last_parent != m_nodes[idx].parent)
                    m_nodes[m_nodes[idx].parent].child[0] = idx;
                else
                    m_nodes[m_nodes[idx].parent].child[1] = idx;
                last_parent = m_nodes[idx].parent;
            }
            if (m_nodes[idx].child[0] != undef) { // if node is not a leaf
                for (uint32_t k=0; k<2; ++k) {       // add children to tree
                    m_nodes[node_cnt] = temp_nodes[ m_nodes[idx].child[k] ];
                    m_nodes[node_cnt].parent = idx;
                    q.push_back(node_cnt);
                    m_nodes[idx].child[k] = node_cnt++;
                }
            }
        }
        // initialize m_c_to_leaf
        for (uint32_t i=0; i<fixed_sigma; ++i)
            m_c_to_leaf[i] = undef; // if c is not in the alphabet m_c_to_leaf[c] = undef
        for (node_type v=0; v < m_nodes.size(); ++v) {
            if (m_nodes[v].child[0] == undef)               // if node is a leaf
                m_c_to_leaf[(uint8_t)m_nodes[v].bv_pos_rank] = v; // calculate value
        }
        // initialize path information
        // Note: In the case of a bfs search order,
        // we can classify nodes as right child and left child with an easy criterion:
        //   node is a left child, if node%2==1
        //   node is a right child, if node%2==0
        for (uint32_t c=0, prev_c=0; c<fixed_sigma; ++c) {
            if (m_c_to_leaf[c] != undef) { // if char exists in the alphabet
                node_type v = m_c_to_leaf[c];
                uint64_t pw = 0; // path
                uint64_t pl = 0; // path len
                while (v != root()) {   // while node is not the root
                    pw <<= 1;
                    if (m_nodes[m_nodes[v].parent].child[1] == v) // if the node is a right child
                        pw |= 1ULL;
                    ++pl;
                    v = m_nodes[v].parent; // go up the tree
                }
                if (pl > 56) {
                    throw std::logic_error("Code depth greater than 56!!!");
                }
                m_path[c] = pw | (pl << 56);
                prev_c = c;
            } else {
                uint64_t pl = 0; // len is  0, good for special case in rank
                m_path[c] = prev_c | (pl << 56);
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

    _byte_tree(const _byte_tree& bt) {
        copy(bt);
    }

    void swap(_byte_tree& bt) {
        std::swap(m_nodes, bt.m_nodes);
        for (uint32_t i=0; i<fixed_sigma; ++i) {
            std::swap(m_c_to_leaf[i], bt.m_c_to_leaf[i]);
            std::swap(m_path[i], bt.m_path[i]);
        }
    }

    _byte_tree& operator=(const _byte_tree& bt) {
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
        written_bytes += write_member(m_nodes_size, out, child, "m_nodes.size()");
        written_bytes += serialize_vector(m_nodes, out, child, "m_nodes");
        out.write((char*) m_c_to_leaf, fixed_sigma*sizeof(m_c_to_leaf[0]));
        written_bytes += fixed_sigma*sizeof(m_c_to_leaf[0]);// bytes from previous loop
        out.write((char*) m_path, fixed_sigma*sizeof(m_path[0]));
        written_bytes += fixed_sigma*sizeof(m_path[0]);// bytes from previous loop
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    //! Loads the data structure from the given istream.
    void load(std::istream& in) {
        uint64_t m_nodes_size = 0;
        read_member(m_nodes_size, in);
        m_nodes = std::vector<data_node>(m_nodes_size);
        load_vector(m_nodes, in);
        in.read((char*) m_c_to_leaf, fixed_sigma*sizeof(m_c_to_leaf[0]));
        in.read((char*) m_path, fixed_sigma*sizeof(m_path[0]));
    }

    //! Get corresponding leaf for symbol c.
    inline node_type c_to_leaf(value_type c)const {
        return m_c_to_leaf[c];
    }
    //! Return the root node of the tree.
    inline static node_type root() {
        return 0;
    }

    //! Return the number of nodes in the tree.
    uint64_t size() const {
        return m_nodes.size();
    }

    //! Return the parent node of v.
    inline node_type parent(node_type v)const {
        return m_nodes[v].parent;
    }
    //! Return left (i=0) or right (i=1) child node of v.
    inline node_type child(node_type v, uint8_t i)const {
        return m_nodes[v].child[i];
    }

    //! Return if v is a leaf node.
    inline bool is_leaf(node_type v)const {
        return m_nodes[v].child[0] == undef;
    }

    //! Return size of an inner node
    inline uint64_t size(node_type v) const {
        auto next_v = t_dfs_shape ? m_nodes[v].child[0] : v+1;
        return bv_pos(next_v) - bv_pos(v);
    }

    //! Return the path as left/right bit sequence in a uint64_t
    inline uint64_t bit_path(value_type c)const {
        return m_path[c];
    }

    //! Return the start of the node in the WT's bit vector
    inline uint64_t bv_pos(node_type v)const {
        return m_nodes[v].bv_pos;
    }

    //! Returns for node v the rank of 1's up to bv_pos(v)
    inline uint64_t bv_pos_rank(node_type v)const {
        return m_nodes[v].bv_pos_rank;
    }

    //! Return if the node is a valid node
    inline bool is_valid(node_type v)const {
        return v != undef;
    }

    //! Return symbol c or the next larger symbol in the wt
    inline std::pair<bool,value_type> symbol_gte(value_type c) const {
        for (uint32_t i=c; i<fixed_sigma; i++) {
            if (m_c_to_leaf[i]!=undef) {
                return {true,i};
            }
        }
        return {false,0};
    }

    //! Return symbol c or the next smaller symbol in the wt
    inline std::pair<bool,value_type> symbol_lte(value_type c) const {
        for (uint32_t i=c; i>0; i--) {
            if (m_c_to_leaf[i]!=undef) {
                return {true,i};
            }
        }
        if (m_c_to_leaf[0]!=undef)
            return {true,0};
        return {false,0};
    }
};

// Strategy class for tree representation of a WT
template<bool t_dfs_shape=false>
struct byte_tree {
    template<class t_wt>
    using type = _byte_tree<t_dfs_shape, t_wt>;
};

// Strategy class for tree representation of a WT
template<bool t_dfs_shape, class t_wt>
struct _int_tree {
    using alphabet_category = int_alphabet_tag;
    using value_type = uint64_t;
    using node_type = uint64_t; // node is represented by index in m_nodes
    using data_node = _node<_int_tree>;
    enum :uint64_t {undef = 0xFFFFFFFFFFFFFFFFULL}; // max uint64_t value
    enum :uint8_t {int_width = 0};                 // width of the input integers is variable



    std::vector<data_node> m_nodes;     // nodes for the prefix code tree structure
    std::vector<node_type> m_c_to_leaf; // map symbol c to a leaf in the tree structure
    // if m_c_to_leaf[c] == undef the char does
    // not exists in the text
    std::vector<uint64_t>  m_path;      // path information for each char; the bits at position
    // 0..55 hold path information; bits 56..63 the length
    // of the path in binary representation

    void copy(const _int_tree& bt) {
        m_nodes     = bt.m_nodes;
        m_c_to_leaf = bt.m_c_to_leaf;
        m_path      = bt.m_path;
    }

    _int_tree() {}

    _int_tree(const std::vector<pc_node>& temp_nodes, uint64_t& bv_size, const t_wt*) {
        m_nodes.resize(temp_nodes.size());
        m_nodes[0] = temp_nodes.back(); // insert root at index 0
        bv_size = 0;
        size_t node_cnt = 1;
        node_type last_parent = undef;
        std::deque<node_type> q;
        q.push_back(0);
        uint64_t max_c = 0;
        while (!q.empty()) {
            node_type idx;
            if (!t_dfs_shape) {
                idx = q.front(); q.pop_front();
            } else {
                idx = q.back(); q.pop_back();
            }
            // frq_sum is store in bv_pos value
            uint64_t frq = m_nodes[idx].bv_pos;
            m_nodes[idx].bv_pos = bv_size;
            if (m_nodes[idx].child[0] != undef) { // if node is not a leaf
                bv_size += frq;               // add frequency
            } else if (max_c < m_nodes[idx].bv_pos_rank) { // node is leaf and contains large symbol
                max_c = m_nodes[idx].bv_pos_rank;
            }
            if (idx > 0) { // node is not the root
                if (last_parent != m_nodes[idx].parent)
                    m_nodes[m_nodes[idx].parent].child[0] = idx;
                else
                    m_nodes[m_nodes[idx].parent].child[1] = idx;
                last_parent = m_nodes[idx].parent;
            }
            if (m_nodes[idx].child[0] != undef) { // if node is not a leaf
                for (uint32_t k=0; k<2; ++k) {       // add children to tree
                    m_nodes[node_cnt] = temp_nodes[ m_nodes[idx].child[k] ];
                    m_nodes[node_cnt].parent = idx;
                    q.push_back(node_cnt);
                    m_nodes[idx].child[k] = node_cnt++;
                }
            }
        }
        // initialize m_c_to_leaf
        // if c is not in the alphabet m_c_to_leaf[c] = undef
        m_c_to_leaf.resize(max_c+1, undef);
        for (node_type v=0; v < m_nodes.size(); ++v) {
            if (m_nodes[v].child[0] == undef) {              // if node is a leaf
                uint64_t c = m_nodes[v].bv_pos_rank;
                m_c_to_leaf[c] = v; // calculate value
                if (c > max_c) max_c = c;
            }
        }
        m_path = std::vector<uint64_t>(m_c_to_leaf.size(), 0);
        // initialize path information
        // Note: In the case of a bfs search order,
        // we can classify nodes as right child and left child with an easy criterion:
        //   node is a left child, if node%2==1
        //   node is a right child, if node%2==0
        for (value_type c=0, prev_c=0; c < m_c_to_leaf.size(); ++c) {
            if (m_c_to_leaf[c] != undef) { // if char exists in the alphabet
                node_type v = m_c_to_leaf[c];
                uint64_t w = 0; // path
                uint64_t l = 0; // path len
                while (v != root()) {   // while node is not the root
                    w <<= 1;
                    if (m_nodes[m_nodes[v].parent].child[1] == v) // if the node is a right child
                        w |= 1ULL;
                    ++l;
                    v = m_nodes[v].parent; // go up the tree
                }
                if (l > 56) {
                    throw std::logic_error("Code depth greater than 56!!!");
                }
                m_path[c] = w | (l << 56);
                prev_c = c;
            } else {
                uint64_t pl = 0; // len is  0, good for special case in rank
                m_path[c] = prev_c | (pl << 56);
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

    _int_tree(const _int_tree& bt) {
        copy(bt);
    }

    void swap(_int_tree& bt) {
        std::swap(m_nodes, bt.m_nodes);
        std::swap(m_c_to_leaf, bt.m_c_to_leaf);
        std::swap(m_path, bt.m_path);
    }

    _int_tree& operator=(const _int_tree& bt) {
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
        written_bytes += write_member(m_nodes_size, out, child, "m_nodes.size()");
        written_bytes += serialize_vector(m_nodes, out, child, "m_nodes");
        uint64_t m_c_to_leaf_size = m_c_to_leaf.size();
        written_bytes += write_member(m_c_to_leaf_size, out, child, "m_c_to_leaf.size()");
        written_bytes += serialize_vector(m_c_to_leaf, out, child, "m_c_to_leaf");
        uint64_t m_path_size = m_path.size();
        written_bytes += write_member(m_path_size, out, child, "m_path.size()");
        written_bytes += serialize_vector(m_path, out, child, "m_path");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    //! Loads the data structure from the given istream.
    void load(std::istream& in) {
        uint64_t m_nodes_size = 0;
        read_member(m_nodes_size, in);
        m_nodes = std::vector<data_node>(m_nodes_size);
        load_vector(m_nodes, in);
        uint64_t m_c_to_leaf_size = 0;
        read_member(m_c_to_leaf_size, in);
        m_c_to_leaf = std::vector<node_type>(m_c_to_leaf_size);
        load_vector(m_c_to_leaf, in);
        uint64_t m_path_size = 0;
        read_member(m_path_size, in);
        m_path = std::vector<uint64_t>(m_path_size);
        load_vector(m_path, in);
    }

    //! Get corresponding leaf for symbol c.
    inline node_type c_to_leaf(value_type c)const {
        if (c >= m_c_to_leaf.size())
            return undef;
        else
            return m_c_to_leaf[c];
    }
    //! Return the root node of the tree.
    inline static node_type root() {
        return 0;
    }

    //! Return the number of nodes in the tree.
    uint64_t size() const {
        return m_nodes.size();
    }

    //! Return the parent node of v.
    inline node_type parent(node_type v)const {
        return m_nodes[v].parent;
    }
    //! Return left (i=0) or right (i=1) child node of v.
    inline node_type child(node_type v, uint8_t i)const {
        return m_nodes[v].child[i];
    }

    //! Return if v is a leaf node.
    inline bool is_leaf(node_type v)const {
        return m_nodes[v].child[0] == undef;
    }

    //! Return size of an inner node
    inline uint64_t size(node_type v) const {
        auto next_v = t_dfs_shape ? m_nodes[v].child[0] : v+1;
        return bv_pos(next_v) - bv_pos(v);
    }

    //! Return the path as left/right bit sequence in a uint64_t
    inline uint64_t bit_path(value_type c)const {
        if (c >= m_path.size()) {
            return m_path.size()-1;
        }
        return m_path[c];
    }

    //! Return the start of the node in the WT's bit vector
    inline uint64_t bv_pos(node_type v)const {
        return m_nodes[v].bv_pos;
    }

    //! Returns for node v the rank of 1's up to bv_pos(v)
    inline uint64_t bv_pos_rank(node_type v)const {
        return m_nodes[v].bv_pos_rank;
    }

    //! Return if the node is a valid node
    inline bool is_valid(node_type v)const {
        return v != undef;
    }

    //! Return symbol c or the next larger symbol in the wt
    inline std::pair<bool,value_type> symbol_gte(value_type c) const {
        if (c >= m_c_to_leaf.size()) {
            return {false,0};
        }
        for (value_type i=c; i<m_c_to_leaf.size(); i++) {
            if (m_c_to_leaf[i]!=undef) {
                return {true,i};
            }
        }
        return {false,0};
    }

    //! Return symbol c or the next smaller symbol in the wt
    inline std::pair<bool,value_type> symbol_lte(value_type c) const {
        if (c >= m_c_to_leaf.size()) {
            // return the largest symbol
            c = m_c_to_leaf.size()-1;
        }
        for (value_type i=c; i>0; i--) {
            if (m_c_to_leaf[i]!=undef) {
                return {true,i};
            }
        }
        if (m_c_to_leaf[0]!=undef)
            return {true,0};
        return {false,0};
    }

};

// Strategy class for tree representation of a WT
template<bool t_dfs_shape=false>
struct int_tree {
    template<class t_wt>
    using type = _int_tree<t_dfs_shape, t_wt>;
};

template<typename t_bv>
class node_bv_container
{
    public:
        typedef typename t_bv::value_type value_type;
        typedef typename t_bv::size_type size_type;
        typedef typename t_bv::difference_type difference_type;
        typedef typename t_bv::const_iterator iterator;
    private:
        iterator m_begin, m_end;
    public:
        node_bv_container(iterator b, iterator e) : m_begin(b), m_end(e) {}
        value_type operator[](size_type i) const { return *(m_begin+i); }
        size_type size() const { return m_end - m_begin; }
        iterator begin() const {
            return m_begin;
        }
        iterator end() const {
            return m_end;
        }
};

template<typename t_bv>
class node_seq_container
{
    public:
        typedef typename t_bv::value_type value_type;
        typedef typename t_bv::size_type size_type;
        typedef typename t_bv::difference_type difference_type;
        typedef typename t_bv::const_iterator iterator;
    private:
        iterator m_begin, m_end;
    public:
        node_seq_container(iterator b, iterator e) : m_begin(b), m_end(e) {}
        value_type operator[](size_type i) const { return *(m_begin+i); }
        size_type size() const { return m_end - m_begin; }
        iterator begin() const {
            return m_begin;
        }
        iterator end() const {
            return m_end;
        }
};


} // end namespace sdsl
#endif
