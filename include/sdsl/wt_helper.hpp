#ifndef INCLUDED_SDSL_WT_HELPER
#define INCLUDED_SDSL_WT_HELPER

#include "int_vector.hpp"

namespace sdsl
{

const uint16_t _undef_node = 65535;

//! Count for each character in [0..255] the number of occurrences in rac[0..size-1]
/*!
 * \return C An array of size 256, which contains for each character the number of occurrences in rac[0..size-1]
 */
void calculate_character_occurences(int_vector_file_buffer<8>& text, const int_vector_size_type size, int_vector_size_type* C);


template<class size_type, class sigma_type>
void calculate_effective_alphabet_size(const size_type* C, sigma_type& sigma)
{
    sigma = 0; // initialize with 0
    for (size_type i=0; i < 256; ++i) // for all possible symbols
        if (C[i] > 0)				  // increase sigma, if it
            ++sigma;				 // exists in the text
}

template<class size_type>
struct _node {
    size_type tree_pos;      // pointer into the bit_vector, which represents the wavelet tree
    size_type tree_pos_rank; // pre-calculated rank for the prefix up to but not including tree_pos
    uint16_t  parent;        // pointer to the parent
    uint16_t  child[2];      // pointer to the children

    _node(size_type tree_pos=0, size_type tree_pos_rank=0, uint16_t parent=_undef_node,
          uint16_t child_left=_undef_node, uint16_t child_right=_undef_node):
        tree_pos(tree_pos), tree_pos_rank(tree_pos_rank), parent(parent) {
        child[0] = child_left;
        child[1] = child_right;
    }

    _node& operator=(const _node& v) {
        if (this != &v) {
            tree_pos      = v.tree_pos;
            tree_pos_rank = v.tree_pos_rank;
            parent        = v.parent;
            child[0]      = v.child[0];
            child[1]      = v.child[1];
        }
        return *this;
    }

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
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




} // end namespace sdsl
#endif
