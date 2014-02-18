/* sdsl - succinct data structures library
    Copyright (C) 2011 Timo Beller, Simon Gog

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
/*! \file nn_dict_dynamic.hpp
    \brief nn_dict_dynamic.hpp contains a class for a dynamic bit vector which also supports the prev and next operations
	\author Timo Beller, Simon Gog
*/

#ifndef INCLUDED_NN_DICT_DYNAMIC
#define INCLUDED_NN_DICT_DYNAMIC

#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>

namespace sdsl
{

class nn_dict_dynamic; // forward declaration

namespace util
{
void set_zero_bits(nn_dict_dynamic& nn);
}


// possible TODO: resize(size_type size)

//! A class for a dynamic bit vector which also supports the prev and next operations
class nn_dict_dynamic
{
    public:
        typedef int_vector<64>::size_type 				size_type;
        class reference; // forward declaration of inner class

        friend class reference;
        friend void util::set_zero_bits(nn_dict_dynamic& nn);
    private:
        uint64_t m_depth;				// Depth of the tree (1 level corresonds to 0, 2 levels corresponds to 1,...)
        uint64_t m_v_begin_leaves;		// Virtual begin of leaves
        size_type m_size;
        int_vector<64> m_offset;	 	// Number of nodes to skip on each level
        int_vector<64> m_tree; 			// Tree

        void copy(const nn_dict_dynamic& nn) {
            m_depth = nn.m_depth;
            m_v_begin_leaves = nn.m_v_begin_leaves;
            m_size = nn.m_size;
            m_offset = nn.m_offset;
            m_tree = nn.m_tree;
        }

    public:

        const uint64_t& depth;

        size_type size()const {
            return m_size;
        }

        //! Constructor
        /*! \param n Number of supported bits
         */
        nn_dict_dynamic(const uint64_t n = 0):depth(m_depth) {
            m_size = n;
            if (n == 0)
                return;
            uint64_t level;				// level indicator
            uint64_t nodes = 1;			// number of nodes (=64 bit integer)
            uint64_t tmp;				// tmp-variable

            /* Calc depth and begin of leaves */
            m_depth = bits::hi(n)/6; // if, n>0 calculate  \f$ \lfloor log_64(n) \rfloor \f$
            m_v_begin_leaves = (1ULL<<(m_depth*6))/63;

            /* Calc how many nodes to skip on each level */
            m_offset = int_vector<64>(m_depth+2, 0);
            level = m_depth;
            tmp = n;
            while (level) {
                tmp = (tmp+63)/64;  // get real number of nodes, of the next higher level
                //                  <number of nodes in the full tree>  - <real number of nodes>
                m_offset[level+1] = (1ULL<<(6*level))                   - tmp;
                nodes += tmp;
                --level;
            }

            /* Calc how many nodes to skip up to each level*/
            for (level = 1; level <= m_depth; ++level) {
                m_offset[level] += m_offset[level-1];
            }

            /* Create Tree incl. leaves */
            m_tree = int_vector<64>(nodes);
        }

        //! Copy constructor
        nn_dict_dynamic(const nn_dict_dynamic& nn):depth(m_depth) {
            copy(nn);
        }

        //! move constructor
        nn_dict_dynamic(nn_dict_dynamic&& nn):depth(m_depth) {
            *this = std::move(nn);
        }

        //! Assignment operator
        nn_dict_dynamic& operator=(const nn_dict_dynamic& nn) {
            if (this != &nn) {
                copy(nn);
            }
            return *this;
        }

        //! Assignment move operator
        nn_dict_dynamic& operator=(nn_dict_dynamic&& nn) {
            if (this != &nn) {
                m_depth = std::move(nn.m_depth);
                m_v_begin_leaves = std::move(nn.m_v_begin_leaves);
                m_size = std::move(nn.m_size);
                m_offset = std::move(nn.m_offset);
                m_tree = std::move(nn.m_tree);
                // set nn to default-constructor state
                nn.m_size = 0;
                nn.m_depth = 0;
                nn.m_v_begin_leaves = 0;
            }
            return *this;
        }

        void swap(nn_dict_dynamic& nn) {
            if (this != &nn) {
                std::swap(m_depth, nn.m_depth);
                std::swap(m_v_begin_leaves, nn.m_v_begin_leaves);
                std::swap(m_size, nn.m_size);
                m_offset.swap(nn.m_offset);
                m_tree.swap(nn.m_tree);
            }
        }

        //! Access the bit at index idx
        /*! \param idx Index
         *  \par Precondition
         *    \f$ 0 \leq  idx < size() \f$
         */
        bool operator[](const size_type& idx)const {
            uint64_t node = m_tree[(m_v_begin_leaves + (idx>>6)) - m_offset[m_depth] ];
            return (node >> (idx&0x3F)) & 1;
        }

        inline reference operator[](const size_type& idx) {
            return reference(this, idx);
        }


        //! Get the leftmost index \f$i\geq idx\f$ where a bit is set.
        /*! \param idx Left border of the search interval. \f$ 0\leq idx < size()\f$
         *
         *  \return If there exists a leftmost index \f$i\geq idx\f$ where a bit is set,
         *          then \f$i\f$ is returned, otherwise size().
         */
        size_type next(const size_type idx)const {
            uint64_t v_node_position; 	// virtual node position
            uint64_t node; 				// current node
            uint64_t depth = m_depth;	// current depth of node
            uint64_t position; 			// position of the 1-bit

            v_node_position = m_v_begin_leaves + (idx>>6);
            uint8_t off = idx & 0x3F; // mod 64

            // Go up until a 1-bit is found
            node = m_tree[ v_node_position-m_offset[depth] ]>>off;
            while (!node or off==64) {
                // Not in the root
                if (v_node_position) {
                    --depth;
                    --v_node_position;
                    off = (v_node_position&0x3F)+1;
                    v_node_position >>= 6;
                    node = m_tree[ v_node_position-m_offset[depth] ]>>off;
                } else {
                    return size();
                }
            }
            // Calculate the position of the 1-bit
            position = bits::lo(node)+off;

            // Go down to the leaf
            while (v_node_position < m_v_begin_leaves) {
                ++depth;
                v_node_position = (v_node_position<<6) + 1 + position;
                node = m_tree[ v_node_position-m_offset[depth] ];

                // Calculate the position of the 1-bit
                position = bits::lo(node);
            }
            return ((v_node_position - m_v_begin_leaves)<<6) + position;
        }

        //! Get the rightmost index \f$i \leq idx\f$ where a bit is set.
        /*! \param idx Right border of the search interval. \f$ 0 \leq idx < size()\f$
         *
         *  \return If there exists a rightmost index \f$i \leq idx\f$ where a bit is set,
         *          then \f$i\f$ is returned, otherwise size().
         */
        size_type prev(const size_type idx)const {
            uint64_t v_node_position; 	// virtual node position
            uint64_t node; 				// current node
            uint64_t depth = m_depth;	// current depth of node
            uint64_t position; 			// position of the 1-bit

            v_node_position = m_v_begin_leaves + (idx>>6);
            uint8_t off = idx & 0x3F; // mod 64

            // Go up until a 1-bit is found
            node = m_tree[ v_node_position-m_offset[depth] ]<<(63-off);
            while (!node or off == (uint8_t)-1) {

                // Not in the root
                if (v_node_position) {
                    --depth;
                    --v_node_position;

                    off = ((uint8_t)(v_node_position&0x3F))-1;
                    v_node_position >>= 6;

                    node = m_tree[ v_node_position-m_offset[depth] ]<<(63-off);
                } else {
                    return size();
                }
            }
            // Calculate the position of the 1-bit
            position = bits::hi(node)-(63-off);

            // Go down to the leaf
            while (v_node_position < m_v_begin_leaves) {
                ++depth;
                v_node_position = (v_node_position<<6) + 1 + position;
                node = m_tree[ v_node_position-m_offset[depth] ];

                // Calculate the position of the 1-bit
                position = bits::hi(node); //-(63-off)
            }
            return ((v_node_position - m_v_begin_leaves)<<6) + position;
        }


        //! Load the data structure
        void load(std::istream& in) {
            read_member(m_depth, in);
            read_member(m_v_begin_leaves, in);
            read_member(m_size, in);
            m_offset.load(in);
            m_tree.load(in);
        }

        //! Serialize the data structure
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_depth, out, child, "depth");
            written_bytes += write_member(m_v_begin_leaves, out, child, "v_begin_leaves");
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += m_offset.serialize(out, child, "offset");
            written_bytes += m_tree.serialize(out, child, "tree");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        class reference
        {
            private:
                nn_dict_dynamic* m_pbv;  // pointer to the bit_vector_nearest_neigbour
                size_type m_idx; 	// virtual node position
            public:
                //! Constructor
                reference(nn_dict_dynamic* pbv,
                          nn_dict_dynamic::size_type idx):m_pbv(pbv),m_idx(idx) {};

                //! Assignment operator for the proxy class
                reference& operator=(bool x) {
                    uint64_t v_node_position; 	// virtual node position
                    uint64_t r_node_position; 	// real    node position
                    uint64_t depth = m_pbv->m_depth;	// current depth of node

                    v_node_position = m_pbv->m_v_begin_leaves + (m_idx>>6);
                    uint8_t offset = m_idx & 0x3F; // pos mod 64
                    if (x) {
                        while (true) {
                            r_node_position = v_node_position - m_pbv->m_offset[depth];
                            uint64_t w = m_pbv->m_tree[r_node_position];
                            if ((w>>offset) & 1) { // if the bit was already set
                                return *this;
                            } else {
                                m_pbv->m_tree[r_node_position] |= (1ULL<<offset); // set bit
                                if (!w and depth) { // go up in the tree
                                    --depth; --v_node_position;
                                    offset = v_node_position&0x3F;
                                    v_node_position >>= 6;
                                } else {
                                    return *this;
                                }
                            }
                        }
                    } else {
                        while (true) {
                            r_node_position = v_node_position - m_pbv->m_offset[depth];
                            uint64_t w = m_pbv->m_tree[r_node_position];
                            if (!((w>>offset) & 1)) { // if the bit is already 0
                                return *this;
                            } else {
                                m_pbv->m_tree[r_node_position] &= (~(1ULL<<offset)); // unset bit
                                if (!m_pbv->m_tree[r_node_position] and depth) { // go up in the tree
                                    --depth; --v_node_position;
                                    offset = v_node_position&0x3F;
                                    v_node_position >>= 6;
                                } else {
                                    return *this;
                                }
                            }
                        }
                    }
                    return *this;
                }

                reference& operator=(const reference& x) {
                    return *this = bool(x);
                }

                //! Cast the reference to a bool
                operator bool()const {
                    uint64_t node = m_pbv->m_tree[(m_pbv->m_v_begin_leaves + (m_idx>>6)) - m_pbv->m_offset[m_pbv->m_depth] ];
                    return (node >> (m_idx&0x3F)) & 1;
                }

                bool operator==(const reference& x)const {
                    return bool(*this) == bool(x);
                }

                bool operator<(const reference& x)const {
                    return !bool(*this) and bool(x);
                }
        };

};



} // end of namespace

#endif // end file
