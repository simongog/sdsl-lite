/* sdsl - succinct data structures library
    Copyright (C) 2014 Simon Gog

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
/*! \file k2_tree_factory.hpp
    \brief k2_tree.hpp contains a compact k^2-tree.
    \author Jan Broß, based on the k2 treap code of Simon Gog
*/
#ifndef INCLUDED_SDSL_K2_TREE_FACTORY
#define INCLUDED_SDSL_K2_TREE_FACTORY

//! Namespace for the succinct data structure library.
#include "k2_tree.hpp"

namespace sdsl {
    template<typename t_vector>
        std::shared_ptr<k2_tree<>> create_k_2_comp_tree(uint8_t access_shortcut_size, t_vector& coords, uint64_t max_hint) ;
    template<typename t_vector>
        std::shared_ptr<k2_tree<>> create_k_2_uncomp_tree(uint8_t access_shortcut_size, t_vector& coords, uint64_t max_hint) ;
    template<typename t_vector>
        std::shared_ptr<k2_tree<>> create_k_4_comp_tree(uint8_t access_shortcut_size, t_vector& coords, uint64_t max_hint) ;
    template<typename t_vector>
    std::shared_ptr<k2_tree<>> create_k_4_uncomp_tree(uint8_t access_shortcut_size, t_vector& coords, uint64_t max_hint) ;
    template<typename t_vector>
    std::shared_ptr<k2_tree<>> create_k_8_comp_tree(uint8_t access_shortcut_size, t_vector& coords, uint64_t max_hint) ;
    template<typename t_vector>
    std::shared_ptr<k2_tree<>> create_k_8_uncomp_tree(uint8_t access_shortcut_size, t_vector& coords, uint64_t max_hint) ;


    template <typename t_vector>
    std::shared_ptr<k2_tree> get_k2_tree(uint8_t k, uint8_t access_shortcut_size, bool compress_leaves, t_vector& coords, uint64_t max_hint = 0) {
        switch (k) {
            case 2:
                if (compress_leaves) {
                    return create_k_2_comp_tree(access_shortcut_size, coords, max_hint);
                }   else {
                    return create_k_2_uncomp_tree(access_shortcut_size, coords, max_hint);
                }
            case 4:
                if (compress_leaves) {
                    return create_k_4_comp_tree(access_shortcut_size, coords, max_hint);
                }   else {
                    return create_k_4_uncomp_tree(access_shortcut_size, coords, max_hint);
                }
            case 8:
                if (compress_leaves) {
                    return create_k_8_comp_tree(access_shortcut_size, coords, max_hint);
                }   else {
                    return create_k_8_uncomp_tree(access_shortcut_size, coords, max_hint);
                }
            default:
                throw std::runtime_error("Currently only K=2, K=4 and K=8 are supported when using the factory");
        }
    }

    /*
    template <typename t_vector>
    std::shared_ptr<k2_tree> get_k2_tree_hybrid(uint8_t t_k_l_1, uint8_t t_k_l_1_size, uint8_t t_k_l_2, uint8_t access_shortcut_size, bool compress_leaves, t_vector& coords, uint64_t max_hint = 0) {
        switch (t_k_l_1) {
            case 2:
                if (compress_leaves) {
                    return create_hybrid_k_2_comp_tree(access_shortcut_size, coords, max_hint);
                }   else {
                    return create_hybrid_k_2_uncomp_tree(access_shortcut_size, coords, max_hint);
                }
            case 4:
                if (compress_leaves) {
                    return create_hybrid_k_4_comp_tree(access_shortcut_size, coords, max_hint);
                }   else {
                    return create_hybrid_k_4_uncomp_tree(access_shortcut_size, coords, max_hint);
                }
            case 8:
                if (compress_leaves) {
                    return create_hybrid_k_8_comp_tree(access_shortcut_size, coords, max_hint);
                }   else {
                    return create_hybrid_k_8_uncomp_tree(access_shortcut_size, coords, max_hint);
                }
            default:
                throw std::runtime_error("Currently only K=2, K=4 and K=8 are supported when using the factory");
        }
    }*/

    template<typename t_vector>
    std::shared_ptr<k2_tree> create_k_2_comp_tree(uint8_t access_shortcut_size, t_vector& coords, uint64_t max_hint) {
        switch (access_shortcut_size){
            case 0: return std::shared_ptr(new k2_tree<2,bit_vector,bit_vector,true,0>("", false, coords, max_hint - 1));
            case 2: return std::shared_ptr(new k2_tree<2,bit_vector,bit_vector,true,2>("", false, coords, max_hint - 1));
            case 4: return std::shared_ptr(new k2_tree<2,bit_vector,bit_vector,true,4>("", false, coords, max_hint - 1));
            case 16: return std::shared_ptr(new k2_tree<2,bit_vector,bit_vector,true,16>("", false, coords, max_hint - 1));
        }
    }

    template<typename t_vector>
    std::shared_ptr<k2_tree> create_k_2_uncomp_tree(uint8_t access_shortcut_size, t_vector& coords, uint64_t max_hint) {
        switch (access_shortcut_size){
            case 0: return std::shared_ptr(new k2_tree<2,bit_vector,bit_vector,false,0>("", false, coords, max_hint - 1));
            case 2: return std::shared_ptr(new k2_tree<2,bit_vector,bit_vector,false,2>("", false, coords, max_hint - 1));
            case 4: return std::shared_ptr(new k2_tree<2,bit_vector,bit_vector,false,4>("", false, coords, max_hint - 1));
            case 16: return std::shared_ptr(new k2_tree<2,bit_vector,bit_vector,false,16>("", false, coords, max_hint - 1));
        }
    }

    template<typename t_vector>
    std::shared_ptr<k2_tree> create_k_4_comp_tree(uint8_t access_shortcut_size, t_vector& coords, uint64_t max_hint) {
        switch (access_shortcut_size){
            case 0: return std::shared_ptr(new k2_tree<4,bit_vector,bit_vector,true,0>("", false, coords, max_hint - 1));
            case 2: return std::shared_ptr(new k2_tree<4,bit_vector,bit_vector,true,2>("", false, coords, max_hint - 1));
            case 4: return std::shared_ptr(new k2_tree<4,bit_vector,bit_vector,true,4>("", false, coords, max_hint - 1));
            case 16: return std::shared_ptr(new k2_tree<4,bit_vector,bit_vector,true,16>("", false, coords, max_hint - 1));
        }
    }

    template<typename t_vector>
    std::shared_ptr<k2_tree> create_k_4_uncomp_tree(uint8_t access_shortcut_size, t_vector& coords, uint64_t max_hint) {
        switch (access_shortcut_size){
            case 0: return std::shared_ptr(new k2_tree<4,bit_vector,bit_vector,false,0>("", false, coords, max_hint - 1));
            case 2: return std::shared_ptr(new k2_tree<4,bit_vector,bit_vector,false,2>("", false, coords, max_hint - 1));
            case 4: return std::shared_ptr(new k2_tree<4,bit_vector,bit_vector,false,4>("", false, coords, max_hint - 1));
            case 16: return std::shared_ptr(new k2_tree<4,bit_vector,bit_vector,false,16>("", false, coords, max_hint - 1));
        }
    }

    template<typename t_vector>
    std::shared_ptr<k2_tree> create_k_8_comp_tree(uint8_t access_shortcut_size, t_vector& coords, uint64_t max_hint) {
        switch (access_shortcut_size){
            case 0: return std::shared_ptr(new k2_tree<8,bit_vector,bit_vector,true,0>("", false, coords, max_hint - 1));
            case 2: return std::shared_ptr(new k2_tree<8,bit_vector,bit_vector,true,2>("", false, coords, max_hint - 1));
            case 4: return std::shared_ptr(new k2_tree<8,bit_vector,bit_vector,true,4>("", false, coords, max_hint - 1));
            case 16: return std::shared_ptr(new k2_tree<8,bit_vector,bit_vector,true,16>("", false, coords, max_hint - 1));
        }
    }

    template<typename t_vector>
    std::shared_ptr<k2_tree> create_k_8_uncomp_tree(uint8_t access_shortcut_size, t_vector& coords, uint64_t max_hint) {
        switch (access_shortcut_size){
            case 0: return std::shared_ptr(new k2_tree<8,bit_vector,bit_vector,false,0>("", false, coords, max_hint - 1));
            case 2: return std::shared_ptr(new k2_tree<8,bit_vector,bit_vector,false,2>("", false, coords, max_hint - 1));
            case 4: return std::shared_ptr(new k2_tree<8,bit_vector,bit_vector,false,4>("", false, coords, max_hint - 1));
            case 16: return std::shared_ptr(new k2_tree<8,bit_vector,bit_vector,false,16>("", false, coords, max_hint - 1));
        }
    }


}
#endif

