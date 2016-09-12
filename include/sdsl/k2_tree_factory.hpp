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
#include "k2_tree_hybrid.hpp"
#include "k2_tree_partitioned.hpp"

/**
 * this is hackyy (!) implementation of a factory, if time permits look at Object Factory in Modern C++ Design
 */

namespace sdsl {

    template<typename t_vector>
    std::shared_ptr<k2_tree> get_k2_tree(uint8_t k, bool use_counting_sort, t_vector &coords, uint64_t max_hint = 0) {
        switch (k) {
            case 2:
                return std::shared_ptr(new k2_tree<2>("", use_counting_sort, coords, max_hint));
            case 3:
                return std::shared_ptr(new k2_tree<3>("", use_counting_sort, coords, max_hint));
            case 4:
                return std::shared_ptr(new k2_tree<4>("", use_counting_sort, coords, max_hint));
            case 5:
                return std::shared_ptr(new k2_tree<5>("", use_counting_sort, coords, max_hint));
            case 6:
                return std::shared_ptr(new k2_tree<6>("", use_counting_sort, coords, max_hint));
            case 7:
                return std::shared_ptr(new k2_tree<7>("", use_counting_sort, coords, max_hint));
            case 8:
                return std::shared_ptr(new k2_tree<8>("", use_counting_sort, coords, max_hint));
            case 9:
                return std::shared_ptr(new k2_tree<9>("", use_counting_sort, coords, max_hint));
            case 10:
                return std::shared_ptr(new k2_tree<10>("", use_counting_sort, coords, max_hint));
            case 11:
                return std::shared_ptr(new k2_tree<11>("", use_counting_sort, coords, max_hint));
            case 12:
                return std::shared_ptr(new k2_tree<12>("", use_counting_sort, coords, max_hint));
            case 13:
                return std::shared_ptr(new k2_tree<13>("", use_counting_sort, coords, max_hint));
            case 14:
                return std::shared_ptr(new k2_tree<14>("", use_counting_sort, coords, max_hint));
            case 15:
                return std::shared_ptr(new k2_tree<15>("", use_counting_sort, coords, max_hint));
            case 16:
                return std::shared_ptr(new k2_tree<16>("", use_counting_sort, coords, max_hint));
            default:
                throw std::runtime_error(
                        "Currently only k = [2,16] is supported when using the factory, z order sort is only supported for k values being a power of two");
        }
    }

    template<typename t_vector>
    std::shared_ptr<k2_tree_hybrid>
    get_k2_tree_hybrid(uint8_t t_k_l_1, uint8_t t_k_l_1_size, uint8_t t_k_l_2, uint8_t t_k_leaves,
                       bool use_counting_sort, t_vector &coords, uint64_t max_hint = 0) {
        switch (t_k_l_1) {
            case 4:
                switch (t_k_l_1_size) {
                    case 2:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 2, 2, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 2, 2, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 2, 2, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 2, 4, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 2, 4, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 2, 4, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 2, 8, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 2, 8, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 2, 8, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 3:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 3, 2, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 3, 2, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 3, 2, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 3, 4, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 3, 4, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 3, 4, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 3, 8, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 3, 8, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 3, 8, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 4:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 4, 2, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 4, 2, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 4, 2, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 4, 4, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 4, 4, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 4, 4, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 4, 8, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 4, 8, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 4, 8, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 5:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 5, 2, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 5, 2, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 5, 2, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 5, 4, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 5, 4, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 5, 4, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 5, 8, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 5, 8, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 5, 8, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 6:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 6, 2, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 6, 2, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 6, 2, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 6, 4, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 6, 4, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 6, 4, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 6, 8, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 6, 8, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 6, 8, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 7:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 7, 2, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 7, 2, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 7, 2, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 7, 4, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 7, 4, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 7, 4, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 7, 8, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 7, 8, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 7, 8, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 8:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 8, 2, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 8, 2, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 8, 2, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 8, 4, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 8, 4, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 8, 4, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 8, 8, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 8, 8, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 8, 8, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 9:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 9, 2, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 9, 2, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 9, 2, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 9, 4, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 9, 4, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 9, 4, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 9, 8, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<4, 9, 8, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<4, 9, 8, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }

                    default:
                        throw std::runtime_error("Currently only t_k_l_1_size = [2,9] is supported in factory");
                }

            case 8:
                switch (t_k_l_1_size) {
                    case 2:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 2, 2, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 2, 2, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 2, 2, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 2, 4, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 2, 4, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 2, 4, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 2, 8, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 2, 8, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 2, 8, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 3:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 3, 2, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 3, 2, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 3, 2, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 3, 4, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 3, 4, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 3, 4, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 3, 8, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 3, 8, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 3, 8, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 4:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 4, 2, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 4, 2, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 4, 2, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 4, 4, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 4, 4, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 4, 4, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 4, 8, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 4, 8, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 4, 8, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 5:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 5, 2, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 5, 2, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 5, 2, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 5, 4, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 5, 4, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 5, 4, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 5, 8, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 5, 8, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 5, 8, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 6:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 6, 2, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 6, 2, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 6, 2, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 6, 4, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 6, 4, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 6, 4, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 6, 8, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 6, 8, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 6, 8, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 7:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 7, 2, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 7, 2, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 7, 2, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 7, 4, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 7, 4, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 7, 4, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 7, 8, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 7, 8, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 7, 8, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 8:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 8, 2, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 8, 2, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 8, 2, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 8, 4, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 8, 4, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 8, 4, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 8, 8, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 8, 8, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 8, 8, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 9:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 9, 2, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 9, 2, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 9, 2, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 9, 4, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 9, 4, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 9, 4, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 9, 8, 4>("", use_counting_sort, coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(
                                                new k2_tree_hybrid<8, 9, 8, 8>("", use_counting_sort, coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<8, 9, 8, 16>("", use_counting_sort,
                                                                                               coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                }

            case 16:
                switch (t_k_l_1_size) {
                    case 2:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 2, 2, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 2, 2, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 2, 2, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 2, 4, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 2, 4, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 2, 4, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 2, 8, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 2, 8, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 2, 8, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 3:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 3, 2, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 3, 2, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 3, 2, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 3, 4, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 3, 4, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 3, 4, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 3, 8, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 3, 8, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 3, 8, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 4:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 4, 2, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 4, 2, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 4, 2, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 4, 4, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 4, 4, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 4, 4, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 4, 8, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 4, 8, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 4, 8, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 5:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 5, 2, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 5, 2, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 5, 2, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 5, 4, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 5, 4, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 5, 4, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 5, 8, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 5, 8, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 5, 8, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 6:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 6, 2, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 6, 2, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 6, 2, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 6, 4, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 6, 4, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 6, 4, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 6, 8, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 6, 8, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 6, 8, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 7:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 7, 2, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 7, 2, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 7, 2, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 7, 4, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 7, 4, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 7, 4, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 7, 8, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 7, 8, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 7, 8, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }
                    case 8:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 8, 2, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 8, 2, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 8, 2, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 8, 4, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 8, 4, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 8, 4, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 8, 8, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 8, 8, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 8, 8, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4,8} is supported in factory");
                        }
                    case 9:
                        switch (t_k_l_2) {
                            case 2:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 9, 2, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 9, 2, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 9, 2, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }
                            case 4:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 9, 4, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 9, 4, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 9, 4, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }

                            case 8:
                                switch (t_k_leaves) {
                                    case 4:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 9, 8, 4>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 8:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 9, 8, 8>("", use_counting_sort,
                                                                                               coords.max_hint));
                                    case 16:
                                        return std::shared_ptr(new k2_tree_hybrid<16, 9, 8, 16>("", use_counting_sort,
                                                                                                coords.max_hint));
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                        }

                    default:
                        throw std::runtime_error(
                                "Currently only K=4, K=8 and K=16 are supported in top level, when using the factory");
                }
        }

    }

    template<typename t_vector>
    std::shared_ptr<k2_tree_partitioned>
    get_k2_tree_partitioned(uint8_t t_k_0, uint8_t t_k, bool use_counting_sort, t_vector &coords, uint64_t max_hint = 0) {
        switch (t_k_0){
            case 2:
                switch (t_k) {
                    case 2:
                        return std::shared_ptr(new k2_tree_partitioned<2, k2_tree<2>>("", use_counting_sort, coords, max_hint));
                    case 3:
                        return std::shared_ptr(new k2_tree_partitioned<2,k2_tree<3>>("", use_counting_sort, coords, max_hint));
                    case 4:
                        return std::shared_ptr(new k2_tree_partitioned<2,k2_tree<4>>("", use_counting_sort, coords, max_hint));
                    case 5:
                        return std::shared_ptr(new k2_tree_partitioned<2,k2_tree<5>>("", use_counting_sort, coords, max_hint));
                    case 6:
                        return std::shared_ptr(new k2_tree_partitioned<2,k2_tree<6>>("", use_counting_sort, coords, max_hint));
                    case 7:
                        return std::shared_ptr(new k2_tree_partitioned<2,k2_tree<7>>("", use_counting_sort, coords, max_hint));
                    case 8:
                        return std::shared_ptr(new k2_tree_partitioned<2,k2_tree<8>>("", use_counting_sort, coords, max_hint));
                    case 9:
                        return std::shared_ptr(new k2_tree_partitioned<2,k2_tree<9>>("", use_counting_sort, coords, max_hint));
                    case 10:
                        return std::shared_ptr(new k2_tree_partitioned<2,k2_tree<10>>("", use_counting_sort, coords, max_hint));
                    case 11:
                        return std::shared_ptr(new k2_tree_partitioned<2,k2_tree<11>>("", use_counting_sort, coords, max_hint));
                    case 12:
                        return std::shared_ptr(new k2_tree_partitioned<2,k2_tree<12>>("", use_counting_sort, coords, max_hint));
                    case 13:
                        return std::shared_ptr(new k2_tree_partitioned<2,k2_tree<13>>("", use_counting_sort, coords, max_hint));
                    case 14:
                        return std::shared_ptr(new k2_tree_partitioned<2,k2_tree<14>>("", use_counting_sort, coords, max_hint));
                    case 15:
                        return std::shared_ptr(new k2_tree_partitioned<2,k2_tree<15>>("", use_counting_sort, coords, max_hint));
                    case 16:
                        return std::shared_ptr(new k2_tree_partitioned<2,k2_tree<16>>("", use_counting_sort, coords, max_hint));
                    default:
                        throw std::runtime_error(
                                "Currently only k = [2,16] is supported when using the factory, z order sort is only supported for k values being a power of two");
                }

            case 4:
                switch (t_k) {
                    case 2:
                        return std::shared_ptr(new k2_tree_partitioned<4, k2_tree<2>>("", use_counting_sort, coords, max_hint));
                    case 3:
                        return std::shared_ptr(new k2_tree_partitioned<4,k2_tree<3>>("", use_counting_sort, coords, max_hint));
                    case 4:
                        return std::shared_ptr(new k2_tree_partitioned<4,k2_tree<4>>("", use_counting_sort, coords, max_hint));
                    case 5:
                        return std::shared_ptr(new k2_tree_partitioned<4,k2_tree<5>>("", use_counting_sort, coords, max_hint));
                    case 6:
                        return std::shared_ptr(new k2_tree_partitioned<4,k2_tree<6>>("", use_counting_sort, coords, max_hint));
                    case 7:
                        return std::shared_ptr(new k2_tree_partitioned<4,k2_tree<7>>("", use_counting_sort, coords, max_hint));
                    case 8:
                        return std::shared_ptr(new k2_tree_partitioned<4,k2_tree<8>>("", use_counting_sort, coords, max_hint));
                    case 9:
                        return std::shared_ptr(new k2_tree_partitioned<4,k2_tree<9>>("", use_counting_sort, coords, max_hint));
                    case 10:
                        return std::shared_ptr(new k2_tree_partitioned<4,k2_tree<10>>("", use_counting_sort, coords, max_hint));
                    case 11:
                        return std::shared_ptr(new k2_tree_partitioned<4,k2_tree<11>>("", use_counting_sort, coords, max_hint));
                    case 12:
                        return std::shared_ptr(new k2_tree_partitioned<4,k2_tree<12>>("", use_counting_sort, coords, max_hint));
                    case 13:
                        return std::shared_ptr(new k2_tree_partitioned<4,k2_tree<13>>("", use_counting_sort, coords, max_hint));
                    case 14:
                        return std::shared_ptr(new k2_tree_partitioned<4,k2_tree<14>>("", use_counting_sort, coords, max_hint));
                    case 15:
                        return std::shared_ptr(new k2_tree_partitioned<4,k2_tree<15>>("", use_counting_sort, coords, max_hint));
                    case 16:
                        return std::shared_ptr(new k2_tree_partitioned<4,k2_tree<16>>("", use_counting_sort, coords, max_hint));
                    default:
                        throw std::runtime_error(
                                "Currently only k = [2,16] is supported when using the factory, z order sort is only supported for k values being a power of two");
                }

            case 8:
                switch (t_k) {
                    case 2:
                        return std::shared_ptr(new k2_tree_partitioned<8, k2_tree<2>>("", use_counting_sort, coords, max_hint));
                    case 3:
                        return std::shared_ptr(new k2_tree_partitioned<8,k2_tree<3>>("", use_counting_sort, coords, max_hint));
                    case 4:
                        return std::shared_ptr(new k2_tree_partitioned<8,k2_tree<4>>("", use_counting_sort, coords, max_hint));
                    case 5:
                        return std::shared_ptr(new k2_tree_partitioned<8,k2_tree<5>>("", use_counting_sort, coords, max_hint));
                    case 6:
                        return std::shared_ptr(new k2_tree_partitioned<8,k2_tree<6>>("", use_counting_sort, coords, max_hint));
                    case 7:
                        return std::shared_ptr(new k2_tree_partitioned<8,k2_tree<7>>("", use_counting_sort, coords, max_hint));
                    case 8:
                        return std::shared_ptr(new k2_tree_partitioned<8,k2_tree<8>>("", use_counting_sort, coords, max_hint));
                    case 9:
                        return std::shared_ptr(new k2_tree_partitioned<8,k2_tree<9>>("", use_counting_sort, coords, max_hint));
                    case 10:
                        return std::shared_ptr(new k2_tree_partitioned<8,k2_tree<10>>("", use_counting_sort, coords, max_hint));
                    case 11:
                        return std::shared_ptr(new k2_tree_partitioned<8,k2_tree<11>>("", use_counting_sort, coords, max_hint));
                    case 12:
                        return std::shared_ptr(new k2_tree_partitioned<8,k2_tree<12>>("", use_counting_sort, coords, max_hint));
                    case 13:
                        return std::shared_ptr(new k2_tree_partitioned<8,k2_tree<13>>("", use_counting_sort, coords, max_hint));
                    case 14:
                        return std::shared_ptr(new k2_tree_partitioned<8,k2_tree<14>>("", use_counting_sort, coords, max_hint));
                    case 15:
                        return std::shared_ptr(new k2_tree_partitioned<8,k2_tree<15>>("", use_counting_sort, coords, max_hint));
                    case 16:
                        return std::shared_ptr(new k2_tree_partitioned<8,k2_tree<16>>("", use_counting_sort, coords, max_hint));
                    default:
                        throw std::runtime_error(
                                "Currently only k = [2,16] is supported when using the factory, z order sort is only supported for k values being a power of two");
                }

            case 16:
                switch (t_k) {
                    case 2:
                        return std::shared_ptr(new k2_tree_partitioned<16, k2_tree<2>>("", use_counting_sort, coords, max_hint));
                    case 3:
                        return std::shared_ptr(new k2_tree_partitioned<16,k2_tree<3>>("", use_counting_sort, coords, max_hint));
                    case 4:
                        return std::shared_ptr(new k2_tree_partitioned<16,k2_tree<4>>("", use_counting_sort, coords, max_hint));
                    case 5:
                        return std::shared_ptr(new k2_tree_partitioned<16,k2_tree<5>>("", use_counting_sort, coords, max_hint));
                    case 6:
                        return std::shared_ptr(new k2_tree_partitioned<16,k2_tree<6>>("", use_counting_sort, coords, max_hint));
                    case 7:
                        return std::shared_ptr(new k2_tree_partitioned<16,k2_tree<7>>("", use_counting_sort, coords, max_hint));
                    case 8:
                        return std::shared_ptr(new k2_tree_partitioned<16,k2_tree<8>>("", use_counting_sort, coords, max_hint));
                    case 9:
                        return std::shared_ptr(new k2_tree_partitioned<16,k2_tree<9>>("", use_counting_sort, coords, max_hint));
                    case 10:
                        return std::shared_ptr(new k2_tree_partitioned<16,k2_tree<10>>("", use_counting_sort, coords, max_hint));
                    case 11:
                        return std::shared_ptr(new k2_tree_partitioned<16,k2_tree<11>>("", use_counting_sort, coords, max_hint));
                    case 12:
                        return std::shared_ptr(new k2_tree_partitioned<16,k2_tree<12>>("", use_counting_sort, coords, max_hint));
                    case 13:
                        return std::shared_ptr(new k2_tree_partitioned<16,k2_tree<13>>("", use_counting_sort, coords, max_hint));
                    case 14:
                        return std::shared_ptr(new k2_tree_partitioned<16,k2_tree<14>>("", use_counting_sort, coords, max_hint));
                    case 15:
                        return std::shared_ptr(new k2_tree_partitioned<16,k2_tree<15>>("", use_counting_sort, coords, max_hint));
                    case 16:
                        return std::shared_ptr(new k2_tree_partitioned<16,k2_tree<16>>("", use_counting_sort, coords, max_hint));
                    default:
                        throw std::runtime_error(
                                "Currently only k = [2,16] is supported when using the factory, z order sort is only supported for k values being a power of two");
                }

                throw std::runtime_error(
                        "Currently only part-factor = {2,4,8,16} is supported when using the factory");
        }
    }

    template<typename t_vector>
    std::shared_ptr<k2_tree_partitioned>
    get_k2_tree_partitioned(uint8_t t_k_0, uint8_t t_k_l_1, uint8_t t_k_l_1_size, uint8_t t_k_l_2, uint8_t t_k_leaves, bool use_counting_sort, t_vector &coords, uint64_t max_hint = 0) {
        switch (t_k_0) {
            case 4:
                switch (t_k_l_1) {
                    case 4:
                        switch (t_k_l_1_size) {
                            case 2:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 2, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 2, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 2, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 2, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 2, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 2, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 2, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 2, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 2, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 3:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 3, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 3, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 3, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 3, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 3, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 3, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 3, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 3, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 3, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 4:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 4, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 4, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 4, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 4, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 4, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 4, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 4, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 4, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 4, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 5:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 5, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 5, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 5, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 5, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 5, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 5, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 5, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 5, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 5, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 6:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 6, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 6, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 6, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 6, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 6, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 6, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 6, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 6, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 6, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 7:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 7, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 7, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 7, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 7, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 7, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 7, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 7, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 7, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 7, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 8:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 8, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 8, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 8, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 8, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 8, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 8, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 8, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 8, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 8, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 9:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 9, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 9, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 9, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 9, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 9, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 9, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 9, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<4, 9, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<4, 9, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_1_size = [2,9] is supported in factory");
                        }

                    case 8:
                        switch (t_k_l_1_size) {
                            case 2:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 2, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 2, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 2, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 2, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 2, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 2, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 2, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 2, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 2, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 3:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 3, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 3, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 3, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 3, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 3, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 3, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 3, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 3, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 3, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 4:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 4, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 4, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 4, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 4, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 4, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 4, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 4, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 4, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 4, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 5:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 5, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 5, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 5, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 5, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 5, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 5, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 5, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 5, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 5, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 6:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 6, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 6, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 6, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 6, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 6, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 6, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 6, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 6, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 6, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 7:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 7, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 7, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 7, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 7, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 7, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 7, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 7, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 7, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 7, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 8:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 8, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 8, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 8, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 8, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 8, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 8, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 8, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 8, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 8, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 9:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 9, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 9, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 9, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 9, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 9, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 9, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 9, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<4, k2_tree_hybrid<8, 9, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<8, 9, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                        }

                    case 16:
                        switch (t_k_l_1_size) {
                            case 2:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 2, 2, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 2, 2, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 2, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 2, 4, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 2, 4, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 2, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 2, 8, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 2, 8, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 2, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 3:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 3, 2, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 3, 2, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 3, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 3, 4, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 3, 4, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 3, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 3, 8, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 3, 8, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 3, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 4:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 4, 2, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 4, 2, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 4, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 4, 4, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 4, 4, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 4, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 4, 8, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 4, 8, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 4, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 5:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 5, 2, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 5, 2, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 5, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 5, 4, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 5, 4, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 5, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 5, 8, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 5, 8, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 5, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 6:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 6, 2, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 6, 2, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 6, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 6, 4, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 6, 4, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 6, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 6, 8, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 6, 8, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 6, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 7:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 7, 2, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 7, 2, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 7, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 7, 4, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 7, 4, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 7, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 7, 8, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 7, 8, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 7, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 8:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 8, 2, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 8, 2, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 8, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 8, 4, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 8, 4, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 8, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 8, 8, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 8, 8, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 8, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4,8} is supported in factory");
                                }
                            case 9:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 9, 2, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 9, 2, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 9, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 9, 4, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 9, 4, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 9, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 9, 8, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 9, 8, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<4, k2_tree_hybrid<16, 9, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }

                            default:
                                throw std::runtime_error(
                                        "Currently only K=4, K=8 and K=16 are supported in top level, when using the factory");
                        }
                }

            case 8:
                switch (t_k_l_1) {
                    case 4:
                        switch (t_k_l_1_size) {
                            case 2:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 2, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 2, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 2, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 2, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 2, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 2, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 2, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 2, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 2, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 3:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 3, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 3, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 3, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 3, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 3, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 3, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 3, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 3, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 3, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 4:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 4, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 4, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 4, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 4, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 4, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 4, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 4, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 4, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 4, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 5:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 5, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 5, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 5, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 5, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 5, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 5, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 5, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 5, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 5, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 6:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 6, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 6, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 6, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 6, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 6, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 6, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 6, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 6, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 6, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 7:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 7, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 7, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 7, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 7, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 7, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 7, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 7, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 7, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 7, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 8:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 8, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 8, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 8, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 8, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 8, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 8, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 8, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 8, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 8, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 9:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 9, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 9, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 9, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 9, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 9, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 9, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 9, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<4, 9, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<4, 9, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_1_size = [2,9] is supported in factory");
                        }

                    case 8:
                        switch (t_k_l_1_size) {
                            case 2:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 2, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 2, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 2, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 2, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 2, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 2, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 2, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 2, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 2, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 3:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 3, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 3, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 3, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 3, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 3, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 3, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 3, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 3, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 3, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 4:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 4, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 4, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 4, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 4, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 4, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 4, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 4, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 4, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 4, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 5:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 5, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 5, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 5, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 5, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 5, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 5, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 5, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 5, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 5, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 6:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 6, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 6, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 6, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 6, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 6, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 6, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 6, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 6, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 6, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 7:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 7, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 7, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 7, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 7, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 7, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 7, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 7, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 7, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 7, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 8:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 8, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 8, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 8, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 8, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 8, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 8, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 8, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 8, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 8, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 9:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 9, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 9, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 9, 2, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 9, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 9, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 9, 4, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 9, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<8, k2_tree_hybrid<8, 9, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<8, 9, 8, 16>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                        }

                    case 16:
                        switch (t_k_l_1_size) {
                            case 2:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 2, 2, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 2, 2, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 2, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 2, 4, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 2, 4, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 2, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 2, 8, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 2, 8, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 2, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 3:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 3, 2, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 3, 2, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 3, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 3, 4, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 3, 4, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 3, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 3, 8, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 3, 8, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 3, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 4:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 4, 2, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 4, 2, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 4, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 4, 4, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 4, 4, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 4, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 4, 8, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 4, 8, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 4, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 5:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 5, 2, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 5, 2, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 5, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 5, 4, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 5, 4, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 5, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 5, 8, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 5, 8, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 5, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 6:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 6, 2, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 6, 2, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 6, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 6, 4, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 6, 4, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 6, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 6, 8, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 6, 8, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 6, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 7:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 7, 2, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 7, 2, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 7, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 7, 4, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 7, 4, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 7, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 7, 8, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 7, 8, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 7, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 8:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 8, 2, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 8, 2, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 8, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 8, 4, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 8, 4, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 8, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 8, 8, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 8, 8, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 8, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4,8} is supported in factory");
                                }
                            case 9:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 9, 2, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 9, 2, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 9, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 9, 4, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 9, 4, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 9, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 9, 8, 4>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 9, 8, 8>>("", use_counting_sort,
                                                                                                                               coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<8, k2_tree_hybrid<16, 9, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }

                            default:
                                throw std::runtime_error(
                                        "Currently only K=4, K=8 and K=16 are supported in top level, when using the factory");
                        }
                }

            case 16:
                switch (t_k_l_1) {
                    case 4:
                        switch (t_k_l_1_size) {
                            case 2:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 2, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 2, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 2, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 2, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 2, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 2, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 2, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 2, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 2, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 3:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 3, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 3, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 3, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 3, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 3, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 3, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 3, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 3, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 3, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 4:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 4, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 4, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 4, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 4, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 4, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 4, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 4, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 4, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 4, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 5:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 5, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 5, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 5, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 5, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 5, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 5, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 5, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 5, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 5, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 6:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 6, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 6, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 6, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 6, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 6, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 6, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 6, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 6, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 6, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 7:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 7, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 7, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 7, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 7, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 7, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 7, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 7, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 7, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 7, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 8:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 8, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 8, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 8, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 8, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 8, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 8, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 8, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 8, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 8, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 9:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 9, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 9, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 9, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 9, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 9, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 9, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 9, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<4, 9, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<4, 9, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }

                            default:
                                throw std::runtime_error("Currently only t_k_l_1_size = [2,9] is supported in factory");
                        }

                    case 8:
                        switch (t_k_l_1_size) {
                            case 2:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 2, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 2, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 2, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 2, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 2, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 2, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 2, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 2, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 2, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 3:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 3, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 3, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 3, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 3, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 3, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 3, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 3, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 3, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 3, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 4:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 4, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 4, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 4, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 4, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 4, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 4, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 4, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 4, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 4, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 5:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 5, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 5, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 5, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 5, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 5, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 5, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 5, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 5, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 5, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 6:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 6, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 6, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 6, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 6, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 6, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 6, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 6, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 6, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 6, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 7:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 7, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 7, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 7, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 7, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 7, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 7, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 7, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 7, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 7, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 8:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 8, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 8, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 8, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 8, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 8, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 8, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 8, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 8, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 8, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 9:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 9, 2, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 9, 2, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 9, 2, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 9, 4, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 9, 4, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 9, 4, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 9, 8, 4>>("", use_counting_sort, coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(
                                                        new k2_tree_partitioned<16, k2_tree_hybrid<8, 9, 8, 8>>("", use_counting_sort, coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<8, 9, 8, 16>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                        }

                    case 16:
                        switch (t_k_l_1_size) {
                            case 2:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 2, 2, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 2, 2, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 2, 2, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 2, 4, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 2, 4, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 2, 4, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 2, 8, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 2, 8, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 2, 8, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 3:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 3, 2, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 3, 2, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 3, 2, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 3, 4, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 3, 4, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 3, 4, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 3, 8, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 3, 8, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 3, 8, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 4:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 4, 2, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 4, 2, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 4, 2, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 4, 4, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 4, 4, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 4, 4, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 4, 8, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 4, 8, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 4, 8, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 5:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 5, 2, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 5, 2, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 5, 2, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 5, 4, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 5, 4, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 5, 4, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 5, 8, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 5, 8, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 5, 8, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 6:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 6, 2, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 6, 2, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 6, 2, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 6, 4, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 6, 4, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 6, 4, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 6, 8, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 6, 8, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 6, 8, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 7:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 7, 2, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 7, 2, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 7, 2, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 7, 4, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 7, 4, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 7, 4, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 7, 8, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 7, 8, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 7, 8, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }
                            case 8:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 8, 2, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 8, 2, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 8, 2, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 8, 4, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 8, 4, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 8, 4, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 8, 8, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 8, 8, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 8, 8, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4,8} is supported in factory");
                                }
                            case 9:
                                switch (t_k_l_2) {
                                    case 2:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 9, 2, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 9, 2, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 9, 2, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }
                                    case 4:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 9, 4, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 9, 4, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 9, 4, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }

                                    case 8:
                                        switch (t_k_leaves) {
                                            case 4:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 9, 8, 4>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 8:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 9, 8, 8>>("", use_counting_sort,
                                                                                                                                coords.max_hint));
                                            case 16:
                                                return std::shared_ptr(new k2_tree_partitioned<16, k2_tree_hybrid<16, 9, 8, 16>>("", use_counting_sort,
                                                                                                                                 coords.max_hint));
                                        }

                                    default:
                                        throw std::runtime_error("Currently only t_k_l_2 = {2,4} is supported in factory");
                                }

                            default:
                                throw std::runtime_error(
                                        "Currently only K=4, K=8 and K=16 are supported in top level, when using the factory");
                        }
                }

            default:
                throw std::runtime_error("For partitioing only k0=4, k=8 and k0=16 are supported when using the factory with hybrid trees");
        }
    }
}
#endif

