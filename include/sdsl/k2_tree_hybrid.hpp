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
/*! \file hybrid_k2_tree.hpp
    \brief hybrid_k2_tree.hpp contains a compact hybrid k^2-tree.
    \author Jan Broß, based on the k2 treap code of Simon Gog
*/
#ifndef INCLUDED_SDSL_HYBRID_K2_TREE
#define INCLUDED_SDSL_HYBRID_K2_TREE

#include "vectors.hpp"
#include "bits.hpp"
#include "k2_tree_helper.hpp"
#include "k2_tree_base.hpp"
#include "k2_tree_hybrid_compressed.hpp"
#include "k2_tree_hash_table.hpp"
#include "k2_tree_vocabulary.h"
#include "../../external/dacs/include/dacs.h"
#include "k2_tree_compressor.hpp"
#include <tuple>
#include <algorithm>
#include <stxxl/vector>
#include <climits>
#include <vector>
#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/rank_support_v.hpp>
#include <gtest/gtest_prod.h>

//! Namespace for the succinct data structure library.
namespace sdsl {

/*! A hybrid k2 tree implementation
 *  \par References
 *       [1] Nieves R. Brisaboa, Susana Ladra, and Gonzalo Navarro:
 *           Compact representation of Web graphs with extended functionality.
 *           Inf. Syst. 39 (January 2014), 152-174. DOI=http://dx.doi.org/10.1016/j.is.2013.08.003
 */
    template<uint8_t t_k_l_1,
            uint8_t t_k_l_1_size,
            uint8_t t_k_l_2,
            uint8_t t_k_leaves,
            typename t_lev=bit_vector,
            typename t_leaf=bit_vector,
            typename t_rank=typename t_lev::rank_1_type>
    class k2_tree_hybrid : public k2_tree_base<t_lev, t_leaf, t_rank> {
        static_assert(t_k_l_1 > 1, "t_k has to be larger than 1.");
        static_assert(t_k_l_1 <= 16, "t_k has to be smaller than 17.");
        static_assert(t_k_l_2 > 1, "t_k has to be larger than 1.");
        static_assert(t_k_l_2 <= 16, "t_k has to be smaller than 17.");
        static_assert(t_k_leaves > 1, "t_k has to be larger than 1.");
        static_assert(t_k_leaves <= 16, "t_k has to be smaller than 17.");

    public:

        typedef stxxl::VECTOR_GENERATOR<std::pair<uint32_t, uint32_t>>::result stxxl_32bit_pair_vector;
        typedef stxxl::VECTOR_GENERATOR<std::pair<uint64_t, uint64_t>>::result stxxl_64bit_pair_vector;
        using k2_tree_base<t_lev, t_leaf, t_rank>::operator=;
        using k2_tree_base<t_lev, t_leaf, t_rank>::operator==;

        k2_tree_hybrid() = default;

        /*k2_tree_hybrid(const k2_tree_hybrid &tr) {
            *this = tr;
        }

        k2_tree_hybrid(k2_tree_hybrid &&tr) {
            *this = std::move(tr);
        }*/

        template<typename t_vector>
        k2_tree_hybrid(std::string temp_file_prefix, bool use_counting_sort, uint access_shortcut_size, t_vector &v,
                       uint64_t max_hint = 0) {
            this->m_access_shortcut_size = access_shortcut_size;

            this->m_tree_height = get_tree_height(v, max_hint);

            if (v.size() > 0) {
                if (use_counting_sort) {
                    k2_tree_base<t_lev, t_leaf, t_rank>::template construct_counting_sort(v, temp_file_prefix);
                    //construct_bottom_up(v, temp_file_prefix);
                } else {
                    construct(v, temp_file_prefix);
                }
            }
        }

        k2_tree_hybrid(int_vector_buffer<> &buf_x,
                       int_vector_buffer<> &buf_y, bool use_counting_sort = false, uint access_shortcut_size = 0,
                       uint64_t max_hint = 0) {
            using namespace k2_treap_ns;
            typedef int_vector_buffer<> *t_buf_p;
            std::vector<t_buf_p> bufs = {&buf_x, &buf_y};

            this->m_access_shortcut_size = access_shortcut_size;

            uint64_t max;
            if (max_hint != 0) {
                max = max_hint;//temporarily set
            } else {
                auto max_element = [](int_vector_buffer<> &buf) {
                    uint64_t max_val = 0;
                    for (auto val : buf) {
                        max_val = std::max((uint64_t) val, max_val);
                    }
                    return max_val;
                };

                auto max_buf_element = [&]() {
                    uint64_t max_v = 0;
                    for (auto buf : bufs) {
                        uint64_t _max_v = max_element(*buf);
                        max_v = std::max(max_v, _max_v);
                    }
                    return max_v;
                };

                max = max_buf_element();
            };


            uint8_t res = 1;
            this->m_max_element = t_k_l_1;
            while (this->m_max_element < max) {
                uint8_t k;
                if (res < t_k_l_1_size) {
                    k = t_k_l_1;
                } else if ((uint) ceil((float) max / this->m_max_element) <= t_k_leaves) {
                    k = t_k_leaves;
                } else {
                    k = t_k_l_2;
                }

                this->m_max_element = this->m_max_element * k;
                ++res;
            }

            this->m_tree_height = res;

            if (res == 65) {
                throw std::logic_error("Maximal element of input is too big.");
            }

            if (this->m_max_element <= std::numeric_limits<uint32_t>::max()) {
                auto v = k2_tree_base<t_lev, t_leaf, t_rank>::template read<uint32_t, uint32_t>(bufs);
                if (use_counting_sort) {
                    k2_tree_base<t_lev, t_leaf, t_rank>::template construct_counting_sort(v, buf_x.filename());
                } else {
                    construct(v, buf_x.filename());
                }

            } else {
                auto v = k2_tree_base<t_lev, t_leaf, t_rank>::template read<uint64_t, uint64_t>(bufs);
                if (use_counting_sort) {
                    k2_tree_base<t_lev, t_leaf, t_rank>::template construct_counting_sort(v, buf_x.filename());
                } else {
                    construct(v, buf_x.filename());
                }
            }

            if (this->m_access_shortcut_size > 0) {
                //construct_access_shortcut();
            }
        }

        inline uint8_t get_k(uint8_t level) const {
            if (level < t_k_l_1_size) {
                return t_k_l_1;
            } else if (level < this->m_tree_height - 1) {
                return t_k_l_2;
            } else {
                return t_k_leaves;
            }
        }

        /**
         * Checks whether link from p = link.first to q = link.second is present i.e. matrix entry a_pq = 1
         */
        template<typename t_x, typename t_y>
        bool check_link(std::pair<t_x, t_y> link) const {

            //Patological case happening e.g. when using k2part
            if (this->m_leaves.size() == 0) {
                return false;
            }

            return check_link_internal(0, this->m_max_element, link.first, link.second, 0);
        }

        std::shared_ptr<k2_tree_hybrid_compressed<t_k_l_1,t_k_l_1_size,t_k_l_2,t_k_leaves,t_lev,t_leaf,t_rank>> compress_leaves() const {
            std::shared_ptr<k2_tree_hybrid_compressed<t_k_l_1,t_k_l_1_size,t_k_l_2,t_k_leaves,t_lev,t_leaf,t_rank>> t;

            FreqVoc(*this, [&](const HashTable &table,
                               std::shared_ptr<Vocabulary> voc) {
                t = compress_leaves(table, voc);
            });
            return t;
        }

        /**
        * Iterates over the words in the leaf level.
        *
        * @param fun Pointer to function, functor or lambda expecting a pointer to
        * each word.
        */
        template<typename Function>
        void words(Function fun) const {
            size_t cnt = words_count();
            uint size = word_size();

            size_t bit = 0;
            for (size_t i = 0; i < cnt; ++i) {
                uchar *word = new uchar[size];
                std::fill(word, word + size, 0);

                uint k_leaf_squared = get_k(this->m_tree_height-1) * get_k(this->m_tree_height-1);
                for (uint j = 0; j < k_leaf_squared ; ++j, ++bit) {
                    if (this->m_leaves[bit])
                        word[j / kUcharBits] |= (uchar) (1 << (j % kUcharBits));
                }
                fun(word);
                delete[] word;
            }
        }

        std::shared_ptr<k2_tree_hybrid_compressed<t_k_l_1,t_k_l_1_size,t_k_l_2,t_k_leaves,t_lev,t_leaf,t_rank>> compress_leaves(
                const HashTable &table,
                std::shared_ptr<Vocabulary> voc) const {
            size_t cnt = words_count();
            uint size = word_size();
            uint *codewords;
            try {
                codewords = new uint[cnt];
            } catch (std::bad_alloc ba) {
                std::cerr << "[HybridK2Tree::CompressLeaves] Error: " << ba.what() << "\n";
                exit(1);
            }

            size_t i = 0;

            words([&] (const uchar *word) {
                size_t addr;
                if (!table.search(word, size, &addr)) {
                    std::cerr << "[HybridK2Tree::CompressLeaves] Error: Word not found\n";
                    exit(1);
                }
                codewords[i++] = table[addr].codeword;
            });

            FTRep *compressL;
            try {
                // TODO Port to 64-bits
                compressL = createFT(codewords, cnt);
            } catch (...) {
                std::cerr << "[HybridK2Tree::CompressLeaves] Error: Could not create DAC\n";
                exit(1);
            }

            delete[] codewords;

            return std::shared_ptr<k2_tree_hybrid_compressed<t_k_l_1,t_k_l_1_size,t_k_l_2,t_k_leaves,t_lev,t_leaf,t_rank>>(
                    new k2_tree_hybrid_compressed<t_k_l_1, t_k_l_1_size, t_k_l_2, t_k_leaves, t_lev, t_leaf, t_rank>(
                            this->m_levels, this->m_levels_rank, compressL, voc,
                            this->m_tree_height, this->m_max_element, this->m_size, this->m_access_shortcut_size,
                            this->m_access_shortcut, this->m_access_shortcut_rank_01_support, this->m_access_shortcut_select_1_support)
            );
        }

        /**
        * Returns the number of words of \f$k_leaves^2\f$ bits in the leaf level.
        *
        * @return Number of words.
        */
        size_t words_count() const {
            return this->m_leaves.size() / get_k(this->m_tree_height-1) / get_k(this->m_tree_height-1);
        }

        /**
        * Return the number of bytes necessary to store a word.
        *
        * @return Size of a word.
        */
        uint word_size() const {
            return div_ceil((uint) get_k(this->m_tree_height-1) * get_k(this->m_tree_height-1), kUcharBits);
        }

    private:

        /**
        * Calculates the smalles integer gretear or equal to x/y
        */
        template<typename T>
        inline T div_ceil(T x, T y) const {
            static_assert(std::is_integral<T>::value, "Parameter is not integral type");
            return (x % y) ? x / y + 1 : x / y;
        }

        /**
         * Checks wether the edge p-->q exists recursively
         * @param level
         *  current level, initialy 0
         * @param p
         *  source_node
         * @param q
         *  target_node
         * @param index
         *  contains the index of the first child of the previous node, initially set to 0
         * @return
         */
        template<typename t_x, typename t_y>
        bool check_link_internal(uint8_t level, uint64_t n, t_x p, t_y q, int64_t index) const {
            using namespace k2_treap_ns;

            const uint8_t k = get_k(level);

            uint64_t current_submatrix_size = n / k;
            int64_t y = index + k * (p / current_submatrix_size) + (q / current_submatrix_size);

            if (this->is_leaf_level(level)) {
                return this->m_leaves[y];
            } else if (this->m_levels[level][y]) {
                return check_link_internal(level + 1, current_submatrix_size, p % current_submatrix_size,
                                           q % current_submatrix_size, this->get_child_index(0, y, level));
            } else {
                return false;
            }
        }

        /**
         * Constructs the tree corresponding to the points in the links vector by partitioning the input multiple times
         * @param links
         * @param temp_file_prefix
         */
        template<typename t_vector>
        void construct(t_vector &links, std::string temp_file_prefix = "") {
            using namespace k2_treap_ns;
            typedef decltype(links[0].first) t_x;
            typedef decltype(links[0].second) t_y;
            using t_e = std::pair<t_x, t_y>;

            this->m_size = links.size();

            if (this->m_tree_height == 0) {//might occur in k2part
                return;
            }

            std::string id_part = util::to_string(util::pid())
                                  + "_" + util::to_string(util::id());

            {
                std::vector<int_vector_buffer<1>> level_buffers = this->create_level_buffers(temp_file_prefix, id_part);

                auto end = std::end(links);
                //uint64_t last_level_bits = 0;
                //uint64_t level_bits = 0;

                //recursively partition that stuff
                uint64_t submatrix_size = this->m_max_element;
                for (int l = this->m_tree_height; l + 1 > 0; --l) {

                    //std::cout << "Processing Level " << l << std::endl;
                    //level_bits = 0;

                    const uint8_t k = get_k(this->m_tree_height - l);
                    uint64_t current_submatrix_size = submatrix_size / k;

                    auto sp = std::begin(links);
                    for (auto ep = sp; ep != end;) {

                        //Iterator which only returns the nodes within a certain subtree
                        ep = std::find_if(sp, end, [&submatrix_size, &k, &sp, &l](const t_e &e) {
                            auto x1 = std::get<0>(*sp);
                            auto y1 = std::get<1>(*sp);
                            auto x2 = std::get<0>(e);
                            auto y2 = std::get<1>(e);

                            bool in_sub_tree = (x1 / submatrix_size) != (x2 / submatrix_size)
                                               or (y1 / submatrix_size) != (y2 / submatrix_size);

                            return in_sub_tree;
                        });

                        if (l > 0) {
                            auto _sp = sp;

                            for (uint8_t i = 0; i < k; ++i) {
                                auto _ep = ep;
                                if (i + 1 < k) {  //partition t_k -1 times vertically (1 in the case of k=2)
                                    _ep = std::partition(_sp, _ep,
                                                         [& current_submatrix_size, &k, &i, &l](const t_e &e) {
                                                             return (std::get<0>(e) / current_submatrix_size) % k <= i;
                                                         });
                                }
                                auto __sp = _sp;

                                for (uint8_t j = 0;
                                     j < k; ++j) { //partition the t_k vertical partitions t_k -1 times horizontally
                                    auto __ep = _ep;
                                    if (j + 1 < k) {
                                        __ep = std::partition(__sp, _ep,
                                                              [& current_submatrix_size, &k, &j, &l](const t_e &e) {
                                                                  return (std::get<1>(e) / current_submatrix_size) %
                                                                         k <= j;
                                                              });
                                    }
                                    bool not_empty = __ep > __sp;
                                    level_buffers[this->m_tree_height - l].push_back(not_empty);
                                    //level_bits++;
                                    __sp = __ep;
                                }
                                _sp = _ep;
                            }
                        }
                        if (ep != end) {
                            //++ep;
                            sp = ep;
                        }
                    }
                    //last_level_bits = level_bits;
                    submatrix_size = current_submatrix_size;
                }
            }

            this->load_vectors_from_file(temp_file_prefix, id_part);

        }

        template<typename t_tv>
        uint8_t get_tree_height(const t_tv &v, uint64_t max_seed = 0) {
            using namespace k2_treap_ns;
            using t_e = typename t_tv::value_type;

            if (v.size() == 0) {
                return 0;
            }

            uint64_t max;
            if (max_seed != 0) {
                max = max_seed;
            } else {
                auto tupmax = [](t_e a) {
                    return std::max(a.first, a.second);
                };
                auto max_it = std::max_element(std::begin(v), std::end(v), [&](t_e a, t_e b) {
                    return tupmax(a) < tupmax(b);
                });

                max = tupmax(*max_it);
            }

            uint8_t res = 1;
            this->m_max_element = t_k_l_1;
            while (this->m_max_element < max) {
                uint8_t k;
                if (res < t_k_l_1_size) {
                    k = t_k_l_1;
                } else if ((uint) ceil((float) max / this->m_max_element) <= t_k_leaves) {
                    k = t_k_leaves;
                } else {
                    k = t_k_l_2;
                }

                this->m_max_element = this->m_max_element * k;
                ++res;
            }

            if (res <= t_k_l_1_size) {
                std::cerr
                        << "The tree height is smaller than t_k_l_1_size, ignoring t_k_l_2 and t_k_leaves"
                        << std::endl;
            } else if (res == t_k_l_1_size + 1) {
                std::cerr << "The tree equals t_k_l_1_size+1, only using t_k_l_1 and t_k_leaves" << std::endl;
            }

            return res;
        }
    };
}
#endif

