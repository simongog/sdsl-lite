#ifndef INCLUDED_SDSL_BASE_K2_TREE
#define INCLUDED_SDSL_BASE_K2_TREE

/*! \file k2_tree_base.hpp
    \brief k2_tree_base.hpp contains the common operations of k2trees and hybrid k2trees
    \author Jan Broß, based on the k2 treap code of Simon Gog, leaf compression is based on the libk2tree implementation which uses DACs impelemented by Brisaboa, Ladra et.al.
*/

#include "k2_tree_vocabulary.hpp"
#include "k2_tree_helper.hpp"
#include "k2_tree_dacs.hpp"
#include "k2_tree_hash_table.hpp"
#include "k2_tree_compressor.hpp"
#include "wt_huff.hpp"
#include "construct.hpp"
#include "wavelet_trees.hpp"
#include <stxxl/vector>
#include <tuple>
#include <algorithm>
#include <climits>
#include <vector>
#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/rank_support_v.hpp>
#include <gtest/gtest_prod.h>

//! Namespace for the succinct data structure library.
namespace sdsl {

/*! A base class for the k2 and hybrid k2 tree implementation
 *  \par References
 *       [1] Nieves R. Brisaboa, Susana Ladra, and Gonzalo Navarro:
 *           Compact representation of Web graphs with extended functionality.
 *           Inf. Syst. 39 (January 2014), 152-174. DOI=http://dx.doi.org/10.1016/j.is.2013.08.003
 *
 *
 */



    template<uint8_t k0,
            typename t_lev,
            typename t_leaf,
            typename t_rank>
    class k2_tree_base {

    public:
        typedef stxxl::VECTOR_GENERATOR<std::pair<uint32_t, uint32_t>>::result stxxl_32bit_pair_vector;
        typedef stxxl::VECTOR_GENERATOR<std::pair<uint64_t, uint64_t>>::result stxxl_64bit_pair_vector;

        uint64_t m_max_element = 0; //FIXME: this is an ugly hack for k2part, could also be figured out with tree height and k!
        uint8_t m_tree_height = 0;

        using node_type = k2_treap_ns::node_type;
        using point_type = k2_treap_ns::point_type;
        using t_p = k2_treap_ns::t_p;

        k2_tree_base() = default;

        k2_tree_base(const k2_tree_base &tr) {
            *this = tr;
        }

        k2_tree_base(k2_tree_base &&tr) {
            *this = std::move(tr);
        }


        k2_tree_base(uint8_t access_shortcut_size, bool dac_compress) : m_access_shortcut_size(access_shortcut_size),
                                                                        m_is_dac_comp(dac_compress) {
        }

        //virtual ~k2_tree_base() = 0;
        virtual ~k2_tree_base() {
            /*if (t_comp){
                m_comp_leaves.destroy();
                m_vocabulary.destroy();
            }*/
        }

    protected:
        std::vector<t_lev> m_levels;
        std::vector<t_rank> m_levels_rank;

        t_leaf m_leaves;
        size_type m_size = 0;
        /** BitArray containing Gog's B vector. */
        bit_vector m_access_shortcut;
        //Rank support for pattern 01 and 1
        rank_support_v<10, 2> m_access_shortcut_rank_10_support;
        bit_vector::select_1_type m_access_shortcut_select_1_support;

        DAC m_comp_leaves;
        uint8_t m_access_shortcut_size = 0;
        bool m_is_dac_comp = false;

        /** For compressed version **/
        bool m_vocabulary_is_shared = false;
        std::shared_ptr<Vocabulary> m_vocabulary;

        bool m_is_wt_comp = false;
        wt_huff<hyb_vector<>> m_leaves_wt;
        static constexpr uint8_t m_wt_word_size = 8;

        virtual uint8_t get_k(uint8_t) const = 0;

        virtual uint8_t get_tree_height(const uint64_t max) = 0;

    private:
        //sl = shortcut_level
        uint64_t m_field_size_on_sl = 0;
        uint8_t m_real_size_on_sl = 0;
        uint64_t m_submatrix_in_row_on_sl = 0;

    public:

        template<typename t_x>
        void direct_links2(t_x source_id, std::vector<t_x> &result) const {
            using namespace k2_treap_ns;
            result.clear();

            //Patological case happening e.g. when using k2part
            if (m_tree_height == 0 || source_id > m_max_element) {
                return;
            }

            //direct_links2_internal(m_max_element, 0, source_id, t_x(0), 0, result);
            if (m_is_dac_comp) {
                direct_links2_internal_queue(source_id, result,
                                             [this](int64_t pos, t_x offset, uint8_t leafK, std::vector<t_x> &result) {
                                                 check_leaf_bits_direct_comp(pos, offset, leafK, result);
                                             });
            } else if (m_is_wt_comp) {
                direct_links2_internal_queue(source_id, result,
                                             [this](int64_t pos, t_x offset, uint8_t leafK, std::vector<t_x> &result) {
                                                 check_leaf_bits_direct_wt(pos, offset, leafK, result);
                                             });
            } else {
                direct_links2_internal_queue(source_id, result,
                                             [this](int64_t pos, t_x offset, uint8_t leafK, std::vector<t_x> &result) {
                                                 check_leaf_bits_direct_uncomp(pos, offset, leafK, result);
                                             });
            }

        }

        template<typename t_x>
        void direct_links_shortcut(t_x source_id, std::vector<t_x> &result) const {
            if (m_is_dac_comp) {
                direct_links_shortcut_internal(source_id, result, [this](int64_t pos, t_x offset, uint8_t leafK,
                                                                         std::vector<t_x> &result) {
                    check_leaf_bits_direct_comp(pos, offset, leafK, result);
                });
            } else {
                direct_links_shortcut_internal(source_id, result, [this](int64_t pos, t_x offset, uint8_t leafK,
                                                                         std::vector<t_x> &result) {
                    check_leaf_bits_direct_uncomp(pos, offset, leafK, result);
                });
            }
        }

        template<typename t_x>
        void direct_links_shortcut_2(t_x source_id, std::vector<t_x> &result) const {
            if (m_is_dac_comp) {
                direct_links_shortcut_internal_2(source_id, result, [this](int64_t pos, t_x offset, uint8_t leafK,
                                                                           std::vector<t_x> &result) {
                    check_leaf_bits_direct_comp(pos, offset, leafK, result);
                });
            } else {
                direct_links_shortcut_internal_2(source_id, result, [this](int64_t pos, t_x offset, uint8_t leafK,
                                                                           std::vector<t_x> &result) {
                    check_leaf_bits_direct_uncomp(pos, offset, leafK, result);
                });
            }
        }


        template<typename t_x, typename Function>
        void direct_links_shortcut_internal(t_x source_id, std::vector<t_x> &result, Function check_leaf_bits) const {
            using namespace k2_treap_ns;
            result.clear();

            if (m_access_shortcut_size == 0) {
                throw std::runtime_error("Cannot use check_link_shortcut if m_access_shortcut_size == 0");
            }

            //Patological case happening e.g. when using k2part
            if (m_tree_height == 0 || source_id > m_max_element) {
                return;
            }

            for (uint j = 0; j < m_submatrix_in_row_on_sl; ++j) {
                t_x column_offset = j * m_field_size_on_sl;
                uint64_t z = access_shortcut_helper<k0>::corresponding_subtree(column_offset, source_id,
                                                                               m_real_size_on_sl,
                                                                               m_access_shortcut_size);
                uint64_t y = this->m_access_shortcut_select_1_support(z + 1);
                //std::cout << "current y " << y <<std::endl;
                //std::cout << "current column offset " << column_offset <<std::endl;

                //check if exists and if B_[y-1] == 0 otherwise no link
                if (!(this->m_access_shortcut[y + 1] == true)) {
                    //rank 01 pattern on B[0,p] to find out how many non-empty trees are there until p
                    //directly get corresponding data from leaf array

                    uint64_t index = this->m_access_shortcut_rank_10_support(y + 1);
                    direct_links2_internal(m_field_size_on_sl, m_access_shortcut_size,
                                           t_x(source_id % m_field_size_on_sl),
                                           column_offset, index, result, check_leaf_bits);
                }
            }
        }

        template<typename t_x, typename Function>
        void direct_links_shortcut_internal_2(t_x source_id, std::vector<t_x> &result, Function check_leaf_bits) const {
            using namespace k2_treap_ns;
            result.clear();

            if (m_access_shortcut_size == 0) {
                throw std::runtime_error("Cannot use check_link_shortcut if m_access_shortcut_size == 0");
            }

            //Patological case happening e.g. when using k2part
            if (m_tree_height == 0 || source_id > m_max_element) {
                return;
            }

            t_x column_offset = 0;
            for (uint j = 0; j < m_submatrix_in_row_on_sl / k0; ++j) {
                uint64_t z = access_shortcut_helper<k0>::corresponding_subtree(column_offset, source_id,
                                                                               m_real_size_on_sl,
                                                                               m_access_shortcut_size);
                uint64_t y = this->m_access_shortcut_select_1_support(z + 1);
                uint64_t index = this->m_access_shortcut_rank_10_support(y + 1);
                //y--;
                for (int i = 0; i < k0; ++i) {
                    //dont use select support, but directly look in access_shortcut once position has been obtained

                    //std::cout << "current y " << y <<std::endl;
                    //std::cout << "current column offset " << column_offset <<std::endl;
                    //check if exists and if B_[y-1] == 0 otherwise no link
                    if (this->m_access_shortcut[y + 1] == false) {
                        //rank 01 pattern on B[0,p] to find out how many non-empty trees are there until p
                        //directly get corresponding data from leaf array

                        direct_links2_internal(m_field_size_on_sl, m_access_shortcut_size,
                                               t_x(source_id % m_field_size_on_sl), column_offset, index, result,
                                               check_leaf_bits);
                        index++;
                        y += 2;
                    } else {
                        y += 1;
                    }

                    column_offset += m_field_size_on_sl;
                }
                z += k0 * k0;
            }
        }

        template<typename t_x>
        void inverse_links_shortcut(t_x target_id, std::vector<t_x> &result) const {
            if (m_is_dac_comp) {
                inverse_links_shortcut_internal(target_id, result, [this](int64_t pos, t_x offset, uint8_t leafK,
                                                                          std::vector<t_x> &result) {
                    check_leaf_bits_inverse_comp(pos, offset, leafK, result);
                });
            } else {
                inverse_links_shortcut_internal(target_id, result, [this](int64_t pos, t_x offset, uint8_t leafK,
                                                                          std::vector<t_x> &result) {
                    check_leaf_bits_inverse_uncomp(pos, offset, leafK, result);
                });
            }
        }

        template<typename t_x, typename Function>
        void inverse_links_shortcut_internal(t_x target_id, std::vector<t_x> &result, Function check_leaf_bits) const {
            using namespace k2_treap_ns;
            if (m_access_shortcut_size == 0) {
                throw std::runtime_error("Cannot use check_link_shortcut if m_access_shortcut_size == 0");
            }

            result.clear();

            //Patological case happening e.g. when using k2part
            if (m_tree_height == 0 || target_id > m_max_element) {
                return;
            }

            for (uint j = 0; j < m_submatrix_in_row_on_sl; ++j) {
                t_x row_offset = j * m_field_size_on_sl;
                uint64_t z = access_shortcut_helper<k0>::corresponding_subtree(target_id, row_offset, m_real_size_on_sl,
                                                                               m_access_shortcut_size);
                uint64_t y = this->m_access_shortcut_select_1_support(z + 1);
                //check if exists and if B_[y-1] == 0 otherwise no link
                if (!(this->m_access_shortcut[y + 1] == true)) {
                    //rank 01 pattern on B[0,p] to find out how many non-empty trees are there until p
                    //directly get corresponding data from leaf array

                    uint64_t index = this->m_access_shortcut_rank_10_support(y + 1);
                    inverse_links2_internal(m_field_size_on_sl, m_access_shortcut_size,
                                            t_x(target_id % m_field_size_on_sl), row_offset, index, result,
                                            check_leaf_bits);
                }
            }
        }

        template<typename t_x>
        void inverse_links2(t_x target_id, std::vector<t_x> &result) const {
            using namespace k2_treap_ns;
            result.clear();

            if (m_tree_height == 0 || target_id > m_max_element) {
                return;
            }

            //inverse_links2_internal(m_max_element, 0, source_id, t_x(0), 0, result);
            if (m_is_dac_comp) {
                inverse_links2_internal_queue(target_id, result,
                                              [this](int64_t pos, t_x offset, uint8_t leafK, std::vector<t_x> &result) {
                                                  check_leaf_bits_inverse_comp(pos, offset, leafK, result);
                                              });
            } else if (m_is_wt_comp) {
                inverse_links2_internal_queue(target_id, result,
                                              [this](int64_t pos, t_x offset, uint8_t leafK, std::vector<t_x> &result) {
                                                  check_leaf_bits_inverse_wt(pos, offset, leafK, result);
                                              });
            } else {
                inverse_links2_internal_queue(target_id, result,
                                              [this](int64_t pos, t_x offset, uint8_t leafK, std::vector<t_x> &result) {
                                                  check_leaf_bits_inverse_uncomp(pos, offset, leafK, result);
                                              });
            }
        }

        //! Move assignment operator
        virtual k2_tree_base &operator=(const k2_tree_base &&tr) {
            if (this != &tr) {
                m_tree_height = tr.m_tree_height;
                m_size = tr.m_size;
                m_max_element = tr.m_max_element;
                m_levels = std::move(tr.m_levels);
                m_levels_rank = std::move(tr.m_levels_rank);
                for (uint i = 0; i < m_levels_rank.size(); ++i) {
                    m_levels_rank[i].set_vector(&m_levels[i]);
                }
                m_leaves = std::move(tr.m_leaves);
                m_access_shortcut_size = tr.m_access_shortcut_size;
                m_is_dac_comp = tr.m_is_dac_comp;
                m_access_shortcut = std::move(tr.m_access_shortcut);
                m_access_shortcut_rank_10_support = std::move(tr.m_access_shortcut_rank_10_support);
                m_access_shortcut_select_1_support = std::move(tr.m_access_shortcut_select_1_support);
                m_comp_leaves = tr.m_comp_leaves;
                m_vocabulary_is_shared = tr.m_vocabulary_is_shared;
                m_vocabulary = tr.m_vocabulary;
                m_field_size_on_sl = tr.m_field_size_on_sl;
                m_real_size_on_sl = tr.m_real_size_on_sl;
                m_submatrix_in_row_on_sl = tr.m_submatrix_in_row_on_sl;
                m_is_wt_comp = tr.m_is_wt_comp;
                m_leaves_wt = tr.m_leaves_wt;
            }
            return *this;
        }

        //! Assignment operator
        virtual k2_tree_base &operator=(const k2_tree_base &tr) {
            if (this != &tr) {
                m_tree_height = tr.m_tree_height;
                m_size = tr.m_size;
                m_max_element = tr.m_max_element;
                m_levels = tr.m_levels;
                m_levels_rank = tr.m_levels_rank;
                for (uint64_t i = 0; i < m_levels_rank.size(); ++i) {
                    m_levels_rank[i].set_vector(&m_levels[i]);
                }
                m_leaves = tr.m_leaves;
                m_access_shortcut_size = tr.m_access_shortcut_size;
                m_is_dac_comp = tr.m_is_dac_comp;
                m_access_shortcut = tr.m_access_shortcut;
                m_access_shortcut_rank_10_support = tr.m_access_shortcut_rank_10_support;
                m_access_shortcut_rank_10_support.set_vector(&m_access_shortcut);
                m_access_shortcut_select_1_support = tr.m_access_shortcut_select_1_support;
                m_access_shortcut_select_1_support.set_vector(&m_access_shortcut);
                m_comp_leaves = tr.m_comp_leaves;
                m_vocabulary_is_shared = tr.m_vocabulary_is_shared;
                m_vocabulary = tr.m_vocabulary;
                m_field_size_on_sl = tr.m_field_size_on_sl;
                m_real_size_on_sl = tr.m_real_size_on_sl;
                m_submatrix_in_row_on_sl = tr.m_submatrix_in_row_on_sl;
                m_is_wt_comp = tr.m_is_wt_comp;
                m_leaves_wt = tr.m_leaves_wt;
            }
            return *this;
        }

        //! Equals operator
        virtual bool operator==(const k2_tree_base &tr) const {
            if (m_tree_height != tr.m_tree_height)
                return false;
            if (m_size != tr.m_size)
                return false;

            if (m_max_element != tr.m_max_element)
                return false;

            if (m_levels.size() != tr.m_levels.size()) {
                std::cout << "m_levels.size() differs" << std::endl;
                return false;
            }
            for (uint i = 0; i < m_levels.size(); ++i) {
                if (m_levels[i].size() != tr.m_levels[i].size()) {
                    std::cout << "m_levels[" << i << "].size() differs" << std::endl;
                    return false;
                }
                for (uint j = 0; j < m_levels[i].size(); ++j) {
                    if (m_levels[i][j] != tr.m_levels[i][j]) {
                        std::cout << "m_levels vectors differ at " << i << std::endl;
                        return false;
                    }
                }
            }

            if (m_is_dac_comp != tr.m_is_dac_comp) {
                std::cout << "one is compressed, the other not" << std::endl;
                return false;
            }

            if (m_is_wt_comp != tr.m_is_wt_comp) {
                std::cout << "one is compressed, the other not" << std::endl;
                return false;
            }

            if (m_tree_height > 0) {
                if (m_is_dac_comp) {
                    if (!(m_comp_leaves == tr.m_comp_leaves)) {
                        std::cout << "comp leaves differ" << std::endl;
                        return false;
                    }

                    if (m_vocabulary_is_shared != tr.m_vocabulary_is_shared){
                        std::cout << "One vocabulary shared, other not"<< std::endl;
                        return false;
                    }


                    if (!m_vocabulary.get()->operator==(*tr.m_vocabulary.get())) {
                        std::cout << "vocabulary differs" << std::endl;
                        return false;
                    }

                } else if (m_is_wt_comp != tr.m_is_wt_comp) {
                    if (m_wt_word_size != tr.m_wt_word_size) {
                        std::cout << "Wavelet Tree word size differs" << std::endl;
                        return false;
                    }

                    /*
                    if (m_leaves_wt != tr.m_leaves_wt) {
                        return false;
                    }*/
                } else {
                    if (m_leaves.size() != tr.m_leaves.size()) {
                        std::cout << "m_leaves.size() differs" << std::endl;
                        return false;
                    }
                    for (uint i = 0; i < m_leaves.size(); ++i) {
                        if (m_leaves[i] != tr.m_leaves[i]) {
                            std::cout << "m_leaves vectors differ at " << i << std::endl;
                            return false;
                        }

                    }
                }
            }

            if (m_access_shortcut_size != tr.m_access_shortcut_size) {
                std::cout << "Shortcut size differs" << std::endl;
                return false;
            }
            //don't compare other access_shortcut vetors as they have to be the same when access_shortcut_size is the same

            return true;
        }

        //! Number of points in the 2^k treap
        size_type
        size() const {
            return m_size;
        }

        //! Swap operator
        virtual void swap(k2_tree_base &tr) {
            if (this != &tr) {
                std::swap(m_tree_height, tr.m_tree_height);
                std::swap(m_size, tr.m_size);
                std::swap(m_max_element, tr.m_max_element);

                uint tr_m_levels_size = tr.m_levels.size();
                uint m_levels_size = m_levels.size();

                uint64_t biggerSize;
                if (tr_m_levels_size > m_levels_size) {
                    biggerSize = tr_m_levels_size;
                    m_levels.resize(biggerSize);
                    m_levels_rank.resize(biggerSize);
                } else {
                    biggerSize = m_levels_size;
                    tr.m_levels.resize(biggerSize);
                    tr.m_levels_rank.resize(biggerSize);
                }

                for (uint64_t j = 0; j < biggerSize; ++j) {
                    m_levels[j].swap(tr.m_levels[j]);
                }
                std::swap(tr_m_levels_size, m_levels_size);

                for (uint64_t j = 0; j < biggerSize; ++j) {
                    util::swap_support(m_levels_rank[j], tr.m_levels_rank[j], &m_levels[j], &(tr.m_levels[j]));
                }

                m_levels.resize(m_levels_size);
                tr.m_levels.resize(tr_m_levels_size);
                m_levels_rank.resize(m_levels_size);
                tr.m_levels_rank.resize(tr_m_levels_size);

                m_leaves.swap(tr.m_leaves);
                std::swap(m_access_shortcut_size, tr.m_access_shortcut_size);
                m_access_shortcut.swap(tr.m_access_shortcut);
                util::swap_support(m_access_shortcut_rank_10_support, tr.m_access_shortcut_rank_10_support,
                                   &m_access_shortcut, &tr.m_access_shortcut);
                util::swap_support(m_access_shortcut_select_1_support, tr.m_access_shortcut_select_1_support,
                                   &m_access_shortcut, &tr.m_access_shortcut);

                std::swap(m_is_dac_comp, tr.m_is_dac_comp);
                std::swap(m_vocabulary_is_shared, tr.m_vocabulary_is_shared);
                std::swap(m_vocabulary, tr.m_vocabulary);
                m_comp_leaves.swap(tr.m_comp_leaves);
                std::swap(m_is_wt_comp, tr.m_is_wt_comp);
                std::swap(m_leaves_wt, tr.m_leaves_wt);
                std::swap(m_field_size_on_sl, tr.m_field_size_on_sl);
                std::swap(m_real_size_on_sl, tr.m_real_size_on_sl);
                std::swap(m_submatrix_in_row_on_sl, tr.m_submatrix_in_row_on_sl);
            }
        }

        //! Serializes the data structure into the given ostream
        virtual size_type serialize(std::ostream &out, structure_tree_node *v = nullptr,
                                    std::string name = "") const {
            structure_tree_node *child = structure_tree::add_child(
                    v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_tree_height, out, child, "t");
            written_bytes += write_member(m_size, out, child, "s");
            written_bytes += write_member(m_max_element, out, child, "max_element");
            if (m_tree_height > 0) {
                for (int i = 0; i < m_tree_height - 1; ++i) {
                    written_bytes += m_levels[i].serialize(out, child, ("level" + std::to_string(i)));
                }
                for (int i = 0; i < m_tree_height - 1; ++i) {
                    written_bytes += m_levels_rank[i].serialize(out, child, "levels_rank");
                }

                written_bytes += write_member(m_is_dac_comp, out, child, "m_is_dac_comp");
                written_bytes += write_member(m_is_wt_comp, out, child, "m_is_wt_comp");
                if (m_is_dac_comp) {
                    written_bytes += write_member(m_vocabulary_is_shared, out, child, "m_voc_is_shared");
                    if (!m_vocabulary_is_shared){
                        written_bytes += m_vocabulary->serialize(out, child, "voc");
                    }
                    written_bytes += m_comp_leaves.serialize(out, child, "comp_leafs");
                } else if (m_is_wt_comp) {
                    written_bytes += m_leaves_wt.serialize(out, child, "wt_huff_int");
                } else {
                    written_bytes += m_leaves.serialize(out, child, "leafv");
                }
            }

            written_bytes += write_member(m_access_shortcut_size, out, child, "access_shortcut_size");
            if (m_access_shortcut_size > 0) {
                written_bytes += m_access_shortcut.serialize(out, child, "access_shortcut");
                written_bytes += m_access_shortcut_rank_10_support.serialize(out, child, "access_rank");
                written_bytes += m_access_shortcut_select_1_support.serialize(out, child, "access_select");
            }
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        virtual void load(std::istream &in) {
            read_member(m_tree_height, in);
            read_member(m_size, in);
            read_member(m_max_element, in);

            if (m_tree_height > 0) {
                m_levels.resize(m_tree_height - 1);
                for (uint64_t i = 0; i < m_levels.size(); ++i) {
                    m_levels[i].load(in);
                }

                m_levels_rank.resize(m_tree_height - 1);
                for (uint64_t i = 0; i < m_levels_rank.size(); ++i) {
                    m_levels_rank[i].load(in);
                    m_levels_rank[i].set_vector(&m_levels[i]);
                }

                read_member(m_is_dac_comp, in);
                read_member(m_is_wt_comp, in);
                if (m_is_dac_comp) {
                    read_member(m_vocabulary_is_shared, in);
                    if (!m_vocabulary_is_shared){
                        m_vocabulary = std::shared_ptr<Vocabulary>(new Vocabulary());
                        m_vocabulary->load(in);
                    }
                    m_comp_leaves.load(in);
                } else if (m_is_wt_comp) {
                    m_leaves_wt.load(in);
                } else {
                    m_leaves.load(in);
                }

            }

            read_member(m_access_shortcut_size, in);
            if (m_access_shortcut_size > 0) {
                m_access_shortcut.load(in);
                m_access_shortcut_rank_10_support.load(in);
                m_access_shortcut_rank_10_support.set_vector(&m_access_shortcut);
                m_access_shortcut_select_1_support.load(in);
                m_access_shortcut_select_1_support.set_vector(&m_access_shortcut);
            }
        }

        void words(std::vector<uchar>& result, bool append = false) const {
            if (m_tree_height == 0) {
                return;
            }


            size_t cnt = words_count();
            uint size = word_size();
            uint64_t offset = 0;
            if (!append){
                result.clear();
                result.resize(cnt*size);
            } else {
                offset = result.size();
                result.resize(result.size()+cnt*size);
            }


            uint k_leaf_squared = get_k(this->m_tree_height - 1) * get_k(this->m_tree_height - 1);

            //this can still be optimized e.g. for 64 bits
            for (uint k = 0; k < cnt; ++k) {
                for (uint i = 0; i < size-1; ++i) {
                    result[offset+k*size + i] = (this->m_leaves.get_int(k*k_leaf_squared+i*kUcharBits, kUcharBits));
                }

                if (k_leaf_squared%kUcharBits){
                    result[offset+k*size + size-1] = (this->m_leaves.get_int(k*k_leaf_squared+(size-1)*kUcharBits,k_leaf_squared%kUcharBits ));
                } else {
                    result[offset+k*size + size-1] = (this->m_leaves.get_int(k*k_leaf_squared+(size-1)*kUcharBits, kUcharBits));
                }
            }
        }
        /**
        * Returns the type as string
        */
        virtual std::string get_type_string() const = 0;

        /**
        * Returns the number of words of \f$k_leaves^2\f$ bits in the leaf level.
        *
        * @return Number of words.
        */
        virtual size_t words_count() const = 0;

        /**
        * Return the number of bytes necessary to store a word.
        *
        * @return Size of a word.
        */
        virtual uint word_size() const = 0;

        template<typename t_x, typename t_y>
        bool check_link_shortcut(std::pair<t_x, t_y> link) const {
            if (m_is_dac_comp) {
                return check_link_shortcut_internal(link, [this](int64_t pos, uint8_t leafK) {
                    return this->is_leaf_bit_set_comp(pos, leafK);
                });
            } else {
                return check_link_shortcut_internal(link, [this](int64_t pos, uint8_t) {
                    return this->is_leaf_bit_set(pos);
                });
            }
        }


        /**
        * Used for accelerating the check whether a certain link exists by skipping m_access_shortcut_size levels
        *
        * @param p Identifier of first object.
        * @param q Identifier of second object.
        *
        * @return Returns true/false depending on wehter link is present or not
        */
        template<typename t_x, typename t_y, typename Function>
        bool check_link_shortcut_internal(std::pair<t_x, t_y> link, Function check_leaf_bits) const {
            using namespace k2_treap_ns;

            //Patological case happening e.g. when using k2part
            if (this->m_tree_height == 0) {
                return false;
            }

            if (m_access_shortcut_size == 0) {
                throw std::runtime_error("Cannot use check_link_shortcut if m_access_shortcut_size == 0");
            }

            //FIXME: height if k_L tree!, it depends as we're not only targeting the last level anymore
            //FIXME: check which points are in the same tree and only fetch once
            //how to get corresponding subtree on level x of a point efficiently? (for k=2^x, interleave x-bitwise the top h bits
            //implement subtree calculation in general and for 2^x special-cases manually, think about precomp in the case of k=3

            auto p = link.first;
            auto q = link.second;

            uint64_t z = access_shortcut_helper<k0>::corresponding_subtree(q, p, m_real_size_on_sl,
                                                                           m_access_shortcut_size);
            //y = zth 1 via rank on B_
            uint64_t y = this->m_access_shortcut_select_1_support(z + 1);
            //check if exists and if B_[y-1] == 0 otherwise no link
            if (this->m_access_shortcut[y + 1] == true) {
                return false;
            }
            //rank 01 pattern on B[0,p] to find out how many non-empty trees are there until p
            //directly get corresponding data from leaf array

            uint64_t number_of_present_trees_searched_value_is_in = this->m_access_shortcut_rank_10_support(y + 1);
            uint64_t index = number_of_present_trees_searched_value_is_in * get_k(m_access_shortcut_size) *
                             get_k(m_access_shortcut_size);

            //std::cout << "For " << p << "," << q << " the index is " << index << "and relative coordinates are " << p%field_size << "," << q%field_size << std::endl;

            return check_link_internal(m_access_shortcut_size, m_field_size_on_sl, p % m_field_size_on_sl,
                                       q % m_field_size_on_sl, index, check_leaf_bits);
        }

        /**
        * Checks whether link from p = link.first to q = link.second is present i.e. matrix entry a_pq = 1
        */
        template<typename t_x, typename t_y>
        bool check_link(std::pair<t_x, t_y> link) const {

            //Patological case happening e.g. when using k2part
            if (this->m_tree_height == 0) {
                return false;
            }
            if (m_is_dac_comp) {
                return check_link_internal(0, m_max_element, link.first, link.second, 0,
                                           [this](int64_t pos, uint8_t leafK) {
                                               return this->is_leaf_bit_set_comp(pos, leafK);
                                           });
            } else if (m_is_wt_comp) {
                return check_link_internal(0, m_max_element, link.first, link.second, 0,
                                           [this](int64_t pos, uint8_t) {
                                               return this->is_leaf_bit_set_wt(pos);
                                           });
            } else {
                return check_link_internal(0, m_max_element, link.first, link.second, 0, [this](int64_t pos, uint8_t) {
                    return this->is_leaf_bit_set(pos);
                });
            }
        }

        //WARNING: only to be used from within k2_tree_partition, FIXME: encapsulate this better
        void clear_leaves() {
            m_leaves = t_leaf();
        }

        void compress_leaves(const HashTable &table, std::shared_ptr<Vocabulary> voc){
            if (m_tree_height == 0){
                return;
            }
            std::vector<uchar> leaf_words;
            words(leaf_words);
            compress_leaves(table, voc, leaf_words);
            m_is_dac_comp = true;
            m_vocabulary_is_shared = true;
        }

        void set_vocabulary(const std::shared_ptr<Vocabulary> vocabulary){
            m_vocabulary_is_shared = true;
            m_vocabulary = vocabulary;
        }

    protected:


        void perform_access_shortcut_precomputations() {
            if (m_access_shortcut_size <= this->m_tree_height - 2) {

                m_field_size_on_sl = m_max_element;//height cannot be used when using hybrid k trees --> buffer this value in this case
                for (int i = 0; i < m_access_shortcut_size; ++i) {
                    m_field_size_on_sl /= get_k(i);
                }

                m_real_size_on_sl = 0;
                //calculate real bits used by m_max_element
                while (1ULL << (m_real_size_on_sl) < m_max_element) { ++m_real_size_on_sl; }

                m_submatrix_in_row_on_sl = /*m_max_element / */precomp<k0>::exp(
                        m_access_shortcut_size);//same as max_element/field_size
            } else {
                std::cerr << "Should crash earlier ;-)" << std::endl;
            }
        }

        void compress_leaves(const HashTable &table, std::shared_ptr<Vocabulary> voc, const std::vector<uchar>& leaf_words) {
            size_t cnt = words_count();
            uint size = word_size();
            uint *codewords;
            try {
                codewords = new uint[cnt];
            } catch (std::bad_alloc ba) {
                std::cerr << "[k2_tree_base::compress_leaves] Error: " << ba.what() << "\n";
                exit(1);
            }

            size_t addr;
            for (size_t i = 0; i < cnt; ++i) {
                if (!table.search(&leaf_words[i*size], size, &addr)) {
                    std::cerr << "[k2_tree_base::compress_leaves] Error: Word not found\n";
                    exit(1);
                } else {
                    //std::cout << "Codeword: " << table[addr].codeword << std::endl;
                    codewords[i] = table[addr].codeword;
                }
            }

            m_leaves = t_leaf();

            try {
                // TODO Port to 64-bits
                /*std::cout << "CodeWords" << std::endl;
                for (int j = 0; j < cnt; ++j) {
                    std::cout << codewords[j] << "\t";
                }
                std::cout << std::endl;

                std::cout << "Count" << cnt << std::endl;*/
                m_comp_leaves = DAC(codewords, cnt);
            } catch (...) {
                std::cerr << "[k2_tree_base::compress_leaves] Error: Could not create DAC\n";
                exit(1);
            }

            delete[] codewords;
            m_vocabulary = voc;
        }

        /**
        * gets the index of the ith child of node x
        */
        inline uint64_t get_child_index(uint i, int64_t x, uint8_t level) const {
            uint rank = m_levels_rank[level](x);
            return rank * get_k(level + 1) * get_k(level + 1) + i;
        }

    public:
        void compress_leaves(uint64_t hash_size = 0) {
            if (m_tree_height == 0){
                return;
            }

            std::cout << "Compressing leaves" << std::endl;
            std::vector<uchar> leaf_words;
            words(leaf_words);

            FreqVoc(leaf_words, word_size(), words_count(), [&](const HashTable &table, std::shared_ptr<Vocabulary> voc, const std::vector<uchar>& leaf_words) {
                compress_leaves(table, voc, leaf_words);
                m_is_dac_comp = true;
            }, hash_size);
        }

        /**
         * Recursive function for getting the successors of a certain node.
         * Detailed in the "Compact representation of Web graphs with extended functionality" Paper
         * @param n
         *      current submatrix size/initially the maximal node id that is theoretically possible for a tree with k=t_k and height m_tree_height
         * @param source_id
         *      starting node for which the successors are searched
         * @param column_offset
         *      of the current submatrix
         * @param index
         *      of the upper left corner of the submatrix in the concatentation of m_levels and m_leafs vector
         * @param result
         */
        template<typename t_x, typename Function>
        void direct_links2_internal(uint64_t n, uint8_t level, t_x source_id, t_x column_offset, int64_t index,
                                    std::vector<t_x> &result, Function check_leaf_bits_direct) const {
            const uint8_t k = get_k(level);
            uint64_t submatrix_size = n / k;
            uint64_t y = index * k * k + k * (source_id / submatrix_size);

            if (this->is_leaf_level(level)) {
                check_leaf_bits_direct(y, column_offset, get_k(level), result);
            } else { //internal node
                for (uint j = 0; j < k; ++j) {
                    if (m_levels[level][y + j] == 1) {
                        direct_links2_internal(submatrix_size, level + 1, t_x(source_id % submatrix_size),
                                               t_x(column_offset + submatrix_size * j), m_levels_rank[level](y + j),
                                               result, check_leaf_bits_direct);
                    }
                }
            }
        }

        void compress_leaves_huf_wt() {
            if (m_tree_height == 0) {
                return;
            }

            m_is_wt_comp = true;

            uint64_t number_of_words = (m_leaves.size() + m_wt_word_size - 1) / m_wt_word_size;
            int_vector<m_wt_word_size> leafs(number_of_words);
            uint64_t ctr = 0;
            for (uint i = 0; i < number_of_words; ++i, ++ctr) {
                leafs[ctr] = m_leaves.get_int(i * m_wt_word_size, m_wt_word_size);
            }

            m_leaves = t_leaf();

            /*
            std::cout << "Leafs2: ";
            for (uint m = 0; m < leafs.size(); ++m) {
                std::cout << std::to_string(leafs[m]) << "\t";
            }
            std::cout << std::endl;
            */
            construct_im(m_leaves_wt, leafs);
            /*
            std::cout << "m_leaves_wt after construction: ";
            for (uint k = 0; k < m_leaves_wt.size(); ++k) {
                std::cout << std::to_string(m_leaves_wt[k]) << "\t";
            }
            std::cout << std::endl;
            */
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
        template<typename t_x, typename t_y, typename Function>
        bool check_link_internal(uint8_t level, uint64_t n, t_x p, t_y q, int64_t index, Function check_leaf) const {
            using namespace k2_treap_ns;

            const uint8_t k = get_k(level);

            uint64_t current_submatrix_size = n / k;
            int64_t y = index + k * (p / current_submatrix_size) + (q / current_submatrix_size);

            if (this->is_leaf_level(level)) {
                return check_leaf(y, k);
            } else if (this->m_levels[level][y]) {
                return check_link_internal(level + 1, current_submatrix_size, p % current_submatrix_size,
                                           q % current_submatrix_size, this->get_child_index(0, y, level), check_leaf);
            } else {
                return false;
            }
        }

        /**
         * Variant of direct_links2_internal using a queue
         * @param source_id
         * @param result
         */
        template<typename t_x, typename Function>
        void
        direct_links2_internal_queue(t_x source_id, std::vector<t_x> &result, Function check_leaf_bits_direct) const {
            using namespace k2_treap_ns;
            //n, level, source_id, column_offset, index
            std::queue<std::tuple<uint64_t, uint8_t, t_x, t_x, int64_t>> queue;
            queue.push(std::make_tuple(m_max_element, 0, source_id, (t_x) 0, 0));

            while (!queue.empty()) {
                auto current_element = queue.front();
                uint64_t n = std::get<0>(current_element);
                uint8_t level = std::get<1>(current_element);
                t_x source_id = std::get<2>(current_element);
                t_x column_offset = std::get<3>(current_element);
                int64_t index = std::get<4>(current_element);
                queue.pop();

                const uint8_t k = get_k(level);

                uint64_t submatrix_size = n / k;
                int64_t y = index * k * k + k * (source_id / submatrix_size);

                if (is_leaf_level(level)) {
                    check_leaf_bits_direct(y, column_offset, get_k(level), result);
                } else { //internal node
                    for (uint j = 0; j < k; ++j) {
                        if (m_levels[level][y + j] == 1) {
                            queue.push(std::make_tuple(submatrix_size, level + 1, t_x(source_id % submatrix_size),
                                                       t_x(column_offset + submatrix_size * j),
                                                       m_levels_rank[level](y + j)));
                        }
                    }
                }
            }
        }

        bool is_leaf_level(int level) const { return level == m_tree_height - 1; }

        /**
         * Recursive function for getting the predecessor of a certain node.
         * Detailed in the "Compact representation of Web graphs with extended functionality" Paper
         * @param n
         *      current submatrix size/initially the maximal node id that is theoretically possible for a tree with k=t_k and height m_tree_height
         * @param source_id
         *      starting node for which the predecessor are searched
         * @param column_offset
         *      of the current submatrix
         * @param index
         *      of the upper left corner of the submatrix in the concatentation of m_levels and m_leafs vector, initially 0
         * @param result
         */
        template<typename t_x, typename Function>
        void inverse_links2_internal(uint64_t n, uint8_t level, t_x source_id, t_x row_offset, int64_t index,
                                     std::vector<t_x> &result, Function check_leaf_bits_inverse) const {

            const uint8_t k = get_k(level);
            uint64_t submatrix_size = n / k;
            int64_t y = index * k * k + (source_id / submatrix_size);

            if (is_leaf_level(level)) {
                check_leaf_bits_inverse(y, row_offset, get_k(level), result);
            } else { //internal node
                for (int j = 0; j < k; ++j) {
                    if (m_levels[level][y + (j * k)]) {
                        inverse_links2_internal(submatrix_size, level + 1, t_x(source_id % submatrix_size),
                                                t_x(row_offset + submatrix_size * j), m_levels_rank[level](y + (j * k)),
                                                result, check_leaf_bits_inverse);
                    }
                }
            }
        }

        /**
         * Variant of inverse_links2_internal using a queue
         * @param source_id
         * @param result
         */
        template<typename t_x, typename Function>
        void
        inverse_links2_internal_queue(t_x source_id, std::vector<t_x> &result, Function check_leaf_bits_inverse) const {
            using namespace k2_treap_ns;
            //n, level, source_id, column_offset, index
            std::queue<std::tuple<uint64_t, uint8_t, t_x, t_x, int64_t>> queue;
            queue.push(std::make_tuple(m_max_element, 0, source_id, (t_x) 0, 0));

            while (!queue.empty()) {
                auto current_element = queue.front();
                t_x n = std::get<0>(current_element);
                uint8_t level = std::get<1>(current_element);
                t_x source_id = std::get<2>(current_element);
                t_x row_offset = std::get<3>(current_element);
                int64_t index = std::get<4>(current_element);
                queue.pop();

                const uint8_t k = get_k(level);

                uint64_t submatrix_size = n / k;
                int64_t y = index * k * k + (source_id / submatrix_size);

                if (is_leaf_level(level)) {
                    check_leaf_bits_inverse(y, row_offset, get_k(level), result);
                } else { //internal node
                    for (int j = 0; j < k; ++j) {
                        if (m_levels[level][y + (j * k)]) {
                            queue.push(std::make_tuple(submatrix_size, level + 1, t_x(source_id % submatrix_size),
                                                       t_x(row_offset + submatrix_size * j),
                                                       m_levels_rank[level](y + (j * k))));
                        }
                    }
                }
            }
        }

    protected:
        template<typename t_x=uint64_t, typename t_y=uint64_t>
        std::vector<std::pair<t_x, t_y>>
        read(std::vector<int_vector_buffer<> *> &bufs) {
            typedef std::vector<std::pair<t_x, t_y>> t_tuple_vec;
            t_tuple_vec v = t_tuple_vec(bufs[0]->size());
            for (uint64_t j = 0; j < v.size(); ++j) {
                std::get<0>(v[j]) = (*(bufs[0]))[j];
            }
            for (uint64_t j = 0; j < v.size(); ++j) {
                std::get<1>(v[j]) = (*(bufs[1]))[j];
            }

            return v;
        }

        //use only for testing purposes (remove and use mock)
        void set_height(uint height) {
            m_tree_height = height;
        }

        void load_vectors_from_file(const std::string &temp_file_prefix, const std::string &id_part) {
            {
                std::vector<bit_vector> levels;
                levels.resize(m_tree_height - 1);
                m_levels.reserve(levels.size());
                for (int i = 0; i < m_tree_height - 1; ++i) {
                    std::string levels_file = temp_file_prefix + "_level_" + std::to_string(i) + "_" + id_part
                                              + ".sdsl";
                    load_from_file(levels[i], levels_file);
                    m_levels.emplace_back(levels[i]);
                    sdsl::remove(levels_file);
                }
            }

            {
                bit_vector leafs;
                std::string leafs_file =
                        temp_file_prefix + "_level_" + std::to_string(m_tree_height - 1) + "_" + id_part
                        + ".sdsl";
                load_from_file(leafs, leafs_file);
                remove(leafs_file);

                m_leaves = t_leaf(leafs);
            }

            m_levels_rank.resize(m_levels.size());
            for (uint64_t i = 0; i < m_levels.size(); ++i) {
                util::init_support(m_levels_rank[i], &m_levels[i]);
            }
        }

        std::vector<int_vector_buffer<1>>
        create_level_buffers(const std::string temp_file_prefix, std::string &id_part) const {
            std::vector<int_vector_buffer<1>> level_buffers;
            level_buffers.reserve(m_tree_height);

            for (uint64_t i = 0; i < m_tree_height; ++i) {
                std::string levels_file = temp_file_prefix + "_level_" + std::to_string(i) + "_" + id_part
                                          + ".sdsl";
                level_buffers.emplace_back(levels_file, std::ios::out);
            }

            return level_buffers;
        }

        /**
        * Constructs the tree corresponding to the points in the links vector by using counting
        * sort with k² buckets and rearranging the input on every level of the tree
        * @param links
        * @param temp_file_prefix
        */
        template<typename t_vector>
        void construct_counting_sort(t_vector &links, std::string temp_file_prefix = "") {
            using namespace k2_treap_ns;
            typedef decltype(links[0].first) t_x;
            typedef decltype(links[0].second) t_y;
            using t_e = std::pair<t_x, t_y>;

            m_size = links.size();

            if (m_size == 0) {
                return;
            }

            //std::cout << "Using counting sort" << std::endl;

            std::string id_part = util::to_string(util::pid())
                                  + "_" + util::to_string(util::id());

            {
                std::vector<int_vector_buffer<1>> level_buffers = create_level_buffers(temp_file_prefix, id_part);
                //                  upper left          lower right                 interval in links   level
                typedef std::tuple<std::pair<t_x, t_y>, std::pair<uint64_t, uint64_t>, t_e, uint8_t> t_queue;
                std::queue<t_queue> queue;

                //partition recursively until reaching the leaves
                queue.push(t_queue(std::make_pair<t_x, t_y>(0, 0), std::make_pair(m_max_element - 1, m_max_element - 1),
                                   t_e(0, links.size()), 0));

                uint64_t number_of_bits = 0; //for speed comparison purposes of different k

                //std::vector<std::vector<double>> utilization(m_tree_height);
                //std::vector<std::vector<int>> item_count(m_tree_height);

                //std::cout << "Setring m_level_begin_idx["<<m_tree_height-1<<"] =" << 0 << std::endl;
                while (!queue.empty()) {
                    auto upper_left = std::get<0>(queue.front());
                    auto lower_right = std::get<1>(queue.front());
                    auto links_interval = std::get<2>(queue.front());
                    auto current_level = std::get<3>(queue.front());

                    const uint8_t k = get_k(current_level);

                    auto submatrix_size =
                            lower_right.first - upper_left.first + 1;
                    std::vector<t_x> intervals(k * k + 1);

                    queue.pop();


                    //do counting sort
                    auto x1 = upper_left.first;
                    auto y1 = upper_left.second;
                    auto subDivK = (submatrix_size / k);
                    for (uint64_t j = links_interval.first; j < links_interval.second; ++j) {
                        auto x = links[j].first;
                        auto y = links[j].second;
                        uint p1 = (x - x1) / subDivK;
                        uint p2 = (y - y1) / subDivK;
                        uint corresponding_matrix = p1 * k + p2;
                        intervals[corresponding_matrix +
                                  1]++;//offset corresponding matrix by one to allow for more efficient in interval comparision
                    }

                    intervals[0] = 0;

                    //append bits to level_vectors[level] based on result
                    for (uint i = 1; i < intervals.size(); ++i) {
//                        item_count[current_level].push_back(intervals[i])
//                        utilization[current_level].push_back(
//                                ((double) intervals[i]) / (submatrix_size * submatrix_size) * 100);

                        if (intervals[i] > 0) {
                            level_buffers[current_level].push_back(1);
                            //std::cout << "1";
                            number_of_bits++;
                        } else {
                            level_buffers[current_level].push_back(0);
                            //std::cout << "0";
                            number_of_bits++;
                        }
                        intervals[i] += intervals[i - 1]; //build prefix sum
                    }

                    //leaves not reached yet --> append to level_vector & reorder
                    if (submatrix_size > k) {
                        std::vector<t_x> offset(k * k);
                        offset[0] = 0;
                        for (size_t l = 1; l < offset.size(); ++l) {
                            offset[l] = intervals[l] + 1;
                        }

                        auto begin = links.begin() + links_interval.first;
                        auto it = begin;

                        //reorder links based on counting sort offsets
                        uint64_t index = 0;
                        while (it != links.begin() + links_interval.second) {
                            uint corresponding_matrix =
                                    (((*it).first - x1) / subDivK) * k + ((*it).second - y1) / subDivK;

                            if (index >= intervals[corresponding_matrix] &&
                                index < intervals[corresponding_matrix + 1]) {
                                //element is at correct position
                                offset[corresponding_matrix]++;

                                it++;
                                index++;
                            } else {
                                //swap to correct position either swapping a match or a non match back
                                //there are at most m swaps if all m elements are at the wrong position
                                //every value is at most read twice
                                std::iter_swap(it, begin + offset[corresponding_matrix] - 1);
                                offset[corresponding_matrix]++;
                            }
                        }

                        //enqueue submatrixes
                        auto new_submatrix_size = submatrix_size / k;
                        for (uint x = 0; x < k; ++x) {
                            for (uint y = 0; y < k; ++y) {
                                auto new_interval = std::make_pair(intervals[x * k + y] + links_interval.first,
                                                                   intervals[x * k + y + 1] + links_interval.first);
                                if (new_interval.first != new_interval.second) {
                                    auto new_upper_left = std::make_pair<t_x, t_y>(
                                            x * new_submatrix_size + upper_left.first,
                                            y * new_submatrix_size + upper_left.second);
                                    auto new_lower_right = std::make_pair<t_x, t_y>(
                                            (x + 1) * new_submatrix_size - 1 + upper_left.first,
                                            (y + 1) * new_submatrix_size - 1 + upper_left.second);
                                    queue.push(std::make_tuple(new_upper_left, new_lower_right, new_interval,
                                                               current_level + 1));
                                }
                            }
                        }
                    }
                }
/*
                uint avg = 0;
                for (uint m = 0; m < item_count.size(); ++m) {
                    std::cout << "Level " << m << "\n";
                    for (auto item: item_count[m]) {
                        //std::cout << item << ",";
                        avg += item;
                    }
                    std::cout << "\n";
                    std::cout << "Level Size " << item_count[m].size() << std::endl;
                    std::cout << "Average " << avg / item_count[m].size() << std::endl;
                    avg = 0;
                }
*/
                /*
                double average_utilization = 0;
                std::cout << "Utilization" << std::endl;
                for (uint m = 0; m < 4; ++m) {
                    //std::cout << "Level "<< m << "\n";
                    for (auto item: utilization[m]) {
                        std::cout << item << ",";
                        average_utilization += item;
                    }
                    std::cout << "\n";
                    //std::cout << "Average " << average_utiliztation/item_count[m].size() << std::endl;
                    average_utilization = 0;
                    //std::sort(utilization[m].begin(), utilization[m].end());

                }*/
            }

            load_vectors_from_file(temp_file_prefix, id_part);
        }

        node_type root() const {
            return node_type(0, t_p(0, 0), 0);
        }


    public:
        /**
        * Hier noch eine Idee um den k^2-tree zu beschleunigen: Um nicht erst durch h Levels zu navigieren kann man sich erst ein bit_vector B bauen, der aus 4^h Einsen und höchstens 4^h Nullen besteht. Für jeden der 4^h Teilbäume schreibt man eine Eins; für nichtleere Teilbäume zusätzlich eine Null vor der entsprechenden Eins.
        *  Also für das Beispiel in Abb.1.3 (in Jans Bericht) mit h=2:
        *
        *       0 12 345678 901 2345
        *  B = 010110111111011101111
        *  P = 3 4  5      8   1
        *                      2
        *  Zum Teilbaum der Koordinate (x,y) kommt man indem man
        *   (1) die oberen h Bits von x und y interleaved; nennen wir das z
        *   (2) die Position p der (z+1)te Eins selektieren
        *   (3) Falls p=0 oder B[p-1]=1 ist, so ist der Teilbaum leer. Andernfalls
        *        ist der Teilbaum nicht leer und durch das Ergebnis r einer
        *       rank Operation auf das Bitpattern ,01' in B[0,p]
        *       adressieren wir einen Array P, der die Präfixsummer der
        *       Subbaumgrößen enthält. Der Eintrag P[r] kann als Pointer auf die
        *       k^2 Repräsentation des nichtleeren Subtrees dienen.
        *
        *  Das ist praktikabel für h=8. Der Bitvektor würde höchstens 16kBytes
        *  benötigen und P im worst case 2MB. Da nur ein rank und ein select
        *  gemacht werden sollte das deutlich schneller sein als die 8 ranks
        *  in der vorherigen Implementierung.
         */
        virtual void construct_access_shortcut(uint8_t access_shortcut_size) {
            using namespace k2_treap_ns;

            m_access_shortcut_size = access_shortcut_size;
            //maximal size of shortcut is tree height
            if (m_access_shortcut_size > this->m_tree_height - 2) {
                std::cerr << "shortcut size must be smaller than tree height -2";
                if (m_tree_height - 2 > 0) {
                    m_access_shortcut_size = m_tree_height - 2;
                } else {
                    m_access_shortcut_size = 0;
                    return;
                }
            }

            //Use 1 to code empty trees in level height-1 and 01 to code non-empty trees, height has to be calculated with kL_ as height of hybrid tree is different
            //coresponds to the amount of non-empty trees in level h-1
            //amount of Zeros = actually existent number of trees --> (level_begin_idx[level+2] - level_begin_idx[level+1])/k² or rank l(evel_begin_idx[level], level_begin_idx[level+1])

            uint64_t amountOfZeros = this->m_levels[m_access_shortcut_size].size() / (get_k(m_access_shortcut_size) *
                                                                                      get_k(m_access_shortcut_size)); //spares rank of comp. level
            //corresponds to the theoretical amount of trees in level m_access_shortcut_size (round up (in case not divisible by k^2)
            uint64_t amountOfOnes = ipow(k0 * k0, m_access_shortcut_size);
            bit_vector access_shortcut(amountOfOnes + amountOfZeros, 1);

            //BitArray<uint> B(amountOfOnes + amountOfZeros);
            uint counter = 0;
            construct_access_shortcut_by_dfs(access_shortcut, this->root(), this->m_max_element, counter);
            this->m_access_shortcut.swap(access_shortcut);

            /*
            std::cout << "Access shortcut: " << std::endl;
            for (int i = 0; i < m_access_shortcut.size(); ++i) {
                std::cout << m_access_shortcut[i];
            }
            std::cout << std::endl;*/

            sdsl::util::init_support(this->m_access_shortcut_rank_10_support, &this->m_access_shortcut);
            sdsl::util::init_support(this->m_access_shortcut_select_1_support, &this->m_access_shortcut);

            perform_access_shortcut_precomputations();
        }

    protected:
        /*##################### Leaf access for uncompressed version#################################################**/
        inline bool is_leaf_bit_set(uint64_t pos) const {
            return m_leaves[pos];
        }

    private:
        template<typename t_x>
        inline void
        check_leaf_bits_direct_uncomp(int64_t pos, t_x result_offset, uint8_t leafK, std::vector<t_x> &result) const {
            for (int j = 0; j < leafK; ++j) {
                if (m_leaves[pos + j] == 1) {
                    result.push_back(result_offset + j);
                }
            }
        }


        template<typename t_x>
        inline void
        check_leaf_bits_inverse_uncomp(int64_t pos, t_x result_offset, uint8_t leafK, std::vector<t_x> &result) const {
            for (int j = 0; j < leafK; ++j) {
                if (m_leaves[pos + j * leafK] == 1) {
                    result.push_back(result_offset + j);
                }
            }
        }
        /*##################### Leaf Access for compressed version###################################################**/
    protected:
        /**
        * Returns word containing the bit at the given position
        * It access the corresponding word in the DAC.
        *
        * @param pos Position in the complete sequence of bit of the last level.
        * @return Pointer to the first position of the word.
        */
        inline bool is_leaf_bit_set_comp(uint64_t pos, uint8_t leafK) const {
            uint64_t subtree_number = pos / (leafK * leafK);
            uint iword = m_comp_leaves.accessFT(subtree_number);
            pos = pos - (subtree_number * leafK * leafK);
            const uchar *word = m_vocabulary->get(iword);
            bool bitSet = ((word[pos / kUcharBits] >> (pos % kUcharBits)) & 1);
            return bitSet;
        }


        void postInit(uint64_t hash_size) {
            if (m_tree_height > 0) {
                if (m_is_dac_comp) {
                    compress_leaves(hash_size);
                }
                if (m_access_shortcut_size > 0) {
                    construct_access_shortcut(m_access_shortcut_size);
                }
            }
        }

    private:

        /** Checks the leaf bits relevant for a direct neighbor query starting from leaf position pos.
         *
         * @param pos
         * @param result_offset
         * @param leafK
         * @param result
         */
        template<typename t_x>
        inline void
        check_leaf_bits_direct_comp(int64_t pos, t_x result_offset, uint8_t leafK, std::vector<t_x> &result) const {
            uint64_t subtree_number = pos / (leafK * leafK);
            uint iword = m_comp_leaves.accessFT(subtree_number);
            const uchar *word = m_vocabulary->get(iword);
            pos = pos - (subtree_number * leafK * leafK);
            for (int i = 0; i < leafK; ++i) {
                if ((word[(pos + i) / kUcharBits] >> ((pos + i) % kUcharBits)) & 1) {
                    result.push_back(i + result_offset);
                }
            }
        }

        /** Checks the leaf bits relevant for a inverse neighbor query starting from leaf position pos.
         *
         * @param pos
         * @param result_offset
         * @param leafK
         * @param result
         */
        template<typename t_x>
        inline void
        check_leaf_bits_inverse_comp(int64_t pos, t_x result_offset, uint8_t leafK, std::vector<t_x> &result) const {
            uint64_t subtree_number = pos / (leafK * leafK);
            uint iword = m_comp_leaves.accessFT(subtree_number);
            const uchar *word = m_vocabulary->get(iword);
            pos = pos - (subtree_number * leafK * leafK);
            for (int i = 0; i < leafK; ++i) {
                if ((word[(pos + i * leafK) / kUcharBits] >> ((pos + i * leafK) % kUcharBits)) & 1) {
                    result.push_back(i + result_offset);
                }
            }
        }


        /** Checks the leaf bits relevant for a direct neighbor query starting from leaf position pos given
         *  in case leaves are compressed using a huffman shaped wavelet tree.
         *
         * @param pos
         * @param result_offset
         * @param leafK
         * @param result
         */
        template<typename t_x>
        inline void
        check_leaf_bits_direct_wt(int64_t pos, t_x result_offset, uint8_t leafK, std::vector<t_x> &result) const {
            auto word_number = pos / m_wt_word_size;
            auto word = m_leaves_wt[word_number];
            auto offset = 0;
            for (int j = 0; j < leafK; ++j) {
                auto current_word_number = (pos + j) / m_wt_word_size;
                if (current_word_number > word_number) {
                    word_number = current_word_number;
                    word = m_leaves_wt[word_number];
                }

                offset = (pos + j) % m_wt_word_size;
                if (word >> (offset) & 1) {
                    result.push_back(j + result_offset);
                }
            }
        }

        /** Checks the leaf bits relevant for a direct neighbor query starting from leaf position pos given
        *  in case leaves are compressed using a huffman shaped wavelet tree.
        *
        * @param pos
        * @param result_offset
        * @param leafK
        * @param result
        */
        template<typename t_x>
        inline void
        check_leaf_bits_inverse_wt(int64_t pos, t_x result_offset, uint8_t leafK, std::vector<t_x> &result) const {
            //std::cout << "Checking posistion" << pos << std::endl;

            auto word_number = pos / m_wt_word_size;
            auto word = m_leaves_wt[word_number];
            auto offset = 0;
            for (int j = 0; j < leafK; ++j) {
                auto current_word_number = (pos + j * leafK) / m_wt_word_size;
                if (current_word_number > word_number) {
                    word_number = current_word_number;
                    word = m_leaves_wt[word_number];
                }

                offset = (pos + j * leafK) % m_wt_word_size;
                if (word >> (offset) & 1) {
                    result.push_back(j + result_offset);
                }
            }
        }

        /*##################### Leaf access for wt compressed version#################################################**/
        inline bool is_leaf_bit_set_wt(uint64_t pos) const {
            auto word_number = pos / m_wt_word_size;
            auto word = m_leaves_wt[word_number];
            auto offset = pos % m_wt_word_size;
            return (word >> (offset) & 1);
        }


        /**
         * Constructs the bitvector m_access_shortcut used for speeding up tree traversal/link checks
         */
        void construct_access_shortcut_by_dfs(bit_vector &access_shortcut, node_type root, uint64_t matrix_size,
                                              uint &counter) {
            using namespace k2_treap_ns;
            uint64_t rank = this->m_levels_rank[root.t](root.idx);
            auto x = std::real(root.p);
            auto y = std::imag(root.p);
            auto k = get_k(0);
            auto submatrix_size = matrix_size / k;

            for (size_t i = 0; i < k; ++i) {
                auto _x = x + i * submatrix_size;
                for (size_t j = 0; j < k; ++j) {
                    // get_int better for compressed bitvectors
                    // or introduce cache for bitvectors
                    if (root.t == (m_access_shortcut_size - 1)) {
                        if (this->m_levels[root.t][root.idx + k * i + j]) { //if subtree present
                            counter++;
                            access_shortcut[counter] = 0;//save 01 at counter position (m_access_shortcut gets initialised with 1s)
                        }
                        counter++;


                    } else { //continue dfs tree traversal
                        if (this->m_levels[root.t][root.idx + k * i + j]) { //if subtree present
                            auto _y = y + j * submatrix_size;
                            node_type subtree_root(root.t + 1, t_p(_x, _y), rank * k * k);
                            ++rank;
                            construct_access_shortcut_by_dfs(access_shortcut, subtree_root, submatrix_size, counter);
                        } else {
                            counter += ipow(k * k, ((m_access_shortcut_size - 1) - root.t));
                        }
                    }
                }
            }
        }
    };
}

#endif

