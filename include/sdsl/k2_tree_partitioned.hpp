#ifndef INCLUDED_SDSL_K2_TREE_PARTITONED
#define INCLUDED_SDSL_K2_TREE_PARTITONED

#include "vectors.hpp"
#include "bits.hpp"
#include <tuple>
#include <algorithm>
#include <climits>
#include <vector>
#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/rank_support_v.hpp>
#include <gtest/gtest_prod.h>

namespace sdsl {

    /**
     * t_comp can be used to compress the leaves of the trees,
     * pay attention: currently for leaf compression of the partitioned tree set the compression
     * for the subtrees to false. In this case a global vocabulary is used.
     */
    template<uint8_t t_k0,
            typename subk2_tree,
            bool t_comp = false>
    class k2_tree_partitioned {
        static_assert(t_k0 > 1, "t_k has to be larger than 1.");
        static_assert(t_k0 <= 16, "t_k has to be smaller than 17.");

    public:

        typedef int_vector<>::size_type size_type;

        enum {
            k0 = t_k0
        };

        size_type m_size = 0;
        uint64_t m_part_matrix_size;

    private:
        //typedef k2_tree<t_k,t_bv,t_rank> k2;
        std::vector<subk2_tree> m_k2trees;
        int_vector<> m_words_prefix_sum; //FIXME: think about using an sdsl compressed vector!!!
        DAC m_comp_leaves;
        /** For compressed version **/
        Vocabulary m_vocabulary;

    public:
        k2_tree_partitioned() = default;

        k2_tree_partitioned(const k2_tree_partitioned &tr) {
            *this = tr;
        }

        k2_tree_partitioned(k2_tree_partitioned &&tr) {
            *this = std::move(tr);
        }

        template<typename t_vector>
        k2_tree_partitioned(std::string temp_file_prefix, bool bottom_up, t_vector &v, uint64_t max_hint) {
            //FIXME: build multiple k trees
            //partition into k02 parts
            using namespace k2_treap_ns;
            uint64_t matrix_size = precomp<t_k0>::exp(get_tree_height(v, max_hint));
            m_part_matrix_size = matrix_size / t_k0;
            build_k2_trees(v, temp_file_prefix, bottom_up);

            if (t_comp){
                compress_leaves();
            }
        }

        k2_tree_partitioned(int_vector_buffer<> &buf_x,
                            int_vector_buffer<> &buf_y, bool bottom_up=false, uint64_t max_hint = 0) {
            using namespace k2_treap_ns;
            typedef int_vector_buffer<> *t_buf_p;
            std::vector<t_buf_p> bufs = {&buf_x, &buf_y};


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

            uint8_t res = 0;
            while (res <= 64 and precomp<t_k0>::exp(res) <= max) { ++res; }
            if (res == 65) {
                throw std::logic_error("Maximal element of input is too big.");
            }
            uint64_t matrix_size = precomp<t_k0>::exp(res);
            m_part_matrix_size = matrix_size / t_k0;

            if (precomp<t_k0>::exp(res) <= std::numeric_limits<uint32_t>::max()) {
                auto v = read < uint32_t, uint32_t>(bufs);
                build_k2_trees(v, buf_x.filename(), bottom_up);
            } else {
                auto v = read < uint64_t, uint64_t>(bufs);
                build_k2_trees(v, buf_x.filename(), bottom_up);
            }

            if (t_comp){
                compress_leaves();
            }
        }

        template<typename t_x=uint64_t, typename t_y=uint64_t>
        std::vector<std::pair<t_x, t_y>> read(std::vector<int_vector_buffer<> * > &bufs) {
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

        template<typename t_x>
        void direct_links(t_x source_id, std::vector<t_x> &result) const {
            result.clear();

            if (t_comp){
                std::cerr << "direct_links access method not implemented for compressed version, use direct_links2" << std::endl;
                return;
            }

            uint y = source_id/m_part_matrix_size;
            //uint submatrix_size = m_matrix_size/k0;

            //TODO: might be slow to use extra vector as reference, maybe it's better to remove clear() in k2treap and add directly to endresult
            std::vector<t_x> tmp_result;
            for (int j = 0; j < t_k0; ++j) {
                    m_k2trees[y*t_k0+j].direct_links( (t_x) (source_id - y * m_part_matrix_size), tmp_result);
                    for (auto item : tmp_result){
                        result.push_back(item + j*m_part_matrix_size);
                    }
            }
        }

        template<typename t_x>
        void direct_links_shortcut(t_x source_id, std::vector<t_x> &result) const {
            result.clear();
            uint y = source_id/m_part_matrix_size;
            //uint submatrix_size = m_matrix_size/k0;

            //TODO: might be slow to use extra vector as reference, maybe it's better to remove clear() in k2treap and add directly to endresult
            std::vector<t_x> tmp_result;
            if (t_comp){
                for (int j = 0; j < t_k0; ++j) {
                    uint index = y*t_k0+j;
                    tmp_result.clear();
                    m_k2trees[index].direct_links_shortcut_internal((t_x) (source_id - y * m_part_matrix_size), tmp_result, [index,this](int64_t pos, t_x offset, uint8_t leafK, std::vector<t_x> & result){
                        check_leaf_bits_direct_comp(index, pos, offset, leafK, result);
                    });
                    for (auto item : tmp_result){
                        result.push_back(item + j*m_part_matrix_size);
                    }
                }
            } else {
                for (int j = 0; j < t_k0; ++j) {
                    m_k2trees[y*t_k0+j].direct_links_shortcut((t_x) (source_id - y * m_part_matrix_size), tmp_result);
                    for (auto item : tmp_result){
                        result.push_back(item + j*m_part_matrix_size);
                    }
                }
            }
        }

        template<typename t_x>
        void direct_links_shortcut_2(t_x source_id, std::vector<t_x> &result) const {
            result.clear();
            uint y = source_id/m_part_matrix_size;
            //uint submatrix_size = m_matrix_size/k0;

            //TODO: might be slow to use extra vector as reference, maybe it's better to remove clear() in k2treap and add directly to endresult
            std::vector<t_x> tmp_result;
            if (t_comp){
                for (int j = 0; j < t_k0; ++j) {
                    uint index = y*t_k0+j;
                    tmp_result.clear();
                    m_k2trees[index].direct_links_shortcut_internal_2((t_x) (source_id - y * m_part_matrix_size), tmp_result, [index,this](int64_t pos, t_x offset, uint8_t leafK, std::vector<t_x> & result){
                        check_leaf_bits_direct_comp(index, pos, offset, leafK, result);
                    });
                    for (auto item : tmp_result){
                        result.push_back(item + j*m_part_matrix_size);
                    }
                }
            } else {
                for (int j = 0; j < t_k0; ++j) {
                    m_k2trees[y*t_k0+j].direct_links_shortcut_2((t_x) (source_id - y * m_part_matrix_size), tmp_result);
                    for (auto item : tmp_result){
                        result.push_back(item + j*m_part_matrix_size);
                    }
                }
            }
        }

        template<typename t_x>
        void direct_links2(t_x source_id, std::vector<t_x> &result) const {
            result.clear();
            t_x y = source_id/m_part_matrix_size;
            t_x y_offsetted = source_id - y * m_part_matrix_size;
            //uint submatrix_size = m_matrix_size/k0;

            //TODO: might be slow to use extra vector as reference, maybe it's better to remove clear() in k2treap and add directly to endresult
            std::vector<t_x> tmp_result;
            if (t_comp){
                for (int j = 0; j < t_k0; ++j) {
                    uint index = y*t_k0+j;
                    if (m_k2trees[index].m_tree_height > 0 && m_k2trees[index].m_max_element >= y_offsetted){
                        tmp_result.clear();
                        m_k2trees[index].direct_links2_internal_queue(y_offsetted, tmp_result, [index,this](int64_t pos, t_x offset, uint8_t leafK, std::vector<t_x> & asd_result){
                            check_leaf_bits_direct_comp(index, pos, offset, leafK, asd_result);
                        });
                        for (auto item : tmp_result){
                            result.push_back(item + j*m_part_matrix_size);
                        }
                    }
                }
            } else {
                for (int j = 0; j < t_k0; ++j) {
                    m_k2trees[y*t_k0+j].direct_links2(y_offsetted, tmp_result);
                    for (auto item : tmp_result){
                        result.push_back(item + j*m_part_matrix_size);
                    }
                }
            }
        }

        template<typename t_x>
        void inverse_links(t_x source_id, std::vector<t_x> &result) const {
            result.clear();

            if (t_comp){
                std::cerr << "inverse_links access method not implemented for compressed version, use inverse_links2" << std::endl;
                return;
            }

            uint x = source_id/m_part_matrix_size;
            //uint submatrix_size = m_matrix_size/k0;

            //TODO: might be slow to use extra vector as reference, maybe it's better to remove clear() in k2treap and add directly to endresult
            std::vector<t_x> tmp_result;
            for (int j = 0; j < t_k0; ++j) {
                m_k2trees[j*t_k0+x].inverse_links( (t_x) (source_id - x * m_part_matrix_size), tmp_result);
                for (auto item : tmp_result){
                    result.push_back(item + j*m_part_matrix_size);
                }
            }
        }

        template<typename t_x>
        void inverse_links2(t_x source_id, std::vector<t_x> &result) const {
            result.clear();
            t_x x = source_id/m_part_matrix_size;
            t_x x_offsetted = source_id - x * m_part_matrix_size;
            //uint submatrix_size = m_matrix_size/k0;

            //TODO: might be slow to use extra vector as reference, maybe it's better to remove clear() in k2treap and add directly to endresult
            std::vector<t_x> tmp_result;

            if (t_comp) {
                for (int j = 0; j < t_k0; ++j) {
                    uint index = j * t_k0 + x;
                    if (m_k2trees[index].m_max_element >= x_offsetted && m_k2trees[index].m_tree_height > 0) {
                        tmp_result.clear();
                        m_k2trees[index].inverse_links2_internal_queue(x_offsetted, tmp_result, [index, this](int64_t pos, t_x offset,
                                                                                     uint8_t leafK,
                                                                                     std::vector<t_x> &result) {
                                                                           check_leaf_bits_inverse_comp(index, pos,
                                                                                                        offset, leafK,
                                                                                                        result);
                                                                       });
                        for (auto item : tmp_result) {
                            result.push_back(item + j * m_part_matrix_size);
                        }
                    }
                }
            } else {
                for (int j = 0; j < t_k0; ++j) {
                    m_k2trees[j * t_k0 + x].inverse_links2(x_offsetted, tmp_result);
                    for (auto item : tmp_result) {
                        result.push_back(item + j * m_part_matrix_size);
                    }
                }
            }
        }

        template<typename t_x>
        inline void  check_leaf_bits_direct_comp(uint index, int64_t pos, t_x result_offset, uint8_t leafK, std::vector<t_x> & result) const {
            uint64_t subtree_number = pos/(leafK*leafK);
            uint64_t global_subtree_number = subtree_number+m_words_prefix_sum[index];
            uint iword = m_comp_leaves.accessFT(global_subtree_number);
            const uchar * word = m_vocabulary.get(iword);
            pos = pos - (subtree_number*leafK*leafK);
            for (int i = 0; i < leafK; ++i) {
                if ((word[(pos+i)/kUcharBits] >> ((pos+i)%kUcharBits)) & 1){
                    result.push_back(i+result_offset);
                }
            }
        }

        template<typename t_x>
        void check_leaf_bits_inverse_comp(uint index, int64_t pos, t_x result_offset, uint8_t leafK, std::vector<t_x> &result) const {
            uint64_t subtree_number = pos/(leafK*leafK);
            uint64_t global_subtree_number = subtree_number+m_words_prefix_sum[index];
            uint iword = m_comp_leaves.accessFT(global_subtree_number);
            const uchar * word = m_vocabulary.get(iword);
            pos = pos - (subtree_number*leafK*leafK);
            for (int i = 0; i < leafK; ++i) {
                if ((word[(pos+i*leafK)/kUcharBits] >> ((pos+i*leafK)%kUcharBits)) & 1){
                    result.push_back(i+result_offset);
                }
            }
        }

        inline bool is_leaf_bit_set_comp(uint index, uint64_t pos, uint8_t leafK) const {
            //std::cout << "Pos " << std::to_string(pos) << std::endl;
            uint64_t subtree_number = pos/(leafK*leafK);
            uint64_t global_subtree_number = subtree_number+m_words_prefix_sum[index];
            //std::cout << "Checking leaf bit of global subtree " << global_subtree_number << std::endl;
            uint iword = m_comp_leaves.accessFT(global_subtree_number);
            pos = pos - (subtree_number*leafK*leafK);
            const uchar * word = m_vocabulary.get(iword);
            //std::cout << "Word " << std::to_string(word[pos/kUcharBits]) << std::endl;
            //std::cout << "RelPos " << std::to_string(pos) << std::endl;
            bool bitSet = ((word[pos/kUcharBits] >> (pos%kUcharBits)) & 1);
            //std::cout << "Resultabab: " << bitSet << std::endl;
            return bitSet;
        }

        template<typename t_x>
        void inverse_links_shortcut(t_x source_id, std::vector<t_x> &result) const {
            result.clear();
            uint x = source_id/m_part_matrix_size;
            //TODO: might be slow to use extra vector as reference, maybe it's better to remove clear() in k2treap and add directly to endresult
            std::vector<t_x> tmp_result;
            if (t_comp) {
                for (int j = 0; j < t_k0; ++j) {
                    uint index = j * t_k0 + x;
                    tmp_result.clear();
                    m_k2trees[index].inverse_links_shortcut_internal((t_x) (source_id - x * m_part_matrix_size), tmp_result,[index,this](int64_t pos, t_x offset, uint8_t leafK, std::vector<t_x> & result){
                        check_leaf_bits_inverse_comp(index, pos, offset, leafK, result);
                    });
                    for (auto item : tmp_result) {
                        result.push_back(item + j * m_part_matrix_size);
                    }
                }
            } else {
                for (int j = 0; j < t_k0; ++j) {
                    m_k2trees[j * t_k0 + x].inverse_links_shortcut((t_x) (source_id - x * m_part_matrix_size), tmp_result);
                    for (auto item : tmp_result) {
                        result.push_back(item + j * m_part_matrix_size);
                    }
                }
            }
        }

        /**
         * Checks whether link from p = link.first to q = link.second is present i.e. matrix entry a_pq = 1
         */
        template<typename t_x, typename t_y>
        bool check_link(std::pair<t_x, t_y> link) const {
            t_x x = link.first/m_part_matrix_size;
            t_y y = link.second/m_part_matrix_size;

            t_x x_offsetted = link.first - x * m_part_matrix_size;
            t_y y_offsetted = link.second -y * m_part_matrix_size;

            uint index = x*t_k0+y;
            if (t_comp){
                if (m_k2trees[index].m_tree_height > 0) {
                    //std::cout << "Checking link in tree " << index << std::endl;
                    return m_k2trees[index].check_link_internal(0, m_k2trees[index].m_max_element, x_offsetted, y_offsetted, 0, [index,this](int64_t pos, uint8_t leafK){
                        return is_leaf_bit_set_comp(index, pos, leafK);;
                    });
                } else {
                    return false;
                }
            } else {
                return m_k2trees[index].check_link(std::make_pair(x_offsetted, y_offsetted));
            }

        }

        /**
         * Checks whether link from p = link.first to q = link.second is present i.e. matrix entry a_pq = 1
         */
        template<typename t_x, typename t_y>
        bool check_link_shortcut(std::pair<t_x, t_y> link) const {
            uint x = link.first/m_part_matrix_size;
            uint y = link.second/m_part_matrix_size;

            t_x x_offsetted = link.first - x * m_part_matrix_size;
            t_y y_offsetted = link.second -y * m_part_matrix_size;

            uint index = x*t_k0+y;
            if (t_comp){
                if (m_k2trees[index].m_tree_height > 0) {
                    //std::cout << "Checking link in tree " << index << std::endl;
                    return m_k2trees[index].check_link_shortcut_internal(std::make_pair(x_offsetted, y_offsetted),[index,this](int64_t pos, uint8_t leafK){
                        return is_leaf_bit_set_comp(index, pos, leafK);;
                    });
                } else {
                    return false;
                }
            } else {
                return m_k2trees[index].check_link_shortcut(std::make_pair(x_offsetted, y_offsetted));
            }
        }

        //! Move assignment operator
        k2_tree_partitioned &operator=(k2_tree_partitioned &&tr) {
            if (this != &tr) {
                m_size = tr.m_size;
                m_k2trees = std::move(tr.m_k2trees);
                m_part_matrix_size = std::move(tr.m_part_matrix_size);
                m_words_prefix_sum = std::move(m_words_prefix_sum);
            }
            return *this;
        }

        //! Assignment operator
        k2_tree_partitioned &operator=(k2_tree_partitioned &tr) {
            if (this != &tr) {
                m_size = tr.m_size;
                m_k2trees = tr.m_k2trees;
                m_part_matrix_size = tr.m_part_matrix_size;
                m_comp_leaves = tr.m_comp_leaves;
                m_vocabulary = tr.m_vocabulary;
                m_words_prefix_sum = m_words_prefix_sum;
            }
            return *this;
        }

        //! Equals operator
        bool operator==(const k2_tree_partitioned &tr) const {
            if (m_size != tr.m_size)
                return false;

            if (m_k2trees != tr.m_k2trees){
                return false;
            }

            if (m_part_matrix_size != tr.m_part_matrix_size){
                return false;
            }

            if (t_comp) {
                if (!(m_comp_leaves == tr.m_comp_leaves)) {
                    std::cout << "comp leaves differ" << std::endl;
                    return false;
                }


                if (!m_vocabulary.operator==(tr.m_vocabulary)) {
                    std::cout << "vocabulary differs" << std::endl;
                    return false;
                }

                if (m_words_prefix_sum != m_words_prefix_sum){
                    std::cout << "Prefix sum differs" << std::endl;
                    return false;
                }
            }

            return true;
        }

        //! Number of points in the 2^k treap
        size_type
        size() const {
            return m_size;
        }

        //! Swap operator
        void swap(k2_tree_partitioned &tr) {
            if (this != &tr) {
                std::swap(m_size, tr.m_size);
                std::swap(m_k2trees, tr.m_k2trees);
                std::swap(m_part_matrix_size, tr.m_part_matrix_size);
                std::swap(m_vocabulary, tr.m_vocabulary);
                m_comp_leaves.swap(tr.m_comp_leaves);
                m_words_prefix_sum.swap(tr.m_words_prefix_sum);
            }
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream &out, structure_tree_node *v = nullptr,
                            std::string name = "") const {

            structure_tree_node *child = structure_tree::add_child(
                    v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "s");
            written_bytes += write_member(m_part_matrix_size, out, child, "matrix_size");

            for (uint j = 0; j < m_k2trees.size(); ++j) {
                written_bytes += m_k2trees[j].serialize(out, child, "k2_tree_"+std::to_string(j));
            }

            if (t_comp){
                written_bytes += m_vocabulary.serialize(out, child, "voc");
                written_bytes += m_comp_leaves.serialize(out, child, "comp_leafs");
                written_bytes += m_words_prefix_sum.serialize(out, child, "words_prefix_sum");
            }

            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream &in) {
            read_member(m_size, in);
            read_member(m_part_matrix_size, in);
            m_k2trees.resize(t_k0*t_k0);
            for (uint j = 0; j < t_k0*t_k0; ++j) {
                m_k2trees[j].load(in);
            }

            if (t_comp){
                m_vocabulary.load(in);
                m_comp_leaves.load(in);
                m_words_prefix_sum.load(in);
            }
        }

        /**
        * Return the number of bytes necessary to store a word.
        *
        * @return Size of a word.
        */
        uint word_size() const {
            return m_k2trees[0].word_size();
        }

        /**
        * Returns the number of words of \f$k_leaves^2\f$ bits in the leaf level.
        *
        * @return Number of words.
        */
        size_t words_count() const{
            size_t words_count = 0;
            for (uint i = 0; i < m_k2trees.size(); ++i){
                words_count += m_k2trees[i].words_count();
                //std::cout << "Words in tree " << i << ": " <<m_k2trees[i].words_count() << std::endl;
            }
            return words_count;
        }
        /**
        * Iterates over the words in the leaf level.
        *
        * @param fun Pointer to function, functor or lambda expecting a pointer to
        * each word.
        */
        template<typename Function>
        void words(Function fun) const {
            for (uint i = 0; i < m_k2trees.size(); ++i){
                m_k2trees[i].words(fun);
            }
        }

    private:
        template<typename t_tv>
        uint8_t get_tree_height(const t_tv &v, uint64_t max_hint) {
            using namespace k2_treap_ns;
            if (v.size() == 0) {
                return 0;
            }

            uint64_t max;
            if (max_hint != 0){
                max = max_hint;
            } else {
                using t_e = typename t_tv::value_type;
                auto tupmax = [](t_e a) {
                    return std::max(a.first, a.second);
                };
                auto max_it = std::max_element(std::begin(v), std::end(v), [&](t_e a, t_e b) {
                    return tupmax(a) < tupmax(b);
                });
                max = tupmax(*max_it);
            }


            uint8_t res = 0;
            while (precomp<t_k0>::exp(res) <= max) { ++res; }
            return res;
        }

        template<typename t_vector>
        void build_k2_trees(t_vector &links, std::string temp_file_prefix = "", bool bottom_up = false) {
            using namespace k2_treap_ns;
            typedef decltype(links[0].first) t_x;
            typedef decltype(links[0].second) t_y;
            using t_e = std::pair<t_x, t_y>;

            m_size = links.size();

            //typedef typename stxxl::VECTOR_GENERATOR<t_e>::result stxxl_pair_vector;

            //FIXME replace with stxxl vector
            std::vector<std::vector<t_e>> buffers;
            buffers.resize(t_k0*t_k0);

            {
                for (uint64_t j = 0; j < links.size(); ++j) {
                    auto x = links[j].first;
                    auto y = links[j].second;
                    uint p1 = x / m_part_matrix_size;
                    uint p2 = y / m_part_matrix_size;
                    uint corresponding_matrix = p1 * t_k0 + p2;
                    buffers[corresponding_matrix].push_back(t_e(links[j].first - p1 * m_part_matrix_size, links[j].second - p2 * m_part_matrix_size));
                }
            }

            m_k2trees.resize(t_k0*t_k0);
            for (int l = 0; l < t_k0*t_k0; ++l) {
                const subk2_tree k2tree(temp_file_prefix, bottom_up, buffers[l]);
                m_k2trees[l] = k2tree;//buffers[l]);//, temp_file_prefix, bottom_up, access_shortcut_size);
            }
        }

        void compress_leaves() {
            //std::cout << "Words count " << words_count() << std::endl;
            //std::cout << "Word size " << word_size() << std::endl;


            /*std::cout << "Words" << std::endl;
            words([&] (const uchar *word) {
                std::cout << std::to_string(*word) << "\t";
            });
            std::cout << std::endl;
            */

            FreqVoc(*this, [&](const HashTable &table, Vocabulary& voc) {
                compress_leaves(table, voc);
            });
        }

        void compress_leaves(const HashTable &table, Vocabulary& voc) {
            size_t cnt = words_count();
            uint size = word_size();
            uint *codewords;
            try {
                codewords = new uint[cnt];
            } catch (std::bad_alloc ba) {
                std::cerr << "[k2_tree_base::compress_leaves] Error: " << ba.what() << "\n";
                exit(1);
            }

            size_t i = 0;

            words([&] (const uchar *word) {
                size_t addr;
                if (!table.search(word, size, &addr)) {
                    std::cerr << "[k2_tree_base::compress_leaves] Error: Word not found\n";
                    exit(1);
                }
                codewords[i++] = table[addr].codeword;
            });


            try {
                // TODO Port to 64-bits
                m_comp_leaves = DAC(codewords, cnt);
            } catch (...) {
                std::cerr << "[HybridK2Tree::CompressLeaves] Error: Could not create DAC\n";
                exit(1);
            }

            delete[] codewords;

            m_vocabulary = voc;

            m_words_prefix_sum.resize(t_k0*t_k0);
            m_words_prefix_sum[0] = 0;
            for (uint i = 1; i < m_k2trees.size(); ++i) {
                size_t count = m_k2trees[i-1].words_count();
                m_words_prefix_sum[i] = m_words_prefix_sum[i-1] + count;
            }
            //util::bit_compress(m_words_prefix_sum);
            for (uint i = 0; i < m_k2trees.size(); ++i){
                m_k2trees[i].clear_leaves();
            }

            /*
            std::cout << "Words 2" << std::endl;
            for(uint i = 0; i < 9; ++i){
                uint iword = m_comp_leaves.accessFT(i);
                const uchar * word = m_vocabulary.get(iword);
                std::cout << std::to_string(*word) << "\t";
            }
            std::cout << std::endl;
            */
            /*
            std::cout << "words_prefix_sum" << std::endl;
            for (auto asd: m_words_prefix_sum){
                std::cout << asd << std::endl;
            }
            */
        }

    };
}

#endif