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
#include "mem_monitor.hpp"

namespace sdsl {

    /**
     * t_comp can be used to compress the leaves of the trees,
     * pay attention: currently for leaf compression of the partitioned tree set the compression
     * for the subtrees to false. In this case a global vocabulary is used.
     */
    template<uint8_t t_k0,
            typename subk2_tree>
    class k2_tree_partitioned {
        static_assert(t_k0 > 1, "t_k has to be larger than 1.");

    public:

        typedef int_vector<>::size_type size_type;

        enum {
            k0 = t_k0
        };

    private:
        //typedef k2_tree<t_k,t_bv,t_rank> k2;
        std::vector<subk2_tree> m_k2trees;
        /** For compressed version **/
        std::shared_ptr<Vocabulary> m_vocabulary;
        size_type m_size = 0;
        //uint64_t m_max_element = 0;
        uint64_t m_matrix_dimension = 0;
        bool m_is_dac_comp = false;
        uint8_t m_access_shortcut_size = 0;

    public:
        k2_tree_partitioned() = default;
        uint64_t m_max_element = 0; //for speed test

        /*
        k2_tree_partitioned(const k2_tree_partitioned &tr) {
            *this = tr;
        }

        k2_tree_partitioned(k2_tree_partitioned &&tr) {
            *this = std::move(tr);
        }*/

        template<typename t_vector>
        k2_tree_partitioned(std::string temp_file_prefix, bool use_counting_sort, t_vector &v, uint64_t max_hint=0, uint8_t access_shortcut_size=0, bool dac_compress = false, uint64_t hash_size = 0) : m_is_dac_comp(dac_compress), m_access_shortcut_size(access_shortcut_size) {
            //FIXME: build multiple k trees
            //partition into k02 parts
            using namespace k2_treap_ns;

            if (max_hint != 0){
                m_max_element = max_hint;
            } else {
                m_max_element = get_maximum(v);
            }


            calculate_matrix_dimension_and_submatrix_count();
            build_k2_trees(v, temp_file_prefix, use_counting_sort);

            if (m_is_dac_comp){
                compress_leaves(hash_size);
            }
        }

        k2_tree_partitioned(int_vector_buffer<> &buf_x,
                            int_vector_buffer<> &buf_y, bool use_counting_sort=false, uint64_t max_hint = 0, uint8_t access_shortcut_size=0, bool dac_compress = false, uint64_t hash_size = 0) : m_is_dac_comp(dac_compress), m_access_shortcut_size(access_shortcut_size){
            using namespace k2_treap_ns;
            typedef int_vector_buffer<> *t_buf_p;
            std::vector<t_buf_p> bufs = {&buf_x, &buf_y};


            if (max_hint != 0) {
                m_max_element = max_hint;//temporarily set
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

                m_max_element = max_buf_element();
            };

            calculate_matrix_dimension_and_submatrix_count();

            if (m_max_element <= std::numeric_limits<uint32_t>::max()) {
                auto v = read < uint32_t, uint32_t>(bufs);
                build_k2_trees(v, buf_x.filename(), use_counting_sort);
            } else {
                auto v = read < uint64_t, uint64_t>(bufs);
                build_k2_trees(v, buf_x.filename(), use_counting_sort);
            }

            if (m_is_dac_comp){
                compress_leaves(hash_size);
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
        void direct_links_shortcut(t_x source_id, std::vector<t_x> &result) const {
            result.clear();

            if (m_access_shortcut_size == 0){
                throw std::runtime_error("Access shortcut not present"); // throw for now
            }

            if (source_id > m_max_element){
                return;
            }
            uint y = source_id/m_matrix_dimension;
            //uint submatrix_size = m_matrix_size/k0;

            //TODO: might be slow to use extra vector as reference, maybe it's better to remove clear() in k2treap and add directly to endresult
            std::vector<t_x> tmp_result;
            for (uint j = 0; j < t_k0; ++j) {
                m_k2trees[y * t_k0 + j].direct_links_shortcut((t_x) (source_id - y * m_matrix_dimension), tmp_result);
                for (auto item : tmp_result) {
                    result.push_back(item + j * m_matrix_dimension);
                }
            }
        }

        template<typename t_x>
        void direct_links_shortcut_2(t_x source_id, std::vector<t_x> &result) const {
            result.clear();

            if (m_access_shortcut_size == 0){
                throw std::runtime_error("Access shortcut not present"); // throw for now
            }

            if (source_id > m_max_element){
                return;
            }

            uint y = source_id/m_matrix_dimension;
            //uint submatrix_size = m_matrix_size/k0;

            //TODO: might be slow to use extra vector as reference, maybe it's better to remove clear() in k2treap and add directly to endresult
            std::vector<t_x> tmp_result;
            for (uint j = 0; j < t_k0; ++j) {
                m_k2trees[y * t_k0 + j].direct_links_shortcut_2((t_x) (source_id - y * m_matrix_dimension), tmp_result);
                for (auto item : tmp_result) {
                    result.push_back(item + j * m_matrix_dimension);
                }
            }
        }

        template<typename t_x>
        void direct_links2(t_x source_id, std::vector<t_x> &result) const {
            result.clear();
            if (source_id > m_max_element){
                return;
            }

            t_x y = source_id/m_matrix_dimension;
            t_x y_offsetted = source_id - y * m_matrix_dimension;
            //uint submatrix_size = m_matrix_size/k0;

            //TODO: might be slow to use extra vector as reference, maybe it's better to remove clear() in k2treap and add directly to endresult
            std::vector<t_x> tmp_result;
            for (uint j = 0; j < t_k0; ++j) {
                    m_k2trees[y*t_k0+j].direct_links2(y_offsetted, tmp_result);
                    for (auto item : tmp_result){
                        result.push_back(item + j*m_matrix_dimension);
                    }
            }
        }

        template<typename t_x>
        void inverse_links2(t_x source_id, std::vector<t_x> &result) const {
            result.clear();

            if (source_id > m_max_element){
                return;
            }
            t_x x = source_id/m_matrix_dimension;
            t_x x_offsetted = source_id - x * m_matrix_dimension;
            //uint submatrix_size = m_matrix_size/k0;

            //TODO: might be slow to use extra vector as reference, maybe it's better to remove clear() in k2treap and add directly to endresult
            std::vector<t_x> tmp_result;
            for (uint j = 0; j < t_k0; ++j) {
                    m_k2trees[j * t_k0 + x].inverse_links2(x_offsetted, tmp_result);
                    for (auto item : tmp_result) {
                        result.push_back(item + j * m_matrix_dimension);
                    }
            }
        }

        template<typename t_x>
        void inverse_links_shortcut(t_x source_id, std::vector<t_x> &result) const {
            result.clear();

            if (source_id > m_max_element){
                return;
            }

            uint x = source_id/m_matrix_dimension;
            //TODO: might be slow to use extra vector as reference, maybe it's better to remove clear() in k2treap and add directly to endresult
            std::vector<t_x> tmp_result;
            for (uint j = 0; j < t_k0; ++j) {
                    m_k2trees[j * t_k0 + x].inverse_links_shortcut((t_x) (source_id - x * m_matrix_dimension), tmp_result);
                    for (auto item : tmp_result) {
                        result.push_back(item + j * m_matrix_dimension);
                    }
            }
        }

        /**
         * Checks whether link from p = link.first to q = link.second is present i.e. matrix entry a_pq = 1
         */
        template<typename t_x, typename t_y>
        bool check_link(std::pair<t_x, t_y> link) const {
            if (link.first > m_max_element || link.second >m_max_element){
                return false;
            }

            t_x x = link.first/m_matrix_dimension;
            t_y y = link.second/m_matrix_dimension;

            t_x x_offsetted = link.first - x * m_matrix_dimension;
            t_y y_offsetted = link.second -y * m_matrix_dimension;

            uint index = x*t_k0+y;
            return m_k2trees[index].check_link(std::make_pair(x_offsetted, y_offsetted));
        }

        /**
         * Checks whether link from p = link.first to q = link.second is present i.e. matrix entry a_pq = 1
         */
        template<typename t_x, typename t_y>
        bool check_link_shortcut(std::pair<t_x, t_y> link) const {
            if (link.first > m_max_element || link.second >m_max_element){
                return false;
            }

            uint x = link.first/m_matrix_dimension;
            uint y = link.second/m_matrix_dimension;

            t_x x_offsetted = link.first - x * m_matrix_dimension;
            t_y y_offsetted = link.second -y * m_matrix_dimension;

            uint index = x*t_k0+y;
            return m_k2trees[index].check_link_shortcut(std::make_pair(x_offsetted, y_offsetted));
        }

        //! Move assignment operator
        k2_tree_partitioned &operator=(k2_tree_partitioned &&tr) {
            if (this != &tr) {
                m_size = tr.m_size;
                m_max_element = tr.m_max_element;
                m_is_dac_comp = tr.m_is_dac_comp;
                m_access_shortcut_size = tr.m_access_shortcut_size;
                m_k2trees = std::move(tr.m_k2trees);
                m_matrix_dimension = std::move(tr.m_matrix_dimension);
            }
            return *this;
        }

        //! Assignment operator
        k2_tree_partitioned &operator=(k2_tree_partitioned &tr) {
            if (this != &tr) {
                m_size = tr.m_size;
                m_max_element = tr.m_max_element;
                m_is_dac_comp = tr.m_is_dac_comp;
                m_access_shortcut_size = tr.m_access_shortcut_size;
                m_k2trees = tr.m_k2trees;
                m_matrix_dimension = tr.m_matrix_dimension;
                m_vocabulary = tr.m_vocabulary;
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

            if (m_matrix_dimension != tr.m_matrix_dimension){
                return false;
            }

            if (m_max_element != tr.m_max_element){
                return false;
            }

            if (m_access_shortcut_size != tr.m_access_shortcut_size){
                std::cout << "Shortcut size differs" << std::endl;
                return false;
            }

            if (m_is_dac_comp != tr.m_is_dac_comp){
                std::cout << "one is compressed, the other not" << std::endl;
            }

            if (m_is_dac_comp) {
                if (!m_vocabulary.get()->operator==(*tr.m_vocabulary.get())) {
                    std::cout << "vocabulary differs" << std::endl;
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
                std::swap(m_max_element, tr.m_max_element);
                std::swap(m_is_dac_comp, tr.m_is_dac_comp);
                std::swap(m_access_shortcut_size, tr.m_access_shortcut_size);
                std::swap(m_k2trees, tr.m_k2trees);
                std::swap(m_matrix_dimension, tr.m_matrix_dimension);
                std::swap(m_vocabulary, tr.m_vocabulary);
            }
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream &out, structure_tree_node *v = nullptr,
                            std::string name = "") const {

            structure_tree_node *child = structure_tree::add_child(
                    v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "s");
            written_bytes += write_member(m_matrix_dimension, out, child, "matrix_size");
            written_bytes += write_member(m_max_element, out, child, "m_max_element");
            written_bytes += write_member(m_is_dac_comp, out, child, "m_is_dac_comp");
            written_bytes += write_member(m_access_shortcut_size, out, child, "m_access_shortcut_size");

            for (uint j = 0; j < m_k2trees.size(); ++j) {
                written_bytes += m_k2trees[j].serialize(out, child, "k2_tree_"+std::to_string(j));
            }

            if (m_is_dac_comp){
                written_bytes += m_vocabulary->serialize(out, child, "voc");
            }

            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream &in) {
            read_member(m_size, in);
            read_member(m_matrix_dimension, in);
            read_member(m_max_element, in);
            read_member(m_is_dac_comp, in);
            read_member(m_access_shortcut_size, in);

            m_k2trees.resize(t_k0*t_k0);
            for (uint j = 0; j < t_k0*t_k0; ++j) {
                m_k2trees[j].load(in);
            }

            if (m_is_dac_comp){
                m_vocabulary = std::shared_ptr<Vocabulary>(new Vocabulary());
                m_vocabulary->load(in);
            }

            for (uint j = 0; j < t_k0*t_k0; ++j) {
                m_k2trees[j].set_vocabulary(m_vocabulary);
            }
        }

        void load_from_ladrabin(std::string fileName, bool use_counting_sort = false, uint8_t access_shortcut_size = 0, bool dac_compress = false,
                                uint64_t hash_size = 0, std::string temp_file_prefix = ""){
            using namespace k2_treap_ns;
            if(!has_ending(fileName, ".ladrabin")){
                fileName.append(".ladrabin");
                std::cout << "Appending .ladrabin to filename as file has to be in .ladrabin format" << std::endl;
            }

            m_is_dac_comp = dac_compress;
            m_access_shortcut_size = access_shortcut_size;

            std::fstream fileStream(fileName, std::ios_base::in);

            if (fileStream.is_open()){
                uint number_of_nodes;
                ulong number_of_edges;

                read_member(number_of_nodes, fileStream);
                read_member(number_of_edges, fileStream);

                m_max_element = number_of_nodes -1;
                m_size = number_of_edges;
                calculate_matrix_dimension_and_submatrix_count();

                uint nodes_read = 0;
                uint source_id;
                uint source_id_offsetted = 0;
                int target_id;

                std::vector<std::vector<std::pair<uint, uint>>> buffers;
                std::vector<uint> maximum_in_buffer(t_k0, 0);
                buffers.resize(t_k0); //build one row at a time
                m_k2trees.resize(t_k0*t_k0);

                uint current_matrix_row = 0;
                for (uint64_t i = 0; i < number_of_nodes + number_of_edges; i++) {
                    read_member(target_id, fileStream);
                    if (target_id < 0) {
                        nodes_read++;
                        source_id = nodes_read - 1;
                        source_id_offsetted = source_id%m_matrix_dimension;
                        uint corresponding_row = (nodes_read -1)/m_matrix_dimension;
                        if (corresponding_row > current_matrix_row){
                            construct_trees_from_buffers(current_matrix_row, use_counting_sort, temp_file_prefix,
                                                         buffers, maximum_in_buffer);

                            //in case of a complete empty row

                            for (uint k = 0; k < (corresponding_row - current_matrix_row - 1); ++k) {
                                current_matrix_row++;
                                std::cout << "Appending completely empty row: " << current_matrix_row << std::endl;
                                construct_trees_from_buffers(current_matrix_row, use_counting_sort, temp_file_prefix,
                                                             buffers, maximum_in_buffer);
                            }

                            current_matrix_row = corresponding_row;

                        }
                    } else {
                        uint column_in_matrix = (target_id)/m_matrix_dimension;
                        auto target_id_offsetted = target_id%m_matrix_dimension;

                        if (source_id_offsetted > maximum_in_buffer[column_in_matrix]){
                            maximum_in_buffer[column_in_matrix] = source_id_offsetted;
                        }
                        if (target_id_offsetted > maximum_in_buffer[column_in_matrix]){
                            maximum_in_buffer[column_in_matrix] = target_id_offsetted;
                        }

                        buffers[column_in_matrix].push_back(std::make_pair(source_id_offsetted, target_id_offsetted));
                    }
                }

                //cover leftovers
                construct_trees_from_buffers(current_matrix_row, use_counting_sort, temp_file_prefix,
                                             buffers, maximum_in_buffer);

                if (m_is_dac_comp){
                    compress_leaves(hash_size);
                }
            } else {
                throw std::runtime_error("Could not open file to load ladrabin graph");
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
        * Returns all the words on the leaf level in a vector, which has size word_size*word_count (word_size is the amount of byte-sized words needed per leaf)
         * FIXME: should use 64 Bit words in future
        *
        */

        void words(std::vector<uchar>& result) const {
            result.clear();
            result.reserve(words_count()*word_size());

            for (uint i = 0; i < m_k2trees.size(); ++i){
                m_k2trees[i].words(result, true);
            }
        }

        void construct_access_shortcut(uint8_t access_shortcut_size) {
            m_access_shortcut_size = access_shortcut_size;
            for (uint i = 0; i < m_k2trees.size(); ++i) {
                m_k2trees[i].construct_access_shortcut(access_shortcut_size);
            }
        }

    private:
        template<typename t_vector>
        void build_k2_trees(t_vector &links, std::string temp_file_prefix = "", bool use_counting_sort = false) {
            using namespace k2_treap_ns;
            typedef decltype(links[0].first) t_x;
            typedef decltype(links[0].second) t_y;
            using t_e = std::pair<t_x, t_y>;

            m_size = links.size();

            //typedef typename stxxl::VECTOR_GENERATOR<t_e>::result stxxl_pair_vector;

            //FIXME replace with stxxl vector
            std::vector<std::vector<t_e>> buffers;
            std::vector<uint64_t> maximum_in_buffer;
            buffers.resize(t_k0*t_k0);
            maximum_in_buffer.resize(t_k0*t_k0);

            {
                for (uint64_t j = 0; j < links.size(); ++j) {
                    auto x = links[j].first;
                    auto y = links[j].second;
                    auto p1 = x / m_matrix_dimension;
                    auto p2 = y / m_matrix_dimension;
                    auto corresponding_matrix = p1 * t_k0 + p2;
                    x = x - p1 * m_matrix_dimension;
                    y = y - p2 * m_matrix_dimension;
                    if (x > maximum_in_buffer[corresponding_matrix]){
                        maximum_in_buffer[corresponding_matrix] = x;
                    }
                    if (y > maximum_in_buffer[corresponding_matrix]){
                        maximum_in_buffer[corresponding_matrix] = y;
                    }
                    buffers[corresponding_matrix].push_back(t_e(x,y));
                }
            }

            uint64_t hash_size = 0;//for now
            m_k2trees.reserve(buffers.size());
            for (uint l = 0; l < buffers.size(); ++l) {
//                const subk2_tree k2tree(temp_file_prefix, use_counting_sort, buffers[l]);
                m_k2trees.emplace_back(temp_file_prefix, use_counting_sort, buffers[l], maximum_in_buffer[l], m_access_shortcut_size, false, hash_size);
                buffers[l].clear();
            }
        }

    public:
        void compress_leaves(uint64_t hash_size = 0) {
            //std::cout << "Words count " << words_count() << std::endl;
            //std::cout << "Word size " << word_size() << std::endl;


            /*std::cout << "Words" << std::endl;
            words([&] (const uchar *word) {
                std::cout << std::to_string(*word) << "\t";
            });
            std::cout << std::endl;
            */

            m_is_dac_comp = true;
            std::cout << "Compressing Leaves" << std::endl;

            std::vector<uchar> leaf_words;
            words(leaf_words);

            FreqVoc(leaf_words, word_size(), words_count(), [&](const HashTable &table, std::shared_ptr<Vocabulary> voc, const std::vector<uchar>&) {
                m_vocabulary = voc;
                compress_leaves(table, voc);
            }, hash_size);
        }

        std::string get_type_string() const {
            return "k2_tree_partitioned<"+std::to_string(t_k0)+","+m_k2trees[0].get_type_string()+">";
        }
    private:
        void compress_leaves(const HashTable table, std::shared_ptr<Vocabulary> voc) {
            std::cout << "After FreqVoc" << std::endl;
            std::cout << "Clearing Leaves" << std::endl;
            //util::bit_compress(m_words_prefix_sum);
            #pragma omp parallel for
            for (uint i = 0; i < m_k2trees.size(); ++i){
                m_k2trees[i].compress_leaves(table, voc, false);
                m_k2trees[i].clear_leaves();
	        }
        }


        void calculate_matrix_dimension_and_submatrix_count() {
            m_matrix_dimension = (m_max_element + t_k0) / t_k0; //round up also in the case divisiable as element ids start at 0
            std::cout << "Matrix dimension: " << m_matrix_dimension << std::endl;
            std::cout << "Submatrix amount per row: " << std::to_string(t_k0) << std::endl;
        }

        inline void
        construct_trees_from_buffers(uint current_matrix_row, bool use_counting_sort, std::string &temp_file_prefix,
                                     std::vector<std::vector<std::pair<uint, uint>>> &buffers, std::vector<uint>& maximum_in_buffer) {
            #pragma omp parallel for
            for (uint j = 0; j < t_k0; ++j) {
                //std::cout << "Constructing tree "<< current_matrix_row * t_k0 + j << std::endl;
                /*if (buffers[j].size() != 0) {
                    std::cout << "Size of " << current_matrix_row * t_k0 + j << ": "
                              << buffers[j].size() * 64 / 8 / 1024 << "kByte" << std::endl;
                }*/
                uint64_t hash_size = 0; //for now
                subk2_tree tree(temp_file_prefix, use_counting_sort, buffers[j], maximum_in_buffer[j], m_access_shortcut_size, false, hash_size);
                m_k2trees[current_matrix_row*t_k0+j].swap(tree);
                buffers[j].clear();
                maximum_in_buffer[j] = 0;
                //std::cout << "Assigning tree " << current_matrix_row * t_k0 + j << std::endl;
            }
        }

    };
}

#endif
