#ifndef INCLUDED_SDSL_K2_TREE_PARTITONED
#define INCLUDED_SDSL_K2_TREE_PARTITONED

#include "vectors.hpp"
#include "k2_tree_utility.hpp"
#include "bits.hpp"
#include <tuple>
#include <algorithm>
#include <climits>
#include <vector>
#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/rank_support_v.hpp>
#include "mem_monitor.hpp"

namespace sdsl {

    uint64_t word_iteration = 0;
    uint64_t frequency_encoding = 0;
    uint64_t dac_compression = 0;

    uint64_t subtree_construction_duration = 0;

    /**
     * t_comp can be used to compress the leaves of the trees,
     * pay attention: currently for leaf compression of the partitioned tree set the compression
     * for the subtrees to false. In this case a global vocabulary is used.
     */
    template<typename subk2_tree>
    class k2_tree_partitioned {

    public:

        typedef int_vector<>::size_type size_type;

    private:
        template <typename T> struct empty { };
        //typedef k2_tree_comp<t_k,t_bv,t_rank> k2;
        std::vector<subk2_tree> m_k2trees;
        /** For compressed version **/
        std::shared_ptr<k2_tree_vocabulary> m_vocabulary;
        /**For newer compression algos **/
        std::shared_ptr<int_vector<>> m_dictionary;
        size_type m_size = 0;
        //uint64_t m_max_element = 0;
        uint64_t m_matrix_dimension = 0;
        uint64_t m_submatrix_per_dim_count = 0; //amount of submatrices per dimension (rows, columns)
        leaf_compression_type  m_used_compression = UNCOMPRESSED;
        uint_fast8_t m_access_shortcut_size = 0;
        uint_fast8_t m_submatrix_shift = 1; //submartix size = 2^m_submatrix_shift

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
        k2_tree_partitioned(t_vector &v, construction_algorithm construction_algo = ZORDER_SORT,
                            uint_fast8_t submatrix_shift = 0, uint64_t max_hint=0, std::string temp_file_prefix = ""){
            //FIXME: build multiple k trees
            //partition into k02 parts
            using namespace k2_tree_ns;

            if (max_hint != 0){
                m_max_element = max_hint;
            } else {
                m_max_element = get_maximum(v);
            }

            m_submatrix_shift = calculateSubmatrixShift(submatrix_shift);

            calculate_matrix_dimension_and_submatrix_count();
            build_k2_trees(v, temp_file_prefix, construction_algo);
        }

        k2_tree_partitioned(int_vector_buffer<> &buf_x,
                            int_vector_buffer<> &buf_y, construction_algorithm construction_algo = ZORDER_SORT,
                            uint_fast8_t submatrix_shift = 0, uint64_t max_hint = 0){
            using namespace k2_tree_ns;
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

            m_submatrix_shift = calculateSubmatrixShift(submatrix_shift);

            calculate_matrix_dimension_and_submatrix_count();

            if (m_max_element <= std::numeric_limits<uint32_t>::max()) {
                auto v = read < uint32_t, uint32_t>(bufs);
                build_k2_trees(v, buf_x.filename(), construction_algo);
            } else {
                auto v = read < uint64_t, uint64_t>(bufs);
                build_k2_trees(v, buf_x.filename(), construction_algo);
            }
        }

        uint_fast8_t calculateSubmatrixShift(uint_fast8_t submatrix_shift) const {
            if (submatrix_shift == 0){
                submatrix_shift = bits::hi(m_max_element >> 3);
            }
            return submatrix_shift;
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
            //same as division by m_matrix_dimension
            uint y = source_id >> m_submatrix_shift;
            t_x y_offsetted = source_id - (y << m_submatrix_shift);

            std::vector<std::vector<t_x>> tmp_results;
            #pragma omp parallel num_threads(m_submatrix_per_dim_count)
            {
                #pragma omp single
                {
                    tmp_results.resize(omp_get_num_threads());
                }

                auto t = omp_get_thread_num();
                m_k2trees[y*m_submatrix_per_dim_count+t].direct_links_shortcut(y_offsetted, tmp_results[t]);
            }

            for (size_t t = 0; t < tmp_results.size(); ++t){
                for (auto item : tmp_results[t]){
                    result.push_back(item + (t<< m_submatrix_shift));
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

            uint y = source_id >> m_submatrix_shift;
            t_x y_offsetted = source_id - (y << m_submatrix_shift);

            std::vector<std::vector<t_x>> tmp_results;
            #pragma omp parallel num_threads(m_submatrix_per_dim_count)
            {
                #pragma omp single
                {
                    tmp_results.resize(omp_get_num_threads());
                }

                auto t = omp_get_thread_num();
                m_k2trees[y*m_submatrix_per_dim_count+t].direct_links_shortcut_2(y_offsetted, tmp_results[t]);
            }

            for (size_t t = 0; t < tmp_results.size(); ++t){
                for (auto item : tmp_results[t]){
                    result.push_back(item + (t<< m_submatrix_shift));
                }
            }



        }

	 template <typename t_x>
	 void direct_links2(t_x source_id, std::vector<t_x> &result) const {
            result.clear();
            if (source_id > m_max_element){
                return;
            }

            t_x y = source_id >> m_submatrix_shift;
            t_x y_offsetted = source_id - (y << m_submatrix_shift);
            //uint submatrix_size = m_matrix_size/k0;

            //TODO: might be slow to use extra vector as reference, maybe it's better to remove clear() in k2treap and add directly to endresult

            std::vector<std::vector<t_x>> tmp_results;
            #pragma omp parallel num_threads(m_submatrix_per_dim_count)
            {
                #pragma omp single 
                {
                        tmp_results.resize(omp_get_num_threads());
                }

                auto t = omp_get_thread_num();
                m_k2trees[y*m_submatrix_per_dim_count+t].direct_links2(y_offsetted, tmp_results[t]);
            }

            for (size_t t = 0; t < tmp_results.size(); ++t){
                for (auto item : tmp_results[t]){
                        result.push_back(item + (t<< m_submatrix_shift));
                }    
            }
        }


        template<typename t_x>
        void inverse_links2(t_x source_id, std::vector<t_x> &result) const {
            result.clear();

            if (source_id > m_max_element){
                return;
            }
            t_x x = source_id >> m_submatrix_shift;
            t_x x_offsetted = source_id - (x << m_submatrix_shift);
            //uint submatrix_size = m_matrix_size/k0;

            std::vector<std::vector<t_x>> tmp_results;
            #pragma omp parallel num_threads(m_submatrix_per_dim_count)
            {
                #pragma omp single 
                {
                        tmp_results.resize(omp_get_num_threads());
                }

                auto t = omp_get_thread_num();
                m_k2trees[t*m_submatrix_per_dim_count+x].inverse_links2(x_offsetted, tmp_results[t]);
            }

            for (size_t t = 0; t < tmp_results.size(); ++t){
                for (auto item : tmp_results[t]){
                        result.push_back(item + (t<< m_submatrix_shift));
                }    
            }
        }

        template<typename t_x>
        void inverse_links_shortcut(t_x source_id, std::vector<t_x> &result) const {
            result.clear();

            if (source_id > m_max_element){
                return;
            }

            t_x x = source_id >> m_submatrix_shift;
            t_x x_offsetted = source_id - (x << m_submatrix_shift);
            std::vector<std::vector<t_x>> tmp_results;
            #pragma omp parallel num_threads(m_submatrix_per_dim_count)
            {
                #pragma omp single
                {
                    tmp_results.resize(omp_get_num_threads());
                }

                auto t = omp_get_thread_num();
                m_k2trees[t*m_submatrix_per_dim_count+x].inverse_links_shortcut(x_offsetted, tmp_results[t]);
            }

            for (size_t t = 0; t < tmp_results.size(); ++t){
                for (auto item : tmp_results[t]){
                    result.push_back(item + (t<< m_submatrix_shift));
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

            t_x x = link.first >> m_submatrix_shift;
            t_y y = link.second >> m_submatrix_shift;

            t_x x_offsetted = link.first - (x << m_submatrix_shift);
            t_y y_offsetted = link.second - (y << m_submatrix_shift);

            uint index = x*m_submatrix_per_dim_count+y;
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

            t_x x = link.first >> m_submatrix_shift;
            t_y y = link.second >> m_submatrix_shift;

            t_x x_offsetted = link.first - (x << m_submatrix_shift);
            t_y y_offsetted = link.second - (y << m_submatrix_shift);

            uint index = x*m_submatrix_per_dim_count+y;
            return m_k2trees[index].check_link_shortcut(std::make_pair(x_offsetted, y_offsetted));
        }

        //! Move assignment operator
        k2_tree_partitioned &operator=(k2_tree_partitioned &&tr) {
            if (this != &tr) {
                m_size = tr.m_size;
                m_max_element = tr.m_max_element;
                m_used_compression = tr.m_used_compression;
                m_dictionary = std::move(tr.m_dictionary);
                m_access_shortcut_size = tr.m_access_shortcut_size;
                m_k2trees = std::move(tr.m_k2trees);
                m_matrix_dimension = std::move(tr.m_matrix_dimension);
                m_submatrix_per_dim_count = tr.m_submatrix_per_dim_count;
                m_submatrix_shift = tr.m_submatrix_shift;
            }
            return *this;
        }

        //! Assignment operator
        k2_tree_partitioned &operator=(k2_tree_partitioned &tr) {
            if (this != &tr) {
                m_size = tr.m_size;
                m_max_element = tr.m_max_element;
                m_used_compression = tr.m_used_compression;
                m_dictionary = tr.m_dictionary;
                m_access_shortcut_size = tr.m_access_shortcut_size;
                m_k2trees = tr.m_k2trees;
                m_matrix_dimension = tr.m_matrix_dimension;
                m_submatrix_per_dim_count = tr.m_submatrix_per_dim_count;
                m_vocabulary = tr.m_vocabulary;
                m_submatrix_shift = tr.m_submatrix_shift;
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

            if (m_submatrix_per_dim_count != tr.m_submatrix_per_dim_count){
                return false;
            }

            if (m_submatrix_shift != tr.m_submatrix_shift){
                return false;
            }

            if (m_max_element != tr.m_max_element){
                return false;
            }

            if (m_access_shortcut_size != tr.m_access_shortcut_size){
                std::cout << "Shortcut size differs" << std::endl;
                return false;
            }

            if (m_used_compression != tr.m_used_compression){
                std::cout << "Used compression differs" << std::endl;
            }

            if (m_used_compression == LEGACY_DAC) {
                if (!m_vocabulary.get()->operator==(*tr.m_vocabulary.get())) {
                    std::cout << "vocabulary differs" << std::endl;
                    return false;
                }
            } else if ((m_used_compression == DAC)|| m_used_compression == WT_INT_DICT){
                if (!m_dictionary.get()->operator==(*tr.m_dictionary.get())){
                    std::cout << "Dictionary differs" << std::endl;
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
                std::swap(m_used_compression, tr.m_used_compression);
                std::swap(m_dictionary, tr.m_dictionary);
                std::swap(m_access_shortcut_size, tr.m_access_shortcut_size);
                std::swap(m_k2trees, tr.m_k2trees);
                std::swap(m_matrix_dimension, tr.m_matrix_dimension);
                std::swap(m_vocabulary, tr.m_vocabulary);
                std::swap(m_submatrix_per_dim_count, tr.m_submatrix_per_dim_count);
                std::swap(m_submatrix_shift, tr.m_submatrix_shift);
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
            written_bytes += write_member(m_access_shortcut_size, out, child, "m_access_shortcut_size");
            written_bytes += write_member(m_submatrix_per_dim_count, out, child, "m_submatrix_per_dim_count");

            for (uint j = 0; j < m_k2trees.size(); ++j) {
                written_bytes += m_k2trees[j].serialize(out, child, "k2_tree_"+std::to_string(j));
            }

            written_bytes += write_member(m_used_compression, out, child, "used_compression");
            if (m_used_compression == LEGACY_DAC){
                written_bytes += m_vocabulary->serialize(out, child, "voc");
            } else if ((m_used_compression == WT_INT_DICT) || (m_used_compression == DAC)){
                written_bytes += m_dictionary->serialize(out, child, "dict");
            }

            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream &in) {
            read_member(m_size, in);
            read_member(m_matrix_dimension, in);
            read_member(m_max_element, in);
            //for backwards compatibility
            m_submatrix_shift = bits::hi(m_matrix_dimension);
            read_member(m_access_shortcut_size, in);
            read_member(m_submatrix_per_dim_count, in);

            m_k2trees.resize(m_submatrix_per_dim_count*m_submatrix_per_dim_count);
            for (uint j = 0; j < m_submatrix_per_dim_count*m_submatrix_per_dim_count; ++j) {
                m_k2trees[j].load(in);
            }

            read_member(m_used_compression, in);
            if (m_used_compression == LEGACY_DAC){
                m_vocabulary = std::shared_ptr<k2_tree_vocabulary>(new k2_tree_vocabulary());
                m_vocabulary->load(in);

                for (uint j = 0; j < m_submatrix_per_dim_count*m_submatrix_per_dim_count; ++j) {
                    m_k2trees[j].set_vocabulary(m_vocabulary);
                }
            } else if ((m_used_compression == WT_INT_DICT) || (m_used_compression == DAC)){
                m_dictionary = std::shared_ptr<int_vector<>>(new int_vector<>());
                m_dictionary->load(in);

                for (uint j = 0; j < m_submatrix_per_dim_count*m_submatrix_per_dim_count; ++j) {
                    m_k2trees[j].set_dictionary(m_dictionary);
                }
            }
        }

        void compress_leaves(leaf_compression_type compression, uint64_t hash_size = 0){
            switch (compression) {
                case UNCOMPRESSED:
                    break;

                case DAC:
                    dac_compress();
                    break;

                case DAC_MMAPS:
                    dac_compress(true);
                    break;

                case LEGACY_DAC:
                    legacy_dac_compress(hash_size);
                    break;

                case WT_INT_DICT:
                    wt_huff_int_dict_compress();
                    break;

                case WT_INT:
                    wt_huff_int_compress();
                    break;

                default:
                    throw new std::runtime_error("invalid value for compression method");
            }
        }

        /**
        * Return the number of bytes necessary to store a word.
        *
        * @return Size of a word.
        */
        virtual inline uint word_size() const {
            return m_k2trees[0].word_size();
        }

        /**
        * Returns the number of words of \f$k_leaves^2\f$ bits in the leaf level.
        *
        * @return Number of words.
        */
        virtual inline size_t words_count() const{
            size_t words_count = 0;
            for (size_t i = 0; i < m_k2trees.size(); ++i){
                words_count += m_k2trees[i].words_count();
                //std::cout << "Words in tree " << i << ": " <<m_k2trees[i].words_count() << std::endl;
            }
            return words_count;
        }

        void dac_compress(bool use_multiple_maps = false){
            if (is_compressed()){
                return;
            }
            using timer = std::chrono::high_resolution_clock;
            //typedef decltype(edges[0].second) t_y;

            auto start = timer::now();
            int_vector<> leaf_words;
            words(leaf_words);
            auto stop = timer::now();
            std::cout << "Word Iteration Finished" << std::endl;
            word_iteration += duration_cast<milliseconds>(stop-start).count();
            start = timer::now();
            std::unordered_map<int_vector<>::value_type, uint64_t> codeword_map; //maps word w to codeword c (the code word is chosen based on the frequency of word w ("huffman"))
            if (use_multiple_maps){
                frequency_encode_using_multiple_maps(leaf_words, m_dictionary, codeword_map);
            } else {
                frequency_encode_using_sort(leaf_words, m_dictionary, codeword_map);
            }
            std::cout << "Frequency Encoding Finished" << std::endl;

            stop = timer::now();
            frequency_encoding += duration_cast<milliseconds>(stop-start).count();

            start = timer::now();
            std::cout << "Performing dac compression of subtrees" << std::endl;
            #pragma omp parallel for
            for (size_t i = 0; i < m_k2trees.size(); ++i){
                m_k2trees[i].dac_compress(codeword_map, m_dictionary);//compress using shared vocabulary
            }
            m_used_compression = DAC;

            stop = timer::now();
            dac_compression += duration_cast<milliseconds>(stop-start).count();
        }

        void wt_huff_int_compress(bool per_tree=true){
            if (is_compressed()){
                return;
            }

            if (per_tree){
                #pragma omp parallel for
                for (size_t i = 0; i < m_k2trees.size(); ++i){
                    m_k2trees[i].wt_huff_int_compress();//compress using shared vocabulary
                }
            } else {

                throw new std::runtime_error("not yet implemented, needs access operations and direct storage here");
                /*const int_vector<> leaf_words;
                words(leaf_words);
                wt_huff_int_compress(leaf_words, m_leaves_wt);
                m_is_wt_direct_comp = true;

                #pragma omp parallel for
                for (uint i = 0; i < m_k2trees.size(); ++i){
                    m_k2trees[i].clear_leaves();//compress using shared vocabulary
                }
                */

            }
            m_used_compression = WT_INT;
        }

        void wt_huff_int_dict_compress(){
            if (is_compressed()){
                return;
            }

            int_vector<> leaf_words;
            words(leaf_words);
            std::unordered_map<int_vector<>::value_type, uint> codeword_map;
            frequency_encode(leaf_words, m_dictionary, codeword_map);

            int_vector<>().swap(leaf_words);//save some memory

            #pragma omp parallel for
            for (size_t i = 0; i < m_k2trees.size(); ++i){
                m_k2trees[i].wt_huff_int_dict_compress(codeword_map, m_dictionary);//compress using shared vocabulary
            }

            m_used_compression = WT_INT_DICT;
        }


        void legacy_dac_compress(uint64_t hash_size = 0) {
            if (is_compressed()){
                return;
            }

            std::cout << "Compressing leaves using legacy dac compression" << std::endl;

            m_used_compression = LEGACY_DAC;
            std::vector<uchar> leaf_words;
            auto start = timer::now();
            words(leaf_words);
            auto stop = timer::now();
            word_iteration += duration_cast<milliseconds>(stop-start).count();

            start = timer::now();
            HashTable codeword_map;
            frequency_encode(leaf_words, word_size(), words_count(), m_vocabulary, codeword_map, hash_size);
            std::cout << "Frequency encoding finished" << std::endl;
            stop = timer::now();
            frequency_encoding += duration_cast<milliseconds>(stop-start).count();

            start = timer::now();
            #pragma omp parallel for
            for (size_t i = 0; i < m_k2trees.size(); ++i){
                m_k2trees[i].legacy_dac_compress(codeword_map, m_vocabulary, false);//compress using shared vocabulary
            }

            stop = timer::now();
            dac_compression += duration_cast<milliseconds>(stop-start).count();
        }

        std::string get_type_string() const {
            return "k2_tree_partitioned<"+std::to_string(m_submatrix_per_dim_count)+","+get_compression_name(m_used_compression)+","+m_k2trees[0].get_type_string()+">";
        }

        /**
   * Returns all the words on the leaf level in a vector, which has size word_size*word_count (word_size is the amount of byte-sized words needed per leaf)
   */
        void words(std::vector<uchar>& result) const {
            result.clear();
            result.resize(words_count()*word_size());

            std::vector<uint64_t> offset(m_k2trees.size());
            offset[0] = 0;
            for (uint i = 1; i < m_k2trees.size(); ++i){
                offset[i] = offset[i-1] + m_k2trees[i-1].words_count()*word_size();
            }

            #pragma omp parallel for
            for (size_t i = 0; i < m_k2trees.size(); ++i) {
                m_k2trees[i].words(result, true, offset[i]);
            }
        }

        void words(int_vector<>& result) const {
            result = int_vector<>(words_count(), 0, word_size()*8);

            std::vector<uint64_t> offset(m_k2trees.size());
            offset[0] = 0;
            for (uint i = 1; i < m_k2trees.size(); ++i){
                offset[i] = offset[i-1] + m_k2trees[i-1].words_count();
            }

            if (word_size() == 8) { //only do it in parallel for 64 bit words for now to avoid writing to the same 64 bit block
                #pragma omp parallel for
                for (size_t i = 0; i < m_k2trees.size(); ++i) {
                    m_k2trees[i].words(result, true, offset[i]);
                }
            } else {
                for (size_t i = 0; i < m_k2trees.size(); ++i) {
                    m_k2trees[i].words(result, true, offset[i]);
                }
            }
        }

        void construct_access_shortcut(uint8_t access_shortcut_size) {
            if (access_shortcut_size == 0) return ;

	    m_access_shortcut_size = access_shortcut_size;
            for (size_t i = 0; i < m_k2trees.size(); ++i) {
                m_k2trees[i].construct_access_shortcut(access_shortcut_size);
            }
        }


        void load_from_ladrabin(std::string fileName, construction_algorithm  construction_algo = COUNTING_SORT,
                                        std::string temp_file_prefix = "",
                                        uint8_t access_shortcut_size = 0, uint_fast8_t submatrix_shift  = 0, bool serialize_during_construction = true){

            load_from_ladrabin_internal(empty<uint32_t>(), fileName, construction_algo, submatrix_shift, temp_file_prefix,
                                        access_shortcut_size, serialize_during_construction);
        }

        void load_from_ladrabin64(std::string fileName, construction_algorithm  construction_algo = COUNTING_SORT,
                                uint_fast8_t submatrix_shift  = 0, std::string temp_file_prefix = "",
                                uint8_t access_shortcut_size = 0, bool serialize_during_construction = true){

            load_from_ladrabin_internal(empty<uint64_t>(), fileName, construction_algo, submatrix_shift, temp_file_prefix,
                                        access_shortcut_size, serialize_during_construction);
        }

    private:



        template <typename T>
        void load_from_ladrabin_internal(empty<T>, std::string fileName,
                                         construction_algorithm construction_algo = COUNTING_SORT,
                                         uint_fast8_t submatrix_shift = 0, std::string temp_file_prefix = "",
                                         uint8_t access_shortcut_size = 0, bool serialize_during_construction = true){
            using namespace k2_tree_ns;
            typedef typename std::make_signed<T>::type S_T;


            if(!has_ending(fileName, ".ladrabin")){
                fileName.append(".ladrabin");
                std::cout << "Appending .ladrabin to filename as file has to be in .ladrabin format" << std::endl;
            }

            m_access_shortcut_size = access_shortcut_size;

            std::fstream fileStream(fileName, std::ios_base::in);

            if (fileStream.is_open()){
                T number_of_nodes;
                ulong number_of_edges;

                read_member(number_of_nodes, fileStream);
                read_member(number_of_edges, fileStream);

                m_max_element = number_of_nodes -1;
                m_size = number_of_edges;

                m_submatrix_shift = calculateSubmatrixShift(submatrix_shift);

                calculate_matrix_dimension_and_submatrix_count();

                T nodes_read = 0;
                T source_id;
                T source_id_offsetted = 0;
                S_T target_id;

                typedef std::pair<uint, uint> t_e;
                //typedef typename stxxl::VECTOR_GENERATOR<t_e>::result stxxl_pair_vector;

                //FIXME replace with stxxl vector
                std::vector<std::vector<t_e>> buffers(m_submatrix_per_dim_count);
                //std::vector<std::vector<t_e>> buffers;
                std::vector<uint> maximum_in_buffer(m_submatrix_per_dim_count, m_matrix_dimension-1);
                m_k2trees.resize(m_submatrix_per_dim_count*m_submatrix_per_dim_count);

                uint current_matrix_row = 0;
                for (uint64_t i = 0; i < number_of_nodes + number_of_edges; i++) {
                    read_member(target_id, fileStream);
                    if (target_id < 0) {

                        source_id = nodes_read;
                        uint corresponding_row = nodes_read >> m_submatrix_shift;
                        source_id_offsetted = source_id - (corresponding_row << m_submatrix_shift); //same as source_id % m_matrix_dimension/ source_id % (1Ull << m_submatrix_shift)
                        if (corresponding_row > current_matrix_row){
                            construct_trees_from_buffers(current_matrix_row, construction_algo, temp_file_prefix,
                                                         buffers, maximum_in_buffer, serialize_during_construction, fileName);

                            //in case of a complete empty row

                            for (uint k = 0; k < (corresponding_row - current_matrix_row - 1); ++k) {
                                current_matrix_row++;
                                construct_trees_from_buffers(current_matrix_row, construction_algo, temp_file_prefix,
                                                             buffers, maximum_in_buffer, serialize_during_construction, fileName);
                            }

                            current_matrix_row = corresponding_row;

                        }
                        nodes_read++;
                    } else {
                        uint column_in_matrix = target_id >> m_submatrix_shift;
                        auto target_id_offsetted = target_id - (column_in_matrix << m_submatrix_shift);
                        buffers[column_in_matrix].push_back(std::make_pair(source_id_offsetted, target_id_offsetted));
                    }
                }

                //cover leftovers
                construct_trees_from_buffers(current_matrix_row, construction_algo, temp_file_prefix,
                                             buffers, maximum_in_buffer, serialize_during_construction, fileName);

                if (serialize_during_construction){
                    for (uint i = 0; i < m_k2trees.size(); i++){
                        load_from_file(m_k2trees[i], fileName+"tr"+std::to_string(i)+".sdsl");
                        sdsl::remove(fileName+"tr"+std::to_string(i)+".sdsl");
                    }
                }

		construct_access_shortcut(m_access_shortcut_size);
            } else {
                throw std::runtime_error("Could not open file to load ladrabin graph");
            }
        }



        //for hyperlinkgraph
        void load_from_binary_edgelist(std::string fileName,
                                         construction_algorithm construction_algo = COUNTING_SORT,
                                         uint_fast8_t submatrix_shift = 0, std::string temp_file_prefix = "",
                                         uint8_t access_shortcut_size = 0, bool serialize_during_construction = true){
            using namespace k2_tree_ns;
            m_access_shortcut_size = access_shortcut_size;

            std::fstream fileStream(fileName, std::ios_base::in);

            if (fileStream.is_open()){
                size_t number_of_edges;
                read_member(number_of_edges, fileStream);
                uint32_t number_of_nodes = 3563666998;
                m_max_element = number_of_nodes -1;
                m_size = number_of_edges;

                m_submatrix_shift = calculateSubmatrixShift(submatrix_shift);

                calculate_matrix_dimension_and_submatrix_count();

                uint64_t edges_read = 0;
                uint32_t source_id;
                uint32_t source_id_offsetted = 0;
                uint32_t target_id;

                typedef std::pair<uint, uint> t_e;
                //typedef typename stxxl::VECTOR_GENERATOR<t_e>::result stxxl_pair_vector;

                //FIXME replace with stxxl vector
                std::vector<std::vector<t_e>> buffers(m_submatrix_per_dim_count);
                //std::vector<std::vector<t_e>> buffers;
                std::vector<uint> maximum_in_buffer(m_submatrix_per_dim_count, m_matrix_dimension - 1);
                m_k2trees.resize(m_submatrix_per_dim_count * m_submatrix_per_dim_count);

                uint current_matrix_row = 0;
                for (uint64_t i = 0; i < number_of_edges; i++) {
                    read_member(source_id, fileStream);
                    read_member(target_id, fileStream);

                    uint corresponding_row = source_id >> m_submatrix_shift;
                    source_id_offsetted = source_id - (corresponding_row << m_submatrix_shift); //same as source_id % m_matrix_dimension/ source_id % (1Ull << m_submatrix_shift)
                    if (corresponding_row > current_matrix_row) {
                        construct_trees_from_buffers(current_matrix_row, construction_algo, temp_file_prefix,
                                                     buffers, maximum_in_buffer, serialize_during_construction,
                                                     fileName);

                        //in case of a complete empty row

                        for (uint k = 0; k < ((int)corresponding_row - current_matrix_row - 1); ++k) {
                            current_matrix_row++;
                            construct_trees_from_buffers(current_matrix_row, construction_algo, temp_file_prefix,
                                                         buffers, maximum_in_buffer, serialize_during_construction,
                                                         fileName);
                        }

                        current_matrix_row = corresponding_row;
                    }

                    uint column_in_matrix = target_id >> m_submatrix_shift;
                    auto target_id_offsetted = target_id - (column_in_matrix << m_submatrix_shift);
                    buffers[column_in_matrix].push_back(std::make_pair(source_id_offsetted, target_id_offsetted));
                }

                //cover leftovers
                construct_trees_from_buffers(current_matrix_row, construction_algo, temp_file_prefix,
                                             buffers, maximum_in_buffer, serialize_during_construction, fileName);

                if (serialize_during_construction) {
                    for (uint i = 0; i < m_k2trees.size(); i++) {
                        load_from_file(m_k2trees[i], fileName + "tr" + std::to_string(i) + ".sdsl");
                        sdsl::remove(fileName + "tr" + std::to_string(i) + ".sdsl");
                    }
                }

                construct_access_shortcut(m_access_shortcut_size);
            } else {
                throw std::runtime_error("Could not open file to load ladrabin graph");
            }
        }


        template<typename t_vector>
        void build_k2_trees(t_vector &links, std::string temp_file_prefix = "", construction_algorithm construction_algo = COUNTING_SORT) {
            using namespace k2_tree_ns;
            typedef decltype(links[0].first) t_x;
            typedef decltype(links[0].second) t_y;
            using t_e = std::pair<t_x, t_y>;

            m_size = links.size();

            //typedef typename stxxl::VECTOR_GENERATOR<t_e>::result stxxl_pair_vector;

            //FIXME replace with stxxl vector
            std::vector<std::vector<t_e>> buffers(m_submatrix_per_dim_count*m_submatrix_per_dim_count);
            std::vector<uint64_t> maximum_in_buffer;
            maximum_in_buffer.resize(m_submatrix_per_dim_count*m_submatrix_per_dim_count);

            {
                for (size_t j = 0; j < links.size(); ++j) {
                    auto x = links[j].first;
                    auto y = links[j].second;
                    auto p1 = x >> m_submatrix_shift;
                    auto p2 = y >> m_submatrix_shift;
                    auto corresponding_matrix = p1 * m_submatrix_per_dim_count + p2;
                    x = x - (p1 << m_submatrix_shift);
                    y = y - (p2 << m_submatrix_shift);
                    if (x > maximum_in_buffer[corresponding_matrix]){
                        maximum_in_buffer[corresponding_matrix] = x;
                    }
                    if (y > maximum_in_buffer[corresponding_matrix]){
                        maximum_in_buffer[corresponding_matrix] = y;
                    }
                    buffers[corresponding_matrix].push_back(t_e(x,y));
                }
            }

            m_k2trees.reserve(buffers.size());
            for (size_t l = 0; l < buffers.size(); ++l) {
//                const subk2_tree k2tree(temp_file_prefix, use_counting_sort, buffers[l]);
                m_k2trees.emplace_back(buffers[l], construction_algo, maximum_in_buffer[l], temp_file_prefix);
                buffers[l].clear();
            }
        }

    private:
        inline bool is_compressed(){
            if (m_used_compression == LEGACY_DAC || m_used_compression == DAC || m_used_compression == WT_INT|| m_used_compression == WT_INT_DICT){
                std::cout << "Already compressed, aborting";
                return true;
            }

            return false;
        }

        void calculate_matrix_dimension_and_submatrix_count() {
            m_matrix_dimension = 1ULL << m_submatrix_shift;
            m_submatrix_per_dim_count = (m_max_element+m_matrix_dimension-1) >> m_submatrix_shift;
            std::cout << "Submatrix shift: " << std::to_string(m_submatrix_shift) << std::endl;
            std::cout << "Matrix dimension: " << m_matrix_dimension << std::endl;
            if (m_submatrix_per_dim_count == 0){
                std::cout << "Please choose a smaller Partition size as " << std::to_string(m_submatrix_shift) << std::endl;
                std::cout << "m_submatrix_shift leads to only one partition as the maximum element is " << m_max_element << std::endl;
                m_submatrix_per_dim_count = 1;
            }
            std::cout << "Submatrix amount per row: " << std::to_string(m_submatrix_per_dim_count ) << std::endl;
        }

        template <typename t_vector, typename t_x>
        inline void
        construct_trees_from_buffers(uint current_matrix_row, construction_algorithm construction_algo, std::string &temp_file_prefix,
                                     std::vector<t_vector> &buffers, std::vector<t_x>& maximum_in_buffer, bool serialize, std::string& filename_prefix) {

            auto start = timer::now();

            //#pragma omp parallel for
            for (uint j = 0; j < m_submatrix_per_dim_count; ++j) {
                //std::cout << "Constructing tree "<< current_matrix_row * m_submatrix_per_dim_count + j << std::endl;
                /*if (buffers[j].size() != 0) {
                    std::cout << "Size of " << current_matrix_row * m_submatrix_per_dim_count + j << ": "
                              << buffers[j].size() * 64 / 8 / 1024 << "kByte" << std::endl;
                }*/
                subk2_tree tree(buffers[j], construction_algo, maximum_in_buffer[j], temp_file_prefix);

                if (serialize){
                    store_to_file(tree, filename_prefix+"tr"+std::to_string(current_matrix_row*m_submatrix_per_dim_count+j)+".sdsl");
                    subk2_tree().swap(tree);
                } else {
                    m_k2trees[current_matrix_row*m_submatrix_per_dim_count+j].swap(tree);
                }

                buffers[j].clear();
                //maximum_in_buffer[j] = 0;
                //std::cout << "Assigning tree " << current_matrix_row * m_submatrix_per_dim_count + j << std::endl;
            }

            //little hack: swap buffers along the diagonal to save memory allocation as most of the data is along the diagonal, it is best to the same vector along the diagonal
            if (current_matrix_row < (m_submatrix_per_dim_count-1)){
                std::swap(buffers[current_matrix_row], buffers[current_matrix_row+1]);
            }

            auto stop = timer::now();
            subtree_construction_duration += duration_cast<milliseconds>(stop - start).count();
        }

    };
}

#endif
