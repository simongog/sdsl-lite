
#ifndef INCLUDE_COMPRESSION_COMPRESSOR_H_
#define INCLUDE_COMPRESSION_COMPRESSOR_H_

#include <sdsl/k2_tree_dacs.hpp>
#include <sdsl/k2_tree_hash_table.hpp>
#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>
#include <climits>
#include <parallel/algorithm>
#include "construct.hpp"
#include "wavelet_trees.hpp"
#include <chrono>

/**
 * This class implements several compression methods for the leaf level of k2 trees.
 * This is only intended to be used directly in one of the k2 tree implementations e.g. k2_tree_base, k2_tree_partitioned
 */

namespace sdsl {

    enum leaf_compression_type {UNCOMPRESSED, LEGACY_DAC, DAC, WT_INT, WT_INT_DICT};

    std::string get_compression_name(leaf_compression_type compression_type) {
        switch (compression_type) {
            case UNCOMPRESSED:
                return "uncompressed";
            case LEGACY_DAC:
                return "legacy_dac";
            case DAC:
                return "dac";
            case WT_INT:
                return "wt_int";
            case WT_INT_DICT:
                return "wt_int_dict";
            default:
                return "unknown";
        }
    }

    template <typename t_map>
    void frequency_encode(const int_vector<> &leaf_words, std::shared_ptr<int_vector<>>& dictionary,
                     t_map &codeword_map);

    template <typename t_map>
    void
    frequency_encode_using_sort( int_vector<> &leaf_words, std::shared_ptr<int_vector<>>& dictionary,
                                t_map &codeword_map);

    template <typename t_map>
    void merge_maps(t_map& first, t_map& second);

    template <typename t_map>
    void
    frequency_encode_using_multiple_maps( int_vector<> &leaf_words, std::shared_ptr<int_vector<>>& dictionary,
                                 t_map &codeword_map);

    void frequency_encode(const std::vector<uchar> &leaf_words, const uint word_size, const size_t word_count, std::shared_ptr<k2_tree_vocabulary>& voc, HashTable& table, uint64_t hash_size);

    template <typename t_map>
    void construct_codewords(const int_vector<> &leaf_vector,
                        const t_map& codeword_map,
                        int_vector<> &codewords);

    void construct_codewords(const std::vector<uchar>& leaf_words, const size_t word_size, const size_t words_count, const HashTable& table, std::vector<uint>& codewords);

    inline void construct_legacy_dac(std::vector<uint>& codewords, size_t word_count,  uint64_t voc_size, k2_tree_dac& comp_leaves);

    void perform_dac_compression(const int_vector<> &leaf_words, std::shared_ptr<int_vector<>>& dictionary,
                                 dac_vector<> &compressed_leaves) {
        std::unordered_map<int_vector<>::value_type, uint> codeword_map;
        frequency_encode(leaf_words, dictionary, codeword_map);
        int_vector<> codewords; //size is known: bits:hi for hashmap.size()? or distinct_values.size()
        construct_codewords(leaf_words, codeword_map, codewords);
        dac_vector<> tmp(codewords);
        tmp.swap(compressed_leaves);
    }

    template <typename t_map>
    void perform_dac_compression_with_shared_vocabulary(const int_vector<> &leaf_words,
                                                        t_map& codeword_map,
                                                        dac_vector<> &compressed_leaves) {
        int_vector<> codewords; //size is known: bits:hi for hashmap.size()? or distinct_values.size()
        construct_codewords(leaf_words, codeword_map, codewords);
        dac_vector<> tmp(codewords);
        tmp.swap(compressed_leaves);
    }

    template<typename wt>
    void perform_wt_huff_int_dict_compression(const int_vector<> &leaf_words, std::shared_ptr<int_vector<>>& dictionary,
                                              wt &compressed_leaves) {
        std::unordered_map<int_vector<>::value_type, uint> codeword_map;
        frequency_encode(leaf_words, dictionary, codeword_map);
        int_vector<> codewords; //size is known: bits:hi for hashmap.size()? or distinct_values.size()
        construct_codewords(leaf_words, codeword_map, codewords);
        wt tmp;
        construct_im(tmp, codewords);
        tmp.swap(compressed_leaves);
    }

    template<typename wt>
    void perform_wt_huff_int_shared_voc_dict_compression(const int_vector<>& leaf_words,
                                                         std::unordered_map<int_vector<>::value_type, uint>& codeword_map,
                                                         wt &compressed_leaves) {
        int_vector<> codewords; //size is known: bits:hi for hashmap.size()? or distinct_values.size()
        construct_codewords(leaf_words, codeword_map, codewords);
        wt tmp;
        construct_im(tmp, codewords);
        tmp.swap(compressed_leaves);
    }

    void perfrom_legacy_dac_compression(const std::vector<uchar> &leaf_words, const uint word_size,
                                        const size_t word_count, std::shared_ptr<k2_tree_vocabulary>& voc,
                                        k2_tree_dac &compressed_leafs, uint64_t hash_size = 0) {
        HashTable hashtable;
        frequency_encode(leaf_words, word_size, word_count, voc, hashtable, hash_size);
        std::vector<uint> codewords;
        construct_codewords(leaf_words, word_size, word_count, hashtable, codewords);
        construct_legacy_dac(codewords, word_count, voc->word_count(), compressed_leafs);
    }

    void perfdorm_legacy_dac_compress_with_shared_vocabulary(const HashTable& hashtable,
                                                             const std::vector<uchar>& leaf_words,
                                                             const uint word_size, const size_t word_count,
                                                             size_t voc_size,
                                                             k2_tree_dac &compressed_leafs) {
        std::vector<uint> codewords;
        construct_codewords(leaf_words, word_size,  word_count, hashtable, codewords);
        construct_legacy_dac(codewords, word_count, voc_size, compressed_leafs);
    }


    inline void construct_legacy_dac(std::vector<uint>& codewords, size_t word_count,  uint64_t voc_size, k2_tree_dac& comp_leaves){
        try {
            k2_tree_dac tmp(codewords, word_count, voc_size);
            tmp.swap(comp_leaves);
        } catch (...) {
            std::cerr << "[k2_tree_compressor::construct_legacy_dac] Error: Could not create k2_tree_dac\n";
            exit(1);
        }
    }

    template <typename t_map>
    void
    frequency_encode_using_multiple_maps(int_vector<> &leaf_words, std::shared_ptr<int_vector<>>& dictionary,
                                          t_map &codeword_map){
        std::vector<t_map> codeword_maps;
        uint num_threads = 0;
        std::string tmp_file = ram_file_name(util::to_string(util::pid()) + "_" + util::to_string(util::id()));
        int_vector_buffer<> dictionary_buffer(tmp_file, std::ios::out);

        #pragma omp parallel shared(leaf_words, num_threads)
        {
            #pragma omp single
            {
                num_threads = omp_get_num_threads();
                codeword_maps.resize(num_threads);
            }

            auto thread_num = omp_get_thread_num();

            #pragma omp for
            for (size_t i = 0; i < leaf_words.size(); ++i){
                auto it = codeword_maps[thread_num].find(leaf_words[i]);
                if (it != codeword_maps[thread_num].end()) {
                    it->second += 1;
                } else {
                    codeword_maps[thread_num].insert(std::make_pair(leaf_words[i], (uint) 1));//consider using sparsehash
                }
            }

            //perform some kind of multiway merge

            uint mod = 2;
            uint index = 1;
            while (mod < num_threads){ //leave 2
                if (thread_num % mod == 0){
                    //if next codeword_map present
                    if ((thread_num + index) < num_threads){
                        merge_maps(codeword_maps[thread_num], codeword_maps[thread_num+index]);
                    }
                }

                mod *= 2;
                index *=2;
                #pragma omp barrier
            }

            #pragma omp single
            {
                merge_maps(codeword_maps[0], codeword_maps[index]);
                codeword_map.swap(codeword_maps[0]);
            }
        }

        int_vector<>* tmp = new int_vector<>();
        dictionary = std::shared_ptr<int_vector<>>(tmp);
        dictionary->resize(codeword_map.size());

        size_t counter = 0;
        for (auto item: codeword_map){
            dictionary->operator[](counter) = item.first;
            counter++;
        }

/*
        using namespace std::chrono;
        using timer = std::chrono::high_resolution_clock;
        //typedef decltype(edges[0].second) t_y;

        auto start = timer::now();
        for (size_t i = 0; i < codeword_maps.size(); i++){
            for (auto& item : codeword_maps[i]){
                auto it = codeword_map.find(item.first);
                if (it != codeword_map.end()) {
                    it->second += item.second;
                } else {
                    codeword_map.insert(std::make_pair(item.first, item.second));//consider using sparsehash
                    dictionary_buffer.push_back(leaf_words[i]);
                }
            }
            t_map().swap(codeword_maps[i]);
        }
        auto stop = timer::now();
        auto duration = duration_cast<milliseconds>(stop-start).count();
        std::cout << "Duration of megrge step: " << duration << std::endl;
*/
        /*
        dictionary_buffer.close();
        int_vector<>* tmp = new int_vector<>();
        load_from_file(*tmp, tmp_file);
        dictionary = std::shared_ptr<int_vector<>>(tmp);
        remove(tmp_file);
         */
        // Sort words by frequency
        __gnu_parallel::sort(dictionary->begin(), dictionary->end(),
                             [&](const int_vector<>::value_type a, const int_vector<>::value_type b) {
                                 return codeword_map.find(a)->second > codeword_map.find(b)->second;
                             });


        if (dictionary->size() > UINT_MAX) {
            std::cerr << "k2_tree_vocabulary size is bigger than 32 Bit, things might break! (have fun anyway"
                      << std::endl;
        }

        #pragma omp parallel for
        for (size_t i = 0; i < dictionary->size(); ++i) {
            auto it = codeword_map.find(dictionary->operator[](i));
            it->second = i;
        }
    }

    template <typename t_map>
    void merge_maps(t_map& first, t_map& second){
        //small possible optimization: swap to insert smaller one
        if (second.size() > first.size()){
            std::swap(first, second);
        }

        for (auto& item : second){
            auto it = first.find(item.first);
            if (it != first.end()) {
                it->second += item.second;
            } else {
                first.insert(std::make_pair(item.first, item.second));//consider using sparsehash
            }
        }
        t_map().swap(second);
    }

    template <typename t_map>
    void
    frequency_encode_using_sort(int_vector<> &leaf_words, std::shared_ptr<int_vector<>>& dictionary,
                     t_map &codeword_map) {

        __gnu_parallel::sort(leaf_words.begin(), leaf_words.end());
        typedef std::vector<std::vector<std::pair<int_vector<>::value_type, uint64_t>>> pair_vec_2_dim;
        pair_vec_2_dim occurrence_count;
        uint num_threads = 0;
        std::vector<uint64_t> intervals;
        std::vector<uint64_t> offsets;
        std::vector<std::pair<int_vector<>::value_type, uint64_t>> merged_occurrences;

        #pragma omp parallel shared(num_threads, leaf_words)
        {
            #pragma omp single
            {
                num_threads = omp_get_num_threads();
                offsets.resize(num_threads);
                occurrence_count.resize(num_threads);
                intervals.resize(num_threads+1);
                intervals[0] = 0;
                for (uint i = 0; i < num_threads; i++){
                    intervals[i] = leaf_words.size()/num_threads * i;
                }
                intervals[num_threads] = leaf_words.size();
            }
            auto thread_num = omp_get_thread_num();

            auto start_index = intervals[thread_num];
            int_vector<>::value_type previous = leaf_words[start_index];
            uint64_t occurence_counter = 0;


            if (!start_index == 0){
                //skip first words (covered by previous thread
                while (leaf_words[start_index] == previous){
                    start_index++;
                }
                previous = leaf_words[start_index];
            }

            uint64_t index;
            for (index  = intervals[thread_num]; index < intervals[thread_num+1]; ++index){
                if (leaf_words[index] == previous){
                    occurence_counter++;
                } else {
                    occurrence_count[thread_num].push_back(std::make_pair(previous, occurence_counter));
                    previous = leaf_words[index];
                    occurence_counter = 0;
                }
            }

            //read first node of next threads memory
            while (index < leaf_words.size() && leaf_words[index] == previous){
                index++;
                occurence_counter++;
            }
            occurrence_count[thread_num].push_back(std::make_pair(previous, occurence_counter));

            if (index < leaf_words.size() && index < intervals[thread_num+1]+1){
               //also cover next word as it directly starts at intervals[thread_num+1] and was therefore skipped by the following thread
                previous = leaf_words[index];
                while (leaf_words[index] == previous){
                    index++;
                    occurence_counter++;
                }
                occurrence_count[thread_num].push_back(std::make_pair(previous, occurence_counter));
            }

            #pragma omp barrier
            #pragma omp single
            {
                offsets[0] = 0;
                for (uint i = 0; i < num_threads-1; ++i){
                    offsets[i+1] = offsets[i] + occurrence_count[i].size();
                }

                merged_occurrences.resize(offsets[num_threads-1]+occurrence_count[num_threads-1].size());
            }

            std::copy(occurrence_count[thread_num].begin(), occurrence_count[thread_num].end(), merged_occurrences.begin()+offsets[thread_num]);
        }
        pair_vec_2_dim().swap(occurrence_count);//free memory

        __gnu_parallel::sort(merged_occurrences.begin(), merged_occurrences.end(), [&](const std::pair<int_vector<>::value_type, uint> a, const std::pair<int_vector<>::value_type, uint> b) {
            return a.second > b.second;
        });

        int_vector<>* tmp = new int_vector<>();
        dictionary = std::shared_ptr<int_vector<>>(tmp);
        dictionary->resize(merged_occurrences.size());

        #pragma omp parallel for
        for (size_t i = 0; i < merged_occurrences.size(); i++){
            dictionary->operator[](i) = merged_occurrences[i].first;
        }

        codeword_map.reserve(merged_occurrences.size());
        for (size_t i = 0; i < merged_occurrences.size(); i++){
            codeword_map.insert(std::make_pair(merged_occurrences[i].first, i));
        }
    }

        /**
        * Creates new k2tree with the leaves compressed.
        * @param tree
        * @param build
        */
    template <typename t_map>
    void
    frequency_encode(const int_vector<> &leaf_words, std::shared_ptr<int_vector<>>& dictionary,
                     t_map &codeword_map) {

        std::string tmp_file = ram_file_name(util::to_string(util::pid()) + "_" + util::to_string(util::id()));
        int_vector_buffer<> dictionary_buffer(tmp_file, std::ios::out);
        for (size_t i = 0; i < leaf_words.size(); ++i) {
            auto it = codeword_map.find(leaf_words[i]);
            if (it != codeword_map.end()) {
                it->second += 1;
            } else {
                codeword_map.insert(std::make_pair(leaf_words[i], (uint) 1));//consider using sparsehash
                dictionary_buffer.push_back(leaf_words[i]);
            }
        }

        dictionary_buffer.close();
        int_vector<>* tmp = new int_vector<>();
        load_from_file(*tmp, tmp_file);
        dictionary = std::shared_ptr<int_vector<>>(tmp);
        remove(tmp_file);

        /* alternative implementation instead of distinct values buffer */
        /*std::vector<std::pair<unsigned long, uint>> pairs;
        for (auto itr = hashmap.begin(); itr != hashmap.end(); ++itr)
            pairs.push_back(*itr);

        sort(pairs.begin(), pairs.end(), [=](std::pair<int, int>& a, std::pair<int, int>& b)
         {
             return a.second < b.second;
         }
        );
         */

        // Sort words by frequency
        __gnu_parallel::sort(dictionary->begin(), dictionary->end(),
                             [&](const int_vector<>::value_type a, const int_vector<>::value_type b) {
                                 return codeword_map.find(a)->second > codeword_map.find(b)->second;
                             });


        if (dictionary->size() > UINT_MAX) {
            std::cerr << "k2_tree_vocabulary size is bigger than 32 Bit, things might break! (have fun anyway)"
                      << std::endl;
        }

        #pragma omp parallel for
        for (size_t i = 0; i < dictionary->size(); ++i) {
            auto it = codeword_map.find(dictionary->operator[](i));
            it->second = i;
        }
    }


    /* The "legacy dac encoding" is an adapted version of nlehmann, which is under the beer-ware license:
    * ----------------------------------------------------------------------------
    * "THE BEER-WARE LICENSE" (Revision 42):
    * <nlehmann@dcc.uchile.cl> wrote this file. As long as you retain this notice you
    * can do whatever you want with this stuff. If we meet some day, and you think
    * this stuff is worth it, you can buy me a beer in return Nicolás Lehmann
    * ----------------------------------------------------------------------------
    */

    void frequency_encode(const std::vector<uchar> &leaf_words, const uint word_size, const size_t word_count, std::shared_ptr<k2_tree_vocabulary>& voc, HashTable& table, uint64_t hash_size = 0){
        try {
            // Insert words in hash
            if (hash_size == 0) {

                std::cerr
                        << "[k2_tree_compressor::FreqVoc] Warning: Hash Size not specified, it will thus be automatically determined by amount of distinct values, which is really slow!"
                        << std::endl;

                k2_tree_vocabulary words(word_count, word_size);

                for (size_t i = 0; i < word_count; ++i) {
                    words.assign(i, leaf_words, i * word_size);
                }

                //FIXME: this seems to be pretty bad, either use a parallel sort or some other technique to find amount of uniqe values, the words dont have to be sorted later on
                //the vocabulary is rebuild anyway
                // Count number of different words
                // We hope there are many repetead words. We need to encode each word in
                // a 32-bit integer.
                size_t res = words.sort();
                if (res > INT_MAX) {
                    std::cerr << "[k2_tree_compressor::FreqVoc]  Too many different words ";
                    std::cerr << "in the vocabulary\n";
                    exit(1);
                }
                hash_size = (uint) res;

            }

            std::cout << "Before putting words in hash" << std::endl;
            HashTable tmp(hash_size);
            table = tmp;
            std::vector<size_t> posInHash;
            posInHash.reserve(hash_size);
            size_t addr;
            for (size_t i = 0; i < word_count; ++i) {
                if (!table.search(&leaf_words[i * word_size], word_size, &addr)) {
                    table.add(&leaf_words[i * word_size], word_size, addr);
                    posInHash.push_back(addr);
                } else {
                    table[addr].weight += 1;
                }
            }


            // Sort words by frequency
            __gnu_parallel::sort(posInHash.begin(), posInHash.end(), [&](const size_t a, const size_t b) {
                return table[a].weight > table[b].weight;
            });

            std::shared_ptr<k2_tree_vocabulary> tmp_voc(new k2_tree_vocabulary(posInHash.size(), word_size));
            tmp_voc.swap(voc);

            #pragma omp parallel for
            for (uint i = 0; i < posInHash.size(); ++i) {
                Nword &w = table[posInHash[i]];
                w.codeword = i;
                voc->assign(i, w.word);
            }

        } catch (std::bad_alloc ba) {
            std::cerr << "[comperssion::FreqVoc] Error:" << ba.what() << "\n";
            exit(1);
        } catch (...) {
            std::cerr << "[comperssion::FreqVoc] Error: unexpected exception\n";
            exit(1);
        }
    }

    template <typename t_map>
    void construct_codewords(const int_vector<>& leaf_vector,
                        const t_map& codeword_map,
                        int_vector<> &codewords) {
        codewords.resize(leaf_vector.size());

        for (size_t i = 0; i < leaf_vector.size(); ++i) {
            auto it = codeword_map.find(leaf_vector[i]);
            if (it != codeword_map.end()) {
                codewords[i] = it->second;
            } else {
                std::cerr << "[k2_tree_compressor::construct_codewords] Error: Word not found\n";
                exit(1);
            }
        }
        //util::bit_compress(codewords);
    }

    void construct_codewords(const std::vector<uchar>& leaf_words, const size_t word_size, const size_t words_count, const HashTable& table, std::vector<uint>& codewords) {
        codewords.resize(words_count);
        size_t addr;
        for (size_t i = 0; i < words_count; ++i) {
            if (!table.search(&leaf_words[i*word_size], word_size, &addr)) {
                std::cerr << "[k2_tree_base::compress_leaves] Error: Word not found\n";
                exit(1);
            } else {
                //std::cout << "Codeword: " << table[addr].codeword << std::endl;
                codewords[i] = table[addr].codeword;
            }
        }
    }
}  // namespace sdsl
#endif  // INCLUDE_COMPRESSION_COMPRESSOR_H_
