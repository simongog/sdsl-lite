/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <nlehmann@dcc.uchile.cl> wrote this file. As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return Nicolás Lehmann
 * ----------------------------------------------------------------------------
 */

#ifndef INCLUDE_COMPRESSION_COMPRESSOR_H_
#define INCLUDE_COMPRESSION_COMPRESSOR_H_

#include <sdsl/k2_tree_hash_table.hpp>
#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>
#include <climits>
#include "k2_tree_helper.hpp"
#include <parallel/algorithm>
#include "k2_tree_vocabulary.hpp"
#include "construct.hpp"
#include "wavelet_trees.hpp"


namespace sdsl {
    using namespace k2_treap_ns;

    template<class Fun>
    void
    frequency_encode(const int_vector<> &leaf_words, int_vector<> &dictionary,
                     std::unordered_map<int_vector<>::reference, uint> &hashmap);


    void
    construct_codewords(const int_vector<> &leaf_vector,
                        std::unordered_map<int_vector<>::reference, uint> codeword_map,
                        int_vector<32> &codewords);



    void dac_compress(const int_vector<> &leaf_words, int_vector<> &dictionary, dac_vector &compressed_leaves) {
        std::unordered_map<int_vector<>::reference, uint> hashmap;
        frequency_encode(leaf_words, dictionary, hashmap);
        int_vector<32> codewords; //size is known: bits:hi for hashmap.size()? or distinct_values.size()
        construct_codewords(leaf_words, hashmap, codewords);
        dac_vector tmp;
        construct_im(tmp, codewords);
        tmp.swap(compressed_leaves);
    }

    template<class Fun>
    void legacy_dac_encode(const std::vector<uchar> &leaf_words, const uint word_size, const size_t word_count,
                           Fun build, uint64_t hash_size = 0) {
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
            HashTable table(hash_size);
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
            std::sort(posInHash.begin(), posInHash.end(), [&](const size_t a, const size_t b) {
                return table[a].weight > table[b].weight;
            });

            std::shared_ptr<k2_tree_vocabulary> voc(new k2_tree_vocabulary(posInHash.size(), word_size));

            #pragma omp parallel for
            for (uint i = 0; i < posInHash.size(); ++i) {
                Nword &w = table[posInHash[i]];
                w.codeword = i;
                voc->assign(i, w.word);
            }


            build(table, voc, leaf_words);
        } catch (std::bad_alloc ba) {
            std::cerr << "[comperssion::FreqVoc] Error:" << ba.what() << "\n";
            exit(1);
        } catch (...) {
            std::cerr << "[comperssion::FreqVoc] Error: unexpected exception\n";
            exit(1);
        }
    }

    void wt_huff_int_dict_compress(const int_vector<> &leaf_words, int_vector<> &dictionary,
                                   wt_huff_int<> &compressed_leaves) {
        std::unordered_map<int_vector<>::reference, uint> hashmap;
        frequency_encode(leaf_words, dictionary, hashmap);
        int_vector<32> codewords; //size is known: bits:hi for hashmap.size()? or distinct_values.size()
        construct_codewords(leaf_words, hashmap, codewords);
        wt_huff_int<> tmp;
        construct_im(tmp, codewords);
        tmp.swap(compressed_leaves);
    }

    void
    wt_huff_int_compress(const int_vector<> &leaf_words, int_vector<> &dictionary,
                         wt_huff_int<> &compressed_leaves) {
        std::unordered_map<int_vector<>::reference, uint> hashmap;
        wt_huff_int<> tmp;
        construct_im(tmp, leaf_words);
        tmp.swap(compressed_leaves);
    }

    /**
    * Creates new k2tree with the leaves compressed.
    * @param tree
    * @param build
    */
    template<class Fun>
    void
    frequency_encode(const int_vector<> &leaf_words, int_vector<> &dictionary,
                     std::unordered_map<int_vector<>::reference, uint> &hashmap) {

        std::cout << "Before putting words in hash" << std::endl;

        std::string tmp_file = ram_file_name(util::to_string(util::pid()) + "_" + util::to_string(util::id()));
        int_vector_buffer<> dictionary_buffer;
        for (size_t i = 0; i < leaf_words.size(); ++i) {
            auto it = hashmap.find(leaf_words[i]);
            if (it != hashmap.end()) {
                it->second += 1;
                dictionary_buffer.push_back(leaf_words[i]);
            } else {
                hashmap.insert(std::make_pair(leaf_words[i], (uint) 1));//consider using sparsehash
            }
        }

        int_vector<> tmp(dictionary_buffer, dictionary_buffer.size());
        tmp.swap(dictionary);

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
        __gnu_parallel::sort(dictionary.begin(), dictionary.end(),
                             [&](const int_vector<>::reference a, const int_vector<>::reference b) {
                                 return hashmap.find(a)->second > hashmap.find(b)->second;
                             });


        if (dictionary.size() > UINT_MAX) {
            std::cerr << "k2_tree_vocabulary size is bigger than 32 Bit, things might break! (have fun anyway"
                      << std::endl;
        }

        #pragma omp parallel for
        for (size_t i = 0; i < dictionary.size(); ++i) {
            auto it = hashmap.find(dictionary[i]);
            it->second = i;
        }
    }


    void
    construct_codewords(const int_vector<> &leaf_vector,
                        std::unordered_map<int_vector<>::reference, uint> codeword_map,
                        int_vector<32> &codewords) {
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
    }

    void construct_codewords(const std::vector<uchar>& leaf_words, const size_t words_count, const size_t word_size, const HashTable& table, std::vector<uint>& codewords) const {
        codewords.reserve(word_count);
        size_t addr;
        for (size_t i = 0; i < word_count; ++i) {
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
