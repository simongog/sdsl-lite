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
#include "k2_tree_vocabulary.hpp"


namespace sdsl {
    using namespace k2_treap_ns;

/**
 * Creates new k2tree with the leaves compressed.
 * @param tree
 * @param build 
 */
template<class Fun>
void FreqVoc(const std::vector<uchar>& leaf_words, const uint word_size, const size_t word_count, Fun build, uint64_t hash_size = 0) {
  try {
      // Insert words in hash
      if (hash_size == 0) {

          std::cerr << "[k2_tree_compressor::FreqVoc] Warning: Hash Size not specified, it will thus be automatically determined by amount of distinct values, which is really slow!" << std::endl;

          Vocabulary words(word_count, word_size);

          for (size_t i = 0; i < word_count; ++i) {
              words.assign(i, leaf_words, i*word_size);
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
          if (!table.search(&leaf_words[i*word_size], word_size, &addr)) {
              table.add(&leaf_words[i*word_size], word_size, addr);
              posInHash.push_back(addr);
          } else {
              table[addr].weight += 1;
          }
      }


    // Sort words by frequency
    std::sort(posInHash.begin(), posInHash.end(), [&](size_t a, size_t b) {
      return table[a].weight > table[b].weight;
    });

    Vocabulary voc(posInHash.size(), word_size);
    for (uint i = 0; i < posInHash.size(); ++i) {
      Nword &w = table[posInHash[i]];
      w.codeword = i;
      voc.assign(i, w.word);
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

}  // namespace sdsl
#endif  // INCLUDE_COMPRESSION_COMPRESSOR_H_
