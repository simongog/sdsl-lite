/* Searches using  Horspool's algorithm adapted to 
search inside a text compressed with End-Tagged Dense Code.
Lightweight Natural Language Text Compression: Information Retrieval 2006

Programmed by Antonio Fari�a.

Author's contact: Antonio Fari�a, Databases Lab, University of
A Coru�a. Campus de Elvi�a s/n. Spain  fari@udc.es

Copyright (C) 2006  Nieves R. Brisaboa, Gonzalo Navarro, Antonio Fari�a and Jos� R. Param�
Author's contact: antonio.fari@gmail.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
aint with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/


/*-----------------------------------------------------------------------
 Hash.h: Definition of a HashTable that uses linear hash
 ------------------------------------------------------------------------*/

#ifndef INCLUDE_COMPRESSION_HASH_H_
#define INCLUDE_COMPRESSION_HASH_H_

#include <cstdio>
#include <vector>
#include "k2_tree_helper.hpp"

// jump done when a collision appears
#define JUMP 101
// a small prime number, used to compute a hash function
#define SMALL_PRIME 1009
#define SEED 1159241
/* Type definitions */


namespace sdsl {

struct Nword {
  Nword() : word(NULL), len(0), weight(0), codeword(0) {}
  const uchar *word;
  uint len;
  uint weight;
  uint codeword;
};

class HashTable {
 public:
    HashTable(size_t sizeVoc, double occup_hash = 1.5) : tam_hash_(NearestPrime((size_t) (occup_hash * sizeVoc))),
                                                         num_elem_(0),
                                                         hash_() {
        if (tam_hash_ <= JUMP)
            tam_hash_ = NearestPrime(JUMP + 1);
        hash_.resize(tam_hash_);
    }
  /**
   * Adds a word to the hash in the specific addres.
   *
   * @param aWord Word
   * @param len Length of the word
   * @param addr Position in the table to store the word.
   * @return Position in the table where the word has stored.
   */
  size_t add(const uchar *aWord, uint len, size_t addr){
      if (addr == tam_hash_) {
          printf("Not enough memory, vocabulary exceeds maximun size !\n");
          exit(1);
      }

      hash_[addr].word = aWord;
      hash_[addr].len = len;
      hash_[addr].weight = 1;
      num_elem_++;

      return addr;
  }
  /**
   * Search a word in the table. After return returnedAddr hold an address in
   * the table. In case the word is found, the address correspond to the word
   * position. In case the word is not present, the address holds the position
   * where to store the word
   *
   * @param aWord Word
   * @param len Length of the word
   * @param returnedAddr Pointer to store the address.
   * @return Ture in case the word is present and false otherwise.
   */
  bool search(const uchar *aWord, uint len, size_t *returnedAddr) const{
      size_t addr;
      addr = hashFunction(aWord, len);

      while ((hash_[addr].word != NULL) &&
             (strcomp(hash_[addr].word, aWord, hash_[addr].len, len) != 0))
          addr = (addr + JUMP) % tam_hash_;
      // position returned
      *returnedAddr = addr;

      // Word was not found
      if (hash_[addr].word  == NULL)
          return false;
      // Word was found
      return true;
  }

  Nword &operator[](size_t i) {
    return hash_[i];
  }

  const Nword &operator[](size_t i) const {
    return hash_[i];
  }

  size_t num_elem() {
    return num_elem_;
  }

 private:
  /**  entries in the hash table.*/
  size_t tam_hash_;
  /** elements already added to the hash table.*/
  size_t num_elem_;
  /** holds a hashTable of words*/
  std::vector<Nword> hash_;

  /*------------------------------------------------------------------
   Modification of Zobel's bitwise function to have into account the 
   lenght of the key explicetely 
   ---------------------------------------------------------------- */
  size_t hashFunction(const uchar *aWord, uint len) const{
      uchar c;
      size_t h;

      h = SEED;

      uint i = 0;
      while (i++ < len) {
          c = *(aWord++);
          h ^= ( (h << 5) + c + (h >> 2) );
      }
      return (h&0x7fffffff) % tam_hash_;
  }

  /*------------------------------------------------------------------
   Modification of Zobel's scmp function compare two strings
   ---------------------------------------------------------------- */
  inline int strcomp(const uchar *s1, const uchar *s2,
                     uint ws1, uint ws2) const{
      if (ws1 !=ws2)
          return -1;

      uint iters;
      iters = 0;
      while (iters < ws1-1 && *s1 == *s2) {
          s1++;
          s2++;
          iters++;
      }
      return( *s1-*s2 );
  }

    /**
    * Find the smallest prime greater or equal to n
    */
    size_t NearestPrime(size_t n) {
      /* the prime number being sought */
      size_t pos;
      size_t i;

      for (pos = n; ; pos++) {
        // checks if those values from 2 to $\sqrt{pos}$ can be factors of $pos$
        for (i = 2; i <= pos/i + 1 && pos % i != 0; i++) ;
        // No factors in that range, therefore a prime number was found
        if (pos % i != 0)
          break;
      }
      return pos;
    }
};

}  // namespace sdsl

#endif  // INCLUDE_COMPRESSION_HASH_H_
