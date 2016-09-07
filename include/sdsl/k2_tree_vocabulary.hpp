/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <nlehmann@dcc.uchile.cl> wrote this file. As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return Nicolás Lehmann
 * ----------------------------------------------------------------------------
 */

#ifndef INCLUDE_COMPRESSION_VOCABULARY_H_
#define INCLUDE_COMPRESSION_VOCABULARY_H_

#include <algorithm>
#include <fstream>
#include "k2_tree_helper.hpp"

namespace sdsl {

    class Vocabulary {

    #ifndef uchar
    #define uchar unsigned char
    #endif

    public:

        Vocabulary() = default;

        /**
         * @param cnt Number of words in the vocabulary
         * @param size Size of each word in bytes
         */
        Vocabulary(size_t cnt, uint size)    : cnt_(cnt),
                                               size_(size),
                                               data_(new uchar[cnt*size_]) {}

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream &out, structure_tree_node *v = nullptr,
                            std::string name = "") const {
            structure_tree_node *child = structure_tree::add_child(
                    v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(cnt_, out, child, "cnt");
            written_bytes += write_member(size_, out, child, "size");

            written_bytes += write_member(data_, cnt_*size_, out, child, "data");

            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream &in) {
            read_member(cnt_, in);
            read_member(size_, in);

            data_ = (uchar*) malloc(sizeof(uchar) *cnt_*size_);
            read_member(&data_, cnt_*size_, in);
        }

        const uchar *operator[](size_t i) const {
            return data_ + i * size_;
        }

        const uchar *get(size_t i) const {
            return data_ + i * size_;
        }

        /**
         * Print the vocabulary
         */
        void print() const{
            for (size_t i = 0; i < cnt_; ++i)
                print(i);
            printf("\n");
        }

        /**
         * Print word on the given position
         * @param i
         */
        void print(size_t i) const{
            const uchar *word = (*this)[i];
            printf("[");
            for (uint bit = 0; bit < size_*sizeof(uchar)*8; ++bit) {
                if ((word[bit/sizeof(uchar)/8] >> (bit%(sizeof(uchar)*8))) & 1)
                    printf("1");
                else
                    printf("0");
            }
            printf("]");
        }

        /**
         * Stores the given word in the specified position.
         * @param p Position to store the word.
         * @param word Word to store.
         */
        void assign(size_t p, const uchar *word){
            std::copy(word, word + size_, data_ + p*size_);
        }

        /**
         * Sort the vocabulary lexicographically and return the number of differents
         * words.
         * @return Number of differents words.
         */
        size_t sort() {
            return quicksort(0, cnt_);
        }

        /**
         * Return the size of the words in the vocabulary.
         *
         * @return Size in bytes of the words.
         */
        uint size() const {
            return size_;
        }

        /**
         * Returns number of words in the vocabulary.
         *
         * @return Number of words in the vocabulary.
         */
        size_t cnt() const {
            return cnt_;
        }

        ~Vocabulary(){
        }

        void destroy(){
            delete [] data_;
        }

        bool operator==(const Vocabulary &rhs) const{
            if (cnt_ != rhs.cnt_ && size_ != rhs.size_)
                return false;
            for (size_t i = 0; i < cnt_*size_; ++i)
                if (data_[i] != rhs.data_[i])
                    return false;
            return true;
        }

        //! Swap operator
        /*void swap(Vocabulary &voc) {
            if (this != &voc) {
                std::swap(cnt_, voc.cnt_);
                std::swap(size_, voc.size_);
                std::swap(data_, voc.data_);
            }
        }*/

    private:
        /**
         * Swap words in the vocabulary.
         * @param a Position of the first word.
         * @param b Position of the second word.
         */
        void swap(size_t a, size_t b){
            for (uint i = 0; i < size_; ++i)
                std::swap(data_[a*size_ + i], data_[b*size_ + i]);
        }

        /**
         * Sort lexicographically the words between positions left and right, and
         * returns the number of differents words.
         * @param left Left position of the subarray.
         * @param right Right position of the subarray.
         * @return Number of differents words.
         */
        size_t quicksort(size_t left, size_t right){
            if (left == right)
                return 0;
            if (left == right - 1)
                return 1;

            size_t i, j, k;
            i = j = k = left;
            while (k < right - 1) {
                // We use the last element as pivot.
                int cmp = strcmp((*this)[k], (*this)[right-1]);

                if (cmp == 0) {
                    swap(j, k);
                    ++j;
                }
                if (cmp < 0) {
                    swap(j, k);
                    swap(j, i);
                    ++j;
                    ++i;
                }
                ++k;
            }
            // Move the pivot to final position
            swap(j, right - 1);
            ++j;

            size_t unique;
            unique = quicksort(left, i);
            unique += quicksort(j, right);
            return unique + 1;
        }

        int strcmp(const uchar *w1, const uchar *w2) {
            uint i = 0;
            while (i < size_ - 1 && *w1 == *w2) {
                ++w1;
                ++w2;
                ++i;
            }
            return (int) *w1 - (int) *w2;
        }

        /** Number of words in the vocabulary*/
        size_t cnt_;
        /** Size in bytes of each word*/
        uint size_;
        /** Array storing words*/
        uchar *data_;
    };

}  // namespace sdsl
#endif  // INCLUDE_COMPRESSION_VOCABULARY_H_
