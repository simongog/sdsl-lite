/* sdsl - succinct data structures library
    Copyright (C) 2012 Simon Gog

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
/*! \file csa_alphabet_strategy.hpp
    \brief csa_alphabet_strategy.hpp includes different strategy classes for representing an alphabet of a CSA.
    \author Simon Gog
*/

#ifndef INCLUDED_CSA_ALPHABET_STRATEGY
#define INCLUDED_CSA_ALPHABET_STRATEGY


// TODO: Strategy with 1-to-1 mapping and C_array type as template parameter
//       This can be also used for a integer based CSA.

/* A alphabet strategy provides the following features:
 *   * Member `sigma` which contains the size (=umber of unique symbols) of the alphabet.
 *   * Method `is_present(char_type c)` which indicates if character c occurs in the text.
 *   * Container `char2comp` which maps a symbol to a number [0..sigma-1]. The alphabetic
 *     order is preserved.
 *   * Container `comp2char` which is the inverse mapping of char2comp.
 *   * Container `C` contains the cumulative counts of occurrences. C[i] is the cumulative
 *     count of occurrences of symbols `comp2char[0]` to `comp2char[i-1]` in the text.
 *   * Typedefs for the four above members:
 *       * char2comp_type
 *       * comp2char_type
 *       * C_type
 *       * sigma_type
 *   * Constructor. Takes a int_vector_buffer<8> for byte-alphabets
 *     and int_vector_buffer<0> for integer-alphabets.
 *
 *    \par Note
 *   sigma_type has to be large enough to represent the alphabet size 2*sigma,
 *   since there is code which will perform a binary search on array `C`.
 */

#include "int_vector.hpp"
#include "sd_vector.hpp"
#include "rank_support.hpp"
#include "select_support.hpp"
#include "sdsl_concepts.hpp"
#include "config.hpp"
#include <string>

namespace sdsl
{

// forward declarations

class byte_alphabet;

template<class bit_vector_type     = bit_vector,
         class rank_support_type   = rank_support_scan<>,   //typename bit_vector_type::rank_1_type,
         class select_support_type = select_support_scan<>, //typename bit_vector_type::select_1_type,
         class C_array_type        = int_vector<>
         >
class succinct_byte_alphabet;

template<class bit_vector_type     = sd_vector<>,
         class rank_support_type   = typename bit_vector_type::rank_1_type,
         class select_support_type = typename bit_vector_type::select_1_type,
         class C_array_type        = int_vector<>
         >
class int_alphabet;

template <uint8_t int_width>
struct key_trait {
    static const char* KEY_BWT;
    static const char* KEY_TEXT;
};

template<>
struct key_trait<8> {
    static const char* KEY_BWT;
    static const char* KEY_TEXT;
};

template<uint8_t int_width>
const char* key_trait<int_width>::KEY_BWT = conf::KEY_BWT_INT;

template<uint8_t int_width>
const char* key_trait<int_width>::KEY_TEXT = conf::KEY_TEXT_INT;

template<class t_alphabet_strategy>
struct alphabet_trait {
    typedef byte_alphabet type;
};

template<>
struct alphabet_trait<int_alphabet_tag> {
    typedef int_alphabet<> type;
};

//! A simple space greedy representation for byte alphabets.
/*!
 *  \par Space consumption:
 *       At least: 2.5 kB
 *       Details:  char2comp + comp2char take  2*256 + 2*8 bytes
 *                 m_C                   takes       257*8 bytes
 *                 m_sigma               takes           2 bytes
 */
class byte_alphabet
{
    public:
        typedef int_vector<>::size_type size_type;
        typedef int_vector<8>           char2comp_type;
        typedef int_vector<8>           comp2char_type;
        typedef int_vector<64>          C_type;
        typedef uint16_t                sigma_type;
        typedef uint8_t                 char_type;
        typedef uint8_t                 comp_char_type;
        typedef std::string             string_type;
        enum { int_width = 8 };

        typedef byte_alphabet_tag       alphabet_category;
    private:
        char2comp_type m_char2comp; // Mapping from a character into the compact alphabet.
        comp2char_type m_comp2char; // Inverse mapping of m_char2comp.
        C_type         m_C;         // Cumulative counts for the compact alphabet [0..sigma].
        sigma_type     m_sigma;     // Effective size of the alphabet.

        void copy(const byte_alphabet&);
    public:

        const char2comp_type& char2comp;
        const comp2char_type& comp2char;
        const C_type&         C;
        const sigma_type&     sigma;

        //! Default constructor
        byte_alphabet();

        //! Construct from a byte-stream
        /*!
         *  \param text_buf Byte stream.
         *  \param len      Length of the byte stream.
         */
        byte_alphabet(int_vector_buffer<8>& text_buf, int_vector_size_type len);

        byte_alphabet(const byte_alphabet&);
        byte_alphabet(byte_alphabet&& b) : byte_alphabet() {
            *this = std::move(b);
        }

        byte_alphabet& operator=(const byte_alphabet&);
        byte_alphabet& operator=(byte_alphabet&&);

        void swap(byte_alphabet&);

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const;

        void load(std::istream& in);
};


//! A space-efficient representation for byte alphabets.
/*!
 *  The mapping `char2comp` and its inverse `comp2char` is realized internally
 *  by a bitvector of size 256 bits and a rank and a select structure. The rank
 *  structure is used to calculate `char2comp`; the select structure is used to
 *  calculate `comp2char`. Array `C` is represented by a bit-compressed
 *  `int_vector` and `sigma` by a uint16_t.
 *  The types to represent `char2comp`, `comp2char`, and `C` can be specified
 *  by template parameters.
 */
template<class bit_vector_type, class rank_support_type, class select_support_type, class C_array_type>
class succinct_byte_alphabet
{
    public:
        class char2comp_wrapper;
        class comp2char_wrapper;
        friend class char2comp_wrapper;
        friend class comp2char_wrapper;

        typedef int_vector<>::size_type size_type;
        typedef char2comp_wrapper       char2comp_type;
        typedef comp2char_wrapper       comp2char_type;
        typedef C_array_type            C_type;
        typedef uint16_t                sigma_type;
        typedef uint8_t                 char_type;
        typedef uint8_t                 comp_char_type;
        typedef std::string             string_type;
        typedef byte_alphabet_tag       alphabet_category;
        enum { int_width = 8 };

        //! Helper class for the char2comp mapping
        class char2comp_wrapper
        {
            private:
                const succinct_byte_alphabet* m_strat;
            public:
                char2comp_wrapper(const succinct_byte_alphabet* strat) : m_strat(strat) {}
                comp_char_type operator[](char_type c) const {
                    return (comp_char_type) m_strat->m_char_rank((size_type)c);
                }
        };

        //! Helper class for the comp2char mapping
        class comp2char_wrapper
        {
            private:
                const succinct_byte_alphabet* m_strat;
            public:
                comp2char_wrapper(const succinct_byte_alphabet* strat) : m_strat(strat) {}
                char_type operator[](comp_char_type c) const {
                    return (char_type) m_strat->m_char_select(((size_type)c)+1);
                }
        };

    private:
        bit_vector_type     m_char;        // `m_char[i]` indicates if character with code i is present or not
        rank_support_type   m_char_rank;   // rank data structure for `m_char` to answer char2comp
        select_support_type m_char_select; // select data structure for `m_char` to answer comp2char
        C_type              m_C;            // cumulative counts for the compact alphabet [0..sigma]
        sigma_type          m_sigma;       // effective size of the alphabet

        void copy(const succinct_byte_alphabet& strat) {
            m_char        = strat.m_char;
            m_char_rank   = strat.m_char_rank;
            m_char_rank.set_vector(&m_char);
            m_char_select = strat.m_char_select;
            m_char_select.set_vector(&m_char);
            m_C           = strat.m_C;
            m_sigma       = strat.m_sigma;
        }
    public:

        const char2comp_type char2comp;
        const comp2char_type comp2char;
        const C_type&        C;
        const sigma_type&    sigma;

        //! Default constructor
        succinct_byte_alphabet() : char2comp(this), comp2char(this), C(m_C), sigma(m_sigma) {
            m_sigma = 0;
        }

        //! Construct from a byte-stream
        /*!
         *  \param text_buf Byte stream.
         *  \param len      Length of the byte stream.
         */
        succinct_byte_alphabet(int_vector_buffer<8>& text_buf, int_vector_size_type len):
            char2comp(this), comp2char(this), C(m_C), sigma(m_sigma) {
            m_sigma = 0;
            if (0 == len or 0 == text_buf.size())
                return;
            assert(len <= text_buf.size());
            // initialize vectors
            int_vector<64> D(257, 0);
            bit_vector tmp_char(256, 0);
            // count occurrences of each symbol
            for (size_type i=0; i < len; ++i) {
                ++D[text_buf[i]];
            }
            assert(1 == D[0]); // null-byte should occur exactly once
            m_sigma = 0;
            for (int i=0; i<256; ++i)
                if (D[i]) {
                    tmp_char[i] = 1;    // mark occurring character
                    D[m_sigma] = D[i];  // compactify m_C
                    ++m_sigma;
                }
            // resize to sigma+1, since CSAs also need the sum of all elements
            m_C = C_type(m_sigma+1, 0, bits::hi(len)+1);

            for (int i=(int)m_sigma; i > 0; --i) m_C[i] = D[i-1];
            m_C[0] = 0;
            for (int i=1; i <= (int)m_sigma; ++i) m_C[i] = m_C[i] + m_C[i-1];
            assert(m_C[sigma]==len);
            m_char = tmp_char;
            util::init_support(m_char_rank, &m_char);
            util::init_support(m_char_select, &m_char);
        }

        //! Copy constructor
        succinct_byte_alphabet(const succinct_byte_alphabet& strat): char2comp(this), comp2char(this), C(m_C), sigma(m_sigma) {
            copy(strat);
        }

        //! Move constructor
        succinct_byte_alphabet(succinct_byte_alphabet&& strat) {
            *this = std::move(strat);
        }

        succinct_byte_alphabet& operator=(const succinct_byte_alphabet& strat) {
            if (this != &strat) {
                copy(strat);
            }
            return *this;
        }

        succinct_byte_alphabet& operator=(succinct_byte_alphabet&& strat) {
            if (this != &strat) {
                m_char        = std::move(strat.m_char);
                m_char_rank   = std::move(strat.m_char_rank);
                m_char_rank.set_vector(&m_char);
                m_char_select = std::move(strat.m_char_select);
                m_char_select.set_vector(&m_char);
                m_C           = std::move(strat.m_C);
                m_sigma       = std::move(strat.m_sigma);
            }
            return *this;
        }

        //! Swap operator
        void swap(succinct_byte_alphabet& strat) {
            m_char.swap(strat.m_char);
            util::swap_support(m_char_rank, strat.m_char_rank, &m_char, &(strat.m_char));
            util::swap_support(m_char_select, strat.m_char_select, &m_char, &(strat.m_char));
            m_C.swap(strat.m_C);
            std::swap(m_sigma,strat.m_sigma);
        }

        //! Serialize method
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_char.serialize(out, child, "m_char");
            written_bytes += m_char_rank.serialize(out, child, "m_char_rank");
            written_bytes += m_char_select.serialize(out, child, "m_char_select");
            written_bytes += m_C.serialize(out, child, "m_C");
            written_bytes += write_member(m_sigma, out, child, "m_sigma");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load method
        void load(std::istream& in) {
            m_char.load(in);
            m_char_rank.load(in);
            m_char_rank.set_vector(&m_char);
            m_char_select.load(in);
            m_char_select.set_vector(&m_char);
            m_C.load(in);
            read_member(m_sigma, in);
        }
};

//! A space-efficient representation for byte alphabets.
/*!
 *  The mapping `char2comp` and its inverse `comp2char` is realized internally
 *  by a bitvector of size sigma bits and a rank and a select structure, if the
 *  alphabet contains not all symbols in the range [0..sigma-1]. If it contains
 *  all symbols, i.e. the alphabet is continuous, then we map the symbols
 *  directly and no extra space is used.
 *
 *  The types to represent `char2comp`, `comp2char`, and `C` can be specified
 *  by template parameters.
 */
template<class bit_vector_type, class rank_support_type, class select_support_type, class C_array_type>
class int_alphabet
{
    public:
        class char2comp_wrapper;
        class comp2char_wrapper;
        friend class char2comp_wrapper;
        friend class comp2char_wrapper;

        typedef int_vector<>::size_type size_type;
        typedef char2comp_wrapper       char2comp_type;
        typedef comp2char_wrapper       comp2char_type;
        typedef C_array_type            C_type;
        typedef uint64_t                sigma_type;
        typedef uint64_t                char_type;
        typedef uint64_t                comp_char_type;
        typedef std::vector<char_type>  string_type;
        typedef int_alphabet_tag        alphabet_category;
        enum { int_width = 0 };

        //! Helper class for the char2comp mapping
        class char2comp_wrapper
        {
            private:
                const int_alphabet* m_strat;
            public:
                char2comp_wrapper(const int_alphabet* strat) : m_strat(strat) {}
                comp_char_type operator[](char_type c) const {
                    if (m_strat->m_char.size() > 0) {  // if alphabet is not continuous
                        return (comp_char_type) m_strat->m_char_rank((size_type)c);
                    } else { // direct map if it is continuous
                        return (comp_char_type) c;
                    }
                }
        };

        //! Helper class for the comp2char mapping
        class comp2char_wrapper
        {
            private:
                const int_alphabet* m_strat;
            public:
                comp2char_wrapper(const int_alphabet* strat) : m_strat(strat) {}
                char_type operator[](comp_char_type c) const {
                    if (m_strat->m_char.size() > 0) {  // if alphabet is not continuous
                        return (char_type) m_strat->m_char_select(((size_type)c)+1);
                    } else { // direct map if it is continuous
                        return (char_type) c;
                    }
                }
        };

    private:
        bit_vector_type     m_char;        // `m_char[i]` indicates if character with code i is present or not
        rank_support_type   m_char_rank;   // rank data structure for `m_char` to answer char2comp
        select_support_type m_char_select; // select data structure for `m_char` to answer comp2char
        C_type              m_C;           // cumulative counts for the compact alphabet [0..sigma]
        sigma_type          m_sigma;       // effective size of the alphabet

        void copy(const int_alphabet& strat) {
            m_char        = strat.m_char;
            m_char_rank   = strat.m_char_rank;
            m_char_rank.set_vector(&m_char);
            m_char_select = strat.m_char_select;
            m_char_select.set_vector(&m_char);
            m_C           = strat.m_C;
            m_sigma       = strat.m_sigma;
        }

        //! Check if the alphabet is continuous.
        bool is_continuous_alphabet(std::map<size_type, size_type>& D) {
            if (D.size() == 0) {  // an empty alphabet is continuous
                return true;
            } else {
                //            max key      + 1  ==  size of map
                return ((--D.end())->first + 1) ==  D.size();
            }
        }
    public:

        const char2comp_type char2comp;
        const comp2char_type comp2char;
        const C_type&        C;
        const sigma_type&    sigma;

        //! Default constructor
        int_alphabet() : char2comp(this), comp2char(this), C(m_C), sigma(m_sigma) {
            m_sigma = 0;
        }

        //! Construct from a byte-stream
        /*!
         *  \param text_buf Byte stream.
         *  \param len      Length of the byte stream.
         */
        int_alphabet(int_vector_buffer<0>& text_buf, int_vector_size_type len):
            char2comp(this), comp2char(this), C(m_C), sigma(m_sigma) {
            m_sigma = 0;
            if (0 == len or 0 == text_buf.size())
                return;
            assert(len <= text_buf.size());
            // initialize vectors
            std::map<size_type, size_type> D;
            // count occurrences of each symbol
            for (size_type i=0; i < len; ++i) {
                D[text_buf[i]]++;
            }
            m_sigma = D.size();
            if (is_continuous_alphabet(D)) {
                // do not initialize m_char, m_char_rank and m_char_select since we can map directly
            } else {
                // note: the alphabet has at least size 1, so the following is safe:
                size_type largest_symbol = (--D.end())->first;
                bit_vector tmp_char(largest_symbol+1, 0);
                for (std::map<size_type, size_type>::const_iterator it = D.begin(), end=D.end(); it != end; ++it) {
                    tmp_char[it->first] = 1;
                }
                m_char = tmp_char;
                util::init_support(m_char_rank, &m_char);
                util::init_support(m_char_select, &m_char);
            }
            assert(D.find(0) != D.end() and 1 == D[0]); // null-byte should occur exactly once

            // resize to sigma+1, since CSAs also need the sum of all elements
            m_C = C_type(m_sigma+1, 0, bits::hi(len)+1);
            size_type sum = 0, idx=0;
            for (std::map<size_type, size_type>::const_iterator it = D.begin(), end=D.end(); it != end; ++it) {
                m_C[idx++] = sum;
                sum += it->second;
            }
            m_C[idx] = sum;  // insert sum of all elements
        }

        //! Copy constructor
        int_alphabet(const int_alphabet& strat): char2comp(this), comp2char(this), C(m_C), sigma(m_sigma) {
            copy(strat);
        }

        //! Copy constructor
        int_alphabet(int_alphabet&& strat) {
            *this = std::move(strat);
        }

        int_alphabet& operator=(const int_alphabet& strat) {
            if (this != &strat) {
                copy(strat);
            }
            return *this;
        }

        int_alphabet& operator=(int_alphabet&& strat) {
            if (this != &strat) {
                m_char        = std::move(strat.m_char);
                m_char_rank   = std::move(strat.m_char_rank);
                m_char_rank.set_vector(&m_char);
                m_char_select = std::move(strat.m_char_select);
                m_char_select.set_vector(&m_char);
                m_C           = std::move(strat.m_C);
                m_sigma       = std::move(strat.m_sigma);
            }
            return *this;
        }

        //! Swap operator
        void swap(int_alphabet& strat) {
            m_char.swap(strat.m_char);
            util::swap_support(m_char_rank, strat.m_char_rank, &m_char, &(strat.m_char));
            util::swap_support(m_char_select, strat.m_char_select, &m_char, &(strat.m_char));
            m_C.swap(strat.m_C);
            std::swap(m_sigma,strat.m_sigma);
        }

        //! Serialize method
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_char.serialize(out, child, "m_char");
            written_bytes += m_char_rank.serialize(out, child, "m_char_rank");
            written_bytes += m_char_select.serialize(out, child, "m_char_select");
            written_bytes += m_C.serialize(out, child, "m_C");
            written_bytes += write_member(m_sigma, out, child, "m_sigma");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load method
        void load(std::istream& in) {
            m_char.load(in);
            m_char_rank.load(in);
            m_char_rank.set_vector(&m_char);
            m_char_select.load(in);
            m_char_select.set_vector(&m_char);
            m_C.load(in);
            read_member(m_sigma, in);
        }
};

} // end namespace sdsl

#endif
