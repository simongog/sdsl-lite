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
 *   * Member `sigma` which contains the size (=number of unique symbols) of the alphabet.
 *   * Container `char2comp` which maps a symbol to a number [0..sigma-1]. The alphabetic
 *     order is preserved.
 *   * Container `comp2char` which is the inverse mapping of char2comp.
 *   * Container `C` contains the cumulative counts of occurrences. C[i] is the cumulative
 *     count of occurrences of symbols `comp2char[0]` to `comp2char[i-1]` in the text.
 *     C is of size `sigma+1`.
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

template<uint8_t t_q                 = 3,
         typename bit_vector_type     = bit_vector,
         typename rank_support_type   = rank_support_scan<>,
         typename select_support_type = select_support_scan<>
         >
class succinct_multibyte_alphabet;

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

// see http://stackoverflow.com/questions/13514587/c-check-for-nested-typedef-of-a-template-parameter-to-get-its-scalar-base-type
// for the next three functions


template<class t_wt, class t_enable = void>
struct wt_alphabet_trait {
    typedef t_enable type;
};

template<class t_wt>
struct wt_alphabet_trait<t_wt, typename enable_if_type<typename t_wt::alphabet_category>::type> {
    using type = typename alphabet_trait<typename t_wt::alphabet_category>::type;
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
        byte_alphabet(byte_alphabet&& b) : byte_alphabet()
        {
            *this = std::move(b);
        }

        byte_alphabet& operator=(const byte_alphabet&);
        byte_alphabet& operator=(byte_alphabet&&);

        void swap(byte_alphabet&);

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const;

        void load(std::istream& in);
};


//! Helper class for the char2comp mapping
template<typename t_alphabet_strat>
class char2comp_wrapper
{
    private:
        const t_alphabet_strat* m_strat;
    public:
        using comp_char_type = typename t_alphabet_strat::comp_char_type;
        using char_type = typename t_alphabet_strat::char_type;
        using size_type = typename t_alphabet_strat::size_type;

        char2comp_wrapper(const t_alphabet_strat* strat) : m_strat(strat) {}

        comp_char_type operator[](char_type c) const // TODO: using a const reference???
        {
            if (c >= m_strat->m_char.size() or !m_strat->m_char[c])
                return (comp_char_type)0;
            return (comp_char_type) m_strat->m_char_rank((size_type)c);
        }

        template<typename t_strat = t_alphabet_strat>
        typename std::enable_if<(t_strat::q>1), typename t_strat::multi_comp_char_type>::type
        operator[](const typename std::enable_if<(t_strat::q>1), typename t_strat::multi_char_type>::type& c) const
        {
            typename t_strat::multi_comp_char_type x {0};
            auto sigma_size =  m_strat->sigma;
            for (size_t i=0; i < t_alphabet_strat::q; ++i) {
                if (c >= m_strat->m_char.size() or !m_strat->m_char[c])
                    return 0ULL;
                x *= m_strat->sigma_q;
                x += m_strat->m_char_rank((size_type)c);
            }
            return x;
        }

};

//! Helper class for the comp2char mapping
template<typename t_alphabet_strat>
class comp2char_wrapper
{
    private:
        const t_alphabet_strat* m_strat;
    public:
        using char_type = typename t_alphabet_strat::char_type;
        using comp_char_type = typename t_alphabet_strat::comp_char_type;
        using size_type = typename t_alphabet_strat::size_type;

        comp2char_wrapper(const t_alphabet_strat* strat) : m_strat(strat) {}

        char_type operator[](comp_char_type c) const // TODO: using a const reference???
        {
            return (char_type) m_strat->m_char_select(((size_type)c)+1);
        }

        template<typename t_strat = t_alphabet_strat>
        typename std::enable_if<(t_strat::q>1), typename t_strat::multi_char_type>::type
        operator[](typename std::enable_if<(t_strat::q>1), typename t_strat::multi_comp_char_type>::type c) const
        {
            std::cout<<"TODO comp2char multi_byte x="<<static_cast<uint64_t>(c)<<" t_alphabet_strat::q="<<(int)t_alphabet_strat::q<<std::endl;
            return 0; // TODO
        }
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
        static constexpr uint8_t q = 1;

        class char2comp_wrapper<succinct_byte_alphabet>;
        class comp2char_wrapper<succinct_byte_alphabet>;
        friend class char2comp_wrapper<succinct_byte_alphabet>;
        friend class comp2char_wrapper<succinct_byte_alphabet>;
        typedef char2comp_wrapper<succinct_byte_alphabet> char2comp_type;
        typedef comp2char_wrapper<succinct_byte_alphabet> comp2char_type;

        typedef int_vector<>::size_type size_type;
        typedef C_array_type            C_type;
        typedef uint16_t                sigma_type;
        typedef sigma_type              multi_sigma_type;
        typedef uint8_t                 char_type;
        typedef uint8_t                 comp_char_type;
        typedef std::string             string_type;
        typedef byte_alphabet_tag       alphabet_category;
        enum { int_width = 8 };

    private:
        bit_vector_type     m_char;        // `m_char[i]` indicates if character with code i is present or not
        rank_support_type   m_char_rank;   // rank data structure for `m_char` to answer char2comp
        select_support_type m_char_select; // select data structure for `m_char` to answer comp2char
        C_type              m_C;            // cumulative counts for the compact alphabet [0..sigma]
        sigma_type          m_sigma;       // effective size of the alphabet

        void copy(const succinct_byte_alphabet& strat)
        {
            m_char        = strat.m_char;
            m_char_rank   = strat.m_char_rank;
            m_char_rank.set_vector(&m_char);
            m_char_select = strat.m_char_select;
            m_char_select.set_vector(&m_char);
            m_C           = strat.m_C;
            m_sigma       = strat.m_sigma;
        }
    public:

        const char2comp_type    char2comp;
        const comp2char_type    comp2char;
        const C_type&           C;
        const sigma_type&       sigma;
        const multi_sigma_type& sigma_q;

        //! Default constructor
        succinct_byte_alphabet() : char2comp(this), comp2char(this), C(m_C), sigma(m_sigma), sigma_q(m_sigma)
        {
            m_sigma = 0;
        }

        //! Construct from a byte-stream
        /*!
         *  \param text_buf Byte stream.
         *  \param len      Length of the byte stream.
         */
        succinct_byte_alphabet(int_vector_buffer<8>& text_buf, int_vector_size_type len):
            char2comp(this), comp2char(this), C(m_C), sigma(m_sigma), sigma_q(m_sigma)
        {
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
        succinct_byte_alphabet(const succinct_byte_alphabet& strat): char2comp(this), comp2char(this), C(m_C), sigma(m_sigma), sigma_q(m_sigma)
        {
            copy(strat);
        }

        //! Move constructor
        succinct_byte_alphabet(succinct_byte_alphabet&& strat)
        {
            *this = std::move(strat);
        }

        succinct_byte_alphabet& operator=(const succinct_byte_alphabet& strat)
        {
            if (this != &strat) {
                copy(strat);
            }
            return *this;
        }

        succinct_byte_alphabet& operator=(succinct_byte_alphabet&& strat)
        {
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
        void swap(succinct_byte_alphabet& strat)
        {
            m_char.swap(strat.m_char);
            util::swap_support(m_char_rank, strat.m_char_rank, &m_char, &(strat.m_char));
            util::swap_support(m_char_select, strat.m_char_select, &m_char, &(strat.m_char));
            m_C.swap(strat.m_C);
            std::swap(m_sigma,strat.m_sigma);
        }

        //! Serialize method
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
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
        void load(std::istream& in)
        {
            m_char.load(in);
            m_char_rank.load(in);
            m_char_rank.set_vector(&m_char);
            m_char_select.load(in);
            m_char_select.set_vector(&m_char);
            m_C.load(in);
            read_member(m_sigma, in);
        }
};

class multibyte_tag;
class multibyte_comp_char;

class multibyte_comp_char
{
    private:
        uint64_t m_x; // value
    public:
        typedef multibyte_tag type;

        template<typename t_alphabet_strat>
        friend bool cyclic_insert_hi(typename t_alphabet_strat::multi_comp_char_type&,
                                     typename t_alphabet_strat::char_type,
                                     const t_alphabet_strat&);

        template<typename t_alphabet_strat>
        friend bool cyclic_insert_lo(typename t_alphabet_strat::multi_comp_char_type&,
                                     typename t_alphabet_strat::char_type,
                                     const t_alphabet_strat&);

        multibyte_comp_char(uint64_t x) : m_x(x) {}

        explicit operator uint64_t() const
        {
            return m_x;
        }

        multibyte_comp_char operator+(uint64_t add)const
        {
            return multibyte_comp_char(m_x+add);
        }
};

template<typename t_alphabet_strat>
bool cyclic_insert_hi(typename t_alphabet_strat::multi_comp_char_type& mc,
                      typename t_alphabet_strat::char_type c,
                      const t_alphabet_strat& alphabet)
{
    auto cc = alphabet.char2comp[c];
    if (cc == 0 and c > 0)
        return false;
//    std::cout<<"mc.mx="<<mc.m_x<<" cc="<<(int)cc<<std::endl;
    mc.m_x /= alphabet.sigma;
    mc.m_x += (cc * alphabet.sigma_q_1);
//    std::cout<<"mc.mx="<<mc.m_x<<" alphabet.sigma_q_1="<<alphabet.sigma_q_1<<std::endl;
    return true;
}

template<typename t_alphabet_strat>
bool cyclic_insert_low(typename t_alphabet_strat::multi_comp_char_type& mc,
                       typename t_alphabet_strat::char_type c,
                       const t_alphabet_strat& alphabet)
{
    auto cc = alphabet.char2comp[c];
    if (cc == 0 and c > 0)
        return false;
    mc.m_x *= alphabet.sigma;
    mc.m_x = (mc.m_x + cc) % alphabet.sigma_q;
}

//! A space-efficient representation for a multi-byte alphabet strategy.
/*!
 *  \tparam t_q
 *  The mapping `char2comp` and its inverse `comp2char` is realized internally
 *  by a bitvector of size 256 bits and a rank and a select structure. The rank
 *  structure is used to calculate `char2comp`; the select structure is used to
 *  calculate `comp2char`. Array `C` is represented by a bit-compressed
 *  `int_vector` and `sigma` by a uint16_t.
 *  The types to represent `char2comp`, `comp2char`, and `C` can be specified
 *  by template parameters.
 */
template<uint8_t t_q,
         typename bit_vector_type,
         typename rank_support_type,
         typename select_support_type
         >
class succinct_multibyte_alphabet
{
    public:
        class char2comp_wrapper<succinct_multibyte_alphabet>;
        class comp2char_wrapper<succinct_multibyte_alphabet>;
        class multibyte_C;
        friend class char2comp_wrapper<succinct_multibyte_alphabet>;
        friend class comp2char_wrapper<succinct_multibyte_alphabet>;
        typedef char2comp_wrapper<succinct_multibyte_alphabet> char2comp_type;
        typedef comp2char_wrapper<succinct_multibyte_alphabet> comp2char_type;
        typedef multibyte_C C_type;
        static constexpr uint8_t q = t_q;

        typedef int_vector<>::size_type size_type;
        typedef uint16_t                sigma_type;
        typedef uint64_t                multi_sigma_type;
        typedef uint8_t                 char_type;
        typedef uint8_t                 comp_char_type;
        typedef std::array<uint8_t,q>   multi_char_type;
        typedef multibyte_comp_char     multi_comp_char_type;
        typedef std::string             string_type;
        typedef byte_alphabet_tag       alphabet_category;
        enum { int_width = 8 };

        struct multibyte_C { // TODO add proper constructor
            int_vector<> C;
            int_vector<> multi_C;

            size_type operator[](comp_char_type c) const
            {
                return C[c];
            }

            typename std::enable_if<
            std::is_same<typename multi_comp_char_type::type,
                multibyte_tag>::value, size_type>::type
                operator[](multi_comp_char_type c) const
            {
                return multi_C[static_cast<uint64_t>(c)];
            }

            void swap(multibyte_C& t)
            {
                C.swap(t.C);
                multi_C.swap(t.multi_C);
            }

            //! Serialize method
            size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
            {
                structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
                size_type written_bytes = 0;
                written_bytes += C.serialize(out, child, "C");
                written_bytes += multi_C.serialize(out, child, "multi_C");
                structure_tree::add_size(child, written_bytes);
                return written_bytes;
            }

            //! Load method
            void load(std::istream& in)
            {
                C.load(in);
                multi_C.load(in);
            }
        };

    private:
        bit_vector_type     m_char;        // `m_char[i]` indicates if character with code i is present or not
        rank_support_type   m_char_rank;   // rank data structure for `m_char` to answer char2comp
        select_support_type m_char_select; // select data structure for `m_char` to answer comp2char
        C_type              m_C;           // cumulative counts for the compact alphabet [0..sigma]
        sigma_type          m_sigma;       // effective size of the alphabet
        multi_sigma_type    m_sigma_q;     // sigma^q
        multi_sigma_type    m_sigma_q_1;   // sigma^{q-1}

        void copy(const succinct_multibyte_alphabet& strat)
        {
            m_char        = strat.m_char;
            m_char_rank   = strat.m_char_rank;
            m_char_rank.set_vector(&m_char);
            m_char_select = strat.m_char_select;
            m_char_select.set_vector(&m_char);
            m_C           = strat.m_C;
            m_sigma       = strat.m_sigma;
        }
    public:

        const char2comp_type    char2comp;
        const comp2char_type    comp2char;
        const C_type&           C;
        const sigma_type&       sigma;
        const multi_sigma_type& sigma_q;
        const multi_sigma_type& sigma_q_1;

        //! Default constructor
        succinct_multibyte_alphabet() : char2comp(this), comp2char(this),
            C(m_C), sigma(m_sigma), sigma_q(m_sigma_q),
            sigma_q_1(m_sigma_q_1)
        {
            m_sigma = 0;
            m_sigma_q = 0;
        }

        //! Construct from a byte-stream
        /*!
         *  \param text_buf Byte stream.
         *  \param len      Length of the byte stream.
         */
        succinct_multibyte_alphabet(int_vector_buffer<8>& text_buf, int_vector_size_type len):
            char2comp(this), comp2char(this), C(m_C), sigma(m_sigma), sigma_q(m_sigma_q),
            sigma_q_1(m_sigma_q_1)
        {
            m_sigma = 0;
            if (0 == len or 0 == text_buf.size())
                return;
            assert(len <= text_buf.size());
            // initialize vectors
            int_vector<64> D(257, 0);
            bit_vector tmp_char(256, 0);
            // count occurrences of each symbol
//            std::cout<<"text=";
            for (size_type i=0; i < len; ++i) {
                ++D[text_buf[i]];
//                std::cout<<(char)text_buf[i];
            }
//            std::cout<<std::endl;
            assert(1 == D[0]); // null-byte should occur exactly once
            m_sigma = 0;
            for (int i=0; i<256; ++i)
                if (D[i]) {
                    tmp_char[i] = 1;    // mark occurring character
                    D[m_sigma] = D[i];  // compactify m_C
                    ++m_sigma;
                }
            m_sigma_q = m_sigma;
            for (uint8_t i=1; i < t_q; ++i) {
                m_sigma_q *= m_sigma;
            }
            m_sigma_q_1 = m_sigma_q/m_sigma;
            // resize to sigma+1, since CSAs also need the sum of all elements
            m_C.C       = int_vector<>(m_sigma+1, 0, bits::hi(len)+1);
            m_C.multi_C = int_vector<>(m_sigma_q+1, 0, bits::hi(len)+1);

            for (int i=(int)m_sigma; i > 0; --i) m_C.C[i] = D[i-1];
            m_C.C[0] = 0;
            for (int i=1; i <= (int)m_sigma; ++i) m_C.C[i] = m_C.C[i] + m_C.C[i-1];
            assert(m_C.C[sigma]==len);
            m_char = tmp_char;
            util::init_support(m_char_rank, &m_char);
            util::init_support(m_char_select, &m_char);
            if (t_q == 1) {
                m_C.multi_C = m_C.C;
            } else if (t_q > 1) {
                int_vector<64> multi_D(m_sigma_q+1, 0);
                // count occurrences of each symbol
                uint64_t x = 0;
                for (size_type i=0; i<q-1; ++i) {
                    x *= m_sigma;
                    x += char2comp[text_buf[i]];
                }
                for (size_type i=q-1; i < len+q-1; ++i) {
                    x *= m_sigma;
                    x += char2comp[text_buf[i%(len)]];
                    x %= m_sigma_q;
                    ++multi_D[x];
//                    std::cout<<"i="<<i<<" x="<<x<<" D[x]="<<multi_D[x]<<std::endl;
                }
                for (size_t i=m_sigma_q; i > 0; --i) {
                    m_C.multi_C[i] = multi_D[i-1];
                }
                m_C.multi_C[0] = 0;
                for (size_t i=1; i <= m_sigma_q; ++i) {
                    m_C.multi_C[i] = m_C.multi_C[i] + m_C.multi_C[i-1];
                }
//                for (size_t i=0; i <= m_sigma_q; ++i) {
//                    std::cout<<"m_C.multi_C["<<i<<"]="<<m_C.multi_C[i]<<std::endl;
//                }
            }
        }

        //! Copy constructor
        succinct_multibyte_alphabet(const succinct_multibyte_alphabet& strat):
            char2comp(this), comp2char(this), C(m_C), sigma(m_sigma), sigma_q(m_sigma_q),
            sigma_q_1(m_sigma_q_1)
        {
            copy(strat);
        }

        //! Move constructor
        succinct_multibyte_alphabet(succinct_multibyte_alphabet&& strat)
        {
            *this = std::move(strat);
        }

        succinct_multibyte_alphabet& operator=(const succinct_multibyte_alphabet& strat)
        {
            if (this != &strat) {
                copy(strat);
            }
            return *this;
        }

        succinct_multibyte_alphabet& operator=(succinct_multibyte_alphabet&& strat)
        {
            if (this != &strat) {
                m_char        = std::move(strat.m_char);
                m_char_rank   = std::move(strat.m_char_rank);
                m_char_rank.set_vector(&m_char);
                m_char_select = std::move(strat.m_char_select);
                m_char_select.set_vector(&m_char);
                m_C           = std::move(strat.m_C);
                m_sigma       = std::move(strat.m_sigma);
                m_sigma_q     = std::move(strat.m_sigma_q);
                m_sigma_q_1   = std::move(strat.m_sigma_q_1);
            }
            return *this;
        }

        //! Swap operator
        void swap(succinct_multibyte_alphabet& strat)
        {
            m_char.swap(strat.m_char);
            util::swap_support(m_char_rank, strat.m_char_rank, &m_char, &(strat.m_char));
            util::swap_support(m_char_select, strat.m_char_select, &m_char, &(strat.m_char));
            m_C.swap(strat.m_C);
            std::swap(m_sigma,strat.m_sigma);
            std::swap(m_sigma_q,strat.m_sigma_q);
            std::swap(m_sigma_q_1,strat.m_sigma_q_1);
        }

        //! Serialize method
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_char.serialize(out, child, "m_char");
            written_bytes += m_char_rank.serialize(out, child, "m_char_rank");
            written_bytes += m_char_select.serialize(out, child, "m_char_select");
            written_bytes += m_C.serialize(out, child, "m_C");
            written_bytes += write_member(m_sigma, out, child, "m_sigma");
            written_bytes += write_member(m_sigma_q, out, child, "m_sigma_q");
            written_bytes += write_member(m_sigma_q_1, out, child, "m_sigma_q_1");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load method
        void load(std::istream& in)
        {
            m_char.load(in);
            m_char_rank.load(in);
            m_char_rank.set_vector(&m_char);
            m_char_select.load(in);
            m_char_select.set_vector(&m_char);
            m_C.load(in);
            read_member(m_sigma, in);
            read_member(m_sigma_q, in);
            read_member(m_sigma_q_1, in);
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
        class char2comp_wrapper_int;
        class comp2char_wrapper_int;
        friend class char2comp_wrapper_int;
        friend class comp2char_wrapper_int;

        typedef int_vector<>::size_type size_type;
        typedef char2comp_wrapper_int       char2comp_type;
        typedef comp2char_wrapper_int       comp2char_type;
        typedef C_array_type            C_type;
        typedef uint64_t                sigma_type;
        typedef uint64_t                char_type;
        typedef uint64_t                comp_char_type;
        typedef std::vector<char_type>  string_type;
        typedef int_alphabet_tag        alphabet_category;
        enum { int_width = 0 };

        //! Helper class for the char2comp mapping
        class char2comp_wrapper_int
        {
            private:
                const int_alphabet* m_strat;
            public:
                char2comp_wrapper_int(const int_alphabet* strat) : m_strat(strat) {}
                comp_char_type operator[](char_type c) const
                {
                    if (m_strat->m_char.size() > 0) {  // if alphabet is not continuous
                        if (c >= m_strat->m_char.size() or !m_strat->m_char[c])
                            return (comp_char_type)0;
                        return (comp_char_type) m_strat->m_char_rank((size_type)c);
                    } else { // direct map if it is continuous
                        if (c >= m_strat->m_sigma)
                            return 0;
                        return (comp_char_type) c;
                    }
                    return 0;
                }
        };

        //! Helper class for the comp2char mapping
        class comp2char_wrapper_int
        {
            private:
                const int_alphabet* m_strat;
            public:
                comp2char_wrapper_int(const int_alphabet* strat) : m_strat(strat) {}
                char_type operator[](comp_char_type c) const
                {
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

        void copy(const int_alphabet& strat)
        {
            m_char        = strat.m_char;
            m_char_rank   = strat.m_char_rank;
            m_char_rank.set_vector(&m_char);
            m_char_select = strat.m_char_select;
            m_char_select.set_vector(&m_char);
            m_C           = strat.m_C;
            m_sigma       = strat.m_sigma;
        }

        //! Check if the alphabet is continuous.
        bool is_continuous_alphabet(std::map<size_type, size_type>& D)
        {
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
        int_alphabet() : char2comp(this), comp2char(this), C(m_C), sigma(m_sigma)
        {
            m_sigma = 0;
        }

        //! Construct from a byte-stream
        /*!
         *  \param text_buf Byte stream.
         *  \param len      Length of the byte stream.
         */
        int_alphabet(int_vector_buffer<0>& text_buf, int_vector_size_type len):
            char2comp(this), comp2char(this), C(m_C), sigma(m_sigma)
        {
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
        int_alphabet(const int_alphabet& strat): char2comp(this), comp2char(this), C(m_C), sigma(m_sigma)
        {
            copy(strat);
        }

        //! Copy constructor
        int_alphabet(int_alphabet&& strat)
        {
            *this = std::move(strat);
        }

        int_alphabet& operator=(const int_alphabet& strat)
        {
            if (this != &strat) {
                copy(strat);
            }
            return *this;
        }

        int_alphabet& operator=(int_alphabet&& strat)
        {
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
        void swap(int_alphabet& strat)
        {
            m_char.swap(strat.m_char);
            util::swap_support(m_char_rank, strat.m_char_rank, &m_char, &(strat.m_char));
            util::swap_support(m_char_select, strat.m_char_select, &m_char, &(strat.m_char));
            m_C.swap(strat.m_C);
            std::swap(m_sigma,strat.m_sigma);
        }

        //! Serialize method
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
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
        void load(std::istream& in)
        {
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
