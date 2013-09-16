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
/*! \file csa_sampling_strategy.hpp
    \brief csa_sampling_strategy.hpp includes different strategy classes for suffix array sampling in the CSAs.
	\author Simon Gog
*/

#ifndef INCLUDED_CSA_SAMPLING_STRATEGY
#define INCLUDED_CSA_SAMPLING_STRATEGY


/*
 *       Text = ABCDEFABCDEF$
 *              0123456789012
 *       sa_sample_dens = 2
 *    *1 SA *2
 *     * 12 *   $
 *       06 *   ABCDEF$
 *     * 00 *   ABCDEFABCDEF$
 *       07     BCDEF$
 *     * 01     BCDEFABCDEF$
 *       08 *   CDEF$
 *     * 02 *   CDEFABCDEF$
 *       09     DEF$
 *     * 03     DEFABCDEF$
 *       10 *   EF$
 *     * 04 *   EFABCDEF$
 *       11     F$
 *     * 05     FABCDEF$
 *
 *    The first sampling (*1) is called suffix order sampling. It has the advantage, that
 *    we don't need to store a bitvector, which marks the sampled suffixes, since a suffix
 *    at index \(i\) in the suffix array is marked if \( 0 \equiv i \mod sa_sample_dens \).
 *
 *   The second sampling (*2) is called text order sampling. It is also called regular in [1].
 *
 * [1] P.Ferragina, J. Siren, R. Venturini: Distribution-Aware Compressed Full-Text Indexes, ESA 2011
 */

#include "int_vector.hpp"
#include "csa_alphabet_strategy.hpp" // for key_trait
#include <set>

namespace sdsl
{

template<class t_csa, uint8_t t_width=0>
class _sa_order_sampling : public int_vector<t_width>
{
    public:
        typedef int_vector<t_width> base_type;
        typedef typename base_type::size_type  size_type;	// make typedefs of base_type visible
        typedef typename base_type::value_type value_type;	//
        enum { sample_dens = t_csa::sa_sample_dens };

        //! Default constructor
        _sa_order_sampling() {}

        //! Constructor
        /*
         * \param cconfig Cache configuration (SA is expected to be cached.).
         * \param csa     Pointer to the corresponding CSA. Not used in this class.
         * \par Time complexity
         *      Linear in the size of the suffix array.
         */
        _sa_order_sampling(const cache_config& cconfig, SDSL_UNUSED const t_csa* csa=nullptr) {
            int_vector_buffer<>  sa_buf(cache_file_name(conf::KEY_SA, cconfig));
            size_type n = sa_buf.size();
            this->width(bits::hi(n)+1);
            this->resize((n+sample_dens-1)/sample_dens);

            for (size_type i=0, cnt_mod=sample_dens, cnt_sum=0; i < n; ++i, ++cnt_mod) {
                size_type sa = sa_buf[i];
                if (sample_dens == cnt_mod) {
                    cnt_mod = 0;
                    (*this)[cnt_sum++] = sa;
                }
            }
        }

        //! Determine if index i is sampled or not
        inline bool is_sampled(size_type i) const {
            return 0 == (i % sample_dens);
        }

        //! Return the suffix array value for the sampled index i
        inline value_type sa_value(size_type i) const {
            return (*this)[i/sample_dens];
        }
};

template<uint8_t t_width=0>
class sa_order_sa_sampling
{
    public:
        template<class t_csa> // template inner class which is used in CSAs to parametrize the
        class type            // sampling strategy class with the Sampling density of the CSA
        {
            public:
                typedef _sa_order_sampling<t_csa, t_width> sample_type;
        };
};


template<class t_csa,
         class bit_vector_type=bit_vector,
         class rank_support_type=typename bit_vector_type::rank_1_type,
         uint8_t t_width=0
         >
class _text_order_sampling : public int_vector<t_width>
{
    private:
        bit_vector_type		m_marked;
        rank_support_type	m_rank_marked;
    public:
        typedef int_vector<t_width> base_type;
        typedef typename base_type::size_type  size_type;	// make typedefs of base_type visible
        typedef typename base_type::value_type value_type;	//
        enum { sample_dens = t_csa::sa_sample_dens };

        //! Default constructor
        _text_order_sampling() {}

        //! Constructor
        /*
         * \param cconfig Cache configuration (SA is expected to be cached.).
         * \param csa    Pointer to the corresponding CSA. Not used in this class.
         * \par Time complexity
         *      Linear in the size of the suffix array.
         */
        _text_order_sampling(const cache_config& cconfig, SDSL_UNUSED const t_csa* csa=nullptr) {
            int_vector_buffer<>  sa_buf(cache_file_name(conf::KEY_SA, cconfig));
            size_type n = sa_buf.size();
            bit_vector marked(n, 0);                // temporary bitvector for the marked text positions
            this->width(bits::hi(n)+1);
            this->resize((n+sample_dens-1)/sample_dens);

            for (size_type i=0, sa_cnt=0; i < n; ++i) {
                size_type sa = sa_buf[i];
                if (0 == (sa % sample_dens)) {
                    marked[i] = 1;
                    (*this)[sa_cnt++] = sa;
                }
            }
            m_marked = std::move(bit_vector_type(marked));
            util::init_support(m_rank_marked, &m_marked);
        }

        //! Copy constructor
        _text_order_sampling(const _text_order_sampling& st) : base_type(st) {
            m_marked = st.m_marked;
            m_rank_marked = st.m_rank_marked;
            m_rank_marked.set_vector(&m_marked);
        }

        //! Determine if index i is sampled or not
        inline bool is_sampled(size_type i) const {
            return m_marked[i];
        }

        //! Return the suffix array value for the sampled index i
        inline value_type sa_value(size_type i) const {
            return (*this)[m_rank_marked(i)];
        }

        //! Assignment operation
        _text_order_sampling& operator=(const _text_order_sampling& st) {
            if (this != &st) {
                base_type::operator=(st);
                m_marked = st.m_marked;
                m_rank_marked = st.m_rank_marked;
                m_rank_marked.set_vector(&m_marked);
            }
            return *this;
        }

        //! Swap operation
        void swap(_text_order_sampling& st) {
            base_type::swap(st);
            m_marked.swap(st.m_marked);
            util::swap_support(m_rank_marked, st.m_rank_marked, &m_marked, &(st.m_marked));
        }

        size_type serialize(std::ostream& out, structure_tree_node* v, std::string name)const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += base_type::serialize(out, child, "samples");
            written_bytes += m_marked.serialize(out, child, "marked");
            written_bytes += m_rank_marked.serialize(out, child, "rank_marked");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            base_type::load(in);
            m_marked.load(in);
            m_rank_marked.load(in);
            m_rank_marked.set_vector(&m_marked);
        }
};

template<class t_bit_vec=bit_vector, class t_rank_sup=typename t_bit_vec::rank_1_type, uint8_t t_width=0>
class text_order_sa_sampling
{
    public:
        template<class t_csa> // template inner class which is used in CSAs to parametrize the
        class type            // sampling strategy class with the Sampling density of the CSA
        {
            public:
                typedef _text_order_sampling<t_csa,
                        t_bit_vec,
                        t_rank_sup,
                        t_width
                        > sample_type;
        };
};

/*
 *       Text = ABCDEFABCDEF$
 *              0123456789012
 *       sa_sample_dens = 4
 *       sa_sample_chars = {B,E}
 *     SA BWT (1)
 *     12  F   * $
 *     06  F     ABCDEF$
 *     00  $   * ABCDEFABCDEF$
 *     07  A     BCDEF$
 *     01  A     BCDEFABCDEF$
 *     08  B   * CDEF$
 *     02  B   * CDEFABCDEF$
 *     09  C     DEF$
 *     03  C     DEFABCDEF$
 *     10  D     EF$
 *     04  D   * EFABCDEF$
 *     11  E   * F$
 *     05  E   * FABCDEF$
 *
 *    In this sampling a suffix x=SA[i] is marked if x \( 0 \equiv x \mod sa_sample_dens \) or
 *    BWT[i] is contained in sa_sample_chars.
 */

template<class t_csa,
         class bit_vector_type=bit_vector,
         class rank_support_type=typename bit_vector_type::rank_1_type,
         uint8_t t_width=0
         >
class _bwt_sampling : public int_vector<t_width>
{
    private:
        bit_vector_type		m_marked;
        rank_support_type	m_rank_marked;
    public:
        typedef int_vector<t_width> base_type;
        typedef typename base_type::size_type  size_type;	// make typedefs of base_type visible
        typedef typename base_type::value_type value_type;	//
        enum { sample_dens = t_csa::sa_sample_dens };

        //! Default constructor
        _bwt_sampling() {}

        //! Constructor
        /*
         * \param cconfig Cache configuration (BWT,SA, and SAMPLE_CHARS are expected to be cached.).
         * \param csa    Pointer to the corresponding CSA. Not used in this class.
         * \par Time complexity
         *      Linear in the size of the suffix array.
         */
        _bwt_sampling(const cache_config& cconfig, SDSL_UNUSED const t_csa* csa=nullptr) {
            int_vector_buffer<>  sa_buf(cache_file_name(conf::KEY_SA, cconfig));
            int_vector_buffer<t_csa::alphabet_type::int_width>
            bwt_buf(cache_file_name(key_trait<t_csa::alphabet_type::int_width>::KEY_BWT,cconfig));
            size_type n = sa_buf.size();
            bit_vector marked(n, 0);                // temporary bitvector for the marked text positions
            this->width(bits::hi(n)+1);
            int_vector<> sample_char;
            typedef typename t_csa::char_type char_type;
            std::set<char_type> char_map;
            if (load_from_cache(sample_char, conf::KEY_SAMPLE_CHAR,cconfig)) {
                for (uint64_t i=0; i<sample_char.size(); ++i) {
                    char_map.insert((char_type)sample_char[i]);
                }
            }
            size_type sa_cnt = 0;
            for (size_type i=0; i < n; ++i) {
                size_type sa  = sa_buf[i];
                char_type bwt = bwt_buf[i];
                if (0 == (sa % sample_dens)) {
                    marked[i] = 1;
                    ++sa_cnt;
                } else if (char_map.find(bwt) != char_map.end()) {
                    marked[i] = 1;
                    ++sa_cnt;
                }
            }
            this->resize(sa_cnt);
            sa_cnt = 0;
            for (size_type i=0; i < n; ++i) {
                size_type sa  = sa_buf[i];
                if (marked[i]) {
                    (*this)[sa_cnt++] = sa;
                }
            }
            util::assign(m_marked, marked);
            util::init_support(m_rank_marked, &m_marked);
        }

        //! Copy constructor
        _bwt_sampling(const _bwt_sampling& st) : base_type(st) {
            m_marked = st.m_marked;
            m_rank_marked = st.m_rank_marked;
            m_rank_marked.set_vector(&m_marked);
        }

        //! Determine if index i is sampled or not
        inline bool is_sampled(size_type i) const {
            return m_marked[i];
        }

        //! Return the suffix array value for the sampled index i
        inline value_type sa_value(size_type i) const {
            return (*this)[m_rank_marked(i)];
        }

        //! Assignment operation
        _bwt_sampling& operator=(const _bwt_sampling& st) {
            if (this != &st) {
                base_type::operator=(st);
                m_marked = st.m_marked;
                m_rank_marked = st.m_rank_marked;
                m_rank_marked.set_vector(&m_marked);
            }
            return *this;
        }

        //! Swap operation
        void swap(_bwt_sampling& st) {
            base_type::swap(st);
            m_marked.swap(st.m_marked);
            util::swap_support(m_rank_marked, st.m_rank_marked, &m_marked, &(st.m_marked));
        }

        size_type serialize(std::ostream& out, structure_tree_node* v, std::string name)const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += base_type::serialize(out, child, "samples");
            written_bytes += m_marked.serialize(out, child, "marked");
            written_bytes += m_rank_marked.serialize(out, child, "rank_marked");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            base_type::load(in);
            m_marked.load(in);
            m_rank_marked.load(in);
            m_rank_marked.set_vector(&m_marked);
        }
};

template<class t_bit_vec=bit_vector, class t_rank_sup=typename t_bit_vec::rank_1_type, uint8_t t_width=0>
class sa_bwt_sampling
{
    public:
        template<class t_csa> // template inner class which is used in CSAs to parametrize the
        class type            // sampling strategy class with the Sampling density of the CSA
        {
            public:
                typedef _bwt_sampling<t_csa,
                        t_bit_vec,
                        t_rank_sup,
                        t_width
                        > sample_type;
        };
};

} // end namespace

#endif
