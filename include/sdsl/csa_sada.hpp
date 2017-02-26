/* sdsl - succinct data structures library
    Copyright (C) 2008-2013 Simon Gog

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
/*! \file csa_sada.hpp
    \brief csa_sada.hpp contains an implementation of the compressed suffix array.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_CSA_SADA
#define INCLUDED_SDSL_CSA_SADA

#include "enc_vector.hpp"
#include "int_vector.hpp"
#include "iterators.hpp"
#include "suffix_array_helper.hpp"
#include "util.hpp"
#include "csa_sampling_strategy.hpp"
#include "csa_alphabet_strategy.hpp"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstring> // for strlen
#include <iomanip>
#include <iterator>

namespace sdsl
{

//! A class for the Compressed Suffix Array (CSA) proposed by Sadakane for practical implementation.
/*!
  *  \tparam t_enc_vec         Space-efficient vector for increasing integer sequences.
  *  \tparam t_dens            Sampling density of SA values
  *  \tparam t_int_dens        Sampling density of ISA values
  *  \tparam t_sa_sample_strat Policy of SA sampling. E.g. sample in SA-order or text-order.
  *  \tparam t_isa             Vector type for ISA sample values.
  *  \tparam t_alphabet_strat  Policy for alphabet representation.
  *
  *  \sa sdsl::csa_wt, sdsl::csa_bitcompressed
  * @ingroup csa
 */
template<class t_enc_vec         = enc_vector<>,          // Vector type used to store the Psi-function
         uint32_t t_dens         = 32,                    // Sample density for suffix array (SA) values
         uint32_t t_inv_dens     = 64,                    // Sample density for inverse suffix array (ISA) values
         class t_sa_sample_strat = sa_order_sa_sampling<>,// Policy class for the SA sampling.
         class t_isa_sample_strat= isa_sampling<>,        // Policy class for ISA sampling.
         class t_alphabet_strat  = byte_alphabet          // Policy class for the representation of the alphabet.
         >
class csa_sada
{
        static_assert(is_enc_vec<t_enc_vec>::value,
                      "First template argument has to be of type env_vector.");
        static_assert(t_dens > 0,
                      "Second template argument has to be greater then 0.");
        static_assert(t_inv_dens > 0,
                      "Third template argument has to be greater then 0.");
        static_assert(std::is_same<typename sampling_tag<t_sa_sample_strat>::type, sa_sampling_tag>::value,
                      "Forth template argument has to be a suffix array sampling strategy.");
        static_assert(std::is_same<typename sampling_tag<t_isa_sample_strat>::type, isa_sampling_tag>::value,
                      "Fifth template argument has to be a inverse suffix array sampling strategy.");
        static_assert(is_alphabet<t_alphabet_strat>::value,
                      "Sixth template argument has to be a alphabet strategy.");

        friend class bwt_of_csa_psi<csa_sada>;
    public:
        enum { sa_sample_dens = t_dens,
               isa_sample_dens = t_inv_dens
             };

        typedef uint64_t                                             value_type;
        typedef random_access_const_iterator<csa_sada>               const_iterator;
        typedef const_iterator                                       iterator;
        typedef const value_type                                     const_reference;
        typedef const_reference                                      reference;
        typedef const_reference*                                     pointer;
        typedef const pointer                                        const_pointer;
        typedef int_vector<>::size_type                              size_type;
        typedef size_type                                            csa_size_type;
        typedef ptrdiff_t                                            difference_type;
        typedef t_enc_vec                                            enc_vector_type;
        typedef enc_vector_type                                      psi_type;
        typedef traverse_csa_psi<csa_sada,false>                     lf_type;
        typedef bwt_of_csa_psi<csa_sada>                             bwt_type;
        typedef isa_of_csa_psi<csa_sada>                             isa_type;
        typedef text_of_csa<csa_sada>                                text_type;
        typedef first_row_of_csa<csa_sada>                           first_row_type;
        typedef typename t_sa_sample_strat::template type<csa_sada>  sa_sample_type;
        typedef typename t_isa_sample_strat::template type<csa_sada> isa_sample_type;
        typedef t_alphabet_strat                                     alphabet_type;
        typedef typename alphabet_type::alphabet_category            alphabet_category;
        typedef typename alphabet_type::comp_char_type               comp_char_type;
        typedef typename alphabet_type::char_type                    char_type; // Note: This is the char type of the CSA not the WT!
        typedef typename alphabet_type::string_type                  string_type;
        typedef csa_sada                                             csa_type;

        typedef csa_tag                                              index_category;
        typedef psi_tag                                              extract_category;

        friend class traverse_csa_psi<csa_sada,true>;
        friend class traverse_csa_psi<csa_sada,false>;

        static const uint32_t linear_decode_limit = 100000;
    private:
        enc_vector_type m_psi;        // psi function
        sa_sample_type  m_sa_sample;  // suffix array samples
        isa_sample_type m_isa_sample; // inverse suffix array samples
        alphabet_type   m_alphabet;   // alphabet component

        mutable std::vector<uint64_t> m_psi_buf; // buffer for decoded psi values

        void copy(const csa_sada& csa)
        {
            m_psi        = csa.m_psi;
            m_sa_sample  = csa.m_sa_sample;
            m_isa_sample = csa.m_isa_sample;
            m_isa_sample.set_vector(&m_sa_sample);
            m_alphabet   = csa.m_alphabet;
        };

        void create_buffer()
        {
            if (enc_vector_type::sample_dens < linear_decode_limit) {
                m_psi_buf = std::vector<uint64_t>(enc_vector_type::sample_dens+1);
            }
        }

    public:
        const typename alphabet_type::char2comp_type& char2comp  = m_alphabet.char2comp;
        const typename alphabet_type::comp2char_type& comp2char  = m_alphabet.comp2char;
        const typename alphabet_type::C_type&         C          = m_alphabet.C;
        const typename alphabet_type::sigma_type&     sigma      = m_alphabet.sigma;
        const psi_type&                               psi        = m_psi;
        const lf_type                                 lf         = lf_type(*this);
        const bwt_type                                bwt        = bwt_type(*this);
        const isa_type                                isa        = isa_type(*this);
        const bwt_type                                L          = bwt_type(*this);
        const first_row_type                          F          = first_row_type(*this);
        const text_type                               text       = text_type(*this);
        const sa_sample_type&                         sa_sample  = m_sa_sample;
        const isa_sample_type&                        isa_sample = m_isa_sample;


        //! Default Constructor
        csa_sada()
        {
            create_buffer();
        }
        //! Default Destructor
        ~csa_sada() { }

        //! Copy constructor
        csa_sada(const csa_sada& csa)
        {
            create_buffer();
            copy(csa);
        }

        //! Move constructor
        csa_sada(csa_sada&& csa)
        {
            *this = std::move(csa);
        }

        csa_sada(cache_config& config);

        //! Number of elements in the \f$\CSA\f$.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         *  \par Time complexity
         *      \f$ \Order{1} \f$
         */
        size_type size()const
        {
            return m_psi.size();
        }

        //! Returns the largest size that csa_sada can ever have.
        /*! Required for the Container Concept of the STL.
         *  \sa size
         */
        static size_type max_size()
        {
            return t_enc_vec::max_size();
        }

        //! Returns if the data structure is empty.
        /*! Required for the Container Concept of the STL.A
         * \sa size
         */
        bool empty()const
        {
            return m_psi.empty();
        }

        //! Swap method for csa_sada
        /*! The swap method can be defined in terms of assignment.
            This requires three assignments, each of which, for a container type, is linear
            in the container's size. In a sense, then, a.swap(b) is redundant.
            This implementation guaranties a run-time complexity that is constant rather than linear.
            \param csa csa_sada to swap.

            Required for the Assignable Concept of the STL.
          */
        void swap(csa_sada& csa);

        //! Returns a const_iterator to the first element.
        /*! Required for the STL Container Concept.
         *  \sa end
         */
        const_iterator begin()const
        {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        /*! Required for the STL Container Concept.
         *  \sa begin.
         */
        const_iterator end()const
        {
            return const_iterator(this, size());
        }

        //! []-operator
        /*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
         * Required for the STL Random Access Container Concept.
         * \par Time complexity
         *      \f$ \Order{s_{SA}\cdot t_{\Psi}} \f$, where every \f$s_{SA}\f$th suffix array entry is sampled and \f$t_{\Psi}\f$
         *           is the access time for an element in the \f$\Psi\f$-function.
         */
        inline value_type operator[](size_type i)const;

        //! Assignment Copy Operator.
        /*!
         *    Required for the Assignable Concept of the STL.
         */
        csa_sada& operator=(const csa_sada& csa)
        {
            if (this != &csa) {
                copy(csa);
            }
            return *this;
        }

        //! Assignment Move Operator.
        /*!
         *    Required for the Assignable Concept of the STL.
         */
        csa_sada& operator=(csa_sada&& csa)
        {
            if (this != &csa) {
                m_psi        = std::move(csa.m_psi);
                m_sa_sample  = std::move(csa.m_sa_sample);
                m_isa_sample = std::move(csa.m_isa_sample);
                m_alphabet   = std::move(csa.m_alphabet);
                m_psi_buf    = std::move(csa.m_psi_buf);
            }
            return *this;
        }

        //! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const;

        //! Load from a stream.
        /*! \param in Input stream to load the data structure from.
         */
        void load(std::istream& in);

        uint32_t get_sample_dens() const
        {
            return t_dens;
        }

    private:

        // Calculates how many symbols c are in the prefix [0..i-1] of the BWT of the original text.
        /*
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in [0..size()]\f$.
         *  \param c The symbol to count the occurrences in the prefix.
         *    \returns The number of occurrences of symbol c in the prefix [0..i-1] of the BWT.
         *  \par Time complexity
         *        \f$ \Order{\log n t_{\Psi}} \f$
         */
        size_type rank_bwt(size_type i, const char_type c)const
        {
            comp_char_type cc = char2comp[c];
            if (cc==0 and c!=0)  // character is not in the text => return 0
                return 0;
            if (i == 0)
                return 0;
            assert(i <= size());

            size_type lower_b, upper_b; // lower_b inclusive, upper_b exclusive

            const size_type sd = m_psi.get_sample_dens();
            size_type lower_sb = (C[cc]+sd-1)/sd; // lower_sb inclusive
            size_type upper_sb = (C[cc+1]+sd-1)/sd; // upper_sb exclusive
            while (lower_sb+1 < upper_sb) {
                size_type mid = (lower_sb+upper_sb)/2;
                if (m_psi.sample(mid) >= i)
                    upper_sb = mid;
                else
                    lower_sb = mid;
            }

            if (lower_sb == upper_sb) { // the interval was smaller than sd
                lower_b = C[cc]; upper_b = C[cc+1];
            } else if (lower_sb > (C[cc]+sd-1)/sd) { // main case
// TODO: don't use get_inter_sampled_values if t_dens is really
//       large
                lower_b = lower_sb*sd;
                if (0 == m_psi_buf.size()) {
                    upper_b = std::min(upper_sb*sd, C[cc+1]);
                    goto finish;
                }
                uint64_t* p = m_psi_buf.data();
                // extract the psi values between two samples
                m_psi.get_inter_sampled_values(lower_sb, p);
                p = m_psi_buf.data();
                uint64_t smpl = m_psi.sample(lower_sb);
                // handle border cases
                if (lower_b + m_psi.get_sample_dens() >= C[cc+1])
                    m_psi_buf[ C[cc+1]-lower_b ] = size()-smpl;
                else
                    m_psi_buf[ m_psi.get_sample_dens() ] = size()-smpl;
                // search the result linear
                while ((*p++)+smpl < i);

                return p-1-m_psi_buf.data() + lower_b - C[cc];
            } else { // lower_b == (m_C[cc]+sd-1)/sd and lower_sb < upper_sb
                if (m_psi.sample(lower_sb) >= i) {
                    lower_b = C[cc];
                    upper_b = lower_sb * sd + 1;
                } else {
                    lower_b = lower_sb * sd;
                    upper_b = std::min(upper_sb*sd, C[cc+1]);
                }
            }
finish:
            // binary search the interval [C[cc]..C[cc+1]-1] for the result
//            size_type lower_b = m_C[cc], upper_b = m_C[cc+1]; // lower_b inclusive, upper_b exclusive
            while (lower_b+1 < upper_b) {
                size_type mid = (lower_b+upper_b)/2;
                if (m_psi[mid] >= i)
                    upper_b = mid;
                else
                    lower_b = mid;
            }
            if (lower_b > C[cc])
                return lower_b - C[cc] + 1;
            else { // lower_b == m_C[cc]
                return m_psi[lower_b] < i;// 1 if m_psi[lower_b]<i, 0 otherwise
            }
        }

        // Calculates the position of the i-th c in the BWT of the original text.
        /*
         *  \param i The i-th occurrence. \f$i\in [1..rank_bwt(size(),c)]\f$.
         *  \param c Symbol c.
         *    \returns The position of the i-th c in the BWT or size() if c does occur less then i times.
         *  \par Time complexity
         *        \f$ \Order{t_{\Psi}} \f$
         */
        size_type select_bwt(size_type i, const char_type c)const
        {
            assert(i > 0);
            comp_char_type cc = char2comp[c];
            if (cc==0 and c!=0)  // character is not in the text => return 0
                return size();
            assert(cc != 255);
            if (C[cc]+i-1 <  C[cc+1]) {
                return m_psi[C[cc]+i-1];
            } else
                return size();
        }
};

// == template functions ==

template<class t_enc_vec, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
csa_sada<t_enc_vec, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::csa_sada(cache_config& config)
{
    create_buffer();
    if (!cache_file_exists(key_trait<alphabet_type::int_width>::KEY_BWT, config)) {
        return;
    }
    int_vector_buffer<alphabet_type::int_width> bwt_buf(cache_file_name(key_trait<alphabet_type::int_width>::KEY_BWT,config));
    size_type n = bwt_buf.size();
    {
        auto event = memory_monitor::event("construct csa-alpbabet");
        alphabet_type tmp_alphabet(bwt_buf, n);
        m_alphabet.swap(tmp_alphabet);
    }

    int_vector<> cnt_chr(sigma, 0, bits::hi(n)+1);
    for (typename alphabet_type::sigma_type i=0; i < sigma; ++i) {
        cnt_chr[i] = C[i];
    }
    // calculate psi
    {
        auto event = memory_monitor::event("construct PSI");
        // TODO: move PSI construct into construct_PSI.hpp
        int_vector<> psi(n, 0, bits::hi(n)+1);
        for (size_type i=0; i < n; ++i) {
            psi[ cnt_chr[ char2comp[bwt_buf[i]] ]++ ] = i;
        }
        std::string psi_file = cache_file_name(conf::KEY_PSI, config);
        if (!store_to_cache(psi, conf::KEY_PSI, config)) {
            return;
        }
    }
    {
        auto event = memory_monitor::event("encode PSI");
        int_vector_buffer<> psi_buf(cache_file_name(conf::KEY_PSI, config));
        t_enc_vec tmp_psi(psi_buf);
        m_psi.swap(tmp_psi);
    }
    {
        auto event = memory_monitor::event("sample SA");
        sa_sample_type tmp_sa_sample(config);
        m_sa_sample.swap(tmp_sa_sample);
    }
    {
        auto event = memory_monitor::event("sample ISA");
        isa_sample_type isa_s(config, &m_sa_sample);
        util::swap_support(m_isa_sample, isa_s, &m_sa_sample, (const sa_sample_type*)nullptr);
    }
}

template<class t_enc_vec, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
inline auto csa_sada<t_enc_vec, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::operator[](size_type i)const -> value_type
{
    size_type off = 0;
    while (!m_sa_sample.is_sampled(i)) {  // while i mod t_dens != 0 (SA[i] is not sampled)
        i = psi[i];                       // go to the position where SA[i]+1 is located
        ++off;                            // add 1 to the offset
    }
    value_type result = m_sa_sample[i];
    if (result < off) {
        return m_psi.size()-(off-result);
    } else
        return result-off;
}


template<class t_enc_vec, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
auto csa_sada<t_enc_vec, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const -> size_type
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_psi.serialize(out, child, "psi");
    written_bytes += m_sa_sample.serialize(out, child, "sa_samples");
    written_bytes += m_isa_sample.serialize(out, child, "isa_samples");
    written_bytes += m_alphabet.serialize(out, child, "alphabet");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<class t_enc_vec, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
void csa_sada<t_enc_vec, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::load(std::istream& in)
{
    m_psi.load(in);
    m_sa_sample.load(in);
    m_isa_sample.load(in, &m_sa_sample);
    m_alphabet.load(in);
}

template<class t_enc_vec, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
void csa_sada<t_enc_vec, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::swap(csa_sada<t_enc_vec, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>& csa)
{
    if (this != &csa) {
        m_psi.swap(csa.m_psi);
        m_sa_sample.swap(csa.m_sa_sample);
        util::swap_support(m_isa_sample, csa.m_isa_sample, &m_sa_sample, &(csa.m_sa_sample));
        m_alphabet.swap(csa.m_alphabet);
    }
}

} // end namespace sdsl
#endif
