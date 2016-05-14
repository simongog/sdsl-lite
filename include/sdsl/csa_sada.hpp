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
#include "enc_vector2.hpp"
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
        const alphabet_type&                          alphabet   = m_alphabet;
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

        //! Returns if the data strucutre is empty.
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

            Required for the Assignable Conecpt of the STL.
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

        // Calculates how many symbols cc are in the prefix [0..i-1] of the BWT of the original text.
        /*
         *  \param i  The exclusive index of the prefix range [0..i-1], so \f$i\in [0..size()]\f$.
         *  \param cc The compactified symbol to count in the prefix.
         *  \returns The number of occurrences of the compactified symbol cc in the prefix [0..i-1].
         *  \par Time complexity
         *        \f$ \Order{\log n t_{\Psi}} \f$
         */
        template<typename t_char>
        size_type rank_comp_bwt(size_type i, const t_char cc)const
        {
//            std::cout<<"rank_comp_bwt("<<i<<" (cc="<<cc<<")"<<std::endl;
            if (i == 0)
                return 0;
            assert(i <= size());
            const auto cc_begin = C[cc];   // begin of interval of context cc (inclusive)
            const auto cc_end   = C[cc+1]; // end of interval of context cc (exclusive)
            const size_type sd  = m_psi.get_sample_dens();
            size_type s_begin   = (cc_begin+sd-1)/sd; // first sample at or after cc_begin
            size_type s_end     = (cc_end+sd-1)/sd;   // first sample at or after cc_end
//            std::cout<<"cc_begin = "<<cc_begin<<" cc_end = "<<cc_end<<" cc_size="<<cc_end-cc_begin<<std::endl;
//            std::cout<<"s_begin = "<<s_begin<<" s_end = "<<s_end<<std::endl;
//            if(s_end - s_begin < 10){
//                std::cout<<"samples in C range: ";
//                for(size_t k=s_begin; k<s_end; ++k){
//                    std::cout<<m_psi.sample(k)<<" (@ "<<s_begin*sd<<") ";
//                }
//                std::cout<<std::endl;
//            }

            if (s_begin == s_end) {
                // Case (1): No sample inside [cc_begin, cc_end)
                //           => search in previous block (s_begin-1)
//                std::cout<<"case (1)"<<std::endl;
            } else if (m_psi.sample(s_begin) >= i) {  // now s_begin < s_end
                // Case (2): Some samples inside [cc_begin, cc_end)
                //           and first sample already larger or equal to i
                //           => search in previous block (s_begin-1)
//                std::cout<<"case (2): "<<m_psi.sample(s_begin)<<" >= " << i << std::endl;
            } else { // still s_begin < s_end
                // Case (3): Some samples inside [cc_begin, cc_end)
                //           and first sample smaller than i
                //           => binary search for first sample >= i
                s_begin = upper_bound(s_begin, s_end, i-1);
                //           => search in previous block (s_begin-1)
//                std::cout<<"case (3): s_begin = " << s_begin << " (s_end=" << s_end <<" )"<< std::endl;
//                std::cout<<">>>>> m_psi.sample(s_begin-1)="<<m_psi.sample(s_begin-1)<<std::endl;
            }
            s_begin -= 1;
            uint64_t smpl = m_psi.sample(s_begin);

            size_t abs_decode_begin = s_begin*sd;
            size_t skip = 0;
            if (abs_decode_begin < cc_begin) {
                skip = cc_begin - abs_decode_begin;
            }
            size_t res = abs_decode_begin + skip - cc_begin;

            if ((s_begin+1)*sd < m_psi.size() and skip == 0 and smpl+sd == m_psi.sample(s_begin+1)) {
//std::cout<<"!!!Special case"<<std::endl;
//std::cout<<"s_begin="<<s_begin<<std::endl;
//std::cout<<"abs_decode_begin="<<abs_decode_begin<<" cc_begin="<<cc_begin<<std::endl;
//std::cout<<"RES="<<res + (i - smpl)<<" res="<<res<<" i="<<i<<" smpl="<<smpl<<std::endl;
                return res + (i - smpl);
            }

            uint64_t* p = m_psi_buf.data();
            // extract the psi values between two samples
            m_psi.get_inter_sampled_values(s_begin, p);
            p = m_psi_buf.data();

            for (auto it = p + skip; (res < cc_end - cc_begin) and it < m_psi_buf.data()+sd; ++it) {
                if ((*it)+smpl >= i)
                    break;
                ++res;
            }
            return res;
        }

        template<typename t_char>
        std::tuple<size_type,size_type> double_rank_comp_bwt(size_type i, size_type j, const t_char cc)const
        {
//            std::cout<<"double_rank_comp_bwt("<<i<<","<<j<<" (cc="<<cc<<")"<<std::endl;
//            return std::make_tuple(rank_comp_bwt(i,cc), rank_comp_bwt(j,cc));
            if (i == 0)
                return std::make_tuple(0, rank_comp_bwt(j,cc));
            assert(i <= size());
            const auto cc_begin = C[cc];   // begin of interval of context cc (inclusive)
            const auto cc_end   = C[cc+1]; // end of interval of context cc (exclusive)
            const size_type sd  = m_psi.get_sample_dens();
            size_type s_begin   = (cc_begin+sd)/sd; // first sample after cc_begin
            size_type s_end     = (cc_end+sd-1)/sd;   // first sample at or after cc_end
            bool answer_j       = false;

            if (s_begin == s_end) {
                // Case (1): No sample inside [cc_begin, cc_end)
                //           => search in previous block (s_begin-1)
                answer_j = true;
            } else if (m_psi.sample(s_begin) >= i) {  // now s_begin < s_end
                // Case (2): Some samples inside [cc_begin, cc_end)
                //           and first sample already larger or equal to i
                //           => search in previous block (s_begin-1)
                answer_j = (m_psi.sample(s_begin) >= j);
            } else { // still s_begin < s_end
                // Case (3): Some samples inside [cc_begin, cc_end)
                //           and first sample smaller than i
                //           => binary search for first sample >= i
                s_begin = upper_bound(s_begin, s_end, i-1);
                //           => search in previous block (s_begin-1)
                answer_j = (s_begin == s_end) or (m_psi.sample(s_begin) >=j);
            }
            s_begin -= 1;
            uint64_t smpl = m_psi.sample(s_begin);

            size_t abs_decode_begin = s_begin*sd;
            size_t skip = 0;
            if (abs_decode_begin < cc_begin) {
                skip = cc_begin - abs_decode_begin;
            }
            size_t res = abs_decode_begin + skip - cc_begin;

            bool uniform_block = (s_begin+1)*sd < m_psi.size() and skip == 0 and smpl+sd == m_psi.sample(s_begin+1);
            if (uniform_block) {
                if (answer_j) {
                    return std::make_tuple(res + (i - smpl), res + (j - smpl));
                } else {
                    return std::make_tuple(res + (i - smpl), rank_comp_bwt(j, cc));
                }
            }

            uint64_t* p = m_psi_buf.data();
            // extract the psi values between two samples
            m_psi.get_inter_sampled_values(s_begin, p);
            p = m_psi_buf.data();

            auto it = p + skip;
            for (; (res < cc_end - cc_begin) and it < m_psi_buf.data()+sd; ++it) {
                if ((*it)+smpl >= i) {
                    break;
                }
                ++res;
            }
            if (answer_j) {
                size_t res2 = res;
                for (; (res2 < cc_end - cc_begin) and it < m_psi_buf.data()+sd; ++it) {
                    if ((*it)+smpl >= j) {
                        break;
                    }
                    ++res2;
                }
                return std::make_tuple(res, res2);
            }
            return std::make_tuple(res, rank_comp_bwt(j, cc));
        }

    private:

        template<typename V>
        size_t upper_bound(size_t first, size_t last, V value) const
        {
            size_t mid;
            size_t count, step;
            count = last-first;

            while (count > 0) {
                mid = first;
                step = count / 2;
                mid += step;
                if (!(value < m_psi.sample(mid))) {
                    first = ++mid;
                    count -= step + 1;
                } else count = step;
            }
            return first;
        }

        // Calculates how many symbols c are in the prefix [0..i-1] of the BWT of the original text.
        /*
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in [0..size()]\f$.
         *  \param c The symbol to count in the prefix.
         *    \returns The number of occurrences of symbol c in the prefix [0..i-1] of the BWT.
         *  \par Time complexity
         *        \f$ \Order{\log n t_{\Psi}} \f$
         */
        // replace const char_type c by const std::array<char_type, alphabet_type::C_depth>& c
        template<typename t_char>
        size_type rank_bwt(size_type i, const t_char c)const
        {
            auto cc = char2comp[c];
            if (cc==0 and c!=0) // character is not in the text => return 0
                return 0;
            if (i == 0)
                return 0;
            return rank_comp_bwt(i, cc);
        }

        template<typename t_char>
        std::array<size_type,2>
        rank_bwt(std::array<size_type,2> ij, const t_char c)const
        {
            return {rank_bwt(ij[0], c), rank_bwt(ij[1],c)};
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
//        alphabet_type tmp_alphabet(bwt_buf, n); // TODO: maybe it is possible to use _buf_buf again for multibyte!!
        int_vector_buffer<alphabet_type::int_width> text_buf(cache_file_name(key_trait<alphabet_type::int_width>::KEY_TEXT,config));
        alphabet_type tmp_alphabet(text_buf, n);
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
        /*
                enc_vector<coder::elias_delta, enc_vector_type::sample_dens> m_psi_check(psi_buf);
                if ( m_psi_check.size() != m_psi.size() ){
                    std::cout<<"m_psi.size()="<<m_psi.size()<<"!="<<m_psi_check.size()<<" m_psi_check.size()"<<std::endl;
                } else {

                    std::vector<uint64_t> buf1 = std::vector<uint64_t>(enc_vector_type::sample_dens+1);
                    std::vector<uint64_t> buf2 = std::vector<uint64_t>(enc_vector_type::sample_dens+1);

                    std::cout<<"m_psi.size()="<<m_psi.size()<<std::endl;
                    for(size_t i=0; i<m_psi.size()/enc_vector_type::sample_dens; ++i){
                        if ( m_psi.sample(i) != m_psi_check.sample(i) ) {
                            std::cout<<"m_psi.sample(i) != m_psi_check.sample(i) for i="<<i<<" "<<m_psi.sample(i)<<"!="<<m_psi_check.sample(i)<<std::endl;
                        }
                        m_psi.get_inter_sampled_values(i, buf1.data());
                        m_psi_check.get_inter_sampled_values(i, buf2.data());
                        bool error = false;
                        for(size_t j=0; j<enc_vector_type::sample_dens and i*enc_vector_type::sample_dens+j<m_psi.size(); ++j) {
                            if ( buf1[j] != buf2[j] ) {
                                std::cout<<"i="<<i<<" j="<<j<<" buf1[j]="<<buf1[j]<<" buf2[j]="<<buf2[j]<<std::endl;
                                error = true;
                            }
                        }
                        if (error) {
                            std::cout<<" m_psi.sample(i)="<<m_psi.sample(i)<<std::endl;
                            std::cout<<" m_psi.sample(i+1)="<<m_psi.sample(i+1)<<" m_psi_check.sample(i+1)="<<m_psi_check.sample(i+1)<<std::endl;
                            throw std::logic_error("error");
                        }
                    }
                }
        */
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
