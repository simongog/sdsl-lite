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
/*! \file csa_sada2.hpp
    \brief csa_sada2.hpp contains an implementation of the compressed suffix array.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_CSA_SADAII
#define INCLUDED_SDSL_CSA_SADAII

#include "bit_vectors.hpp"
#include "int_vector.hpp"
#include "iterators.hpp"
#include "suffix_array_helper.hpp"
#include "util.hpp"
#include "io.hpp"
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

template<typename t_hyb_vec,
         typename t_csa
         >
class uef_psi_support
{
    public:
        typedef typename bit_vector::size_type                size_type;
        typedef size_type                                     value_type;
        typedef typename t_csa::alphabet_type                 alphabet_type;
        typedef typename alphabet_type::comp_char_type        comp_char_type;
        typedef typename alphabet_type::C_type                C_type;
        typedef random_access_const_iterator<uef_psi_support> iterator;
        typedef iterator                                      const_iterator;
        typedef const value_type                              reference;
        typedef const value_type                              const_reference;
        typedef const value_type*                             const_pointer;
        typedef ptrdiff_t                                     difference_type;
        typedef csa_member_tag                                category;
        typedef int_alphabet_tag                              alphabet_category;
        typedef wt_huff_int<bit_vector,
                rank_support_v<>,
                select_support_scan<1>,
                select_support_scan<0>>                sml_wt_type;

    private:
        std::vector<t_hyb_vec>                         m_inc_seq;
        std::vector<typename t_hyb_vec::rank_1_type>   m_inc_seq_rank;
        std::vector<typename t_hyb_vec::select_1_type> m_inc_seq_sel;
        bit_vector                                     m_sml;        // indicates if a context is small or large
        rank_support_v5<>                              m_sml_rank;   // rank for m_sml
        sml_wt_type                                    m_sml_wt;     // wt to get rank to index into
        std::vector<int_vector<>>                      m_sml_inc_seq;// small sequences

        const t_csa*                                   m_csa;

        void set_vector()
        {
            for (size_t i=0; i<m_inc_seq_rank.size(); ++i) {
                m_inc_seq_rank[i].set_vector(&(m_inc_seq[i]));
                m_inc_seq_sel[i].set_vector(&(m_inc_seq[i]));
            }
        }
    public:

        uef_psi_support(const t_csa* csa=nullptr)
        {
            set_vector(csa);
        }

        uef_psi_support(int_vector_buffer<>& psi_buf, const t_csa* csa)
        {
//            std::cout<<"Hello!!!!"<<std::endl;
            g_saved_bits=0;
            set_vector(csa);
            const auto& C = m_csa->C;

            m_sml = bit_vector(C.size()-1,0);
            const auto threshold = t_hyb_vec::block_size;
// (1)      Determine the number of small blocks
            for (size_t i=0; i<C.size()-1; ++i) {
                m_sml[i] = (C[i+1]-C[i]) < threshold;
            }
            m_sml_rank = decltype(m_sml_rank)(&m_sml);
            size_t sigma_small = m_sml_rank(C.size()-1);
            size_t sigma_large = C.size()-1-sigma_small;
            {
                int_vector<> sml(sigma_small, 0, bits::hi(threshold)+1);

// (2)          Create a vector containing only the small context sizes
                for (size_t i=0, ii=0; i<C.size()-1; ++i) {
                    if (m_sml[i] == 1) {
                        sml[ii++] = C[i+1]-C[i];
                    }
                }
// (3)          Greate WT over sml
                construct_im(m_sml_wt, sml, 0);
            }
// (4)      Initialize m_sml_inc_seq
            m_sml_inc_seq.resize(threshold);
            for (uint64_t cs=1; cs<threshold; ++cs) {
                auto size = cs * m_sml_wt.rank(m_sml_wt.size(), cs);
                m_sml_inc_seq[cs] = int_vector<>(size, 0, bits::hi(m_csa->size())+1);
            }

// (5)      Initialize m_inc_seq (to store the larger contexts
            m_inc_seq.resize(sigma_large);
            m_inc_seq_rank.resize(sigma_large);
            m_inc_seq_sel.resize(sigma_large);
            for (size_t i=0,i0=0,i1=0; i<C.size()-1; ++i) {
                int_vector<> v(C[i+1]-C[i]);
                for (size_t j=C[i]; j<C[i+1]; ++j) {
                    v[j-C[i]] = psi_buf[j];
                }
                if (m_sml[i]) {
                    auto rank = m_sml_wt.rank(i1++, v.size());
                    auto start_pos = rank * v.size();
                    for (size_t j=0; j<v.size(); ++j) {
                        m_sml_inc_seq[v.size()][j+start_pos] = v[j];
                    }
                } else {
                    t_hyb_vec tmp(v.begin(), v.end());
                    m_inc_seq[i0++].swap(tmp);
                }
            }
            set_vector();
        }

        uef_psi_support& operator=(const uef_psi_support& psi)
        {
            if (this != &psi) {
                m_inc_seq      = psi.m_inc_seq;
                m_inc_seq_rank = psi.m_inc_seq_rank;
                m_inc_seq_sel  = psi.m_inc_seq_sel;
                m_sml          = psi.m_sml;
                m_sml_rank     = psi.m_sml_rank;
                m_sml_rank.set_vector(&m_sml);
                m_sml_wt       = psi.m_sml_wt;
                m_sml_inc_seq  = psi.m_sml_inc_seq;
                set_vector();
                set_vector(psi.m_csa);
            }
            return *this;
        }

        uef_psi_support& operator=(uef_psi_support&& psi)
        {
            if (this != &psi) {
                m_inc_seq      = std::move(psi.m_inc_seq);
                m_inc_seq_rank = std::move(psi.m_inc_seq_rank);
                m_inc_seq_sel  = std::move(psi.m_inc_seq_sel);
                m_sml          = std::move(psi.m_sml);
                m_sml_rank     = std::move(psi.m_sml_rank);
                m_sml_rank.set_vector(&m_sml);
                m_sml_wt       = std::move(psi.m_sml_wt);
                m_sml_inc_seq  = std::move(psi.m_sml_inc_seq);
                set_vector();
                set_vector(psi.m_csa);
            }
            return *this;
        }

        void set_vector(const t_csa* csa)
        {
            m_csa = csa;
        }

        uint64_t rank(uint64_t i, comp_char_type cc) const
        {
            if (m_sml[cc]) {
                auto cc_sml  = m_sml_rank(cc);
                size_type cs = m_csa->C[cc+1] - m_csa->C[cc]; // context size
                auto rank = m_sml_wt.rank(cc_sml, cs);
                size_type begin = rank*cs;
                for (size_t j=0; j<cs; ++j) {
                    if (m_sml_inc_seq[cs][begin+j] >= i)
                        return j;
                }
                return cs;
            } else {
//                std::cout<<"single_rank: for i="<<i<<std::endl;
                size_type cc_large  = cc - m_sml_rank(cc);
                return m_inc_seq_rank[cc_large](i);
            }
        }

        std::array<uint64_t,2> rank(std::array<uint64_t,2> ij, comp_char_type cc) const
        {
            if (m_sml[cc]) {
                auto cc_sml  = m_sml_rank(cc);
                size_type cs = m_csa->C[cc+1] - m_csa->C[cc]; // context size
                auto rnk = m_sml_wt.rank(cc_sml, cs);
                size_type begin = rnk*cs;
                std::array<uint64_t,2> res = {0,0};
                size_t j=0;
                for (size_t k=0; k<2; ++k) {
                    while (j < cs and  m_sml_inc_seq[cs][begin+j] < ij[k]) {
                        ++j;
                    }
                    res[k] = j;
                }
//                std::array<uint64_t,2> res2 = {rank(ij[0],cc),rank(ij[1],cc)};
//                if ( res != res2 ){
//                    std::cout<<"double rank: res=["<<res[0]<<","<<res[1]<<"] != ";
//                    std::cout<<"["<<res2[0]<<","<<res2[1]<<"] for"<<
//                        ij[0]<<" and "<<ij[1]<<std::endl;
//                }
                return res;
            } else {
                size_type cc_large  = cc - m_sml_rank(cc);
                auto res = m_inc_seq_rank[cc_large](ij);
//                std::array<uint64_t,2> res2 = {rank(ij[0],cc),rank(ij[1],cc)};
//                std::cout<<"_double rank: res=["<<res[0]<<","<<res[1]<<"] != ";
//                std::cout<<"["<<res2[0]<<","<<res2[1]<<"] for "<<
//                    ij[0]<<" and "<<ij[1]<<std::endl;
                return res;
            }
        }

        uint64_t select(uint64_t i, comp_char_type cc) const
        {
            if (m_sml[cc]) {
                auto cc_sml  = m_sml_rank(cc);
                size_type cs = m_csa->C[cc+1] - m_csa->C[cc]; // context size
                auto rank = m_sml_wt.rank(cc_sml, cs);
                return m_sml_inc_seq[cs][rank*cs+(i-1)];
            } else {
                size_type cc_large  = cc - m_sml_rank(cc);
                return m_inc_seq_sel[cc_large](i);
            }
        }

        value_type operator[](const size_type i) const
        {
            size_t cc = std::upper_bound(m_csa->C.begin(), m_csa->C.end(),i) - m_csa->C.begin() - 1;
            size_t cum_sum = m_csa->C[cc];
            if (m_sml[cc]) {
                auto cc_sml  = m_sml_rank(cc);
                size_type cs = m_csa->C[cc+1] - cum_sum; // context size
                auto rank = m_sml_wt.rank(cc_sml, cs);
                return m_sml_inc_seq[cs][rank*cs+(i-cum_sum)];
            } else {
                size_type cc_large  = cc - m_sml_rank(cc);
                return m_inc_seq_sel[cc_large](i-cum_sum+1);
            }
        }

        size_type size() const
        {
            return m_csa->size();
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += sdsl::serialize(m_inc_seq, out, child, "inc_seq");
            written_bytes += sdsl::serialize(m_inc_seq_rank, out, child, "inc_seq_rank");
            written_bytes += sdsl::serialize(m_inc_seq_sel, out, child, "inc_seq_rank");
            written_bytes += sdsl::serialize(m_sml, out, child, "sml");
            written_bytes += sdsl::serialize(m_sml_rank, out, child, "sml_rank");
            written_bytes += sdsl::serialize(m_sml_wt, out, child, "sml_wt");
            written_bytes += sdsl::serialize(m_sml_inc_seq, out, child, "sml_inc_seq");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in, const t_csa* csa = nullptr)
        {
            sdsl::load(m_inc_seq, in);
            sdsl::load(m_inc_seq_rank, in);
            sdsl::load(m_inc_seq_sel, in);
            sdsl::load(m_sml, in);
            sdsl::load(m_sml_rank, in);
            m_sml_rank.set_vector(&m_sml);
            sdsl::load(m_sml_wt, in);
            sdsl::load(m_sml_inc_seq, in);
            set_vector();
            set_vector(csa);
        }


        void swap(uef_psi_support& v)
        {
            m_inc_seq.swap(v.m_inc_seq);
            m_inc_seq_rank.swap(v.m_inc_seq_rank);
            m_inc_seq_sel.swap(v.m_inc_seq_sel);
            m_sml.swap(v.m_sml);
            util::swap_support(m_sml_rank, v.m_sml_rank, &m_sml, &(v.m_sml));
            m_sml_wt.swap(v.m_sml_wt);
            m_sml_inc_seq.swap(v.m_sml_inc_seq);
            set_vector();
        }

        const const_iterator begin()const
        {
            return const_iterator(this, 0);
        }

        const const_iterator end()const
        {
            return const_iterator(this, size());
        }

};

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
template<class t_hyb_sd          = hyb_sd_vector<>,       // Vector type used to store the Psi-function
         uint32_t t_dens         = 32,                    // Sample density for suffix array (SA) values
         uint32_t t_inv_dens     = 64,                    // Sample density for inverse suffix array (ISA) values
         class t_sa_sample_strat = sa_order_sa_sampling<>,// Policy class for the SA sampling.
         class t_isa_sample_strat= isa_sampling<>,        // Policy class for ISA sampling.
         class t_alphabet_strat  = byte_alphabet          // Policy class for the representation of the alphabet.
         >
class csa_sada2
{
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

        friend class bwt_of_csa_psi<csa_sada2>;
    public:
        enum { sa_sample_dens = t_dens,
               isa_sample_dens = t_inv_dens
             };

        typedef uint64_t                                             value_type;
        typedef random_access_const_iterator<csa_sada2>               const_iterator;
        typedef const_iterator                                       iterator;
        typedef const value_type                                     const_reference;
        typedef const_reference                                      reference;
        typedef const_reference*                                     pointer;
        typedef const pointer                                        const_pointer;
        typedef int_vector<>::size_type                              size_type;
        typedef size_type                                            csa_size_type;
        typedef ptrdiff_t                                            difference_type;
        typedef traverse_csa_psi<csa_sada2,false>                     lf_type;
        typedef bwt_of_csa_psi<csa_sada2>                             bwt_type;
        typedef isa_of_csa_psi<csa_sada2>                             isa_type;
        typedef text_of_csa<csa_sada2>                                text_type;
        typedef first_row_of_csa<csa_sada2>                           first_row_type;
        typedef typename t_sa_sample_strat::template type<csa_sada2>  sa_sample_type;
        typedef typename t_isa_sample_strat::template type<csa_sada2> isa_sample_type;
        typedef t_alphabet_strat                                     alphabet_type;
        typedef typename alphabet_type::alphabet_category            alphabet_category;
        typedef typename alphabet_type::comp_char_type               comp_char_type;
        typedef typename alphabet_type::char_type                    char_type; // Note: This is the char type of the CSA not the WT!
        typedef typename alphabet_type::string_type                  string_type;
        typedef csa_sada2                                            csa_type;

        typedef csa_tag                                              index_category;
        typedef psi_tag                                              extract_category;
        typedef uef_psi_support<t_hyb_sd, csa_sada2>                 psi_type;

        friend class traverse_csa_psi<csa_sada2,true>;
        friend class traverse_csa_psi<csa_sada2,false>;

    private:
        psi_type m_psi_support;        // psi function
        sa_sample_type  m_sa_sample;  // suffix array samples
        isa_sample_type m_isa_sample; // inverse suffix array samples
        alphabet_type   m_alphabet;   // alphabet component

        void copy(const csa_sada2& csa)
        {
            m_psi_support = csa.m_psi_support;
            m_psi_support.set_vector(this);
            m_sa_sample  = csa.m_sa_sample;
            m_isa_sample = csa.m_isa_sample;
            m_isa_sample.set_vector(&m_sa_sample);
            m_alphabet   = csa.m_alphabet;
        };

    public:
        const typename alphabet_type::char2comp_type& char2comp  = m_alphabet.char2comp;
        const typename alphabet_type::comp2char_type& comp2char  = m_alphabet.comp2char;
        const typename alphabet_type::C_type&         C          = m_alphabet.C;
        const typename alphabet_type::sigma_type&     sigma      = m_alphabet.sigma;
        const alphabet_type&                          alphabet   = m_alphabet;
        const psi_type&                               psi        = m_psi_support;
        const lf_type                                 lf         = lf_type(*this);
        const bwt_type                                bwt        = bwt_type(*this);
        const isa_type                                isa        = isa_type(*this);
        const bwt_type                                L          = bwt_type(*this);
        const first_row_type                          F          = first_row_type(*this);
        const text_type                               text       = text_type(*this);
        const sa_sample_type&                         sa_sample  = m_sa_sample;
        const isa_sample_type&                        isa_sample = m_isa_sample;


        //! Default Constructor
        csa_sada2() { }
        //! Default Destructor
        ~csa_sada2() { }

        //! Copy constructor
        csa_sada2(const csa_sada2& csa)
        {
            copy(csa);
        }

        //! Move constructor
        csa_sada2(csa_sada2&& csa)
        {
            *this = std::move(csa);
        }

        csa_sada2(cache_config& config);

        //! Number of elements in the \f$\CSA\f$.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         *  \par Time complexity
         *      \f$ \Order{1} \f$
         */
        size_type size()const
        {
            return C[C.size()-1];
        }

        //! Returns the largest size that csa_sada2 can ever have.
        /*! Required for the Container Concept of the STL.
         *  \sa size
         */
        static size_type max_size()
        {
            return int_vector<>::max_size();
        }

        //! Returns if the data strucutre is empty.
        /*! Required for the Container Concept of the STL.A
         * \sa size
         */
        bool empty()const
        {
            return 0==size();
        }

        //! Swap method for csa_sada2
        /*! The swap method can be defined in terms of assignment.
            This requires three assignments, each of which, for a container type, is linear
            in the container's size. In a sense, then, a.swap(b) is redundant.
            This implementation guaranties a run-time complexity that is constant rather than linear.
            \param csa csa_sada2 to swap.

            Required for the Assignable Conecpt of the STL.
          */
        void swap(csa_sada2& csa)
        {
            if (this != &csa) {
                util::swap_support(m_psi_support, csa.m_psi_support, this, &csa);
                m_sa_sample.swap(csa.m_sa_sample);
                util::swap_support(m_isa_sample, csa.m_isa_sample, &m_sa_sample, &(csa.m_sa_sample));
                m_alphabet.swap(csa.m_alphabet);
            }
        }


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
        value_type operator[](size_type i)const
        {
            size_type off = 0;
            while (!m_sa_sample.is_sampled(i)) {  // while i mod t_dens != 0 (SA[i] is not sampled)
                i = psi[i];                       // go to the position where SA[i]+1 is located
                ++off;                            // add 1 to the offset
            }
            value_type result = m_sa_sample[i];
            if (result < off) {
                return psi.size()-(off-result);
            } else
                return result-off;
        }


        //! Assignment Copy Operator.
        /*!
         *    Required for the Assignable Concept of the STL.
         */
        csa_sada2& operator=(const csa_sada2& csa)
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
        csa_sada2& operator=(csa_sada2&& csa)
        {
            if (this != &csa) {
                m_psi_support = std::move(csa.m_psi_support);
                m_psi_support.set_vector(this);
                m_sa_sample   = std::move(csa.m_sa_sample);
                m_isa_sample  = std::move(csa.m_isa_sample);
                m_alphabet    = std::move(csa.m_alphabet);
            }
            return *this;
        }

        //! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_psi_support.serialize(out, child, "psi");
            written_bytes += m_sa_sample.serialize(out, child, "sa_samples");
            written_bytes += m_isa_sample.serialize(out, child, "isa_samples");
            written_bytes += m_alphabet.serialize(out, child, "alphabet");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load from a stream.
        /*! \param in Input stream to load the data structure from.
         */
        void load(std::istream& in)
        {
            m_psi_support.load(in);
            m_psi_support.set_vector(this);
            m_sa_sample.load(in);
            m_isa_sample.load(in, &m_sa_sample);
            m_alphabet.load(in);
        }

    private:

        // Calculates how many symbols c are in the prefix [0..i-1] of the BWT of the original text.
        /*
         *  \tpara Type of index. Should either be an unsigned integer or and std::array<,2> of unsigned integers
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in [0..size()]\f$.
         *  \param c The symbol to count in the prefix.
         *    \returns The number of occurrences of symbol c in the prefix [0..i-1] of the BWT.
         *  \par Time complexity
         *        \f$ \Order{\log n t_{\Psi}} \f$
         */
        // replace const char_type c by const std::array<char_type, alphabet_type::C_depth>& c
        template<typename t_pos , typename t_char>
        t_pos rank_bwt(t_pos i, const t_char c)const
        {
            auto cc = char2comp[c];
            if (cc==0 and c!=0) // character is not in the text => return 0
                return t_pos {0};
            if (i == t_pos {0})
                return t_pos {0};
            return m_psi_support.rank(i, cc);
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
            return m_psi_support.select(i, cc);
        }
};

// == template functions ==

template<class t_enc_vec, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
csa_sada2<t_enc_vec, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::csa_sada2(cache_config& config)
{
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
        m_psi_support = psi_type(psi_buf, this);
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

} // end namespace sdsl
#endif
