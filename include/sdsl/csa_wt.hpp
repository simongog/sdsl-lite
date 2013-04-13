/* sdsl - succinct data structures library
    Copyright (C) 2009-2013 Simon Gog

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
/*! \file csa_wt.hpp
    \brief csa_wt.hpp contains an implementation of the compressed suffix array based on a wavelet tree.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_CSA_WT
#define INCLUDED_SDSL_CSA_WT

#include "wavelet_trees.hpp"
#include "suffix_array_helper.hpp"
#include "algorithms.hpp"
#include "iterators.hpp"
#include "util.hpp"
#include "fast_cache.hpp"
#include "csa_sampling_strategy.hpp"
#include "csa_alphabet_strategy.hpp"
#include <iostream>
#include <algorithm> // for std::swap
#include <cassert>
#include <cstring> // for strlen
#include <iomanip>
#include <iterator>

namespace sdsl
{

template<class t_csa>
class psi_of_csa_wt;  // forward declaration of PSI-array class

template<class t_csa>
class bwt_of_csa_wt;   // forward declaration of BWT-array class


//! A class for the Compressed Suffix Array (CSA) based on a Wavelet Tree (WT) of the Burrow Wheeler Transform of the original text.
/*!
  *  \tparam t_wt              Wavelet tree
  *  \tparam t_dens            Sampling density of SA values
  *  \tparam t_int_dens        Sampling density of ISA values
  *  \tparam t_sa_sample_strat Policy of SA sampling. E.g. sample in SA-order or text-order.
  *  \tparam t_isa             Vector type for ISA sample values.
  *  \tparam t_alphabet_strat  Policy for alphabet representation.
  *
  *  \sa sdsl::csa_sada, sdsl::csa_bitcompressed
  * @ingroup csa
 */
template<class t_wt              = wt_huff<>,              // Wavelet tree type
         uint32_t t_dens         = 32,                     // Sample density for suffix array (SA) values
         uint32_t t_inv_dens     = 64,                     // Sample density for inverse suffix array (ISA) values
         class t_sa_sample_strat = sa_order_sa_sampling<>, // Policy class for the SA sampling. Alternative text_order_sa_sampling.
         class t_isa             = int_vector<>,           // Container for the ISA samples.
         class t_alphabet_strat  = byte_alphabet           // Policy class for the representation of the alphabet.
         >
class csa_wt
{
    public:
        enum { sa_sample_dens = t_dens,
               isa_sample_dens = t_inv_dens
             };

        typedef uint64_t                                                       value_type;
        typedef random_access_const_iterator<csa_wt>                           const_iterator;
        typedef const_iterator                                                 iterator;
        typedef const value_type                                               const_reference;
        typedef const_reference                                                reference;
        typedef const_reference*                                               pointer;
        typedef const pointer                                                  const_pointer;
        typedef int_vector<>::size_type                                        size_type;
        typedef size_type                                                      csa_size_type;
        typedef ptrdiff_t                                                      difference_type;
        typedef psi_of_csa_wt<csa_wt>                                          psi_type;
        typedef bwt_of_csa_wt<csa_wt>                                          bwt_type;
        typedef text_of_csa<csa_wt>                                            text_type;
        typedef t_wt                                                           wavelet_tree_type;
        typedef typename t_sa_sample_strat::template type<csa_wt>::sample_type sa_sample_type;
        typedef t_isa                                                          isa_sample_type;
        typedef t_alphabet_strat                                               alphabet_type;
        typedef typename alphabet_type::char_type                              char_type; // Note: This is the char type of the CSA not the WT!
        typedef typename alphabet_type::comp_char_type                         comp_char_type;
        typedef const char_type*                                               pattern_type;
        typedef csa_wt                                                         csa_type;

        typedef csa_tag                                                        index_category;
        typedef lf_tag                                                         extract_category;
        typedef typename alphabet_type::alphabet_category                      alphabet_category;

    private:
        t_wt            m_wavelet_tree; // the wavelet tree
        sa_sample_type  m_sa_sample; // suffix array samples
        isa_sample_type m_isa_sample; // inverse suffix array samples
        alphabet_type   m_alphabet;
//#define USE_CSA_CACHE
#ifdef USE_CSA_CACHE
        mutable fast_cache csa_cache;
#endif

        void copy(const csa_wt& csa) {
            m_wavelet_tree = csa.m_wavelet_tree;
            m_sa_sample    = csa.m_sa_sample;
            m_isa_sample   = csa.m_isa_sample;
            m_alphabet     = csa.m_alphabet;
        }

    public:
        const typename alphabet_type::char2comp_type& char2comp;
        const typename alphabet_type::comp2char_type& comp2char;
        const typename alphabet_type::C_type&         C;
        const typename alphabet_type::sigma_type&     sigma;
        const psi_type                                psi;
        const bwt_type                                bwt;
        const text_type                               text;
        const sa_sample_type&                         sa_sample;
        const isa_sample_type&                        isa_sample;
        const wavelet_tree_type&                      wavelet_tree;

        //! Default Constructor
        csa_wt(): char2comp(m_alphabet.char2comp), comp2char(m_alphabet.comp2char), C(m_alphabet.C), sigma(m_alphabet.sigma),
            psi(this), bwt(this), text(this), sa_sample(m_sa_sample), isa_sample(m_isa_sample), wavelet_tree(m_wavelet_tree) {}

        //! Copy constructor
        csa_wt(const csa_wt& csa): char2comp(m_alphabet.char2comp), comp2char(m_alphabet.comp2char), C(m_alphabet.C), sigma(m_alphabet.sigma),
            psi(this), bwt(this), text(this), sa_sample(m_sa_sample), isa_sample(m_isa_sample), wavelet_tree(m_wavelet_tree) {
            copy(csa);
        }

        //! Constructor taking a cache_config
        csa_wt(cache_config& config);

        //! Number of elements in the \f$\CSA\f$.
        /*! Required for the Container Concept of the STL.
         *  \sa max_size, empty
         *  \par Time complexity
         *      \f$ \Order{1} \f$
         */
        size_type size()const {
            return m_wavelet_tree.size();
        }

        //! Returns the largest size that csa_wt can ever have.
        /*! Required for the Container Concept of the STL.
         *  \sa size
         */
        static size_type max_size() {
            return bit_vector::max_size();
        }

        //! Returns if the data strucutre is empty.
        /*! Required for the Container Concept of the STL.
         * \sa size
         */
        bool empty()const {
            return m_wavelet_tree.empty();
        }

        //! Swap method for csa_wt
        /*! The swap method can be defined in terms of assignment.
            This requires three assignments, each of which, for a container type, is linear
            in the container's size. In a sense, then, a.swap(b) is redundant.
            This implementation guaranties a run-time complexity that is constant rather than linear.
            \param csa csa_wt to swap.

            Required for the Assignable Conecpt of the STL.
          */
        void swap(csa_wt& csa);

        //! Returns a const_iterator to the first element.
        /*! Required for the STL Container Concept.
         *  \sa end
         */
        const_iterator begin()const {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        /*! Required for the STL Container Concept.
         *  \sa begin.
         */
        const_iterator end()const {
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

        //! ()-operator return inverse suffix array values
        /*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
         *   \par Time complexity
         *      \f$ \Order{s_{SA^{-1}}\cdot t_{\Psi}} \f$, where every \f$s_{SA^{-1}}\f$th suffix array entry is sampled and \f$t_{\Psi}\f$
         *           is the access time for an element in the \f$\Psi\f$-function.
         */
        inline value_type operator()(size_type i)const;

        //! Assignment Operator.
        /*!
         *    Required for the Assignable Concept of the STL.
         */
        csa_wt& operator=(const csa_wt& csa);

        //! Serialize to a stream.
        /*! \param out Output stream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const;

        //! Load from a stream.
        /*! \param in Input stream to load the data structure from.
         */
        void load(std::istream& in);

        uint32_t get_psi_sample_dens() const {
            return psi.get_sample_dens();
        }


        //! Calculates how many symbols c are in the prefix [0..i-1] of the BWT of the original text.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in [0..size()]\f$.
         *  \param c The symbol to count the occurrences in the prefix.
         *    \returns The number of occurrences of symbol c in the prefix [0..i-1] of the BWT.
         *  \par Time complexity
         *        \f$ \Order{\log |\Sigma|} \f$
         */
        size_type rank_bwt(size_type i, const char_type c)const {
            return m_wavelet_tree.rank(i, c);
        }

        //! Calculates the position of the i-th c in the BWT of the original text.
        /*!
         *  \param i The i-th occurrence. \f$i\in [1..rank(size(),c)]\f$.
         *  \param c Symbol c.
         *    \returns The position of the i-th c in the BWT or size() if c does occur less then i times.
         *  \par Time complexity
         *        \f$ \Order{t_{\Psi}} \f$
         */
        size_type select_bwt(size_type i, const char_type c)const {
            assert(i > 0);
            char_type cc = char2comp[c];
            if (cc==0 and c!=0)  // character is not in the text => return size()
                return size();
            assert(cc != 255);
            if (C[cc]+i-1 <  C[cc+1]) {
                return m_wavelet_tree.select(i, c);
            } else
                return size();
        }
};

// == template functions ==

template<class t_wt, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::csa_wt(cache_config& config) :
    char2comp(m_alphabet.char2comp), comp2char(m_alphabet.comp2char), C(m_alphabet.C), sigma(m_alphabet.sigma),
    psi(this), bwt(this), text(this), sa_sample(m_sa_sample),isa_sample(m_isa_sample), wavelet_tree(m_wavelet_tree)
{
    if (!cache_file_exists(key_trait<alphabet_type::int_width>::KEY_BWT, config)) {
        return;
    }
    int_vector_file_buffer<alphabet_type::int_width> bwt_buf(cache_file_name(key_trait<alphabet_type::int_width>::KEY_BWT,config));
    int_vector_file_buffer<>  sa_buf(cache_file_name(constants::KEY_SA, config));
    size_type n = bwt_buf.int_vector_size;
    mm::log("csa-alphabet-construct-begin");
    {
        alphabet_type tmp_alphabet(bwt_buf, n);
        m_alphabet.swap(tmp_alphabet);
    }
    mm::log("csa-alphabet-construct-end");

    mm::log("wt-begin");
    {
        wavelet_tree_type tmp_wt(bwt_buf, n);
        m_wavelet_tree.swap(tmp_wt);
    }
    mm::log("wt-end");

    mm::log("sa-sample-begin");
    {
        sa_sample_type tmp_sa_sample(sa_buf);
        m_sa_sample.swap(tmp_sa_sample);
    }
    mm::log("sa-sample-end");

    mm::log("isa-sample-begin");
    algorithm::set_isa_samples<csa_wt>(sa_buf, m_isa_sample);
    mm::log("isa-sample-end");
}


template<class t_wt, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
inline typename csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::value_type csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::operator[](size_type i)const
{
    size_type off = 0;
    while (!m_sa_sample.is_sampled(i)) {
        i = psi(i);
        ++off;
    }
    value_type result = m_sa_sample.sa_value(i);
    if (result + off < size()) {
        return result + off;
    } else {
        return result + off - size();
    }
    /*
    #ifdef USE_CSA_CACHE
        size_type r = 0;
        if( csa_cache.exists(i, r) ){
            return r;
        }
    #endif
        size_type off = 0;
        while( i % t_dens ){// while i mod t_dens != 0 (SA[i] is not sampled)
            i = psi[i];       // go to the position where SA[i]+1 is located
            ++off;              // add 1 to the offset
        }
        value_type result = m_sa_sample[i/t_dens];
        // TODO: try LF-function for iteration
        if( result < off ){
    #ifdef USE_CSA_CACHE
            r = psi.size()-(off-result);
            return r;
    #endif
            return psi.size()-(off-result);
        }
        else{
    #ifdef USE_CSA_CACHE
            r = result-off;
            return r;
    #endif
            return result-off;
        }
    */
}

template<class t_wt, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
inline typename csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::value_type csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::operator()(size_type i)const
{
    size_type ii;
    value_type result = m_isa_sample[ ii = ((i+t_inv_dens-1)/t_inv_dens) ]; // get the leftmost sampled isa value to the right of i
    ii *= t_inv_dens;
    if (ii >= size()) {
        i = size() - 1 - i;
    } else {
        i = ii - i;
    }
    while (i--) {
        result = psi(result);
    }
    return result;
}

template<class t_wt, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
csa_wt<t_wt,t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>& csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::operator=(const csa_wt<t_wt,t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>& csa)
{
    if (this != &csa) {
        copy(csa);
    }
    return *this;
}


template<class t_wt, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
typename csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::size_type csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_wavelet_tree.serialize(out, child, "wavelet_tree");
    written_bytes += m_sa_sample.serialize(out, child, "sa_samples");
    written_bytes += m_isa_sample.serialize(out, child, "isa_samples");
    written_bytes += m_alphabet.serialize(out, child, "alphabet");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<class t_wt, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
void csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::load(std::istream& in)
{
    m_wavelet_tree.load(in);
    m_sa_sample.load(in);
    m_isa_sample.load(in);
    m_alphabet.load(in);
}

template<class t_wt, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
void csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::swap(csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>& csa)
{
    if (this != &csa) {
        m_wavelet_tree.swap(csa.m_wavelet_tree);
        m_sa_sample.swap(csa.m_sa_sample);
        m_isa_sample.swap(csa.m_isa_sample);
        m_alphabet.swap(csa.m_alphabet);
    }
}

} // end namespace sdsl
#endif
