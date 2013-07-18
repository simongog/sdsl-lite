/* sdsl - succinct data structures library
    Copyright (C) 2010-2013 Simon Gog

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
/*! \file suffix_array_algorithm.hpp
    \brief suffix_array_algorithm.hpp contains algorithms on CSAs
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_SUFFIX_ARRAY_ALGORITHM
#define INCLUDED_SDSL_SUFFIX_ARRAY_ALGORITHM

#include <iterator>
#include "suffix_array_helper.hpp"

namespace sdsl
{

//! Backward search for a character c in an \f$\omega\f$-interval \f$[\ell..r]\f$ in the CSA.
/*!
 * \tparam t_csa CSA type.
 *
 * \param csa    The CSA object.
 * \param l      Left border of the interval \f$ [\ell..r]\f$.
 * \param r      Right border of the interval \f$ [\ell..r]\f$.
 * \param c      Character to be prepended to \f$\omega\f$.
 * \param l_res  New left border.
 * \param r_res  Right border.
 * \return The size of the new interval [\ell_{new}..r_{new}].
 *         Equals zero, if no match is found.
 *
 * \pre \f$ 0 \leq \ell \leq r < csa.size() \f$
 *
 * \par Time complexity
 *         \f$ \Order{ t_{rank\_bwt} } \f$
 * \par Reference
 *         Paolo Ferragina, Giovanni Manzini:
 *         Opportunistic Data Structures with Applications.
 *         FOCS 2000: 390-398
 */
template<class t_csa>
typename t_csa::size_type backward_search(
    const t_csa& csa,
    typename t_csa::size_type l,
    typename t_csa::size_type r,
    typename t_csa::char_type c,
    typename t_csa::size_type& l_res,
    typename t_csa::size_type& r_res,
    SDSL_UNUSED typename std::enable_if<std::is_same<csa_tag, typename t_csa::index_category>::value, csa_tag>::type x = csa_tag()
)
{
    assert(l <= r); assert(r < csa.size());
    typename t_csa::size_type c_begin = csa.C[csa.char2comp[c]];
    l_res = c_begin + csa.rank_bwt(l, c); // count c in bwt[0..l-1]
    r_res = c_begin + csa.rank_bwt(r+1, c) - 1; // count c in bwt[0..r]
    assert(r_res+1-l_res >= 0);
    return r_res+1-l_res;
}


//! Backward search for a pattern in an \f$\omega\f$-interval \f$[\ell..r]\f$ in the CSA.
/*!
 * \tparam t_csa      A CSA type.
 * \tparam t_pat_iter Pattern iterator type.
 *
 * \param csa   The CSA object.
 * \param l     Left border of the lcp-interval \f$ [\ell..r]\f$.
 * \param r     Right border of the lcp-interval \f$ [\ell..r]\f$.
 * \param begin Iterator to the begin of the pattern (inclusive).
 * \param end   Iterator to the end of the pattern (exclusive).
 * \param l_res New left border.
 * \param r_res New right border.
 * \return The size of the new interval [\ell_{new}..r_{new}].
 *         Equals zero, if no match is found.
 *
 * \pre \f$ 0 \leq \ell \leq r < csa.size() \f$
 *
 * \par Time complexity
 *       \f$ \Order{ len \cdot t_{rank\_bwt} } \f$
 * \par Reference
 *         Paolo Ferragina, Giovanni Manzini:
 *         Opportunistic Data Structures with Applications.
 *         FOCS 2000: 390-398
 */
template<class t_csa, class t_pat_iter>
typename t_csa::size_type
backward_search(
    const t_csa& csa,
    typename t_csa::size_type l,
    typename t_csa::size_type r,
    t_pat_iter begin,
    t_pat_iter end,
    typename t_csa::size_type& l_res,
    typename t_csa::size_type& r_res,
    SDSL_UNUSED typename std::enable_if<std::is_same<csa_tag, typename t_csa::index_category>::value, csa_tag>::type x = csa_tag()
)
{
    t_pat_iter it = end;
    while (begin < it and r+1-l > 0) {
        --it;
        backward_search(csa, l, r, (typename t_csa::char_type)*it, l, r);
    }
    l_res = l;
    r_res = r;
    return r+1-l;
}

// Maybe one should replace these traits with a member variable in each wavelet tree
template<class t_T>
struct wt_has_lex_count_trait {
    enum {value = false};
};

template<class t_rac, class t_bv, class t_rs, class t_ss1, class t_ss0>
struct wt_has_lex_count_trait<wt<t_rac, t_bv, t_rs, t_ss1, t_ss0> > {
    enum {value = true};
};

template<class t_bv, class t_rs, class t_ss1, class t_ss0, class t_tree_strat>
struct wt_has_lex_count_trait<wt_hutu<t_bv, t_rs, t_ss1, t_ss0, t_tree_strat> > {
    enum {value = true};
};

//! Bidirectional search for a character c on an interval \f$[l_fwd..r_fwd]\f$ of the suffix array.
/*!
 * \param csa_fwd   The CSA object of the forward text in which the backward_search should be done.
 * \param l_fwd     Left border of the lcp-interval \f$ [l_fwd..r_fwd]\f$ in suffix array of the forward text.
 * \param r_fwd     Right border of the lcp-interval \f$ [l_fwd..r_fwd]\f$ in suffix array of the forward text.
 * \param l_bwd     Left border of the lcp-interval \f$ [l_bwd..r_bwd]\f$ in suffix array of the backward text.
 * \param r_bwd     Right border of the lcp-interval \f$ [l_bwd..r_bwd]\f$ in suffix array of the backward text.
 * \param c         The character c which is the starting character of the suffixes in the resulting interval \f$ [l_fwd_res..r_fwd_res] \f$ .
 * \param l_fwd_res Reference to the resulting left border in suffix array of the forward text.
 * \param r_fwd_res Reference to the resulting right border in suffix array of the forward text.
 * \param l_bwd_res Reference to the resulting left border in suffix array of the backward text.
 * \param r_bwd_res Reference to the resulting right border in suffix array of the backward text.
 * \return The size of the new interval [l_fwd_res..r_fwd_res].
 * \pre \f$ 0 \leq \ell \leq r_fwd < csa_fwd.size() \f$
 * \par Reference
 *         Thomas Schnattinger, Enno Ohlebusch, Simon Gog:
 *         Bidirectional search in a string with wavelet trees and bidirectional matching statistics.
 *         Inf. Comput. 213: 13-22
 */
template<class t_wt, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
typename csa_wt<t_wt>::size_type bidirectional_search(
    const csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>& csa_fwd,
    typename csa_wt<>::size_type l_fwd,
    typename csa_wt<>::size_type r_fwd,
    typename csa_wt<>::size_type l_bwd,
    typename csa_wt<>::size_type r_bwd,
    typename csa_wt<>::char_type c,
    typename csa_wt<>::size_type& l_fwd_res,
    typename csa_wt<>::size_type& r_fwd_res,
    typename csa_wt<>::size_type& l_bwd_res,
    typename csa_wt<>::size_type& r_bwd_res,
    SDSL_UNUSED typename std::enable_if< wt_has_lex_count_trait<t_wt>::value, csa_tag>::type x = csa_tag()
)
{
    assert(l_fwd <= r_fwd); assert(r_fwd < csa_fwd.size());
    typedef typename csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::size_type size_type;
    size_type c_begin = csa_fwd.C[csa_fwd.char2comp[c]];
    size_type s, b;
    size_type rank_l = csa_fwd.wavelet_tree.lex_count(l_fwd, r_fwd+1, c, s, b);
    size_type rank_r = r_fwd - l_fwd - s - b + rank_l;
    l_fwd_res = c_begin + rank_l;
    r_fwd_res = c_begin + rank_r;
    assert(r_fwd_res+1 >= l_fwd_res);
    l_bwd_res = l_bwd + s;
    r_bwd_res = r_bwd - b;
    assert(r_bwd_res-l_bwd_res == r_fwd_res-l_fwd_res);
    return r_fwd_res+1-l_fwd_res;
}

//! Bidirectional search in backward direction.
/*!
 * The function requires a pattern \f$p\f$, an \f$\omega\f$-interval \f$[l_fwd..r_fwd]\f$ in the CSA object
 * of the forward text and an \f$\omega^{rev}\f$-interval \f$[l_bwd..r_bwd]\f$ in the CSA object of the backward text.
 * The function returns the \f$p\omega\f$-interval in the CSA object of the forward text and
 * the \f$\omega^{rev}p^{rev}\f$-interval in the CSA object of the backward text.
 *
 * \tparam t_pat_iter Pattern iterator type.
 *
 * \param csa_fwd   The CSA object of the forward text.
 * \param csa_bwd   The CSA object of the backward text.
 * \param l_fwd     Left border of the lcp-interval \f$ [l_fwd..r_fwd]\f$ in suffix array of the forward text.
 * \param r_fwd     Right border of the lcp-interval \f$ [l_fwd..r_fwd]\f$ in suffix array of the forward text.
 * \param l_bwd     Left border of the lcp-interval \f$ [l_bwd..r_bwd]\f$ in suffix array of the backward text.
 * \param r_bwd     Right border of the lcp-interval \f$ [l_bwd..r_bwd]\f$ in suffix array of the backward text.
 * \param begin     Iterator to the begin of the pattern (inclusive).
 * \param end       Iterator to the end of the pattern (exclusive).
 * \param l_fwd_res Reference to the resulting left border in suffix array of the forward text.
 * \param r_fwd_res Reference to the resulting right border in suffix array of the forward text.
 * \param l_bwd_res Reference to the resulting left border in suffix array of the backward text.
 * \param r_bwd_res Reference to the resulting right border in suffix array of the backward text.
 * \return The size of the new interval [l_fwd_res..r_fwd_res].
 *         Equals zero, if no match is found.
 *
 * \pre \f$ 0 \leq \ell \leq r_fwd < csa_fwd.size() \f$
 * \par Reference
 *         Thomas Schnattinger, Enno Ohlebusch, Simon Gog:
 *         Bidirectional search in a string with wavelet trees and bidirectional matching statistics.
 *         Inf. Comput. 213: 13-22
 */
template<class t_pat_iter, class t_wt, uint32_t t_dens, uint32_t t_inv_dens, class t_sa_sample_strat, class t_isa, class t_alphabet_strat>
typename csa_wt<>::size_type bidirectional_search_backward(
    const csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>& csa_fwd,
    SDSL_UNUSED const csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>& csa_bwd,
    typename csa_wt<>::size_type l_fwd,
    typename csa_wt<>::size_type r_fwd,
    typename csa_wt<>::size_type l_bwd,
    typename csa_wt<>::size_type r_bwd,
    t_pat_iter begin,
    t_pat_iter end,
    typename csa_wt<>::size_type& l_fwd_res,
    typename csa_wt<>::size_type& r_fwd_res,
    typename csa_wt<>::size_type& l_bwd_res,
    typename csa_wt<>::size_type& r_bwd_res,
    SDSL_UNUSED typename std::enable_if< wt_has_lex_count_trait<t_wt>::value, csa_tag>::type x = csa_tag()
)
{
    t_pat_iter it = end;
    while (begin < it and r_fwd+1-l_fwd > 0) {
        --it;
        bidirectional_search(csa_fwd, l_fwd, r_fwd, l_bwd, r_bwd, (typename csa_wt<>::char_type)*it, l_fwd, r_fwd, l_bwd, r_bwd);
    }
    l_fwd_res = l_fwd;
    r_fwd_res = r_fwd;
    l_bwd_res = l_bwd;
    r_bwd_res = r_bwd;
    return r_fwd+1-l_fwd;
}

//! Bidirectional search in forward direction.
/*!
 * The function requires a pattern \f$p\f$, an \f$\omega\f$-interval \f$[l_fwd..r_fwd]\f$ in the CSA object
 * of the forward text and an \f$\omega^{rev}\f$-interval \f$[l_bwd..r_bwd]\f$ in the CSA object of the backward text.
 * The function returns the \f$\omega p\f$-interval in the CSA object of the forward text and
 * the \f$\p^{rev}omega^{rev}\f$-interval in the CSA object of the backward text.
 *
 * \tparam t_pat_iter Pattern iterator type.
 *
 * \param csa_fwd   The CSA object of the forward text.
 * \param csa_bwd   The CSA object of the backward text.
 * \param l_fwd     Left border of the lcp-interval \f$ [l_fwd..r_fwd]\f$ in suffix array of the forward text.
 * \param r_fwd     Right border of the lcp-interval \f$ [l_fwd..r_fwd]\f$ in suffix array of the forward text.
 * \param l_bwd     Left border of the lcp-interval \f$ [l_bwd..r_bwd]\f$ in suffix array of the backward text.
 * \param r_bwd     Right border of the lcp-interval \f$ [l_bwd..r_bwd]\f$ in suffix array of the backward text.
 * \param begin     Iterator to the begin of the pattern (inclusive).
 * \param end       Iterator to the end of the pattern (exclusive).
 * \param l_fwd_res Reference to the resulting left border in suffix array of the forward text.
 * \param r_fwd_res Reference to the resulting right border in suffix array of the forward text.
 * \param l_bwd_res Reference to the resulting left border in suffix array of the backward text.
 * \param r_bwd_res Reference to the resulting right border in suffix array of the backward text.
 * \return The size of the new interval [l_fwd_res..r_fwd_res].
 *         Equals zero, if no match is found.
 *
 * \pre \f$ 0 \leq \ell \leq r_fwd < csa_fwd.size() \f$
 * \par Reference
 *         Thomas Schnattinger, Enno Ohlebusch, Simon Gog:
 *         Bidirectional search in a string with wavelet trees and bidirectional matching statistics.
 *         Inf. Comput. 213: 13-22
 */
template<class t_pat_iter,
         class t_wt,
         uint32_t t_dens,
         uint32_t t_inv_dens,
         class t_sa_sample_strat,
         class t_isa,
         class t_alphabet_strat>
typename csa_wt<t_wt>::size_type
bidirectional_search_forward(
    SDSL_UNUSED const csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>& csa_fwd,
    const csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>& csa_bwd,
    typename csa_wt<>::size_type l_fwd,
    typename csa_wt<>::size_type r_fwd,
    typename csa_wt<>::size_type l_bwd,
    typename csa_wt<>::size_type r_bwd,
    t_pat_iter begin,
    t_pat_iter end,
    typename csa_wt<>::size_type& l_fwd_res,
    typename csa_wt<>::size_type& r_fwd_res,
    typename csa_wt<>::size_type& l_bwd_res,
    typename csa_wt<>::size_type& r_bwd_res,
    SDSL_UNUSED typename std::enable_if< wt_has_lex_count_trait<t_wt>::value, csa_tag>::type x = csa_tag()
)
{
    t_pat_iter it = begin;
    while (it < end and r_fwd+1-l_fwd > 0) {
        bidirectional_search(csa_bwd, l_bwd, r_bwd, l_fwd, r_fwd, (typename csa_wt<>::char_type)*it, l_bwd, r_bwd, l_fwd, r_fwd);
        ++it;
    }
    l_fwd_res = l_fwd;
    r_fwd_res = r_fwd;
    l_bwd_res = l_bwd;
    r_bwd_res = r_bwd;
    return r_fwd+1-l_fwd;
}

//! Counts the number of occurrences of a pattern in a CSA.
/*!
 * \tparam t_csa      CSA type.
 * \tparam t_pat_iter Pattern iterator type.
 *
 * \param csa   The CSA object.
 * \param begin Iterator to the begin of the pattern (inclusive).
 * \param end   Iterator to the end of the pattern (exclusive).
 * \return The number of occurrences of the pattern in the CSA.
 *
 * \par Time complexity
 *        \f$ \Order{ t_{backward\_search} } \f$
 */
template<class t_csa, class t_pat_iter>
typename t_csa::size_type count(
    const t_csa& csa,
    t_pat_iter begin,
    t_pat_iter end,
    csa_tag
)
{
    if (end - begin > (typename std::iterator_traits<t_pat_iter>::difference_type)csa.size())
        return 0;
    typename t_csa::size_type t=0; // dummy variable for the backward_search call
    typename t_csa::size_type result = backward_search(csa, 0, csa.size()-1, begin, end, t, t);
    return result;
}


template<class t_csx, class t_pat_iter>
typename t_csx::size_type count(
    const t_csx& csx,
    t_pat_iter begin,
    t_pat_iter end
)
{
    typename t_csx::index_category tag;
    return count(csx, begin, end, tag);
}

//! Calculates all occurrences of a pattern pat in a CSA.
/*!
 * \tparam t_csa      CSA type.
 * \tparam t_pat_iter Pattern iterator type.
 * \tparam t_rac      Resizeable random access container.
 *
 * \param csa   The CSA object.
 * \param begin Iterator to the begin of the pattern (inclusive).
 * \param end   Iterator to the end of the pattern (exclusive).
 * \param occ   Container object in which the occurrences are stored.
 * \return The number of occurrences of the pattern  in the CSA.
 *
 * \par Time complexity
 *        \f$ \Order{ t_{backward\_search} + z \cdot t_{SA} } \f$, where \f$z\f$ is the number of
 *         occurrences of pattern in the CSA.
 */
template<class t_csa, class t_pat_iter, class t_rac>
typename t_csa::size_type locate(
    const t_csa&  csa,
    t_pat_iter begin,
    t_pat_iter end,
    t_rac& occ,
    SDSL_UNUSED typename std::enable_if<std::is_same<csa_tag, typename t_csa::index_category>::value, csa_tag>::type x = csa_tag()
)
{
    typename t_csa::size_type occ_begin, occ_end, occs;
    occs = backward_search(csa, 0, csa.size()-1, begin, end, occ_begin, occ_end);
    occ.resize(occs);
    for (typename t_csa::size_type i=0; i < occs; ++i) {
        occ[i] = csa[occ_begin+i];
    }
    return occs;
}


//! Writes the substring T[begin..end] of the original text T to text[0..end-begin+1].
/*!
 * \tparam t_csa       CSA type.
 * \tparam t_text_iter Random access iterator type.
 *
 * \param csa   The CSA object.
 * \param begin Position of the first character which should be extracted (inclusive).
 * \param end   Position of the last character which should be extracted (inclusive).
 * \param text  Random access iterator pointing to the start of an container, which can hold at least (end-begin+1) character.
 * \returns The length of the extracted text.
 * \pre \f$begin <= end\f$ and \f$ end < csa.size() \f$
 * \par Time complexity
 *        \f$ \Order{ (end-begin+1) \cdot t_{\Psi} + t_{SA^{-1}} } \f$
 */
template<class t_csa, class t_text_iter>
typename t_csa::size_type extract(
    const t_csa& csa,
    typename t_csa::size_type begin,
    typename t_csa::size_type end,
    t_text_iter text,
    SDSL_UNUSED typename std::enable_if<std::is_same<csa_tag, typename t_csa::index_category>::value, csa_tag>::type x = csa_tag()
)
{
    typename t_csa::extract_category extract_tag;
    return extract(csa, begin, end, text, extract_tag);
}

//! Specialization of extract for \f$\Psi\f$-function based CSAs
template<class t_csa, class t_text_iter>
typename t_csa::size_type extract(
    const t_csa& csa,
    typename t_csa::size_type begin,
    typename t_csa::size_type end,
    t_text_iter text,
    lf_tag
)
{
    assert(end < csa.size());
    assert(begin <= end);
    typename t_csa::size_type steps = end-begin+1;
    for (typename t_csa::size_type order = csa(end);  steps != 0; --steps) {
        text[steps-1] = first_row_symbol(order, csa);
        if (steps != 0) order = csa.psi(order);
    }
    return end-begin+1;
}

//! Specialization of extract for \f$\Psi\f$-function based CSAs
template<class t_csa, class t_text_iter>
typename t_csa::size_type extract(
    const t_csa& csa,
    typename t_csa::size_type begin,
    typename t_csa::size_type end,
    t_text_iter text,
    psi_tag
)
{
    assert(end < csa.size());
    assert(begin <= end);
    typename t_csa::size_type steps = end-begin+1;
    for (typename t_csa::size_type i=0, order = csa(begin); steps != 0; --steps, ++i) {
        text[i] = first_row_symbol(order, csa);
        if (steps != 0) order = csa.psi[order];
    }
    return end-begin+1;
}

//! Reconstructs the substring T[begin..end] of the original text T to text[0..end-begin+1].
/*!
 * \tparam t_rac Random access container which should hold the result.
 * \tparam t_csa CSA type.
 *
 * \param csa   The CSA object.
 * \param begin Position of the first character which should be extracted (inclusive).
 * \param end   Position of the last character which should be extracted (inclusive).
 * \return A t_rac object holding the extracted text.
 * \pre \f$begin <= end\f$ and \f$ end < csa.size() \f$
 * \par Time complexity
 *        \f$ \Order{ (end-begin+1) \cdot t_{\Psi} + t_{SA^{-1}} } \f$
 */
template<class t_rac, class t_csa>
t_rac extract(
    const t_csa& csa,
    typename t_csa::size_type begin,
    typename t_csa::size_type end,
    SDSL_UNUSED typename std::enable_if<std::is_same<csa_tag, typename t_csa::index_category>::value, csa_tag>::type x = csa_tag()
)
{
    assert(end <= csa.size());
    assert(begin <= end);
    t_rac result(end-begin+1, (typename t_csa::char_type)0);
    extract(csa, begin, end, result.begin());
    return result;
}


} // end namespace
#endif
