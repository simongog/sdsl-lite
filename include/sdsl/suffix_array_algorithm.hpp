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

//! Forward search for a pattern in an \f$\omega\f$-interval \f$[\ell..r]\f$ in the CSA.
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
 */
template<class t_csa, class t_pat_iter>
typename t_csa::size_type
forward_search(
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
    assert(l <= r); assert(r < csa.size());

    auto size = csa.size();

    l_res = l;
    r_res = l - 1;
    auto l_res_upper = r + 1;
    auto r_res_upper = r + 1;

    // shortcut for too long patterns
    if ((typename t_csa::size_type)(end - begin) >= size)
        return 0;

    // compares the pattern with CSA-prefix i (truncated to length $|pattern|$).
    auto compare = [&](typename t_csa::size_type i) -> int {
        for (auto current = begin; current != end; current++)
        {
            auto index = csa.char2comp[*current];
            if (index == 0) return -1;
            if (csa.C[index + 1] - 1 < i) return -1;
            if (csa.C[index] > i) return 1;
            i = csa.psi[i];
        }
        return 0;
    };

    // binary search (on min)
    while (l_res < l_res_upper) {
        typename t_csa::size_type sample = l_res + (l_res_upper - l_res) / 2;
        int result = compare(sample);
        if (result == 1)
            l_res = sample + 1;
        else if (result == -1)
            l_res_upper = sample;
        else
            l_res_upper = sample;
    }

    // binary search (on max)
    while (r_res + 1 < r_res_upper) {
        typename t_csa::size_type sample = r_res + (r_res_upper - r_res) / 2;
        int result = compare(sample);
        if (result == 1)
            r_res = sample;
        else if (result == -1)
            r_res_upper = sample;
        else
            r_res = sample;
    }

    return r_res - l_res + 1;
}

//! Forward search for a character c in an \f$\omega\f$-interval \f$[\ell..r]\f$ in the CSA.
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
 */
template<class t_csa>
typename t_csa::size_type forward_search(
    const t_csa& csa,
    typename t_csa::size_type l,
    typename t_csa::size_type r,
    typename t_csa::char_type c,
    typename t_csa::size_type& l_res,
    typename t_csa::size_type& r_res,
    SDSL_UNUSED typename std::enable_if<std::is_same<csa_tag, typename t_csa::index_category>::value, csa_tag>::type x = csa_tag()
)
{
    auto c_ptr = &c;
    return forward_search(csa, l, r, c_ptr, c_ptr + 1, l_res, r_res);
}

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
    typename t_csa::size_type cc = csa.char2comp[c];
    if (cc == 0 and c > 0) {
        l_res = 1;
        r_res = 0;
    } else {
        typename t_csa::size_type c_begin = csa.C[cc];
        if (l == 0 and r+1 == csa.size()) {
            l_res = c_begin;
            r_res = csa.C[cc+1] - 1;
        } else {
            l_res = c_begin + csa.bwt.rank(l, c); // count c in bwt[0..l-1]
            r_res = c_begin + csa.bwt.rank(r+1, c) - 1; // count c in bwt[0..r]
        }
    }
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
    SDSL_UNUSED typename std::enable_if< t_wt::lex_ordered, csa_tag>::type x = csa_tag()
)
{
    assert(l_fwd <= r_fwd); assert(r_fwd < csa_fwd.size());
    typedef typename csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, t_alphabet_strat>::size_type size_type;
    size_type c_begin = csa_fwd.C[csa_fwd.char2comp[c]];
    auto r_s_b =  csa_fwd.wavelet_tree.lex_count(l_fwd, r_fwd+1, c);
    size_type rank_l = std::get<0>(r_s_b);
    size_type s = std::get<1>(r_s_b), b = std::get<2>(r_s_b);
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
    SDSL_UNUSED typename std::enable_if< t_wt::lex_ordered, csa_tag>::type x = csa_tag()
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
 * the \f$p^{rev}omega^{rev}\f$-interval in the CSA object of the backward text.
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
    SDSL_UNUSED typename std::enable_if< t_wt::lex_ordered, csa_tag>::type x = csa_tag()
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

//! Counts the number of occurrences of a pattern in a CSA.
/*!
 * \tparam t_csa      CSA type.
 *
 * \param csa The CSA object.
 * \param pat The pattern.
 * \return The number of occurrences of the pattern in the CSA.
 *
 * \par Time complexity
 *        \f$ \Order{ t_{backward\_search} } \f$
 */

template<class t_csx>
typename t_csx::size_type count(
    const t_csx& csx,
    const typename t_csx::string_type& pat
)
{
    typename t_csx::index_category tag;
    return count(csx, pat.begin(), pat.end(), tag);
}

//! Returns the lexicographic interval of a pattern in the SA
/*!
 * \tparam t_csx      Type of CSA/CST.
 * \tparam t_pat_iter Type of pattern iterator.
 *
 * \param csx         CSA or CST object.
 * \param begin       Iterator to the begin of the pattern.
 * \param end         Iterator to the end (exlusive) of the pattern.
 *
 * \return The interval in the SA in which all suffixes are prefixed
 *         by the given pattern.
 */
template<typename t_csx, typename t_pat_iter>
auto lex_interval(
    const t_csx& csx,
    t_pat_iter begin,
    t_pat_iter end
) -> std::array<typename t_csx::size_type, 2> {
    std::array<typename t_csx::size_type, 2> res;
    backward_search(csx, 0, csx.size()-1, begin, end, res[0], res[1]);
    return res;
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
 * \return A vector containing the occurrences of the pattern  in the CSA.
 *
 * \par Time complexity
 *        \f$ \Order{ t_{backward\_search} + z \cdot t_{SA} } \f$, where \f$z\f$ is the number of
 *         occurrences of pattern in the CSA.
 */
template<class t_csa, class t_pat_iter, class t_rac=int_vector<64>>
t_rac locate(
    const t_csa&  csa,
    t_pat_iter begin,
    t_pat_iter end,
    SDSL_UNUSED typename std::enable_if<std::is_same<csa_tag, typename t_csa::index_category>::value, csa_tag>::type x = csa_tag()
)
{
    typename t_csa::size_type occ_begin, occ_end, occs;
    occs = backward_search(csa, 0, csa.size()-1, begin, end, occ_begin, occ_end);
    t_rac occ(occs);
    for (typename t_csa::size_type i=0; i < occs; ++i) {
        occ[i] = csa[occ_begin+i];
    }
    return occ;
}

//! Calculates all occurrences of a pattern pat in a CSA/CST.
/*!
 * \tparam t_csa      CSA/CST type.
 * \tparam t_rac      Resizeable random access container.
 *
 * \param csa  The CSA/CST object.
 * \param pat  The pattern.
 * \return A vector containing the occurrences of the pattern  in the CSA.
 *
 * \par Time complexity
 *        \f$ \Order{ t_{backward\_search} + z \cdot t_{SA} } \f$, where \f$z\f$ is the number of
 *         occurrences of pattern in the CSA.
 */
template<class t_csx, class t_rac=int_vector<64>>
t_rac locate(
    const t_csx&  csx,
    const typename t_csx::string_type& pat
)
{
    typename t_csx::index_category tag;
    return locate<t_csx, decltype(pat.begin()), t_rac>(csx, pat.begin(), pat.end(), tag);
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

//! Specialization of extract for LF-function based CSAs
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
    auto steps = end-begin+1;
    if (steps > 0) {
        auto order = csa.isa[end];
        text[--steps] = first_row_symbol(order, csa);
        while (steps != 0) {
            auto rc = csa.wavelet_tree.inverse_select(order);
            auto j = rc.first;
            auto c = rc.second;
            order = csa.C[ csa.char2comp[c] ] + j;
            text[--steps] = c;
        }
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
    for (typename t_csa::size_type i=0, order = csa.isa[begin]; steps != 0; --steps, ++i) {
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
template<class t_csa>
typename t_csa::string_type extract(
    const t_csa& csa,
    typename t_csa::size_type begin,
    typename t_csa::size_type end,
    SDSL_UNUSED typename std::enable_if<std::is_same<csa_tag, typename t_csa::index_category>::value, csa_tag>::type x = csa_tag()
)
{
    assert(end <= csa.size());
    assert(begin <= end);
    typedef typename t_csa::string_type string_type;
    string_type result(end-begin+1, (typename string_type::value_type)0);
    extract(csa, begin, end, result.begin());
    return result;
}


} // end namespace
#endif
