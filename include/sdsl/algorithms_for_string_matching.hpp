/* sdsl - succinct data structures library
    Copyright (C) 2010 Simon Gog

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
/*! \file algorithms_for_string_matching.hpp
    \brief algorithms_for_string_matching.hpp contains algorithms for string matching like backward_search, ...
	\author Simon Gog
*/

#ifndef INCLUDED_SDSL_ALGORITHMS_FOR_STRING_MATCHING
#define INCLUDED_SDSL_ALGORITHMS_FOR_STRING_MATCHING

#include "int_vector.hpp"

#include <stdexcept> // for exceptions
#include <iostream>
#include <cassert>
#include <stack>
#include <utility>


namespace sdsl
{

/*!
	\author Simon Gog
 */
namespace algorithm
{

//! Backward search for a character c on an interval \f$[\ell..r]\f$ of the suffix array.
/*!
 * \param csa The csa in which the backward_search should be done.
 * \param l Left border of the lcp-interval \f$ [\ell..r]\f$.
 * \param r Right border of the lcp-interval \f$ [\ell..r]\f$.
 * \param c The character c which is the starting character of the suffixes in the resulting interval \f$ [\ell_{new}..r_{new}] \f$ .
 * \param l_res Reference to the resulting left border.
 * \param r_res Reference to the resulting right border.
 * \return The size of the new interval [\ell_{new}..r_{new}].
 * \pre \f$ 0 \leq \ell \leq r < csa.size() \f$
 *
 * \par Time complexity
 *		 \f$ \Order{ t_{rank\_bwt} } \f$
 * \par References
 * 		Ferragina, P. and Manzini, G. 2000. Opportunistic data structures with applications. FOCS 2000.
 */
template<class Csa>
static typename Csa::csa_size_type backward_search(
    const Csa& csa,
    typename Csa::size_type l,
    typename Csa::size_type r,
    typename Csa::char_type c,
    typename Csa::size_type& l_res,
    typename Csa::size_type& r_res
)
{
    assert(l <= r); assert(r < csa.size());
    typename Csa::size_type c_begin = csa.C[csa.char2comp[c]];
    l_res = c_begin + csa.rank_bwt(l, c); // count c in bwt[0..l-1]
    r_res = c_begin + csa.rank_bwt(r+1, c) - 1; // count c in bwt[0..r]
    assert(r_res+1-l_res >= 0);
    return r_res+1-l_res;
}

//! Backward search for a pattern pat on an interval \f$[\ell..r]\f$ of the suffix array.
/*!
 * \param csa The csa in which the backward_search should be done.
 * \param l Left border of the lcp-interval \f$ [\ell..r]\f$.
 * \param r Right border of the lcp-interval \f$ [\ell..r]\f$.
 * \param pat The string which is the prefix of the suffixes in the resulting interval \f$ [\ell_{new}..r_{new}] \f$ .
 * \param len The length of pat.
 * \param l_res Reference to the resulting left border.
 * \param r_res Reference to the resulting right border.
 * \return The size of the new interval [\ell_{new}..r_{new}].
 * \pre \f$ 0 \leq \ell \leq r < csa.size() \f$
 *
 * \par Time complexity
 * 		\f$ \Order{ len \cdot t_{rank\_bwt} } \f$
 * \par References
 * 		Ferragina, P. and Manzini, G. 2000. Opportunistic data structures with applications. FOCS 2000.
 */
template<class Csa>
static typename Csa::csa_size_type backward_search(const Csa& csa, typename Csa::size_type l,
        typename Csa::size_type r, typename Csa::pattern_type pat,
        typename Csa::size_type len, typename Csa::size_type& l_res, typename Csa::size_type& r_res)
{
    typename Csa::size_type i = 0;
    l_res = l; r_res = r;
    while (i < len and backward_search(csa, l_res, r_res, pat[len-i-1], l_res, r_res)) {
        ++i;
    }
    return r_res+1-l_res;
}

// TODO: forward search. Include original text???

//! Counts the number of occurrences of pattern pat in the string of the compressed suffix array csa.
/*!
 * \param csa The compressed suffix array.
 * \param pat The pattern for which we count the occurences in the string of the compressed suffix array.
 * \param len The length of the pattern.
 * \return The number of occurences of pattern pat in the string of the compressed suffix array.
 *
 * \par Time complexity
 *		\f$ \Order{ t_{backward\_search} } \f$
 */
template<class Csa>
static typename Csa::csa_size_type count(const Csa& csa, typename Csa::pattern_type pat, typename Csa::size_type len)
{
    if (len > csa.size())
        return 0;
    typename Csa::size_type t1,t2; t1=t2=0; // dummy variable for the backward_search call
    return backward_search(csa, 0, csa.size()-1, pat, len, t1, t2);
}

//! Calculates all occurences of pattern pat in the string of the compressed suffix array csa.
/*!
 * \param csa The compressed suffix array.
 * \param pat The pattern for which we get the occurences in the string of the compressed suffix array.
 * \param len The length of the pattern.
 * \param occ A resizable random access container in which the occurences are stored.
 * \return The number of occurences of pattern pat of lenght len in the string of the compressed suffix array.
 *
 * \par Time complexity
 *		\f$ \Order{ t_{backward\_search} + z \cdot t_{SA} } \f$, where \f$z\f$ is the number of
 *							occurences of pat in the text of the compressed suffix array.
 */
template<class Csa, class RandomAccessContainer>
static typename Csa::csa_size_type locate(const Csa&  csa, typename Csa::pattern_type pattern,
        typename Csa::size_type pattern_len, RandomAccessContainer& occ)
{
    typename Csa::size_type occ_begin, occ_end, occs;
    occs = backward_search(csa, 0, csa.size()-1, pattern, pattern_len, occ_begin, occ_end);
    occ.resize(occs);
    for (typename Csa::size_type i=0; i < occs; ++i) {
        occ[i] = csa[occ_begin+i];
    }
    return occs;
}

//! Returns the substring T[begin..end] of the original text T from the corresponding compressed suffix array.
/*!
 * \param csa The compressed suffix array.
 * \param begin Index of the starting position (inclusive) of the substring in the original text.
 * \param end   Index of the end position (inclusive) of the substring in the original text.
 * \param text	A pointer to the extracted text. The memory has to be initialized before the call of the function!
 * \pre text has to be initialized with enough memory (end-begin+2 bytes) to hold the extracted text.
 * \pre \f$begin <= end\f$ and \f$ end < csa.size() \f$
 * \par Time complexity
 *		\f$ \Order{ (end-begin+1) \cdot t_{\Psi} + t_{SA^{-1}} } \f$
 */
// Is it cheaper to call T[i] = BWT[iSA[i+1]]??? Additional ranks but H_0 average access
// TODO: extract backward!!! is faster in most cases!
template<class Csa>
static void extract(const Csa& csa, typename Csa::size_type begin, typename Csa::size_type end, unsigned char* text)
{
    assert(end <= csa.size());
    assert(begin <= end);
    for (typename Csa::size_type i=begin, order = csa(begin); i<=end; ++i, order =  csa.psi[order]) {
        uint16_t c_begin = 1, c_end = 257, mid;
        while (c_begin < c_end) {
            mid = (c_begin+c_end)>>1;
            if (csa.C[mid] <= order) {
                c_begin = mid+1;
            } else {
                c_end = mid;
            }
        }
        text[i-begin] = csa.comp2char[c_begin-1];
    }
    if (text[end-begin]!=0)
        text[end-begin+1] = 0; // set terminal character
}

//! Reconstructs the text from position \f$begin\f$ to position \f$end\f$ (inclusive) from the compressed suffix array.
/*!
 * \param csa The compressed suffix array.
 * \param begin Starting position (inclusive) of the text to extract.
 * \param end   End position (inclusive) of the text to extract.
 * \return A std::string holding the extracted text.
 * \pre \f$begin <= end\f$ and \f$ end < csa.size() \f$
 * \par Time complexity
 *		\f$ \Order{ (end-begin+1) \cdot t_{\Psi} + t_{SA^{-1}} } \f$
 */
// TODO: use extract with four parameters to implement this method
template<class Csa>
static std::string extract(const Csa& csa, typename Csa::size_type begin, typename Csa::size_type end)
{
    assert(end <= csa.size());
    assert(begin <= end);
    std::string result(end-begin+1,' ');
    for (typename Csa::size_type i=begin, order = csa(begin); i<=end; ++i, order =  csa.psi[order]) {
        uint16_t c_begin = 1, c_end = 257, mid;
        while (c_begin < c_end) {
            mid = (c_begin+c_end)>>1;
            if (csa.C[mid] <= order) {
                c_begin = mid+1;
            } else {
                c_end = mid;
            }
        }
        result[i-begin] = csa.comp2char[c_begin-1];
    }
    return result;
}

//! Forward search for a character c on the path on depth \f$d\f$ to node \f$v\f$.
/*!
 *	\param cst 		The compressed suffix tree.
 *	\param v   		The node at the endpoint of the current edge.
 *	\param d   		The current depth of the path. 0 = first character on each edge of the root node.
 *	\param c   		The character c which should be matched at the path on depth \f$d\f$ to node \f$v\f$.
 *	\param char_pos One position in the text, which corresponds to the text that is already matched. If v=cst.root() and d=0 => char_pos=0.
 *
 *	\par Time complexity
 *		\f$ \Order{ t_{\Psi} } \f$ or \f$ \Order{t_{cst.child}} \f$
 */
template<class Cst>
typename Cst::cst_size_type forward_search(const Cst& cst, typename Cst::node_type& v,
        const typename Cst::size_type d, const typename Cst::char_type c,
        typename Cst::size_type& char_pos)
{
    unsigned char cc = cst.csa.char2comp[c]; // check if c occures in the text of the csa
    if (cc==0 and cc!=c)					 	 //   "    " "    "     "  "    "   "  "   "
        return 0;
    typename Cst::size_type depth_node = cst.depth(v);
    if (d < depth_node) {  //
        char_pos = cst.csa.psi[char_pos];
        if (char_pos < cst.csa.C[cc] or char_pos >= cst.csa.C[cc+1])
            return 0;
        return cst.leaves_in_the_subtree(v);
    } else if (d == depth_node) { //
        v = cst.child(v, c, char_pos);
        if (v == cst.root())
            return 0;
        else
            return cst.leaves_in_the_subtree(v);
    } else {
        throw std::invalid_argument("depth d should be smaller or equal than node depth of v!");
        return 0;
    }
}

//! Forward search for a pattern pat on the path on depth \f$d\f$ to node \f$v\f$.
/*!
 *	\param cst 		The compressed suffix tree.
 *	\param v   		The node at the endpoint of the current edge.
 *	\param d   		The current depth of the path. 0 = first character on each edge of the root node.
 *	\param pat   	The character c which should be matched at the path on depth \f$d\f$ to node \f$v\f$.
 *	\param len 		The length of the pattern.
 *	\param char_pos One position in the text, which corresponds to the text that is already matched. If v=cst.root() and d=0 => char_pos=0.
 *
 *	\par Time complexity
 *		\f$ \Order{ t_{\Psi} } \f$ or \f$ \Order{t_{cst.child}} \f$
 */
template<class Cst>
typename Cst::cst_size_type forward_search(const Cst& cst, typename Cst::node_type& v,
        typename Cst::size_type d, typename Cst::pattern_type pat, typename Cst::size_type len,
        typename Cst::size_type& char_pos
                                          )
{
    if (len==0)
        return cst.leaves_in_the_subtree(v);
    typename Cst::size_type i=0, size=0;
    while (i < len and (size=forward_search(cst, v, d, pat[i], char_pos))) {
        ++d;
        ++i;
    }
    return size;
}

//! Calculates the count method for a (compressed) suffix tree of type Cst.
/*! \param cst A const reference to the (compressed) suffix tree.
 *	\param pat The pattern we seach the number of occurences.
 *	\param len The length of the pattern.
 *	\return Number of occurences of the pattern in the text of the (compressed) suffix tree.
 */
template<class Cst>
typename Cst::cst_size_type count(const Cst& cst, typename Cst::pattern_type pat, typename Cst::size_type len)
{
    typename Cst::node_type v = cst.root();
    typename Cst::size_type char_pos = 0;
    return forward_search(cst, v, 0, pat, len, char_pos);
}

//! Calculates the locate method for a (compressed) suffix tree of type Cst.
/*! \param cst A const reference to the (compressed) suffix tree.
 *	\param pat The pattern for which we seach the occurences for.
 *	\param len The length of the pattern.
 *  \param occ A reference to a random access container in which we store the occurences of pattern in the text of the (compressed suffix array).
 *	\return Number of occurences of the pattern in the text of the (compressed) suffix tree.
 */
template<class Cst, class RandomAccessContainer>
typename Cst::cst_size_type locate(const Cst& cst, typename Cst::pattern_type pat,
                                   typename Cst::size_type len, RandomAccessContainer& occ)
{
    typedef typename Cst::size_type size_type;
    typename Cst::node_type v = cst.root();
    size_type char_pos = 0;
    typename Cst::cst_size_type occs = forward_search(cst, v, 0, pat, len, char_pos);
    occ.resize(occs);
    if (occs == 0)
        return 0; // because v equals cst.root()
    size_type left 	= cst.lb(v);   // get range in the
    size_type right = cst.rb(v);
    for (size_type i=left; i <= right; ++i)
        occ[i-left] = cst.csa[i];
    return occs;
}

//! Calculate the concatenation of edge labels from the root to the node v of the (compressed) suffix tree of type Cst.
/*!
 *	\param cst A const reference to the compressed suffix tree.
 *	\param v The node where the concatenation of the edge labels ends.
 *	\param text A pointer in which the string representing the concatenation of edge labels form the root to the node v will be stored.
 *   \pre text has to be initialized with enough memory (\f$ cst.depth(v)+1\f$ bytes) to hold the extracted text.
 */
template<class Cst>
void extract(const Cst& cst, const typename Cst::node_type& v, unsigned char* text)
{
    if (v == cst.root()) {
        text[0] = 0;
        return;
    }
    // first get the suffix array entry of the leftmost leaf in the subtree rooted at v
    typename Cst::size_type begin = cst.csa[cst.lb(v)];
    // then call the extract method on the compressed suffix array
    extract(cst.csa, begin, begin + cst.depth(v) - 1, text);
}

//! Calculate the concatenation of edge labels from the root to the node v of the (compressed) suffix tree of type Cst.
/*!
 *	\param cst A const reference to the compressed suffix tree.
 *	\param v The node where the concatenation of the edge labels ends.
 *	\return The string of the concatenated edge labels from the root to the node v.
 */
template<class Cst>
std::string extract(const Cst& cst, const typename Cst::node_type& v)
{
    if (v==cst.root()) {
        return std::string();
    }
    // first get the suffix array entry of the leftmost leaf in the subtree rooted at v
    typename Cst::size_type begin = cst.csa[cst.lb(v)];
    // then call the extract method on the compressed suffix array
    return extract(cst.csa, begin, begin + cst.depth(v) - 1);
}


/*
template<class Cst>
typename Cst::size_type count(
			const Cst &cst,
			typename Cst::pattern_type pattern,
		    typename Cst::size_type pattern_len){
	if(pattern_len==0){
		return 0;
	}
	typedef typename Cst::size_type size_type;
	typedef typename Cst::node_type node_type;
	node_type node = cst.root();
	for(size_type i=0, char_pos=0; cst.depth(node) < pattern_len; ++i){
		node_type newnode = cst.child(node, (unsigned char)pattern[cst.depth(node)], char_pos);
		if( newnode == cst.root() )// root node, no match found
			return 0;
		// else the first character of the newnode matches the pattern at position depth(node)
		for(size_type j=cst.depth(node)+1; j < cst.depth(newnode) and j < pattern_len; ++j){
			char_pos = cst.csa.psi[char_pos];
			size_type cc = cst.csa.char2comp[pattern[j]];
			if(char_pos < cst.csa.C[cc] or char_pos >= cst.csa.C[cc+1] )
				return 0;
		}
		node = newnode;
	}
	return cst.leaves_in_the_subtree(node);
}
*/
/*
template<class Cst, class RandomAccessContainer>
typename Cst::size_type locate(
		const Cst &cst,
		typename Cst::pattern_type pattern,
		typename Cst::size_type pattern_len,
		RandomAccessContainer &occ){
	occ.resize(0);
	typedef typename Cst::size_type size_type;
	typedef typename Cst::node_type node_type;
	node_type node = cst.root();
	for(size_type i=0, char_pos=0; cst.depth(node) < pattern_len; ++i){
		node_type newnode = cst.child(node, (unsigned char)pattern[cst.depth(node)], char_pos);
		if( newnode == cst.root() )// root node, no match found
			return 0;
		// else the first character of the newnode matches the pattern at position depth(node)
		for(size_type j=cst.depth(node)+1; j < cst.depth(newnode) and j < pattern_len; ++j){
			char_pos = cst.csa.psi[char_pos];
			size_type cc = cst.csa.char2comp[pattern[j]];
			if(char_pos < cst.csa.C[cc] or char_pos >= cst.csa.C[cc+1] )
				return 0;
		}
		node = newnode;
	}
	size_type occs = cst.leaves_in_the_subtree(node);
	occ.resize(occs);
	size_type left 	= cst.leftmost_suffix_array_index_in_the_subtree(node);
	size_type right = cst.rightmost_suffix_array_index_in_the_subtree(node);
	for(size_type i=left; i <= right; ++i)
		occ[i-left] = cst.csa[i];
	return occs;
}
*/


/*
template<class Csa, class Lcp, class Bp_support>
typename cst_sct<Csa, Lcp, Bp_support>::size_type count(
			const cst_sct<Csa, Lcp, Bp_support> &cst,
			typename cst_sct<Csa, Lcp, Bp_support>::pattern_type pattern,
		    typename cst_sct<Csa, Lcp, Bp_support>::size_type pattern_len){
	if(pattern_len==0){
		return 0;
	}
	typedef typename cst_sct<Csa, Lcp, Bp_support>::size_type size_type;
	typedef typename cst_sct<Csa, Lcp, Bp_support>::node_type node_type;
	node_type node = cst.root();
	for(size_type i=0, char_pos=0; node.l < pattern_len; ++i){
		node_type newnode = cst.child(node, (unsigned char)pattern[node.l], char_pos);
		if( newnode.l == 0 )// root node, no match found
			return 0;
		// else the first character of the newnode matches the pattern at position node.l
		for(size_type j=node.l+1; j < newnode.l and j< pattern_len; ++j){
			char_pos = cst.csa.psi[char_pos];
			size_type cc = cst.csa.char2comp[pattern[j]];
			if(char_pos < cst.csa.C[cc] or char_pos >= cst.csa.C[cc+1] )
				return 0;
		}
		node = newnode;
	}
	return cst.leaves_in_the_subtree(node);
}
*/

} // end namespace algorithm

} // end namespace sdsl

#endif
