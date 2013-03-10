/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

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
/*! \file algorithms.hpp
    \brief algorithms.hpp contains algorithms for suffixarrays.
	\author Simon Gog
*/

#ifndef INCLUDED_SDSL_ALGORITHMS
#define INCLUDED_SDSL_ALGORITHMS

#include "int_vector.hpp"

#include <stdexcept> // for exceptions
#include <iostream>
#include <cassert>
#include <stack>
#include <utility>


namespace sdsl
{

//! A helper class containing algorithms for succinct data structures.
/*!
	\author Simon Gog
 */
namespace algorithm
{
//	private:
//		algorithm(); // This helper class can not be instantiated

// Returns if the pair (a1,a2) is lex. less or equal than (b1,b2)
template<typename I>
static bool leq(I a1, I a2, I b1, I b2)
{
    return a1 < b1 || (a1 == b1 && a2 <= b2);
}

// Returns if the triple (a1,a2,a3) is lex. less or equal than (b1,b2,b3)
template<typename I>
static bool leq(I a1, I a2, I a3, I b1, I b2, I b3)
{
    return a1 < b1 || (a1 == b1 && leq(a2,a3,b2,b3));
}


//	public:
//! Calculate the zero-order entropy for a text T
/*!
 *  \param c Pointer to a 0-terminated string.
 *  \return The zero-order entropy of the text.
 */
double H_0(const unsigned char* c);

// Claculate the star entropy of T, see Manzini 2001 for details
double H_0s(const unsigned char* c);


//! Calculate the Inverse Suffix Array from a Suffix Array SA.
/*!
 *    - Time requirement: \f$ O( sa.size() ) \f$ i.e. linear.
 *    - Space requirement: No additional space needed.
 *  \param sa Suffix Array.
 *  \param isa Container to store the resulting inverse Suffix Array.
 */
template<class RandomAccessContainer1, class RandomAccessContainer2>
static void sa2isa(const RandomAccessContainer1& sa, RandomAccessContainer2& isa);

template<class RandomAccessContainer>
static void sa2isa(const RandomAccessContainer& sa, int_vector<>& isa);

//! Calculate the Inverse Permutation of a Permutation from \f$0..sa.size()-1\f$ in-place.
/*!
 * \param sa A reference to the permutation of the numbers \f$0..sa.size()-1\f$ stored in a random access container.
 * \par Time complexity
 *		\f$ \Order{sa.size()}, i.e. linear time complexity \f$
 * \par Note
 *  If there is space for the visited bits in the random access container the procedure is really implemented inplace.
 *  Otherwise we use a additional bit_vector of size sa.size() bits to store the indicator bits for the procedure.
 */
template<class RandomAccessContainer>
static void inverse_permutation_inplace(RandomAccessContainer& sa);

//! Calculate the previous smaller value (psv) array for a random access container a.
/*!
 * \param a Container to calculate the psv array.
 * \param psv Container that contains the result after the calculation.
 * \pre The array \e a contains only non negative values and a.size() == psv.size().
 * \post \f[ psv[i] = \begin{array}{rl} a.size() &\mbox{ if} \min\{ a[j] \mid j<i \} \geq a[i] \\
                                         max\{j\mid a[j] < a[i] \wedge j<i\} &\mbox{ otherwise.}\end{array} \f]
 */
template<class RandomAccessContainer1, class RandomAccessContainer2>
void calculate_psv(const RandomAccessContainer1& a, RandomAccessContainer2& psv);


//! Verify the result of the method caculate_psv
/*! \return True if the RandomAccessContainer psv is the previous smaller array of a.
  */
template<class RandomAccessContainer1, class RandomAccessContainer2>
static bool verify_psv(const RandomAccessContainer1& a, RandomAccessContainer2& psv);

template<class RandomAccessContainer1, class RandomAccessContainer2>
static void calculate_psv2(const RandomAccessContainer1& a, RandomAccessContainer2& psv);

//! Calculate the next smaller value (nsv) array for a random access container a.
/*!
 * \param a Container to calculate the nsv array.
 * \param nsv Container that contains the result after the calculation.
 * \pre The array \e a contains only non negative values and a.size() == nsv.size().
 * \post \f[ nsv[i] = \begin{array}{rl} 0 &\mbox{ if} \min\{ a[j] \mid j>i \} \geq a[i] \\
                                         min\{j\mid a[j] < a[i] \wedge j>i\} &\mbox{ otherwise.}\end{array} \f]
 */
template<class RandomAccessContainer1, class RandomAccessContainer2>
static void calculate_nsv(const RandomAccessContainer1& a, RandomAccessContainer2& nsv);

//! Verify the result of the method of calculate_nsv
template<class RandomAccessContainer1, class RandomAccessContainer2>
static bool verify_nsv(const RandomAccessContainer1& a, RandomAccessContainer2& nsv);

//! TODO: Impelement
template<class RandomAccessContainer, class Text>
static void calculate_lcp(const RandomAccessContainer& sa, const Text& text, RandomAccessContainer& lcp);

//! Verify that a Suffix Array sa is correct for a text c
/*!
 * \param c Text (c-string) to check the suffix array for.
   \param len Length of the text c.
   \param sa Suffix array to check.
   \return If the suffix array is correct for the text c.
   The suffix array sa is correct for c if
     - the values of sa lie in the range [0..len)
	 - all values are different
	 - c+sa[i] < c+sa[i+1] for all i<len-1
 */
template<class RandomAccessContainer>
static bool verify_sa(const unsigned char* c, typename RandomAccessContainer::size_type len, const RandomAccessContainer& sa);

//! Verify that a SelectSupport rs is correct for a int_vector<1>.
template<class SelectSupport>
inline static bool verify_select_support(const SelectSupport& ss, const int_vector<1>& b);

//! Verify that two Containers have the same number of elements and the all corresponding (i-th) elements (0<=i<c1.size()) of both containers are equal.
template<class Container1, class Container2>
static bool equal_container_values(const Container1& c1, Container2& c2);

//! Calculate the Inverse Suffix Array inplace.
/*! \param sa RandomAccessContainer that contains the Suffix Array and is replaced by the Inverse Suffix Array.
 * \par Time complexity
 * 		\f$ \Order{ 2*SA.size() }\f$, i.e. linear.
 * \par Space complexity
 *		Additional sa.size() bits.
 */
template<class RandomAccessContainer>
static void sa2isa_inplace(RandomAccessContainer& sa);

//! Calculate the \f$\Psi\f$-function for a given Burrows and Wheeler Transformation.
/* \param bwt Burrows and Wheeler Transformation.
* \param len Length of the bwt.
* \param psi Container of size len for the result.
* \par Time complexity
	*	\f$\Order{2n}\f$ i.e. linear.
* \par Space complexity
*	Space of bwt (\f$\Order{n}\f$ bits) + space of uncompressed \f$\Psi\f$-function (\f$\Order{4n}\f$ bits).
*/
template<class RandomAccessContainer>
static void bwt2psi(const unsigned char* bwt, typename RandomAccessContainer::size_type len, RandomAccessContainer& psi);

//! Calculate the \f$\Psi\f$-function for a given suffix array.
/*! \param sa Suffix Array to calculate the \f$\Psi\f$-function for.
 * \param psi RandomAccessContainer that will contain the resulting \f$\Psi\f$-function.
 * \par Time complexity
 *		\f$ \Order{3*SA.size() }\f$, i.e. linear.
 * \par Space complexity
 * 		 Additional \f$ sa.size() \cdot \log(RandomAccessContainer::size_type)\f$ bits
 */
template<class RandomAccessContainer1, class RandomAccessContainer2>
static void sa2psi(const RandomAccessContainer1& sa, RandomAccessContainer2& psi);

template<class RandomAccessContainer>
static void sa2psi(const RandomAccessContainer& sa, int_vector<>& psi);

//! Calculate the Longest Common Prefix Table (lcptab).
/*! Algorithm from Kasai et al. "Linear-Time Longest-Common-Prefix Computation in Suffix Arrays and Its Applications"
 *  \param sa Suffix Array to calculate the lcptab for.
 *  \param text Text to calculate the lcptab for.
 *  \param lcp RandomAccessContainer that will contain the resulting lcptab.
 *  \par Time complexity
 *		\f$ \Order{ SA.size() } \f$, i.e. linear.
 */
template<class RandomAccessContainer, class Text>
static void	calculate_lcp12(const RandomAccessContainer& sa, const Text& text, RandomAccessContainer& lcp);

//! Calculate the suffix array SA out of the \f$\Psi\f$-function and \f$ SA^{-1}[0]\f$.
/*! \param psi A \f$\Psi-\f$ function.
 * \param isa_0 \f$ SA^{-1}[0] \f$. If SA[0]=n \f$ SA^{-1}[0]=\Psi(0) \f$.
 * \param sa A RandomAccessContainer that will contain the resulting suffix array.
 * \par Time complexity
 *       \f$\Order{psi.size()}\f$, i.e. linear.
 */
template<class RandomAccessContainer1, class RandomAccessContainer2>
static void psi2sa(const RandomAccessContainer1& psi, const typename RandomAccessContainer1::size_type isa_0,  RandomAccessContainer2& sa);

template<class RandomAccessContainer>
static void psi2sa(const RandomAccessContainer& psi, const typename RandomAccessContainer::size_type isa_0,  int_vector<>& sa);
//! Calculate the inverse suffix array SA out of the \f$\Psi\f$-function and \f$ SA^{-1}[0]\f$.
/*! \param psi A \f$\Psi-\f$ function.
 * \param isa_0 \f$ SA^{-1}[0] \f$. If SA[0]=n \f$ SA^{-1}[0]=\Psi(0) \f$.
 * \param isa A RandomAccessContainer that will contain the resulting inverse suffix array.
 * \par Time complexity
 * 	    \f$\Order{psi.size()}\f$, i.e. linear.
 */
template<class RandomAccessContainer>
static void psi2isa(const RandomAccessContainer& psi, const typename RandomAccessContainer::size_type isa_0, RandomAccessContainer& isa);


template<class RandomAccessContainer1, class RandomAccessContainer2>
void calculate_psv(const RandomAccessContainer1& a, RandomAccessContainer2& psv)
{
    assert(psv.size() == a.size());
    if (a.empty())
        return;
    psv[0] = psv.size();
    assert(psv[0] == psv.size());
    std::stack<typename RandomAccessContainer1::size_type> psv_index;
    typename RandomAccessContainer1::value_type min_element = a[0];
    for (typename RandomAccessContainer1::size_type i=0; i < a.size(); ++i) {
        if (a[i] <= min_element) {
            while (!psv_index.empty())
                psv_index.pop();
            min_element = a[i];
            psv[i] = a.size();
            psv_index.push(i);
        } else { // a[i] > min_element => stack will not be empty
            while (a[psv_index.top()] >= a[i])
                psv_index.pop();
            psv[i] = psv_index.top();
            psv_index.push(i);
        }
    }
}

template<class RandomAccessContainer1, class RandomAccessContainer2>
bool verify_psv(const RandomAccessContainer1& a, RandomAccessContainer2& psv)
{
    if (a.size()!=psv.size())
        return false;
    typename RandomAccessContainer1::value_type min_element = a[0];
    for (typename RandomAccessContainer1::size_type i=0; i<a.size(); ++i) {
        if (a[i] <= min_element) {
            min_element = a[i];
            if (psv[i] != a.size()) // see definition of calculate_psv
                return false;
        } else {
            if (psv[i]>=i)
                return false;
            if (a[psv[i]] >= a[i])
                return false;
            for (typename RandomAccessContainer1::size_type j=psv[i]+1; j<i; ++j)
                if (a[j]<a[i])
                    return false;
        }
    }
    return true;
}


template<class RandomAccessContainer1, class RandomAccessContainer2>
void calculate_psv2(const RandomAccessContainer1& a, RandomAccessContainer2& psv)
{
    assert(psv.size() == a.size());
    if (a.empty())
        return;
    psv[0] = psv.size();
    assert(psv[0] == psv.size());
    // TODO implementing the algorithm with use of a stack
    psv[0] = psv.size();
    typedef std::pair<typename RandomAccessContainer1::value_type, typename RandomAccessContainer1::size_type> tPII;
    std::stack<tPII> psv_stack;
    typename RandomAccessContainer1::value_type min_element = a[0], ai;
    for (typename RandomAccessContainer1::size_type i=0; i < a.size(); ++i) {
        if ((ai=a[i]) <= min_element) {
            while (!psv_stack.empty())
                psv_stack.pop();
            min_element = ai;
            psv[i] = a.size();
            psv_stack.push(tPII(ai, i));
        } else { // a[i] > min_element => stack will not be empty
            while (psv_stack.top().first >= ai)
                psv_stack.pop();
            psv[i] = psv_stack.top().second;
            psv_stack.push(tPII(ai, i));
        }
    }
}

template<class RandomAccessContainer1, class RandomAccessContainer2>
void calculate_nsv(const RandomAccessContainer1& a, RandomAccessContainer2& nsv)
{
    assert(nsv.size() == a.size());
    if (a.empty())
        return;
    nsv[nsv.size()-1] = 0;
    std::stack<typename RandomAccessContainer1::size_type> nsv_index;
    typename RandomAccessContainer1::value_type min_element = a[nsv.size()-1];
    for (typename RandomAccessContainer1::size_type i=nsv.size(); i > 0; --i) {
        if (a[i-1] <= min_element) {
            while (!nsv_index.empty())
                nsv_index.pop();
            min_element = a[i-1];
            nsv[i-1] = 0;
            nsv_index.push(i-1);
        } else { // a[i] > min_element => stack will not be empty
            while (a[nsv_index.top()] >= a[i-1])
                nsv_index.pop();
            nsv[i-1] = nsv_index.top();
            nsv_index.push(i-1);
        }
    }
}


template<class RandomAccessContainer1, class RandomAccessContainer2>
bool verify_nsv(const RandomAccessContainer1& a, RandomAccessContainer2& nsv)
{
    if (a.size() != nsv.size())
        return false;
    typename RandomAccessContainer1::value_type min_element = a[a.size()-1];
    for (typename RandomAccessContainer1::size_type i=a.size(); i>0; --i) {
        if (a[i-1] <= min_element) {
            min_element = a[i-1];
            if (nsv[i-1] != 0) // see definition of calculate_nsv
                return false;
        } else {
            if (nsv[i-1] <= i-1)
                return false;
            if (a[nsv[i-1]] >= a[i-1])
                return false;
            for (typename RandomAccessContainer1::size_type j=i; j<nsv[i-1]; ++j)
                if (a[j]<a[i-1])
                    return false;
        }
    }
    return true;
}





template<class RandomAccessContainer>
void bwt2psi(const unsigned char* bwt, typename RandomAccessContainer::size_type len, RandomAccessContainer& psi)
{
    if (psi.size() != len)
        psi.resize(len);
    typename RandomAccessContainer::size_type C[256] = {0}, index_of_dollar = 0;
    for (typename RandomAccessContainer::size_type i=0; i<len; ++i) {
        ++C[bwt[i]];
        if (bwt[i]=='\0')
            index_of_dollar = i;
    }
    //	std::cerr<<"index of . = "<<index_of_dollar<<std::endl;
    for (uint16_t i=255; i!=0; --i)
        C[i] = C[i-1];
    C[0]=0;
    for (uint16_t i=1; i<256; ++i) {
        C[i] += C[i-1];
    }
    //	assert(C[bwt[0]]==0);
    //	psi[C[bwt[0]]] = index_of_dollar;
    //	++C[bwt[0]];
    for (typename RandomAccessContainer::size_type i=0; i<len; ++i) {
        psi[C[bwt[i]]] = i;
        ++C[bwt[i]];
    }
    /*	for(typename RandomAccessContainer::size_type i=index_of_dollar+1; i<len; ++i){
    		psi[C[bwt[i]]] = i;
    		++C[bwt[i]];
    	}
    */
    /*
    	typename RandomAccessContainer::size_type C[256] = {0}, index_of_dollar = 0;
    	for(typename RandomAccessContainer::size_type i=0; i<len+1; ++i){
    		++C[bwt[i]];
    		if(bwt[i]=='\0')
    			index_of_dollar = i;
    	}
    //	std::cerr<<"index of $ = "<<index_of_dollar<<std::endl;
    	for(uint16_t i=255;i!=0;--i)
    		C[i] = C[i-1];
    	C[0]=0;
    	for(uint16_t i=1;i<256;++i){
    		C[i] += C[i-1];
    	}
    	psi[C[bwt[0]]-1] = index_of_dollar-1;
    	++C[bwt[0]];
    	for(typename RandomAccessContainer::size_type i=1; i<index_of_dollar; ++i){
    		psi[C[bwt[i]]-1] = i-1;
    		++C[bwt[i]];
    	}
    	for(typename RandomAccessContainer::size_type i=index_of_dollar+1; i<len+1; ++i){
    		psi[C[bwt[i]]-1] = i-1;
    		++C[bwt[i]];
    	}
    */
}

template<class RandomAccessContainer1, class RandomAccessContainer2>
void sa2isa(const RandomAccessContainer1& sa, RandomAccessContainer2& isa)
{
    isa.resize(sa.size()); // init isa
    typename RandomAccessContainer1::size_type i = 0;
    for (typename RandomAccessContainer1::const_iterator sa_it = sa.begin(), end = sa.end(); sa_it != end; ++sa_it, ++i) {
        isa[ *sa_it ] = i;
    }
}

template<class RandomAccessContainer>
void sa2isa(const RandomAccessContainer& sa, int_vector<>& isa)
{
    isa.set_int_width(bit_magic::l1BP(sa.size())+1);
    isa.resize(sa.size()); // init isa
    typename RandomAccessContainer::size_type i = 0;
    for (typename RandomAccessContainer::const_iterator sa_it = sa.begin(), end = sa.end(); sa_it != end; ++sa_it, ++i) {
        isa[ *sa_it ] = i;
    }
}

template<class RandomAccessContainer>
void inverse_permutation_inplace(RandomAccessContainer& perm)
{
    if (perm.size() == 0)
        return;
    typedef bit_vector::size_type size_type;
#ifdef SDSL_DEBUG
    {
        // check if input permutation perm is really a permutation from \f$0..perm.size()-1\f$
        bit_vector visited(perm.size(),0);
        for (size_type i=0; i<perm.size(); ++i) {
            if (!(perm[i]>=0 and perm[i]<perm.size() and 0==visited[perm[i]]))
                throw std::invalid_argument("Input for inverse_permutation_inplace is not a permutation!");
            visited[perm[i]]=1;
        }
    }
#endif
    uint8_t  w 		= bit_magic::l1BP(perm.size()-1)+1;
    uint64_t b 		= 1ULL << ((w+1)%64);
    uint64_t biggest = b | (perm.size()-1);
    uint64_t perm_0 = perm[0];
    perm[0] = biggest;
    if (w != 63 and perm[0] == biggest) { // space left for real inplace calculation
#ifdef SDSL_DEBUG
        std::cout<<"using real in-place procedure for inverse_permutation_inplace"<<std::endl;
#endif
        perm[0] = perm_0;
        for (size_type i=0, j,jj,t; i < perm.size(); ++i) {
            if ((perm[i] & b) == 0) {
                j = perm[i]; jj = i;
                while ((perm[j] & b) == 0) {
                    t = perm[j];
                    perm[j] = jj | b; // mark perm[j] as visited
                    jj=j; j=t;
                }
            }
        }
        for (size_type i=0; i < perm.size(); ++i)
            perm[i] = perm[i] & (b-1);
    } else { // not enough space for real inplace calculation =>  use additional
        perm[0] = perm_0;
        bit_vector is_inverse(perm.size(), 0);				// indicator bit_vector!
        for (size_type i=0, j,jj,t; i<perm.size(); ++i) {
            if (!is_inverse[i]) {
                j = perm[i]; jj = i;
                while (!is_inverse[j]) {
                    is_inverse[j] = 1;  // mark perm[j] as visited
                    t = perm[j]; perm[j] = jj; jj=j; j=t;
                }
            }
        }
    }
#ifdef SDSL_DEBUG
    {
        // check if the result is really a permutation form \f$0..perm.size()-1\f$
        bit_vector visited(perm.size(),0);
        for (size_type i=0; i<perm.size(); ++i) {
            assert(perm[i]>=0 and perm[i]<perm.size() and 0==visited[perm[i]]);
            visited[perm[i]]=1;
        }
    }
#endif
}

template<class RandomAccessContainer1, class RandomAccessContainer2>
void sa2psi(const RandomAccessContainer1& sa, RandomAccessContainer2& psi)
{
    RandomAccessContainer2 isa; // temporary array for the inverse suffix array
    sa2isa(sa, isa);
    psi.resize(sa.size());
    typename RandomAccessContainer1::value_type tmp; //
    typename RandomAccessContainer2::iterator psi_it = psi.begin();
    for (typename RandomAccessContainer1::const_iterator sa_it = sa.begin(), end = sa.end(); sa_it != end; ++sa_it, ++psi_it) {
        if ((tmp = *sa_it+1) != sa.size())
            *psi_it = isa[tmp];
        else
            *psi_it = isa[0];
    }
}

template<class RandomAccessContainer>
void sa2psi(const RandomAccessContainer& sa, int_vector<>& psi)
{
    int_vector<> isa; // temporary array for the inverse suffix array
    sa2isa(sa, isa);
    psi.set_int_width(bit_magic::l1BP(sa.size())+1);
    psi.resize(sa.size());
    typename RandomAccessContainer::value_type tmp; //
    int_vector<>::iterator psi_it = psi.begin();
    for (typename RandomAccessContainer::const_iterator sa_it = sa.begin(), end = sa.end(); sa_it != end; ++sa_it, ++psi_it) {
        if ((tmp = *sa_it+1) != sa.size())
            *psi_it = isa[tmp];
        else
            *psi_it = isa[0];
    }
}

template<class RandomAccessContainer, class Text>
void calculate_lcp(const RandomAccessContainer& sa, const Text& text, RandomAccessContainer& lcp)
{
    lcp = sa;
    RandomAccessContainer isa;
    sa2isa(sa, isa);

    lcp[0] = 0;
    typename RandomAccessContainer::size_type i=0,j,k,l=0;
    for (typename RandomAccessContainer::const_iterator isa_it = isa.begin(), end = isa.end(); isa_it != end; ++isa_it, ++i) {
        if ((j = *isa_it)) {
            k = sa[j-1];
            while (text[k+l]==text[i+l])
                ++l;
            lcp[j] = l;
            l = (l==0)?0:l-1;
        }
    }
}

/*
TODO: add implementation and definition
template<class RandomAccessContainer>
void algorithm::calculate_lps(){

}
*/

template<class RandomAccessContainer1, class RandomAccessContainer2>
void psi2sa(const RandomAccessContainer1& psi, const typename RandomAccessContainer1::size_type isa_0, RandomAccessContainer2& sa)
{
    sa.resize(psi.size());
    if (psi.empty())
        return;
    typename RandomAccessContainer1::value_type isa_k = isa_0;
    for (typename RandomAccessContainer1::size_type k = 0, size=psi.size(); k < size; ++k, isa_k = psi[isa_k]) {
        sa[isa_k]	= k;
    }
}

template<class RandomAccessContainer>
void psi2sa(const RandomAccessContainer& psi, const typename RandomAccessContainer::size_type isa_0, int_vector<>& sa)
{
    sa.set_int_width(bit_magic::l1BP(psi.size())+1);
    sa.resize(psi.size());
    if (psi.empty())
        return;
    typename RandomAccessContainer::value_type isa_k = isa_0;
    for (typename RandomAccessContainer::size_type k = 0, size=psi.size(); k < size; ++k, isa_k = psi[isa_k]) {
        sa[isa_k]	= k;
    }
}

template<class RandomAccessContainer>
void psi2isa(const RandomAccessContainer& psi, const typename RandomAccessContainer::size_type isa_0, RandomAccessContainer& isa)
{
    isa = psi;
    if (psi.empty())
        return;
    typename RandomAccessContainer::value_type isa_k = isa_0;
    for (typename RandomAccessContainer::size_type k=0, size=psi.size(); k < size; ++k, isa_k = psi[isa_k]) {
        isa[k]	= isa_k;
    }
}

template<class RandomAccessContainer>
bool verify_sa(const unsigned char* c, typename RandomAccessContainer::size_type len, const RandomAccessContainer& sa)
{
    typedef typename RandomAccessContainer::size_type size_type;
    if (sa.size() != len) { // check length
        std::cerr<<"sa.size()!=len"<<std::endl;
        std::cerr<<sa.size()<<"!="<<len<<std::endl;
        return false;
    }
    {
        // check if values are in the range [0..len) and all are different
        int_vector<> occ(len);
        util::set_zero_bits(occ);
        for (typename RandomAccessContainer::const_iterator it=sa.begin(), end=sa.end(); it!=end; ++it) {
            size_type value = *it;
            if (value < len and !occ[value])
                occ[value] = 1;
            else {
                std::cerr<<"sa is not a permutation"<<std::endl;
                return false;
            }
        }
    }
    // check if the lexicographic order is correct
    if (sa.size()<2)
        return true;
    size_type v1,v2;
    v1 = v2 = sa[(size_t)0];
    for (typename RandomAccessContainer::const_iterator it=sa.begin()+1, end = sa.end(); it!=end; ++it) {
        v1 = v2;
        v2 = *it;
        size_type i=v1, j=v2;
        for (; i!=len and j!=len; ++i, ++j) {
            if (c[i] < c[j])
                break;
            else if (c[i]>c[j]) { // lex order is wrong!
                std::cerr<<"lex order is wrong"<<std::endl;
                std::cerr<<"v1="<<v1<<" v2="<<v2<<std::endl;
                //				std::cerr<<c+v1<<std::endl;
                //				std::cerr<<c+v2<<std::endl;
                return false;
            }
        }
    }
    return true;
}

template<class SelectSupport>
bool verify_select_support(const SelectSupport& ss, const int_vector<1>& b)
{
    uint64_t i=0,j=0;
    for (int_vector<1>::const_iterator it = b.begin(), end = b.end(); it != end; ++it, ++i) {
        if (*it) {
            ++j;// it's the j-th 1 detected
            if (ss.select(j)!=i) return false;
        }
    }
    return true;
}

template<class Container1, class Container2>
bool equal_container_values(const Container1& c1, Container2& c2)
{
    if (c1.size() != c2.size())
        return false;
    typename Container2::const_iterator c2_it = c2.begin();
    for (typename Container1::const_iterator c1_it = c1.begin(), c1_end = c1.end(); c1_it!=c1.end(); ++c1_it, ++c2_it)
        if (*c1_it != *c2_it)
            return false;
    return true;
}


//	template<uint8_t w>
//	void inplace_radix_sort(int_vector<w> &v);
//	template<class RandomAccessContainer, class Text>
//	static void calulate_lcp9(const RandomAccessContainer &sa, const Text &text, RandomAccessContainer& lcp);


} // end namespace algorithm

} // end namespace sdsl

#include "construct_sa.hpp"
#include "algorithms_for_balanced_parentheses.hpp"
#include "algorithms_for_compressed_suffix_arrays.hpp"
#include "algorithms_for_compressed_suffix_trees.hpp"
#include "algorithms_for_string_matching.hpp"



#endif
