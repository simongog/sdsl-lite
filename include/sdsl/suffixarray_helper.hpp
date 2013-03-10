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
/*! \file suffixarray_helper.hpp
    \brief suffixarray_helper.hpp contains some helper classes for (compressed suffix arrays)
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_SUFFIXARRAY_HELPER
#define INCLUDED_SDSL_SUFFIXARRAY_HELPER

#include <stdint.h>
#include <cstdlib>
#include <cassert>
#include "algorithms_for_compressed_suffix_arrays.hpp"
#include "iterators.hpp"

namespace sdsl
{

//! A helper class for the \f$\Psi\f$ and LF function for (compressed) suffix arrays that are based on the compressed \f$\Psi\f$ function.
template<class CsaPsi>
class psi_of_csa_psi
{
    public:
        typedef typename CsaPsi::value_type value_type;
        typedef typename CsaPsi::size_type size_type;
        typedef typename CsaPsi::difference_type difference_type;
        typedef typename CsaPsi::enc_vector_type::const_iterator const_iterator;
    private:
        const CsaPsi* m_csa; //<-	pointer to the (compressed) suffix array that is based on the compressed psi function
    public:

        //! Constructor
        psi_of_csa_psi(const CsaPsi* csa_psi) : m_csa(csa_psi) { }

        //! Copy constructor
        psi_of_csa_psi(const psi_of_csa_psi& psi_of_csa) : m_csa(psi_of_csa.m_csa){ }

        //! Calculate the \f$\Psi\f$ value at position i.
        /*!	\param i The index for which the \f$\Psi\f$ value should be calculated, \f$i\in [0..size()-1]\f$.
         *	\par Time complexity
         *		\f$ \Order{1} \f$, if the enc_vector_type of the csa returns the entry in constant time.
         */
        value_type operator[](size_type i)const {
            assert(m_csa != NULL);
            assert(i < size());
            return m_csa->m_psi[i];
        }

        //! Calculate the LF mapping at position i.
        /*! The LF mapping function is the inverse to the \f$\Psi\f$ function. That is \f$LF[\Psi(i)]=i\f$.
         *  \param i The index for which the LF value should be calculated, \f$i\in [0..size()-1]\f$.
         */
        value_type operator()(size_type i)const {
            // TODO: in case of a very sparse sampling of SA it may be faster to
		    //  use \sigma binary searches on PSI function to determine the 
			// LF values.
            return (*m_csa)(((*m_csa)[i]+size()-1)%size());
        }

        //! Returns the size of the \f$\Psi\f$ function.
        size_type size()const {
            return m_csa->size();
        }

        //! Returns if the \f$\Psi\f$ function is empty.
        size_type empty()const {
            return m_csa->empty();
        }

        //! Swap operation require by the STL.
        void swap(psi_of_csa_psi& psi_of_csa) {
            if (this != &psi_of_csa) {
                ;// std::swap(m_csa, psi.m_csa_wt); // do not exchange the supported structure!
            }
        }

        //! Returns a const_iterator to the first element.
        /*! Required for the STL Container Concept.
         *  \sa end
         */
        const_iterator begin()const {
            return m_csa->m_psi.begin();//const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        /*! Required for the STL Container Concept.
         *  \sa begin.
         */
        const_iterator end()const {
            return m_csa->m_psi.end();//const_iterator(this, size());
        }

        // Get the number of sampled \f$\Psi\f$ values.
        uint32_t get_sample_dens()const {
            return m_csa->m_psi.get_sample_dens();
        }
};

//! A helper class for the \f$\Psi\f$ function for (compressed) suffix arrays which provide also the inverse suffix array values (like sdsl::csa_bitcompressed).
template<class Csa>
class psi_of_sa_and_isa
{
    public:
        typedef typename Csa::value_type value_type;
        typedef typename Csa::size_type size_type;
        typedef typename Csa::difference_type	difference_type;
        typedef random_access_const_iterator<psi_of_sa_and_isa> 		 const_iterator;// STL Container requirement
    private:
        Csa* m_csa; //<- pointer to the (full text index) suffix array
        value_type m_size_m1; // size minus 1
        value_type to_add1[2], to_add2[2];
    public:

        //! Constructor
        psi_of_sa_and_isa(Csa* csa=NULL) {
            m_csa = csa;
            if (csa != NULL)
                m_size_m1 = m_csa->size()-1;
            else
                m_size_m1 = (size_type)-1;
            to_add1[0] = 1;
            to_add1[1] = -m_size_m1;
            to_add2[0] = m_size_m1;
            to_add2[1] = -1;
        }

        // Copy constructor
        psi_of_sa_and_isa(const psi_of_sa_and_isa& csa) {
            m_csa = csa.m_csa;
            m_size_m1 = csa.m_size_m1;
            for (int i=0; i<2; ++i) {
                to_add1[i] = csa.to_add1[i];
                to_add2[i] = csa.to_add2[i];
            }
        }

        //! Calculate the \f$\Psi\f$ value at position i.
        /*!	\param i The index for which the \f$\Psi\f$ value should be calculated, \f$i\in [0..size()-1]\f$.
         *	\par Time complexity
         *		\f$ \Order{\saaccess+\isaaccess} \f$
        */
        value_type operator[](size_type i)const {
            assert(m_csa!=NULL);
            assert(i<size());
            // \f$\Psi[i] = SA^{-1}[SA[i]+1 \mod n]\f$, where \f$n\f$ is the length of the suffix array SA
            value_type sai = (*m_csa)[i];
            return (*m_csa)(sai + to_add1[ sai == m_size_m1 ]);
        }

        //! Calculate the LF mapping at position i.
        /*! The LF mapping function is the inverse to the \f$\Psi\f$ function. That is \f$LF[\Psi(i)]=i\f$.
         *  \param i The index for which the LF value should be calculated, \f$i\in [0..size()-1]\f$.
         */
        value_type operator()(size_type i)const {
            assert(m_csa!=NULL);
            assert(i<size());
            value_type sai = (*m_csa)[i];
            return (*m_csa)(sai + to_add2[(bool)sai ]);
        }

        //! Assignment operator
        psi_of_sa_and_isa& operator=(const psi_of_sa_and_isa& csa) {
            if (this != &csa) {
                m_csa = csa.m_csa;
                m_size_m1 = csa.m_size_m1;
                for (int i=0; i<2; ++i) {
                    to_add1[i] = csa.to_add1[i];
                    to_add2[i] = csa.to_add2[i];
                }
            }
            return *this;
        }

        //! Returns the size of the \f$\Psi\f$ function.
        size_type size()const {
            return m_csa->size();
        }

        //! Returns if the \f$\Psi\f$ function is empty.
        size_type empty()const {
            return m_csa->empty();
        }
        //! Swap operation require by the STL.
        void swap(psi_of_sa_and_isa& psi) {
            if (this != &psi) {
                ;// std::swap(m_csa, psi.m_csa); // do not exchange the supported structure!
                std::swap(m_size_m1, psi.m_size_m1);
                for (int i=0; i<2; ++i) {
                    std::swap(to_add1[i], psi.to_add1[i]);
                    std::swap(to_add2[i], psi.to_add2[i]);
                }
            }
        }

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

        // Get the number of sampled \f$\Psi\f$ values. In this wrapper the number of sampled (i.e. explicitly stored) \f$\Psi\f$ values is zero.
        uint32_t get_sample_dens()const {
            return 0;
        }
};

//! A wrapper for the bwt of a compressed suffix array that is based on the \f$\psi\f$ function.
template<class CsaPsi>
class bwt_of_csa_psi
{
    public:
        typedef typename CsaPsi::char_type value_type;
        typedef typename CsaPsi::size_type size_type;
        typedef typename CsaPsi::difference_type difference_type;
        typedef random_access_const_iterator<bwt_of_csa_psi> const_iterator;// STL Container requirement
    private:
        const CsaPsi* m_csa; //<- pointer to the (compressed) suffix array that is based on the \f$\Psi\f$ function.
        bwt_of_csa_psi() {}
    public:

        //! Constructor
        bwt_of_csa_psi(const CsaPsi* csa) : m_csa(csa) { }

        //! Calculate the Burrows Wheeler Transform (BWT) at position i.
        /*!	\param i The index for which the BWT value should be calculated, \f$i\in [0..size()-1]\f$.
         *	\par Time complexity
         *		\f$ \Order{\log |\Sigma|} \f$
         */
        value_type operator[](size_type i)const {
            assert(m_csa != NULL);
            assert(i < size());
            size_type pos = m_csa->psi(i);
            return algorithm::get_ith_character_of_the_first_row(pos, *m_csa);
        }

        //! Returns the size of the \f$\Psi\f$ function.
        size_type size()const {
            return m_csa->size();
        }

        //! Returns if the bwt is empty.
        size_type empty()const {
            return m_csa->empty();
        }

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
};

//! A wrapper class for the \f$\Psi\f$ and LF function for (compressed) suffix arrays that are based on a wavelet tree (like sdsl::csa_wt).
template<class CsaWT>
class psi_of_csa_wt {
    public:
        typedef typename CsaWT::value_type value_type;
        typedef typename CsaWT::size_type size_type;
        typedef typename CsaWT::char_type char_type;
        typedef typename CsaWT::difference_type difference_type;
        typedef random_access_const_iterator<psi_of_csa_wt> const_iterator;// STL Container requirement
    private:
        const CsaWT* m_csa_wt; //<- pointer to the (compressed) suffix array that is based on a wavelet tree
        psi_of_csa_wt() {};    // disable default constructor
    public:
        //! Constructor
        psi_of_csa_wt(CsaWT* csa_wt) {
            m_csa_wt = csa_wt;
        }
        //! Calculate the \f$\Psi\f$ value at position i.
        /*!	\param i The index for which the \f$\Psi\f$ value should be calculated, \f$i\in [0..size()-1]\f$.
         *	\par Time complexity
         *		\f$ \Order{\log |\Sigma|} \f$
         */
        value_type operator[](size_type i)const {
            assert(m_csa_wt != NULL);
			assert(i < m_csa_wt->size() );
            char_type c = algorithm::get_ith_character_of_the_first_row(i, *m_csa_wt);
            return m_csa_wt->wavelet_tree.select(i - m_csa_wt->C[m_csa_wt->char2comp[c]] + 1 , c);
        }
        //! Apply \f$\Psi\f$ k times to the value at position i.
        /*!	\param i The index for which \f$\Psi\f$ should be applied k times, \f$i\in [0..size()-1]\f$.
         *  \param k Number of times \f$\Psi\f$ should be applied
         *	\par Time complexity
         *		\f$ \Order{\min\{k\cdot \log |\Sigma|, (s_{\SUF}+s_{\ISA})\log|\Sigma|\}} \f$
         */
        value_type psi_k(size_type i, size_type k)const {
            assert(m_csa_wt != NULL);
            if (k < m_csa_wt->sa_sample_dens + m_csa_wt->isa_sample_dens) {
                for (size_type j=0; j<k; ++j) {
                    i = (*this)[i];
                }
                return i;
            } else {
                size_type x = (*m_csa_wt)[i];
                x += k;
                if (x >= m_csa_wt->size()) {
                    x -= m_csa_wt->size();
                }
                return (*m_csa_wt)(x);
            }
        }
        //! Calculate the LF mapping at position i.
        /*! The LF mapping function is the inverse to the \f$\Psi\f$ function. That is \f$LF[\Psi(i)]=i\f$.
         *  \param i The index for which the LF value should be calculated, \f$i\in [0..size()-1]\f$.
         *	\par Time complexity
         *		\f$ \Order{\log |\Sigma|} \f$
         */
        value_type operator()(size_type i)const {
            assert(m_csa_wt != NULL);
            assert(i < size());
            typename CsaWT::char_type c;
            size_type j = m_csa_wt->wavelet_tree.inverse_select(i,c); // see documentation of inverse_select in wt_huff
            return m_csa_wt->C[ m_csa_wt->char2comp[c] ] + j;
        }
        //! Apply LF k times to the value at position i.
        /*!	\param i The index for which LF should be applied k times, \f$i\in [0..size()-1]\f$.
         *  \param k Number of times LF should be applied
         *	\par Time complexity
         *		\f$ \Order{\min\{k\cdot \log |\Sigma|, (s_{\SUF}+s_{\ISA})\log|\Sigma|\}} \f$
         */
        value_type lf_k(size_type i, size_type k)const {
            assert(m_csa_wt != NULL);
            if (k < m_csa_wt->sa_sample_dens + m_csa_wt->isa_sample_dens) {
                for (size_type j=0; j<k; ++j) {
                    i = (*this)(i);
                }
                return i;
            } else {
                size_type x = (*m_csa_wt)[i];
                if (x < k) {
                    x += m_csa_wt->size();
                }
                x -= k;
                return (*m_csa_wt)(x);
            }
        }
        //! Returns the size of the \f$\Psi\f$ function.
        size_type size()const {
            return m_csa_wt->size();
        }
        //! Returns if the \f$\Psi\f$ function is empty.
        size_type empty()const {
            return m_csa_wt->empty();
        }
        //! Returns a const_iterator to the first element.
        const_iterator begin()const {
            return const_iterator(this, 0);
        }
        //! Returns a const_iterator to the element after the last element.
        const_iterator end()const {
            return const_iterator(this, size());
        }
        // Get the number of sampled \f$\Psi\f$ values. In the wavelet tree approach the number of sampled (i.e. explicitly stored) \f$\Psi\f$ values is zero.
        uint32_t get_sample_dens()const {
            return 0;
        }
};

template<class CsaWT>
class bwt_of_csa_wt {
    public:
        typedef const typename CsaWT::char_type value_type;
        typedef typename CsaWT::size_type size_type;
        typedef typename CsaWT::difference_type difference_type;
        typedef random_access_const_iterator<bwt_of_csa_wt> const_iterator;// STL Container requirement
    private:
        const CsaWT* m_csa_wt; //<- pointer to the (compressed) suffix array that is based on a wavelet tree
        bwt_of_csa_wt(){};     // disable default constructor
    public:
        //! Constructor
        bwt_of_csa_wt(CsaWT* csa_wt) {
            m_csa_wt = csa_wt;
        }
        //! Calculate the Burrows Wheeler Transform (BWT) at position i.
        /*!	\param i The index for which the \f$\Psi\f$ value should be calculated, \f$i\in [0..size()-1]\f$.
         *	\par Time complexity
         *		\f$ \Order{\log |\Sigma|} \f$
         */
        value_type operator[](size_type i)const {
            assert(m_csa_wt != NULL);
            assert(i < size());
            return m_csa_wt->wavelet_tree[i];
        }
        //! Returns the size of the BWT function.
        size_type size()const {
            return m_csa_wt->size();
        }
        //! Returns if the BWT function is empty.
        size_type empty()const {
            return m_csa_wt->empty();
        }
        //! Returns a const_iterator to the first element.
        const_iterator begin()const {
            return const_iterator(this, 0);
        }
        //! Returns a const_iterator to the element after the last element.
        const_iterator end()const {
            return const_iterator(this, size());
        }
};



template<class Csa>
class text_of_csa
{
    public:
        typedef typename Csa::char_type value_type;
        typedef typename Csa::size_type size_type;
        typedef typename Csa::difference_type difference_type;
        typedef random_access_const_iterator<text_of_csa> const_iterator;// STL Container requirement
    private:
        const Csa* m_csa; //<- pointer to the (compressed) suffix array 
        text_of_csa() {}
    public:

        //! Constructor
        text_of_csa(const Csa* csa):m_csa(csa) { }

        //! Character at index \f$i\f$ of the original text.
        /*!	\param i Text position , \f$i\in [0..size()-1]\f$.
         *	\par Time complexity
         *		\f$ t_{ISA} \log\sigma \f$
         */
        value_type operator[](size_type i)const {
            assert(m_csa != NULL);
            assert(i < size());
            size_type pos = (*m_csa)(i);
            return algorithm::get_ith_character_of_the_first_row(pos, *m_csa);
        }

        //! Returns the size of the original text.
        size_type size()const {
            return m_csa->size();
        }

        //! Returns if text text has size 0.
        size_type empty()const {
            return m_csa->empty();
        }

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
};




}

#endif
