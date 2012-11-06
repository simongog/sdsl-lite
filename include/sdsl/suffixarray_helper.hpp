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
        psi_of_csa_psi(const CsaPsi* csa_psi=NULL) {
            m_csa = csa_psi;
        }

        //! Copy constructor
        psi_of_csa_psi(const psi_of_csa_psi& psi_of_csa) {
            m_csa = psi_of_csa.m_csa;
        }

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
            // TODO: replace by \sigma binary searches on PSI function??
            return (*m_csa)(((*m_csa)[i]+size()-1)%size());   // TODO:replace % by to_add table
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
        bwt_of_csa_psi(const CsaPsi* csa) {
            m_csa = csa;
        }

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




}

#endif
