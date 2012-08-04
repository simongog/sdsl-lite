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
/*! \file rank_support.hpp
    \brief rank_support.hpp contains classes that support a sdsl::bit_vector with constant time rank information.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_RANK_SUPPORT
#define INCLUDED_SDSL_RANK_SUPPORT

/** \defgroup rank_support_group Rank Support (RS)
 * This group contains data structures which support an sdsl::bit_vector with the rank method.
 */

#include "int_vector.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! The base class of classes supporting rank_queries for a sdsl::bit_vector in constant time.
/*!
*/
class rank_support
{
    protected:
        const int_vector<1>* m_v; //!< Pointer to the rank supported bit_vector
    public:
        typedef int_vector<1>::size_type size_type;

        //! Constructor
        /*! \param v The supported bit_vector.
         */
        rank_support(const int_vector<1>* v = NULL);
        //! Copy constructor
        rank_support(const rank_support& rs);
        //! Destructor
        virtual ~rank_support() {}
        //! Initializes the data structure.
        /*! \param v The supported bit_vector. If v equals NULL the previous set
                   bit_vector is supported. Otherwise v will be supported.
        	\note Call this function before the first call of rank.
            \sa rank
         */
        virtual void init(const int_vector<1>* v=NULL) = 0;

        //! Answers rank queries for the supported bit_vector if init() was called before.
        /*!	\param i Argument for the length of the prefix v[0..i-1].
        	\returns Number of 1-bits in the prefix [0..i-1] of the supported bit_vector.
        	\note Method init has to be called before the first call of rank.
        	\sa init
         */
        virtual const size_type rank(size_type i) const = 0;
        //! Alias for rank(i)
        virtual const size_type operator()(size_type idx) const = 0;
        //! Serializes rank_support.
        /*! \param out Out-Stream to serialize the data to.
        */
        virtual size_type serialize(std::ostream& out, structure_tree_node* v, std::string name)const = 0;
        //! Loads the rank_support.
        /*! \param in In-Stream to load the rank_support data from.
            \param v The supported bit_vector.
         */
        virtual void load(std::istream& in, const int_vector<1>* v=NULL) = 0;
        //! Sets the supported bit_vector to the given pointer.
        /*! \param v The new bit_vector to support.
         *  \note Method init has to be called before the next call of rank.
         *  \sa init, rank
         */
        virtual void set_vector(const int_vector<1>* v=NULL) = 0;
};

inline rank_support::rank_support(const int_vector<1>* v)
{
    m_v = v;
}

inline rank_support::rank_support(const rank_support& rs)
{
    m_v = rs.m_v;
}

}// end namespace sdsl

#include "rank_support_v.hpp"
#include "rank_support_v5.hpp"
#include "rank_support_jmc.hpp"

#endif // end file 
