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
/*! \file select_support.hpp
    \brief select_support.hpp contains classes that support a sdsl::bit_vector with constant time select information.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_SELECT_SUPPORT
#define INCLUDED_SDSL_SELECT_SUPPORT

/** \defgroup select_support_group Select Support (SCS)
 * This group contains data structures which support an sdsl::bit_vector with the select method.
 */

#include "int_vector.hpp"
#include "rank_support.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{
//! The base class of classes supporting select queries for a sdsl::bit_vector in constant time.
/*! Abstract base class for classes supporting select queries.
 */
class select_support
{
    protected:
        const int_vector<1>* m_v; //!< Pointer to the select supported sdsl::bit_vector.
    public:
        typedef int_vector<1>::size_type size_type;
        const bit_vector* v;

        //! Constructor of select_support.
        /*! \param v The bit_vector to support rank queries.
         */
        select_support(const int_vector<1>* f_v=NULL):v(f_v) {
            m_v = f_v;
        }
        //! Copy constructor
        /*! Copy the whole select_support including the  pointer
         *  to the supported bit_vector.
         */
        select_support(const select_support& f_v);
        //! Destructor of select_support.
        virtual ~select_support() {};
        //! Initalization method for select_support.
        /*! Init takes no arguments and should be called
        	before the first call to the select method if not
        	  - load is called to initialize the select_support or
        	  - the constructor is called with the pointer to the supported bit_vector.
        	\sa select, load.
         */
        virtual void init(const int_vector<1>* v=NULL) = 0;
        //! Select returns the index of the i-th 1-bit in the supported bit_vector.
        /*!	\param i Argument to calculate the index of the i-th 1-bit in the supported bit_vector.
        	\return The index \f$\in [0..v.size()-1]\f$ of the i-th 1-bit in the supported bit_vector.
        	Call init or load to initialize the data structure before the first call of this method.
         	\sa init, load.
         */
        virtual const size_type select(size_type i) const = 0;

        //! Alias for select
        virtual const size_type operator()(size_type i) const = 0;
        //! Serialize the select_support to an out file stream.
        virtual size_type serialize(std::ostream& out)const = 0;
        //! Load the select_support from an in file stream.
        /*!	Load an previously serialized select_support from a std::istream.
            This method could replace the call of init before
        	the first call of the select method.
        	\param in The std::istream to load the select_support.
        	\param v The bit_vector to be supported.
        	\sa init, select.
         */
        virtual void load(std::istream& in, const int_vector<1>* v=NULL) = 0;

        //! This method sets the supported bit_vector
        /*! \note Call the init function before you call select
         *   the first time after you changed the supported bit_vector.
         */
        virtual void set_vector(const int_vector<1>* v=NULL) = 0;
};


} // end namespace sdsl

#include "select_support_bs.hpp"
#include "select_support_mcl.hpp"
#include "select_support_dummy.hpp"

#endif
