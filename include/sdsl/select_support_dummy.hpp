/* sdsl - succinct data structures library
    Copyright (C) 2011 Simon Gog

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
/*! \file select_support_dummy.hpp
    \brief select_support_dummy.hpp contains classes that support a sdsl::bit_vector with constant time select information.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_SELECT_SUPPORT_DUMMY
#define INCLUDED_SDSL_SELECT_SUPPORT_DUMMY

#include "int_vector.hpp"
#include "select_support.hpp"
#include <string>
#include <iostream>

//#define SDSL_DEBUG_SELECT_SUPPORT_JMC

#ifdef SDSL_DEBUG_SELECT_SUPPORT_DUMMY
#include "testutils.hpp"
#endif

//! Namespace for the succinct data structure library.
namespace sdsl
{


//! A dummy class for select
/*!
 * @ingroup select_support_group
 */
class select_support_dummy : public select_support
{
    public:
        select_support_dummy(const int_vector<1>* v=NULL);
        select_support_dummy(const select_support_dummy& ss);
        ~select_support_dummy();
        void init(const int_vector<1>* v=NULL);
        //! Select function
        /*! \sa select_support.select
         */
        inline const size_type select(size_type i) const;
        //! Alias for select(i).
        inline const size_type operator()(size_type i)const;
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const;
        void load(std::istream& in, const int_vector<1>* v=NULL);
        void set_vector(const int_vector<1>* v=NULL);
        select_support_dummy& operator=(const select_support_dummy& ss);

        //! Swap operator
        void swap(select_support_dummy& ss);
        //! Equality Operator
        /*! Two select_support_dummys are equal if all member variables are equal.
         * Required for the Equality Comparable Concept of the STL.
         * \sa operator!=
         */
        bool operator==(const select_support_dummy& ss)const;
        //! Unequality Operator
        bool operator!=(const select_support_dummy& ss)const;
};



}

#endif
