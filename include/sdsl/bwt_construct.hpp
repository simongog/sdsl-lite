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
/*! \file bwt_construct.hpp
    \brief bwt_construct.hpp contains a space and time efficient construction method for the Burrows and Wheeler Transform (BWT).
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_BWT_CONSTRUCT
#define INCLUDED_SDSL_BWT_CONSTRUCT

#include "typedefs.hpp"
#include "int_vector.hpp"
#include "util.hpp"
#include "testutils.hpp"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <list>

namespace sdsl
{

/*! Constructs the Burrows and Wheeler Transform (BWT) from text and suffix array
 * \param file_map A map, which contains the paths of the precalculated files like suffix array or text
 * \param dir	   Directory in which the result should be written on disk.
 * \param id	   Id which should be used to build a file name for the calculated BWT.
 * \par Space complexity:
 *        \f$n\f$ bytes
 */
bool construct_bwt(tMSS& file_map, const std::string& dir, const std::string& id);

bool construct_int_bwt(tMSS& file_map, const char *file);

}// end namespace

#endif
