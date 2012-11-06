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
#include "config.hpp" // for cache_config

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <list>

namespace sdsl
{

/*! Constructs the Burrows and Wheeler Transform (BWT) from text and suffix array store on disk
 * \param bwt 		Will contain the resulting BWT after the call of the method.
 * \param config	Config object for location of suffix array and text on disk.
 * \par Space complexity:
 *        \f$n\f$ bytes
 */
void construct_bwt(int_vector<8>& bwt, const cache_config& conf);

void construct_int_bwt(int_vector<> &bwt, const cache_config& conf);

}// end namespace

#endif
