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
/*! \file construct_isa.hpp
    \brief construct_isa.hpp contains a space and time efficient construction method for the inverse suffix array
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_CONSTRUCT_ISA
#define INCLUDED_SDSL_CONSTRUCT_ISA

#include "int_vector.hpp"
#include "util.hpp"

#include <iostream>
#include <stdexcept>
#include <list>

namespace sdsl
{

void construct_isa(cache_config& config);

}// end namespace

#endif
