/* sdsl - succinct data structures library
    Copyright (C) 2012 Simon Gog

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
/*! \file uint256_t.hpp
   \brief uint256_t.hpp contains a convenientunion for YMM registers (256-bits).
   \author Diego Havenstein
*/
#ifndef INCLUDED_SDSL_YMMUNION
#define INCLUDED_SDSL_YMMUNION

namespace sdsl
{

#ifdef __AVX2__
template<typename T>
union YMM_union {
  __m256i ymm;
  T values[32/sizeof(T)];
};
#endif

} // end namespace

#endif
