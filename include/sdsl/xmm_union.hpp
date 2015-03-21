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
/*! \file xmm_union.hpp
   \brief xmm_union.hpp contains a convenientunion for XMM registers (128-bits).
   \author Diego Havenstein
*/
#ifndef INCLUDED_SDSL_XMMUNION
#define INCLUDED_SDSL_XMMUNION

namespace sdsl
{

#ifdef __SSE4_2__
template<typename T>
union XMM_union {
  __m128i xmm;
  T values[16/sizeof(T)];
};
#endif

} // end namespace

#endif
