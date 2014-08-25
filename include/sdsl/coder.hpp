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
/*! \file coder.hpp
    \brief coder.hpp contains the coder namespace and includes the header files of sdsl::coder::fibonacci, sdsl::coder::elias_delta, and sdsl::coder::run_length
	\author Simon Gog
 */
#ifndef SDSL_CODER
#define SDSL_CODER

#include "int_vector.hpp"
#include "coder_fibonacci.hpp"
#include "coder_elias_delta.hpp"
#include "coder_elias_gamma.hpp"
#include "coder_comma.hpp"

namespace sdsl
{

//! Namespace for the different coder of the sdsl.
namespace coder
{

template<class Coder>
class run_length
{
    public:
        typedef uint64_t size_type;
        static void encode(uint64_t x, uint64_t*& z, uint8_t offset);
        static uint64_t encoding_length(const uint64_t* s, uint8_t s_offset, size_type bit_length);
};

template<class Coder>
typename run_length<Coder>::size_type run_length<Coder>::encoding_length(const uint64_t* s, uint8_t s_offset, size_type bit_length)
{
    assert(s_offset < 64);
    size_type i=0;
    uint64_t w = (*s >> s_offset);
    uint8_t last_bit = w&1;
    size_type result = 0;
    while (i < bit_length) {
        size_type len = 0;
        while (last_bit == (w&1) and  i < bit_length) {
//			std::cout<<w<<" "<<i<<std::endl;
            ++len; ++i; ++s_offset;
            w >>= 1;
            if (s_offset == 64) {
                s_offset = 0;
                w = *(++s);
            }
        }
//		std::cout<<"len="<<Coder::encoding_length(len)<<std::endl;
        last_bit = (w&1);
        result += Coder::encoding_length(len);
    }
    return result;
}


} // end namespace coder

} // end namespace sdsl

#endif
