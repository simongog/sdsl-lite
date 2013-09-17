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
/*! \file construct_bwt.hpp
    \brief construct_bwt.hpp contains a space and time efficient construction method for the Burrows and Wheeler Transform (BWT).
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_CONSTRUCT_BWT
#define INCLUDED_SDSL_CONSTRUCT_BWT

#include "int_vector.hpp"
#include "sfstream.hpp"
#include "util.hpp"
#include "config.hpp" // for cache_config

#include <iostream>
#include <stdexcept>
#include <list>

namespace sdsl
{

//! Constructs the Burrows and Wheeler Transform (BWT) from text over byte- or integer-alphabet and suffix array.
/*!	The algorithm constructs the BWT and stores it to disk.
 *  \tparam t_width Width of the text. 0==integer alphabet, 8=byte alphabet.
 *  \param config	Reference to cache configuration
 *  \par Space complexity
 *		\f$ n \log \sigma \f$ bits
 *  \pre Text and Suffix array exist in the cache. Keys:
 *         * conf::KEY_TEXT for t_width=8 or conf::KEY_TEXT_INT for t_width=0
 *         * conf::KEY_SA
 *  \post BWT exist in the cache. Key
 *         * conf::KEY_BWT for t_width=8 or conf::KEY_BWT_INT for t_width=0
 */
template<uint8_t t_width>
void construct_bwt(cache_config& config)
{
    static_assert(t_width == 0 or t_width == 8 , "construct_bwt: width must be `0` for integer alphabet and `8` for byte alphabet");

    typedef int_vector<>::size_type size_type;
    typedef int_vector<t_width> text_type;
    typedef int_vector_buffer<t_width> bwt_type;
    const char* KEY_TEXT = key_text_trait<t_width>::KEY_TEXT;
    const char* KEY_BWT = key_bwt_trait<t_width>::KEY_BWT;

    //  (1) Load text from disk
    text_type text;
    load_from_cache(text, KEY_TEXT, config);
    size_type n = text.size();
    uint8_t bwt_width = text.width();

    //  (2) Prepare to stream SA from disc and BWT to disc
    size_type buffer_size = 1000000; // buffer_size is a multiple of 8!, TODO: still true?
    int_vector_buffer<> sa_buf(cache_file_name(conf::KEY_SA, config), std::ios::in, buffer_size);
    std::string bwt_file = cache_file_name(KEY_BWT, config);
    bwt_type bwt_buf(bwt_file, std::ios::out, buffer_size, bwt_width);

    //  (3) Construct BWT sequentially by streaming SA and random access to text
    size_type to_add[2] = {(size_type)-1,n-1};
    for (size_type i=0; i < n; ++i) {
        bwt_buf[i] = text[ sa_buf[i]+to_add[sa_buf[i]==0] ];
    }
    bwt_buf.close();
    register_cache_file(KEY_BWT, config);
}

}// end namespace

#endif
