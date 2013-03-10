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

//! Constructs the Burrows and Wheeler Transform (BWT) from text over byte- or integer-alphabet and suffix array.
/*!	The algorithm constructs the BWT and stores it to disk.
 *  \tparam t_width Width of the text. 0==integer alphabet, 8=byte alphabet.
 *  \param config	Reference to cache configuration
 *  \par Space complexity
 *		\f$ n \log \sigma \f$ bits
 *  \pre Text and Suffix array exist in the cache. Keys:
 *         * constants::KEY_TEXT for t_width=8 or constants::KEY_TEXT_INT for t_width=0
 *         * constants::KEY_SA
 *  \post BWT exist in the cache. Key
 *         * constants::KEY_BWT for t_width=8 or constants::KEY_BWT_INT for t_width=0
 */
template<uint8_t t_width>
void construct_bwt(cache_config& config){
    typedef int_vector<>::size_type size_type;
    typedef int_vector<t_width> text_type;
    typedef int_vector<t_width> bwt_type;
    const char * KEY_TEXT = key_text_trait<t_width>::KEY_TEXT;
    const char * KEY_BWT = key_bwt_trait<t_width>::KEY_BWT;

    //  (1) Load text from disk
    write_R_output("bwt", "load text", "begin", 1, 0);
    text_type text;
    util::load_from_cache(text, KEY_TEXT, config);
    size_type n = text.size(); 
    uint8_t bwt_width = text.get_int_width();
    write_R_output("bwt", "load text", "end", 1, 0);

    //  (2) Prepare to stream SA from disc and BWT to disc
    write_R_output("bwt", "prepare io", "begin", 1, 0);
    size_type buffer_size = 1000000; // buffer_size is a multiple of 8!
    int_vector_file_buffer<> sa_buf(util::cache_file_name(constants::KEY_SA, config));
    sa_buf.reset(buffer_size);

    bwt_type bwt_buf(buffer_size, 0, bwt_width);

    std::string bwt_file = util::cache_file_name(KEY_BWT, config);
    std::ofstream bwt_out_buf(bwt_file.c_str(), std::ios::binary | std::ios::app | std::ios::out);   // open buffer for bwt
    size_type bit_size = n*bwt_width;
    bwt_out_buf.write((char*) &(bit_size), sizeof(bit_size));	// write size of vector
    if ( t_width != 8) {
        bwt_out_buf.write((char*) &(bwt_width),sizeof(bwt_width));  // write int_width of vector
    }
    write_R_output("bwt", "prepare io", "end", 1, 0);

    //  (3) Construct BWT sequentially by streaming SA and random access to text
    write_R_output("bwt", "construct BWT", "begin", 1, 0);
    size_type wb = 0;  // bytes written into bwt int_vector
    size_type to_add[2] = {(size_type)-1,n-1};
    for (size_type i=0, r_sum=0, r=0; r_sum < n;) {
        for (; i < r_sum+r; ++i) {
            bwt_buf[i-r_sum] = text[ sa_buf[i-r_sum]+to_add[sa_buf[i-r_sum]==0] ]; 
        }
        if (r > 0) {
            size_type cur_wb = (r*bwt_buf.get_int_width()+7)/8;
            bwt_out_buf.write((const char*)bwt_buf.data(), cur_wb);
            wb += cur_wb;
        }
        r_sum += r;
        r = sa_buf.load_next_block();
    }
    if (wb%8) {
        bwt_out_buf.write("\0\0\0\0\0\0\0\0", 8-wb%8);
    }
    bwt_out_buf.close();
    util::register_cache_file(KEY_BWT, config);
    write_R_output("bwt", "construct BWT", "end", 1, 0);
}


}// end namespace

#endif
