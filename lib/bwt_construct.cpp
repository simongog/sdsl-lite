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
#include "sdsl/bwt_construct.hpp"
#include <string>
#include <iostream>
#include <iomanip>

namespace sdsl
{

/*! Constructs the Burrows and Wheeler Transform (BWT) from text and suffix array
 * \param file_map A map, which contains the paths of the pre-calculated files like suffix array or text
 * \param dir	   Directory in which the result should be written on disk.
 * \param id	   Id which should be used to build a file name for the calculated BWT.
 * \par Space complexity:
 *        \f$n\f$ bytes
 */
bool construct_bwt(tMSS& file_map, const std::string& dir, const std::string& id)
{
    typedef int_vector<>::size_type size_type;
    if (file_map.find(constants::KEY_BWT) == file_map.end()) { // if bwt is not already registered in file_map
        std::string bwt_file_name = dir+constants::KEY_BWT+"_"+id;
        std::ifstream bwt_in(bwt_file_name.c_str());
        // check if bwt is already on disk => register it
        if (bwt_in) {
            file_map[constants::KEY_BWT] = bwt_file_name;
            bwt_in.close();
            return true;
        }
        write_R_output("csa", "construct BWT", "begin", 1, 0);
        const size_type buffer_size = 1000000;

		int_vector<8> text;
		util::load_vector_from_file(text, file_map[constants::KEY_TEXT].c_str(), 0 );
        size_type n = text.size();

        int_vector_file_buffer<> sa_buf(file_map[constants::KEY_SA].c_str(), buffer_size);

        std::ofstream bwt_out_buf(bwt_file_name.c_str(), std::ios::binary | std::ios::trunc | std::ios::out); // open out file stream
        file_map[constants::KEY_BWT] = bwt_file_name;														  // and save result to disk

        int_vector<8> bwt_buf(buffer_size);
        size_type bit_size = n*8;
        bwt_out_buf.write((char*) &(bit_size), sizeof(int_vector_size_type));   // write size
        size_type wb = 0; // written bytes
        size_type to_add[2] = {-1,n-1};
        for (size_type i=0, r_sum=0, r = 0; r_sum < n;) {
            for (; i<r_sum+r; ++i) {
                bwt_buf[i-r_sum] = text[ sa_buf[i-r_sum]+to_add[sa_buf[i-r_sum]==0] ]; 
            }
            if (r > 0) {
                bwt_out_buf.write((const char*)bwt_buf.data(), r);
                wb += r;
            }
            r_sum += r; r = sa_buf.load_next_block();
        }
        if (wb%8) {
            bwt_out_buf.write("\0\0\0\0\0\0\0\0", 8-wb%8);
        }
        bwt_out_buf.close();
        write_R_output("csa", "construct BWT", "end", 1, 0);
    }
    return true;
}

bool construct_int_bwt(tMSS& file_map, const char *file){
	int_vector<> text;
	util::load_from_file(text, file_map[constants::KEY_TEXT].c_str());
	int_vector_file_buffer<> sa_buf(file_map[constants::KEY_SA].c_str());
	int_vector<> bwt(text.size(), 0, text.get_int_width());
    int_vector_size_type to_add[2] = {-1,text.size()-1};
	for (int_vector_size_type i=0, r_sum=0, r = 0; r_sum < text.size();) {
		for (; i<r_sum+r; ++i) {
			bwt[i] = text[ sa_buf[i-r_sum]+to_add[sa_buf[i-r_sum]==0] ]; 
		}
		r_sum += r; r = sa_buf.load_next_block();
	}
	util::store_to_file_map(bwt, constants::KEY_BWT, file, file_map);	
	return true;	
}

}// end namespace
