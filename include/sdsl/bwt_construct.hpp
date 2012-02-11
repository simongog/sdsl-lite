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

namespace sdsl{
	
	/*! Constructs the Burrows and Wheeler Transform (BWT) from text and suffix array
	 * \param file_map A map, which contains the paths of the precalculated files like suffix array or text
	 * \param dir	   Directory in which the result should be written on disk. 
	 * \param id	   Id which should be used to build a file name for the calculated BWT.
	 * \par Space complexity: 
	 *        \f$n\f$ bytes
	 */
	bool construct_bwt( tMSS &file_map, const std::string &dir, const std::string &id);

	/* Constructs the Burrows and Wheeler Transform (BWT) from text and suffix array
	 * \param file_map A map, which contains the paths of the precalculated files like suffix array or text
	 * \param dir	   Directory in which the result should be written on disk. 
	 * \param id	   Id which should be used to build a file name for the calculated BWT.
	 * \par Space complexity: 
	 *        \f$2n\f$ bytes
	 */
	bool construct_bwt2( tMSS &file_map, const std::string &dir, const std::string &id);
	
	
/*
	bool construct_bwt( tMSS &file_map, const std::string &dir, const std::string &id){
		typedef int_vector<>::size_type size_type;
		write_R_output("csa", "construct BWT", "begin", 1, 0);
		if( file_map.find("bwt") == file_map.end() ){ // if bwt is not already on disk => calculate it
			int_vector_file_buffer<8> text_buf( file_map["text"].c_str() );
			size_type n = text_buf.int_vector_size;
			int_vector_file_buffer<> sa_buf(file_map["sa"].c_str());
			unsigned char *bwt = new unsigned char[n+1];
			unsigned char *text = NULL;
			util::load_from_int_vector_buffer(text, text_buf);		
			for(size_type i=0, r_sum=0, r = sa_buf.load_next_block(); r_sum < n; ){
				for(; i<r_sum+r; ++i){
					bwt[i] = text[ (sa_buf[i-r_sum]+n-1)%n ];
				}
				r_sum += r; r = sa_buf.load_next_block();
			}
			delete [] text; text = NULL;
			write_R_output("csa", "store BWT", "begin", 1, 0);

			if( !util::store_to_file(char_array_serialize_wrapper<>(bwt,n), (dir+"bwt_"+id).c_str() ) ){
				throw std::ios_base::failure( "#csa_construct: Cannot store bwt to file system!" );
				return false;
			}else{
				file_map["bwt"] = dir+"bwt_"+id;
			}
			delete [] bwt; bwt = NULL;
			write_R_output("csa", "store BWT", "end", 1, 0);
		}
		write_R_output("csa", "construct BWT", "end", 1, 0);
		return true;
	}	
*/	

}// end namespace

#endif
