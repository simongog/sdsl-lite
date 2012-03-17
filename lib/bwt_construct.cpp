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

namespace sdsl{
	
	/*! Constructs the Burrows and Wheeler Transform (BWT) from text and suffix array
	 * \param file_map A map, which contains the paths of the precalculated files like suffix array or text
	 * \param dir	   Directory in which the result should be written on disk. 
	 * \param id	   Id which should be used to build a file name for the calculated BWT.
	 * \par Space complexity: 
	 *        \f$n\f$ bytes
	 */
	bool construct_bwt( tMSS &file_map, const std::string &dir, const std::string &id){
		typedef int_vector<>::size_type size_type;
		if( file_map.find("bwt") == file_map.end() ){ // if bwt is not already registered in file_map 
			std::string bwt_file_name = dir+"bwt_"+id;
			std::ifstream bwt_in( bwt_file_name.c_str() );
			 // check if bwt is already on disk => register it
			if ( bwt_in ){
				file_map["bwt"] = bwt_file_name;
				bwt_in.close();
				return true;
			}
			write_R_output("csa", "construct BWT", "begin", 1, 0);
			const size_type buffer_size = 1000000;
			int_vector_file_buffer<8> text_buf( file_map["text"].c_str(), buffer_size );
			size_type n = text_buf.int_vector_size;
			int_vector_file_buffer<> sa_buf(file_map["sa"].c_str(), buffer_size);
			unsigned char *text = NULL;
			util::load_from_int_vector_buffer(text, text_buf);		

			std::ofstream bwt_out_buf( bwt_file_name.c_str(), std::ios::binary | std::ios::trunc | std::ios::out ); // open out file stream
			file_map["bwt"] = bwt_file_name;																		  // and save result to disk

			int_vector<8> bwt_buf(buffer_size);
			size_type bit_size = n*8;
			bwt_out_buf.write((char *) &(bit_size), sizeof(text_buf.int_vector_size) ); // write size
			size_type wb = 0; // written bytes

			size_type to_add[2] = {-1,n-1};

			for(size_type i=0, r_sum=0, r = 0; r_sum < n; ){
				for(; i<r_sum+r; ++i){
//					Variant (a) uses modulo and is very slow (double time of variant (e))
//					bwt_buf[i-r_sum] = text[ (sa_buf[i-r_sum]+n-1)%n ];
//					Variant (b) uses the if-condition (and is about 10-20 percent slower than variant (e))				
//					if( sa_buf[i-r_sum] )
//						bwt_buf[i-r_sum] = text[ sa_buf[i-r_sum]-1 ];
//					else
//						bwt_buf[i-r_sum] = text[ n-1 ];
//					Variant (c)
//					bwt_buf[i-r_sum] = sa_buf[i-r_sum] ? text[ sa_buf[i-r_sum]-1 ] : 0; // 
//					Variant (d)					
//					bwt_buf[i-r_sum] = text[ sa_buf[i-r_sum]+(sa_buf[i-r_sum]==0)*n-1 ]; // 
//					Variant (e) is the fastest and uses a very small lookup table 2010-11-03
					bwt_buf[i-r_sum] = text[ sa_buf[i-r_sum]+to_add[sa_buf[i-r_sum]==0] ]; // this solution is almost twice as fast
				}
				if(r > 0){
					bwt_out_buf.write((const char*)bwt_buf.data(), r);
					wb += r;
				}
				r_sum += r; r = sa_buf.load_next_block();
			}
			if(wb%8){
				bwt_out_buf.write("\0\0\0\0\0\0\0\0", 8-wb%8);
			}
			delete [] text; text = NULL;
			bwt_out_buf.close();
			write_R_output("csa", "construct BWT", "end", 1, 0);
		}
		return true;
	}	

	/* Constructs the Burrows and Wheeler Transform (BWT) from text and suffix array
	 * \param file_map A map, which contains the paths of the precalculated files like suffix array or text
	 * \param dir	   Directory in which the result should be written on disk. 
	 * \param id	   Id which should be used to build a file name for the calculated BWT.
	 * \par Space complexity: 
	 *        \f$2n\f$ bytes
	 */
	bool construct_bwt2( tMSS &file_map, const std::string &dir, const std::string &id){
		typedef int_vector<>::size_type size_type;
		write_R_output("csa", "construct BWT", "begin", 1, 0);
//		if( file_map.find("bwt") == file_map.end() ){ // if bwt is not already on disk => calculate it
			int_vector_file_buffer<> sa_buf(file_map["sa"].c_str());
			const size_type n = sa_buf.int_vector_size;

			if(n < 3)
				return construct_bwt(file_map, dir, id);

			unsigned char *text = NULL;
			int_vector_file_buffer<8> text_buf( file_map["text"].c_str() );
			util::load_from_int_vector_buffer(text, text_buf);		

			size_type cnt_c[257] = {0};   // counter for each character in the text
			size_type cnt_cc[257] = {0};  // prefix sum of the counter cnt_c
//			unsigned char alphabet[257] = {0};
//			uint8_t sigma = 0;

			write_R_output("csa", "construct C", "begin", 1, 0);
			for(size_type i=0; i<n; ++i){ // initialize cnt_c
				++cnt_c[text[i]+1];
			}
			write_R_output("csa", "construct C", "end", 1, 0);
			for(int i=1; i<257; ++i){  // calculate sigma and initailize cnt_cc
//				if( cnt_c[i] > 0 ){ alphabet[sigma++] = (unsigned char)(i-1); }
				cnt_cc[i] = cnt_c[i] + cnt_cc[i-1];
			}
//			alphabet[sigma] = '\0'; 
			
			size_type to_add[2] = {-1,n-1};

			int_vector<8> bwt(n,0);
			sa_buf.reset();
			for(size_type i=0, sai, r=0, r_sum=0 ; r_sum < n; ){
				for(; i < r_sum+r; ++i ){
					uint8_t bwti;
					sai = sa_buf[i-r_sum];	
					if( bwt[i] ){ // if the current BWT entry is already done
						bwti = bwt[i];
					}else{
						bwti = bwt[i] = text[sai+to_add[sai==0]];	
						size_type lf = cnt_cc[bwti];
						if( lf > i and sai > 1 ){
							bwt[lf] = text[sai-2];
						}
					}
					++cnt_cc[bwti];					 // update counter and therefore the LF information 
				}
				r_sum += r; r = sa_buf.load_next_block();
			}
//		}
		write_R_output("csa", "construct BWT", "end", 1, 0);
		return true;
	}
	
	
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

