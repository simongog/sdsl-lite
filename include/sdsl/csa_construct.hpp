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
/*! \file csa_construct.hpp
    \brief csa_construct.hpp contains a space and time efficient construction method for csa_wt, csa_sada, csa_uncompressed.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_CSA_CONSTRUCT
#define INCLUDED_SDSL_CSA_CONSTRUCT

#include "typedefs.hpp"
#include "int_vector.hpp"
#include "algorithms_for_suffix_array_construction.hpp"
#include "algorithms_for_compressed_suffix_arrays.hpp"
#include "util.hpp"
#include "testutils.hpp"
#include "suffixarrays.hpp"
#include "bwt_construct.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

namespace sdsl{

	template<class Csa>
	static bool construct_csa(std::string file_name, Csa &csa){
		tMSS file_map;
		return construct_csa(file_name, csa, file_map, true, "./","");
	}	

	template<class Csa>
	static bool construct_csa(std::string file_name, Csa &csa, tMSS &file_map, bool delete_files=true, std::string dir="./", std::string id=""){
		uint64_t fs = 0;
		char *ccc = NULL; 	
		if( (fs=file::read_text((const char*)file_name.c_str(), ccc))>0 ){
			typename Csa::size_type n = strlen((const char*)ccc);
			if(id == "")
				id =  util::to_string(util::get_pid())+"_"+util::to_string(util::get_id()).c_str();
			if(fs != n + 1){
				std::cerr << "# WARNING: file \"" << file_name << "\" contains 0-bytes." << std::endl;
				algorithm::shift_text((char*)ccc, fs, true);
				n = fs-1;
			}

			if( !util::store_to_file(char_array_serialize_wrapper<>((unsigned char*)ccc,n+1), (dir+"text_"+id).c_str() ) ){
				throw std::ios_base::failure( "#csa_construct: Cannot store text to file system!" );
				delete [] ccc;
				return false;
			}
			else{
				file_map["text"] = (dir+"text_"+id).c_str();
			}
			delete [] ccc;
			return construct_csa(csa, file_map, delete_files, dir, id);
		}
		return false;
	}	

	template<class Csa>
	static bool construct_csa(Csa &csa, tMSS &file_map, bool delete_files=true, std::string dir="./", std::string id=""){
		write_R_output("csa", "construct CSA", "begin", 1, 0);
		int_vector_file_buffer<8> text_buf( file_map["text"].c_str() );
		typedef int_vector<>::size_type size_type;
		if(id=="")
			id =  util::to_string(util::get_pid())+"_"+util::to_string(util::get_id()).c_str();
		text_buf.reset();
		size_type n = text_buf.int_vector_size;

		// if sa file already exists
		std::ifstream in((dir+"sa_"+id).c_str());
		if( !in ){
			unsigned char *text = NULL;
			util::load_from_int_vector_buffer(text, text_buf);		

			typename Csa::size_type nn = strlen((const char*) text);
			if( nn+1 != n ){
				std::cerr << "# WARNING: text contains 0-bytes. nn=" << nn << " n=" << n << std::endl;
				algorithm::shift_text((char*)text, n);
			}

			write_R_output("csa", "construct SA", "begin", 1, 0);
			int_vector<> sa = int_vector<>(n, 0, bit_magic::l1BP(n+1)+1);
			algorithm::calculate_sa(text,n, sa);	 // calculate the suffix array sa of str
			delete [] text;

			assert(sa.size() == n);

			write_R_output("csa", "construct SA", "end", 1, 0);
			write_R_output("csa", "store SA", "begin", 1, 0);

			if( !util::store_to_file(sa, (dir+"sa_"+id).c_str() ) ){
				throw std::ios_base::failure( "#csa_construct: Cannot store SA to file system!" );
				return false;
			}else{
				file_map["sa"] = dir+"sa_"+id;
			}
			write_R_output("csa", "store SA", "end", 1, 0);
			{
				sa.resize(0);
				int_vector<> temp;
				temp.swap(sa);
			}
		}else{
			file_map["sa"] = dir+"sa_"+id;
		}

		
		write_R_output("csa", "encode CSA", "begin", 1, 0);
		csa.construct(file_map, dir, id); // TODO: for all three choices
// TODO replace line above by swap operation of csa
//		csa = Csa(file_map, dir, id);
		write_R_output("csa", "encode CSA", "end", 1, 0);

		if( delete_files ){
			util::delete_all_files(file_map);
		}
		write_R_output("csa", "construct CSA", "end", 1, 0);
		return true;
	}	


	template<class Csa>
	static bool construct_csa_of_reversed_text(std::string file_name, Csa &csa){
		typedef int_vector<>::size_type size_type;
		std::string tmp_rev_file_name = "./text_rev_"+util::to_string(util::get_pid())+"_"+util::to_string(util::get_id());
		char *text = NULL;
		size_type n = 0;
		if ( (n=file::read_text((const char*)file_name.c_str(), text)) > 0 ){
			--n; // since read_text appends a 0 byte
			for( size_type i=0; i < n/2; ++i ){
				std::swap(text[i], text[n-1-i] );
			}
			file::write_text((const char*)tmp_rev_file_name.c_str(),text, n);
			delete [] text;
		}else{
			std::cout<<"ERROR: text cannot be read from file "<<file_name<<std::endl;
			return false;
		}
		bool res = construct_csa(tmp_rev_file_name, csa);
		std::remove(tmp_rev_file_name.c_str());
		return res;
	}	


}// end namespace

#endif
