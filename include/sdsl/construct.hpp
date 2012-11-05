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
/*! \file construct.hpp
    \brief construct.hpp contains methods to construct indexes (compressed suffix arrays and trees).
	\author Simon Gog
*/

#ifndef SDSL_INCLUDED_CONSTRUCT
#define SDSL_INCLUDED_CONSTRUCT

#include "int_vector.hpp"
#include "sdsl_concepts.hpp"
#include "wavelet_trees.hpp"
#include "suffixarrays.hpp"
#include "suffixtrees.hpp"
#include "qsufsort.hpp"
#include <string>

namespace sdsl{

bool contains_no_zero_symbol(const int_vector<>& text, const char* file);
void append_zero_symbol(int_vector<>& text);



template<class Index>
void construct(Index& idx, const char* file, uint8_t num_bytes=0){
	tMSS file_map;
	construct(idx, file, cache_config(), num_bytes);
}
	
//! Constructs an index object of class Index for a text stored on disk.
/*!
 * \param idx       	Index object.  Any sdsl suffix array of suffix tree.
 * \param file      	Name of the text file. The representation of the file
 *                  	is dependent on the next parameter.
 * \
 * \param num_bytes 	If `num_bytes` equals 0, the file format is a serialized 
 *				    	int_vector<>. Otherwise the file is interpreted as sequence
 *                  	of `num_bytes`-byte integer stored in big endian order.
 */
template<class Index>
void construct(Index& idx, const char* file, cache_config& config, uint8_t num_bytes=0){
	// delegate to CSA or CST construction	
	typename Index::index_category 		index_tag;
	typename Index::alphabet_category 	alphabet_tag;
	construct(idx, file, config, num_bytes, index_tag, alphabet_tag);
}

template<class Index>
void construct(Index& idx, const char* file, cache_config& config, uint8_t num_bytes, csa_tag, int_alphabet_tag){
	std::cout<<"config.dir="<<config.dir<<std::endl;
	std::cout<<"config.id="<<config.id<<std::endl;
	{// (1) check, if the text is cached
		int_vector<> text;
		if ( !util::load_from_cache(text, constants::KEY_TEXT_INT, config) ){
			util::load_vector_from_file(text, file, num_bytes);
			if ( contains_no_zero_symbol(text, file) ){
				append_zero_symbol(text);
				util::store_to_cache(text, constants::KEY_TEXT_INT, config);
			}
		}
	}
	{// (2) check, if the suffix array is cached 
		int_vector<> sa; 
		if ( !util::load_from_cache(sa, constants::KEY_SA, config) ){
			sdsl::qsufsort::construct_sa(sa, config.file_map[constants::KEY_TEXT_INT].c_str(), 0);
			util::store_to_cache(sa, constants::KEY_SA, config);
		}
	//  (3) construct BWT
		int_vector<> bwt;
		if ( !util::load_from_cache(bwt, constants::KEY_BWT_INT, config) ){
			construct_int_bwt(bwt, config);
			util::store_to_cache(bwt, constants::KEY_BWT_INT, config);
		}
	}
	util::assign(idx, Index(config.file_map, config.dir,config.id));
	if ( config.delete_files ){
		util::delete_all_files(config.file_map);	
	}
}


template<class Index>
void construct(Index& idx, const char* file, cache_config& config, uint8_t num_bytes, cst_tag, int_alphabet_tag){
	csa_tag csa_t;
	typename Index::csa_type::alphabet_category alph_t;
	// TODO lookup if csa is cached
	{// (1) check, if the compressed suffix array is cached
		typename Index::csa_type csa;
		if ( !util::load_from_cache(csa, util::class_to_hash(csa).c_str(), config) ){
			construct(csa, file, config, num_bytes, csa_t, alph_t);
		}
		util::store_to_cache(csa, util::class_to_hash(csa).c_str(), config); 
	}
	{// (2) check, if the longest common prefix array is cached
		int_vector<> lcp;
		if ( !util::load_from_cache(lcp, constants::KEY_LCP, config) ){
			construct_int_lcp_kasai(lcp, config);
			util::store_to_cache(lcp, constants::KEY_LCP, config);
		}
	}
	util::assign(idx, Index(config.file_map, config.dir, config.id));
	if ( config.delete_files ){
		util::delete_all_files(config.file_map);	
	}
}


} // end namespace sdsl
#endif
