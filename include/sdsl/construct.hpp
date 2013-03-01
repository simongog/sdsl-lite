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

#ifndef INCLUDED_SDSL_CONSTRUCT
#define INCLUDED_SDSL_CONSTRUCT

#include "sdsl_concepts.hpp"
#include "int_vector.hpp"
#include "lcp_construct.hpp"
#include "bwt_construct.hpp"
#include "qsufsort.hpp"
#include <string>

namespace sdsl{

template<class int_vector>
bool contains_no_zero_symbol(const int_vector& text, const char* file){
	for (int_vector_size_type i=0; i < text.size(); ++i){
		if ( (uint64_t)0 == text[i] ){
			throw std::logic_error((std::string("Error: File \"")+std::string(file)+std::string("\" contains zero symbol.")).c_str());
			return false;
		}
	}
	return true;
}

template<class int_vector>
void append_zero_symbol(int_vector& text){
	text.resize(text.size()+1);
	text[text.size()-1] = 0;
}


template<class Index>
void construct(Index& idx, const char* file, uint8_t num_bytes=0){
	tMSS file_map;
	cache_config config;
	construct(idx, file, config, num_bytes);
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

// Specialization for WTs and integer-alphabet
template<class Index>
void construct(Index& idx, const char* file, cache_config& config, uint8_t num_bytes, wt_tag, int_alphabet_tag){
	int_vector<> text;	
	util::load_vector_from_file(text, file, num_bytes);
	std::string tmp_key = util::to_string(util::get_pid())+"_"+util::to_string(util::get_id());
	std::string tmp_file_name = util::cache_file_name(tmp_key.c_str(), config);
	util::store_to_file(text, tmp_file_name.c_str());
	util::clear(text);
	int_vector_file_buffer<> text_buf(tmp_file_name.c_str());
	{
		Index tmp(text_buf, text_buf.int_vector_size);
		idx.swap(tmp);
	}
	std::remove(tmp_file_name.c_str());
}

// Specialization for WTs and byte-alphabet
template<class Index>
void construct(Index& idx, const char* file, cache_config& config, uint8_t num_bytes, wt_tag, byte_alphabet_tag){
	int_vector<8> text;	
	util::load_vector_from_file(text, file, num_bytes);
	std::string tmp_key = util::to_string(util::get_pid())+"_"+util::to_string(util::get_id());
	std::string tmp_file_name = util::cache_file_name(tmp_key.c_str(), config);
	util::store_to_file(text, tmp_file_name.c_str());
	util::clear(text);
	int_vector_file_buffer<8> text_buf(tmp_file_name.c_str());
	{
		Index tmp(text_buf, text_buf.int_vector_size);
		idx.swap(tmp);
	}
	std::remove(tmp_file_name.c_str());
}

// Specialization for CSAs and integer-alphabet
template<class Index>
void construct(Index& idx, const char* file, cache_config& config, uint8_t num_bytes, csa_tag, int_alphabet_tag){
	{// (1) check, if the text is cached
		if ( !util::cache_file_exists(constants::KEY_TEXT_INT, config) ){
			int_vector<> text;
			util::load_vector_from_file(text, file, num_bytes);
			if ( contains_no_zero_symbol(text, file) ){
				append_zero_symbol(text);
				util::store_to_cache(text, constants::KEY_TEXT_INT, config);
			}
		}
		util::register_cache_file(constants::KEY_TEXT_INT, config);
	}
	{// (2) check, if the suffix array is cached 
		if ( !util::cache_file_exists(constants::KEY_SA, config) ){
			int_vector<> sa; 
			sdsl::qsufsort::construct_sa(sa, config.file_map[constants::KEY_TEXT_INT].c_str(), 0);
			util::store_to_cache(sa, constants::KEY_SA, config);
		}
		util::register_cache_file(constants::KEY_SA, config);
	}
	{//  (3) construct BWT
		if ( !util::cache_file_exists(constants::KEY_BWT_INT, config) ){
			int_vector<> bwt;
			construct_int_bwt(bwt, config);
			util::store_to_cache(bwt, constants::KEY_BWT_INT, config);
		}
		util::register_cache_file(constants::KEY_BWT_INT, config);
	}
	{
		Index tmp(config.file_map, config.dir, config.id);
		idx.swap(tmp);
	}
	if ( config.delete_files ){
		util::delete_all_files(config.file_map);	
	}
}

// Specialization for CSAs and byte-alphabet
template<class Index>
void construct(Index& idx, const char* file, cache_config& config, uint8_t num_bytes, csa_tag, byte_alphabet_tag){
	{// (1) check, if the text is cached
		if ( !util::cache_file_exists(constants::KEY_TEXT, config) ){
			int_vector<8> text;
			util::load_vector_from_file(text, file, num_bytes);
			if ( contains_no_zero_symbol(text, file) ){
				append_zero_symbol(text);
				util::store_to_cache(text, constants::KEY_TEXT, config);
			}
		}
		util::register_cache_file(constants::KEY_TEXT, config);
	}
	{// (2) check, if the suffix array is cached 
		if ( !util::cache_file_exists(constants::KEY_SA, config) ){
			int_vector<8> text;
			util::load_from_cache(text, constants::KEY_TEXT, config);
			// call divsufsort
			int_vector<> sa(text.size(), 0, bit_magic::l1BP(text.size())+1);
			algorithm::calculate_sa((const unsigned char*)text.data(), text.size(), sa);
			util::store_to_cache(sa, constants::KEY_SA, config);
		}
		util::register_cache_file(constants::KEY_SA, config);
	}
	{//  (3) construct BWT
		if ( !util::cache_file_exists(constants::KEY_BWT, config) ){
			int_vector<8> bwt;
			construct_bwt(bwt, config);
			util::store_to_cache(bwt, constants::KEY_BWT, config);
		}
		util::register_cache_file(constants::KEY_BWT, config);
	}
	{
		Index tmp(config.file_map, config.dir, config.id);
		idx.swap(tmp);
	}
	if ( config.delete_files ){
		util::delete_all_files(config.file_map);
	}
}

// Specialization for CSTs and integer-alphabet
template<class Index>
void construct(Index& idx, const char* file, cache_config& config, uint8_t num_bytes, cst_tag, int_alphabet_tag){
	csa_tag csa_t;
	typename Index::csa_type::alphabet_category alph_t;
	{// (1) check, if the compressed suffix array is cached
		typename Index::csa_type csa;
		if ( !util::cache_file_exists(util::class_to_hash(csa).c_str(), config) ){
			cache_config csa_config(false, config.dir, config.id, config.file_map);
			construct(csa, file, csa_config, num_bytes, csa_t, alph_t);
			config.file_map = csa_config.file_map;
			util::store_to_cache(csa, util::class_to_hash(csa).c_str(), config); 
		}
		util::register_cache_file(util::class_to_hash(csa).c_str(), config);
	}
	{// (2) check, if the longest common prefix array is cached
		util::register_cache_file(constants::KEY_TEXT_INT, config);
		util::register_cache_file(constants::KEY_BWT_INT, config);
		util::register_cache_file(constants::KEY_SA, config);
		if ( !util::cache_file_exists(constants::KEY_LCP, config) ){
			int_vector<> lcp;
			construct_int_lcp_kasai(lcp, config);
			util::store_to_cache(lcp, constants::KEY_LCP, config);
		}
		util::register_cache_file(constants::KEY_LCP, config);
	}
	util::assign(idx, Index(config.file_map, config.dir, config.id));
	if ( config.delete_files ){
		util::delete_all_files(config.file_map);	
	}
}

// Specialization for CSTs and byte-alphabet
template<class Index>
void construct(Index& idx, const char* file, cache_config& config, uint8_t num_bytes, cst_tag, byte_alphabet_tag){
	csa_tag csa_t;
	typename Index::csa_type::alphabet_category alph_t;
	{// (1) check, if the compressed suffix array is cached
		typename Index::csa_type csa;
		if ( !util::cache_file_exists(util::class_to_hash(csa).c_str(), config) ){
			cache_config csa_config(false, config.dir, config.id, config.file_map);
			construct(csa, file, csa_config, num_bytes, csa_t, alph_t);
			config.file_map = csa_config.file_map;
			util::store_to_cache(csa, util::class_to_hash(csa).c_str(), config); 
		}
		util::register_cache_file(util::class_to_hash(csa).c_str(), config); 
	}
	{// (2) check, if the longest common prefix array is cached
		util::register_cache_file(constants::KEY_TEXT, config);
		util::register_cache_file(constants::KEY_BWT, config);
		util::register_cache_file(constants::KEY_SA, config);
		if ( !util::cache_file_exists(constants::KEY_LCP, config) ){
			int_vector<> lcp;
			construct_lcp_kasai(lcp, config);
			util::store_to_cache(lcp, constants::KEY_LCP, config);
		}
		util::register_cache_file(constants::KEY_LCP, config);
	}
	util::assign(idx, Index(config.file_map, config.dir, config.id));
	if ( config.delete_files ){
		util::delete_all_files(config.file_map);	
	}
}


} // end namespace sdsl
#endif
