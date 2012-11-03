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

namespace sdsl{

bool contains_no_zero_symbol(const int_vector<>& text, const char* file);
void append_zero_symbol(int_vector<>& text);

template<class Index>
void construct(Index& idx, const char* file, uint8_t num_bytes=0){
	tMSS file_map;
	construct(idx, file, file_map, false, num_bytes);
}
	
//! Constructs an index object of class Index for a text stored as serialized int_vector<>. 
/*!
 * \param idx       Index object.  Any sdsl suffix array of suffix tree.
 * \param file      Name of the text file. The representation of the file
 *                  is dependent on the next parameter.
 * \param num_bytes If `num_bytes` equals 0, the file format is a serialized 
 *				    int_vector<>. Otherwise the file is interpreted as sequence
 *                  of `num_bytes`-byte integer stored in big endian order.
 */
template<class Index>
void construct(Index& idx, const char* file, tMSS& file_map, bool cache=false, uint8_t num_bytes=0){
	// delegate to CSA or CST construction	
	typename Index::index_category 		index_tag;
	typename Index::alphabet_category 	alphabet_tag;
	construct(idx, file, file_map, cache, num_bytes, index_tag, alphabet_tag);
}

template<class Index>
void construct(Index& idx, const char* file, tMSS& file_map, bool cache, uint8_t num_bytes, csa_tag, int_alphabet_tag){
	int_vector<> text;
	util::load_vector_from_file(text, file, num_bytes);   // load text
//	std::cout << "text = " << text << std::endl;
	if ( contains_no_zero_symbol(text, file) ){
		append_zero_symbol(text);
		if ( util::store_to_file_map(text, constants::KEY_TEXT, file, file_map) ){ // TODO lookup if text is cached
			int_vector<> sa; // TODO lookup if sa is cached
			if ( !util::load_from_file_map(sa, constants::KEY_SA, file, file_map) ){
				sdsl::qsufsort::construct_sa(sa, file_map[constants::KEY_TEXT].c_str(), 0);
//			std::cout<<"---------\n sa: ";
//			for(int_vector_size_type i=0; i<sa.size(); ++i){std::cout<<" "<<sa[i]; }
//			std::cout<<"\n---------\n ";
				util::store_to_file_map(sa, constants::KEY_SA, file, file_map);
			}
			util::clear(sa);
			int_vector<> bwt;
			if( !util::load_from_file_map(bwt, constants::KEY_BWT, file, file_map)){
				construct_int_bwt(file_map, file); // TODO lookup if bwt is cached
			}
			util::clear(bwt);
			util::assign(idx, Index(file_map, "",""));
			// TODO: save CSA so that this can be used for constructing the CST
		}
	}
}

template<class Index>
void construct(Index& idx, const char* file, tMSS& file_map, bool cache, uint8_t num_bytes, cst_tag, int_alphabet_tag){
	typename Index::csa_type csa;
	csa_tag csa_t;
	typename Index::csa_type::alphabet_category alph_t;
	// TODO lookup if csa is cached
	construct(csa, file, file_map, cache, num_bytes, csa_t, alph_t);
	util::store_to_file_map(csa, util::class_to_hash(csa).c_str(), file, file_map); 
	util::clear(csa);
	construct_int_lcp_kasai(file_map, file, "./", "");
	util::assign(idx, Index(file_map, "./",""));
}


} // end namespace sdsl
#endif
