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
    \brief csa_construct.hpp contains a space and time efficient construction method for csa_wt, csa_sada, csa_bitcompressed.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_CSA_CONSTRUCT
#define INCLUDED_SDSL_CSA_CONSTRUCT

#include "typedefs.hpp"
#include "sdsl_concepts.hpp"
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

namespace sdsl
{

//! Construct the compressed suffix array of a file.
/*!
 *  \param file_name	Name of the file.
 *  \param csa			CSA object which should be used to store the
 *                      compressed suffix array. 
 */	
template<class Csa>
static bool construct_csa(std::string file_name, Csa& csa) {
	typename Csa::alphabet_category tag;
	return construct_csa(file_name, csa, tag);
}

template<class Csa>
static bool construct_csa(std::string file_name, Csa& csa, byte_alphabet_tag) {
    tMSS file_map;
  //  return construct_csa(file_name, csa, file_map, true, "./","", byte_alphabet_tag);
    return construct_csa(file_name, csa, file_map, true, "./","");
}

template<class Csa>
static bool construct_csa(std::string file_name, Csa& csa, int_alphabet_tag) {
	std::logic_error("Not implemented yet");
	return false;
}

// TODO: this is construct_csa_from_ascii_file 
template<class Csa>
static bool construct_csa(std::string file_name, Csa& csa, tMSS& file_map, bool delete_files=true, std::string dir="./", std::string id="") {
    uint64_t fs = 0;
    char* ccc = NULL;
    if ((fs=file::read_text((const char*)file_name.c_str(), ccc))>0) {
        typename Csa::size_type n = strlen((const char*)ccc);
        if (id == "")
            id =  util::to_string(util::get_pid())+"_"+util::to_string(util::get_id()).c_str();
        if (fs != n + 1) {
			throw std::logic_error("ERROR: text contains 0-bytes.");
			return false;
        }
		std::string text_file = dir+constants::KEY_TEXT+"_"+id;
        if (!util::store_to_file(char_array_serialize_wrapper((unsigned char*)ccc,n+1), text_file.c_str())) {
            throw std::ios_base::failure("#csa_construct: Cannot store text to file system!");
            delete [] ccc;
            return false;
        } else {
            file_map[constants::KEY_TEXT] = text_file;
        }
        delete [] ccc;
        return construct_csa(csa, file_map, delete_files, dir, id);
    }
    return false;
}

template<class Csa>
static bool construct_csa(Csa& csa, tMSS& file_map, bool delete_files=true, std::string dir="./", std::string id="")
{
    write_R_output("csa", "construct CSA", "begin", 1, 0);
    typedef int_vector<>::size_type size_type;
    if (id=="")
        id =  util::to_string(util::get_pid())+"_"+util::to_string(util::get_id()).c_str();

    std::string sa_file_name = dir+constants::KEY_SA+"_"+id;
    // if sa file already exists
    std::ifstream in(sa_file_name.c_str());
    if (!in) {
        int_vector<8> text;
        util::load_vector_from_file(text, file_map[constants::KEY_TEXT].c_str(), 0 );
    	size_type n = text.size();

        typename Csa::size_type nn = strlen((const char*) text.data());
        if (nn+1 != n) {
			throw std::logic_error("ERROR: text contains 0-bytes.");
			return false;
        }

        write_R_output("csa", "construct SA", "begin", 1, 0);
        int_vector<> sa = int_vector<>(n, 0, bit_magic::l1BP(n+1)+1);
        algorithm::calculate_sa((const unsigned char*)text.data(),n, sa);	 // calculate the suffix array sa of str
		util::clear(text);

        assert(sa.size() == n);

        write_R_output("csa", "construct SA", "end", 1, 0);
        write_R_output("csa", "store SA", "begin", 1, 0);

        if (!util::store_to_file(sa, sa_file_name.c_str())) {
            throw std::ios_base::failure("#csa_construct: Cannot store SA to file system!");
            return false;
        } else {
            file_map[constants::KEY_SA] = sa_file_name;
        }
        write_R_output("csa", "store SA", "end", 1, 0);
        {
            sa.resize(0);
            int_vector<> temp;
            temp.swap(sa);
        }
    } else {
        file_map[constants::KEY_SA] = sa_file_name;
        in.close();
    }


    write_R_output("csa", "encode CSA", "begin", 1, 0);
    util::assign(csa, Csa(file_map, dir, id));
    write_R_output("csa", "encode CSA", "end", 1, 0);

    if (delete_files) {
        util::delete_all_files(file_map);
    }
    write_R_output("csa", "construct CSA", "end", 1, 0);
    return true;
}

template<class Csa>
static bool construct_csa_of_reversed_text(std::string file_name, Csa& csa)
{
    tMSS file_map;
    return construct_csa_of_reversed_text(file_name, csa, file_map, true, "./","");
}


/*
template<class Csa>
static bool construct_csa_of_reversed_text(std::string file_name, Csa& csa, tMSS& file_map, bool delete_files=true,
        std::string dir="./", std::string id="")
{
    typedef int_vector<>::size_type size_type;

    if (id=="")
        id =  util::to_string(util::get_pid())+"_"+util::to_string(util::get_id()).c_str();

    std::string tmp_rev_file_name = dir+"text_rev_"+id;
    std::ifstream in(tmp_rev_file_name.c_str());
    if (!in) {
        char* text = NULL;
        size_type n = 0;
        if ((n=file::read_text((const char*)file_name.c_str(), text)) > 0) {
            --n; // since read_text appends a 0 byte
            for (size_type i=0; i < n/2; ++i) {
                std::swap(text[i], text[n-1-i]);
            }
            file::write_text((const char*)tmp_rev_file_name.c_str(),text, n);
            file_map["text_rev"] = tmp_rev_file_name;
            delete [] text;
        } else {
            std::cout<<"ERROR: text cannot be read from file "<<file_name<<std::endl;
            return false;
        }
    } else {
        file_map["text_rev"] = tmp_rev_file_name;
        in.close();
    }
    bool res = construct_csa(tmp_rev_file_name, csa, file_map, delete_files, dir, id);
    if (delete_files) {
        util::delete_all_files(file_map);
    }
    return res;
}
*/


}// end namespace

#endif
