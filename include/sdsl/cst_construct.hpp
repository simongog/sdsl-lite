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
/*! \file cst_construct.hpp
    \brief cst_construct.hpp contains a space and time efficient construction method for compressed suffix trees (cst).
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_CST_CONSTRUCT
#define INCLUDED_SDSL_CST_CONSTRUCT

#include "int_vector.hpp"
#include "typedefs.hpp"
#include "algorithms_for_suffix_array_construction.hpp"
#include "algorithms_for_compressed_suffix_arrays.hpp"
#include "util.hpp"
#include "testutils.hpp"
#include "lcp_construct.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

namespace sdsl
{

//! Constructs a compressed suffix tree (cst) by a semi-external algorithm
/*!
 *  \param file_name	File name of the text, for which the cst should be build
 *  \param cst			A reference to the cst object.
 *
 *  \par Details
 *       The semi-external construction algorithm creates some temporary files in the
 *       current directory during the calculation.
 *
 *  \par Space complexity
 *		\f$5n\f$ bytes
 *	\par Time complexity
 *		\f$\Order{n}\f$
 */
template<class Cst>
bool construct_cst(std::string file_name, Cst& cst)
{
    tMSS file_map;
    return construct_cst(file_name, cst, file_map, true, "./", false, "", "any");
}

//! Constructs a compressed suffix tree (cst) semi-external
/*!
 *  \param file_name 		File name of the text, for which the cst should be build
 *  \param cst		 		A reference to the cst object. cst will hold the result after the execution of the method.
 *  \param file_map  		A map which will contain the file names of structures which are calculated and stored during the construction.
 *  \param delete_files		Boolean flag. If it is true, all files stored in file_map will be removed after the construction.
 *  \param dir				A directory path which points to the directory where the calculated files should be stored during the construction.
 *	\param build_only_bps   Boolean flag. If it is true, only the navigation part of the cst is build. I.e. compressed suffix array and other stuff is ommitted...
 *  \param lcp_method		Specify, which lcp construction algorithm should be used. Possible options are
 *								kasai, PHI, go, go2, any
 *
 *	\par Space complexity
 *		\f$5n\f$ bytes
 *	\par Time complexity
 *		\f$\Order{n}\f$
 */
template<class Cst>
bool construct_cst(std::string file_name, Cst& cst, tMSS& file_map, bool delete_files=true, std::string dir="./", bool build_only_bps=false, std::string id="", std::string lcp_method="any")
{
    uint64_t fs = 0;
    char* ccc = NULL;
    if ((fs=file::read_text((const char*)file_name.c_str(), ccc))>0) {
        typename Cst::size_type n = strlen((const char*)ccc);
        if (id=="")
            id =  util::to_string(util::get_pid())+"_"+util::to_string(util::get_id()).c_str();
        if (fs != n+1) {
			throw std::logic_error( std::string("ERROR: file \""+ file_name+ "\" contains 0-bytes. (cst_construct) ").c_str());
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
        return construct_cst(cst, file_map, delete_files, dir, build_only_bps, id, lcp_method);
    }
    return false;
}

template<class Cst>
static bool construct_cst(Cst& cst, tMSS& file_map, bool delete_files=true, std::string dir="./", bool build_only_bps=false, std::string id="", std::string lcp_method="any")
{
    typedef typename Cst::csa_type csa_type;
    typedef int_vector<>::size_type size_type;
    if (id == "")
        id =  util::to_string(util::get_pid())+"_"+util::to_string(util::get_id()).c_str();
    write_R_output("cst","construct CST", "begin", 1, 0);
    {
        csa_type csa;
        construct_csa(csa, file_map, false, dir, id);
		std::string csa_file = dir+util::class_to_hash(cst.csa)+"_"+id;
        if (!util::store_to_file(csa, csa_file.c_str())) {
            throw std::ios_base::failure("cst_construct: Cannot store CSA to file system!");
            return false;
        } else {
            file_map[util::class_to_hash(cst.csa)] = csa_file;
        }
    }

    {
        // test if lcp is already calculated
		std::string lcp_file = dir+constants::KEY_LCP+"_"+id;
        std::ifstream in(lcp_file.c_str());
        if (in) {
            file_map[constants::KEY_LCP] = lcp_file;
        }
    }

    if (file_map.find("lcp") == file_map.end()) {
        if (lcp_method=="any") {
            if (file_map.find("isa") == file_map.end()) { // if isa is not already on disk => use PHI algorithm for lcp
                if (!construct_lcp_PHI(file_map, dir, id))
                    return false;
            } else {
                if (!construct_lcp_kasai(file_map, dir, id))
                    return false;
            }
        } else {
            if (lcp_method=="go")
                construct_lcp_go(file_map, dir, id);
            else if (lcp_method=="go2") {
                construct_lcp_go2(file_map, dir, id);
            } else if (lcp_method=="goPHI") {
                construct_lcp_goPHI(file_map, dir, id);
            } else if (lcp_method=="PHI")
                construct_lcp_PHI(file_map, dir, id);
            else if (lcp_method=="PHIse")
                construct_lcp_PHI(file_map, dir, id, true);
            else if (lcp_method=="simple_n5")
                construct_lcp_simple_5n(file_map, dir, id);
            else if (lcp_method=="sparse_phi")
                construct_lcp_semi_extern_PHI(file_map, dir, id);
            else
                construct_lcp_kasai(file_map, dir, id);
        }
    }
	util::assign( cst, Cst(file_map, dir, id, build_only_bps) );
    if (delete_files) {
        util::delete_all_files(file_map);
    }
    write_R_output("cst","construct CST", "end", 1, 0);
    return true;
}

}// end namespace

#endif
