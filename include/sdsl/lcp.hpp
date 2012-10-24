/* sdsl - succinct data structures library
    Copyright (C) 2009 Simon Gog

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
/*! \file lcp.hpp
    \brief lcp.hpp contains classes for lcp information.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_LCP
#define INCLUDED_SDSL_LCP

#include "sdsl_concepts.hpp"

#include "int_vector.hpp"
#include "select_support_mcl.hpp"
#include "isa_construct.hpp"
#include <istream>

//! Namespace for the succinct data structure library.
namespace sdsl
{

// construct lcp arrays
template<class Lcp, class Cst, uint8_t int_width, class size_type_class, uint8_t int_width1, class size_type_class1>
void construct_lcp(Lcp& lcp, const Cst& cst,
                   int_vector_file_buffer<int_width, size_type_class>& lcp_buf,
                   int_vector_file_buffer<int_width1, size_type_class1>& isa_buf)
{
    typename Lcp::lcp_category tag;
    construct_lcp(lcp, cst, lcp_buf, isa_buf, tag);
}

template<class Lcp, class Cst, uint8_t int_width, class size_type_class, uint8_t int_width1, class size_type_class1>
void construct_lcp(Lcp& lcp, const Cst& cst,
                   int_vector_file_buffer<int_width, size_type_class>& lcp_buf,
                   int_vector_file_buffer<int_width1, size_type_class1>& isa_buf,
                   lcp_plain_tag)
{
	util::assign( lcp, Lcp(lcp_buf) );
}

template<class Lcp, class Cst, uint8_t int_width, class size_type_class, uint8_t int_width1, class size_type_class1>
void construct_lcp(Lcp& lcp, const Cst& cst,
                   int_vector_file_buffer<int_width, size_type_class>& lcp_buf,
                   int_vector_file_buffer<int_width1, size_type_class1>& isa_buf,
                   lcp_permuted_tag)
{
	util::assign( lcp, Lcp(lcp_buf, isa_buf, &(cst.csa)) );
}

template<class Lcp, class Cst, uint8_t int_width, class size_type_class, uint8_t int_width1, class size_type_class1>
void construct_lcp(Lcp& lcp, const Cst& cst,
                   int_vector_file_buffer<int_width, size_type_class>& lcp_buf,
                   int_vector_file_buffer<int_width1, size_type_class1>& isa_buf,
                   lcp_tree_compressed_tag)
{
	util::assign( lcp, Lcp( lcp_buf, &(cst.bp_support), &(cst.first_child_rank) ) );
}

// construct lcp arrays
template<class Lcp, class Cst>
void construct_lcp(Lcp& lcp, const Cst& cst, tMSS& file_map, const std::string dir, const std::string id)
{
    typename Lcp::lcp_category tag;
    construct_lcp(lcp, cst, file_map, dir, id, tag);
}

template<class Lcp, class Cst>
void construct_lcp(Lcp& lcp, const Cst& cst, tMSS& file_map, const std::string dir, const std::string id, lcp_plain_tag)
{
    int_vector_file_buffer<> lcp_buf(file_map["lcp"].c_str());
	util::assign( lcp, Lcp(lcp_buf) );
}

template<class Lcp, class Cst>
void construct_lcp(Lcp& lcp, const Cst& cst, tMSS& file_map, const std::string dir, const std::string id, lcp_permuted_tag)
{
    int_vector_file_buffer<> lcp_buf(file_map["lcp"].c_str());
    if (file_map.find("isa") == file_map.end()) {
        construct_isa(file_map, dir, id);
    }
    int_vector_file_buffer<> isa_buf(file_map["isa"].c_str());
	util::assign( lcp, Lcp(lcp_buf, isa_buf, &(cst.csa)) );
}

template<class Lcp, class Cst>
void construct_lcp(Lcp& lcp, const Cst& cst, tMSS& file_map, const std::string dir, const std::string id, lcp_tree_compressed_tag)
{
    int_vector_file_buffer<> lcp_buf(file_map["lcp"].c_str());
	util::assign( lcp, Lcp(lcp_buf, &cst) );
}

template<class Lcp, class Cst>
void construct_lcp(Lcp& lcp, const Cst& cst, tMSS& file_map, const std::string dir, const std::string id, lcp_tree_and_lf_compressed_tag)
{
    int_vector_file_buffer<> lcp_buf(file_map["lcp"].c_str());
    int_vector_file_buffer<8> bwt_buf(file_map["bwt"].c_str());
	util::assign( lcp, Lcp(lcp_buf, bwt_buf, &cst) );
}

// copy lcp arrays
template<class Lcp, class Cst>
void copy_lcp(Lcp& lcp, const Lcp& lcp_c, const Cst& cst)
{
    typename Lcp::lcp_category tag;
    copy_lcp(lcp, lcp_c, cst, tag);
}

template<class Lcp, class Cst>
void copy_lcp(Lcp& lcp, const Lcp& lcp_c, const Cst& cst, lcp_plain_tag)
{
    lcp = lcp_c;
}

template<class Lcp, class Cst>
void copy_lcp(Lcp& lcp, const Lcp& lcp_c, const Cst& cst, lcp_permuted_tag)
{
    lcp = lcp_c;
    lcp.set_csa(&(cst.csa));
}

template<class Lcp, class Cst>
void copy_lcp(Lcp& lcp, const Lcp& lcp_c, const Cst& cst, lcp_tree_compressed_tag)
{
    lcp = lcp_c;
    lcp.set_cst(&cst);
}

template<class Lcp, class Cst>
void copy_lcp(Lcp& lcp, const Lcp& lcp_c, const Cst& cst, lcp_tree_and_lf_compressed_tag)
{
    lcp = lcp_c;
    lcp.set_cst(&cst);
}

// swap lcp arrays
template<class Lcp, class Cst>
void swap_lcp(Lcp& lcp1, Lcp& lcp2, const Cst& cst1, const Cst& cst2)
{
    typename Lcp::lcp_category tag;
    swap_lcp(lcp1, lcp2, cst1, cst2, tag);
}

template<class Lcp, class Cst>
void swap_lcp(Lcp& lcp1, Lcp& lcp2, const Cst& cst1, const Cst& cst2, lcp_plain_tag)
{
    lcp1.swap(lcp2);
}

template<class Lcp, class Cst>
void swap_lcp(Lcp& lcp1, Lcp& lcp2, const Cst& cst1, const Cst& cst2, lcp_permuted_tag)
{
    lcp1.swap(lcp2);
    lcp1.set_csa(&(cst1.csa));
    lcp2.set_csa(&(cst2.csa));
}

template<class Lcp, class Cst>
void swap_lcp(Lcp& lcp1, Lcp& lcp2, const Cst& cst1, const Cst& cst2, lcp_tree_compressed_tag)
{
    lcp1.swap(lcp2);
    lcp1.set_cst(&cst1);
    lcp2.set_cst(&cst2);
}

template<class Lcp, class Cst>
void swap_lcp(Lcp& lcp1, Lcp& lcp2, const Cst& cst1, const Cst& cst2, lcp_tree_and_lf_compressed_tag)
{
    lcp1.swap(lcp2);
    lcp1.set_cst(&cst1);
    lcp2.set_cst(&cst2);
}

// load lcp arrays
template<class Lcp, class Cst>
void load_lcp(Lcp& lcp, std::istream& in, const Cst& cst)
{
    typename Lcp::lcp_category tag;
    load_lcp(lcp, in, cst, tag);
}

template<class Lcp, class Cst>
void load_lcp(Lcp& lcp, std::istream& in, const Cst& cst, lcp_plain_tag)
{
    lcp.load(in);
}

template<class Lcp, class Cst>
void load_lcp(Lcp& lcp, std::istream& in, const Cst& cst, lcp_permuted_tag)
{
    lcp.load(in, &(cst.csa));
}

template<class Lcp, class Cst>
void load_lcp(Lcp& lcp, std::istream& in, const Cst& cst, lcp_tree_compressed_tag)
{
    lcp.load(in, &cst);
}

template<class Lcp, class Cst>
void load_lcp(Lcp& lcp, std::istream& in, const Cst& cst, lcp_tree_and_lf_compressed_tag)
{
    lcp.load(in, &cst);
}

} // end namespace sdsl

#include "lcp_support_sada.hpp"     // type (b)
#include "lcp_kurtz.hpp"            // type (a)
#include "lcp_wt.hpp"               // type (a)
#include "lcp_dac.hpp"              // type (a)
#include "lcp_vlc.hpp"              // type (a)
#include "lcp_bitcompressed.hpp"    // type (a)
#include "lcp_support_tree.hpp"     // type (c)
#include "lcp_support_tree2.hpp"    // type (c)


#endif
