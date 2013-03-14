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
#include "csa_alphabet_strategy.hpp" // for key_trait
#include "select_support_mcl.hpp"
#include "construct_isa.hpp"
#include <istream>

//! Namespace for the succinct data structure library.
namespace sdsl
{

// construct lcp arrays
template<class Lcp, class Cst>
void construct_lcp(Lcp& lcp, const Cst& cst, cache_config &config) {
    typename Lcp::lcp_category tag;
    construct_lcp(lcp, cst, config, tag);
}

template<class Lcp, class Cst>
void construct_lcp(Lcp& lcp, const Cst&, cache_config &config, lcp_plain_tag) {
    int_vector_file_buffer<> lcp_buf(config.file_map[constants::KEY_LCP]);
	Lcp tmp_lcp(lcp_buf);
	lcp.swap(tmp_lcp);
}

template<class Lcp, class Cst>
void construct_lcp(Lcp& lcp, const Cst& cst, cache_config &config, lcp_permuted_tag) {
    int_vector_file_buffer<> lcp_buf(config.file_map[constants::KEY_LCP]);
	tMSS::const_iterator key = config.file_map.find(constants::KEY_ISA);
    if ( config.file_map.end() == key ) {
        construct_isa(config);
    }
    int_vector_file_buffer<> isa_buf(config.file_map[constants::KEY_ISA]);
	Lcp tmp_lcp(lcp_buf, isa_buf, &(cst.csa));
	lcp.swap(tmp_lcp);
}

template<class Lcp, class Cst>
void construct_lcp(Lcp& lcp, const Cst& cst, cache_config &config, lcp_tree_compressed_tag) {
    int_vector_file_buffer<> lcp_buf(config.file_map[constants::KEY_LCP]);
	Lcp tmp_lcp(lcp_buf, &cst);
	lcp.swap(tmp_lcp);
}

template<class Lcp, class Cst>
void construct_lcp(Lcp& lcp, const Cst& cst, cache_config &config, lcp_tree_and_lf_compressed_tag) {
    int_vector_file_buffer<> lcp_buf(config.file_map[constants::KEY_LCP]);
    int_vector_file_buffer<Cst::csa_type::alphabet_type::int_width> bwt_buf( config.file_map[key_trait<Cst::csa_type::alphabet_type::int_width>::KEY_BWT] ); 
	Lcp tmp_lcp(lcp_buf, bwt_buf, &cst);
	lcp.swap(tmp_lcp);
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
void swap_lcp(Lcp& lcp1, Lcp& lcp2, const Cst&, const Cst&, lcp_plain_tag)
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
void load_lcp(Lcp& lcp, std::istream& in, const Cst&, lcp_plain_tag)
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
#include "lcp_byte.hpp"            // type (a)
#include "lcp_wt.hpp"               // type (a)
#include "lcp_dac.hpp"              // type (a)
#include "lcp_vlc.hpp"              // type (a)
#include "lcp_bitcompressed.hpp"    // type (a)
#include "lcp_support_tree.hpp"     // type (c)
#include "lcp_support_tree2.hpp"    // type (c)


#endif
