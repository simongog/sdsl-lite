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
template<class t_lcp, class t_cst>
void construct_lcp(t_lcp& lcp, const t_cst& cst, cache_config& config)
{
    typename t_lcp::lcp_category tag;
    construct_lcp(lcp, cst, config, tag);
}

template<class t_lcp, class t_cst>
void construct_lcp(t_lcp& lcp, const t_cst&, cache_config& config, lcp_plain_tag)
{
    t_lcp tmp_lcp(config);
    lcp.swap(tmp_lcp);
}

template<class t_lcp, class t_cst>
void construct_lcp(t_lcp& lcp, const t_cst& cst, cache_config& config, lcp_permuted_tag)
{
    t_lcp tmp_lcp(config, &(cst.csa));
    lcp.swap(tmp_lcp);
}

template<class t_lcp, class t_cst>
void construct_lcp(t_lcp& lcp, const t_cst& cst, cache_config& config, lcp_tree_compressed_tag)
{
    t_lcp tmp_lcp(config, &cst);
    lcp.swap(tmp_lcp);
}

template<class t_lcp, class t_cst>
void construct_lcp(t_lcp& lcp, const t_cst& cst, cache_config& config, lcp_tree_and_lf_compressed_tag)
{
    t_lcp tmp_lcp(config, &cst);
    lcp.swap(tmp_lcp);
}

// copy lcp arrays
template<class t_lcp, class t_cst>
void copy_lcp(t_lcp& lcp, const t_lcp& lcp_c, const t_cst& cst)
{
    typename t_lcp::lcp_category tag;
    copy_lcp(lcp, lcp_c, cst, tag);
}

template<class t_lcp, class t_cst>
void copy_lcp(t_lcp& lcp, const t_lcp& lcp_c, const t_cst&, lcp_plain_tag)
{
    lcp = lcp_c;
}

template<class t_lcp, class t_cst>
void copy_lcp(t_lcp& lcp, const t_lcp& lcp_c, const t_cst& cst, lcp_permuted_tag)
{
    lcp = lcp_c;
    lcp.set_csa(&(cst.csa));
}

template<class t_lcp, class t_cst>
void copy_lcp(t_lcp& lcp, const t_lcp& lcp_c, const t_cst& cst, lcp_tree_compressed_tag)
{
    lcp = lcp_c;
    lcp.set_cst(&cst);
}

template<class t_lcp, class t_cst>
void copy_lcp(t_lcp& lcp, const t_lcp& lcp_c, const t_cst& cst, lcp_tree_and_lf_compressed_tag)
{
    lcp = lcp_c;
    lcp.set_cst(&cst);
}

// move lcp arrays
template<class t_lcp, class t_cst>
void move_lcp(t_lcp& lcp, t_lcp& lcp_c, const t_cst& cst)
{
    typename t_lcp::lcp_category tag;
    move_lcp(lcp, lcp_c, cst, tag);
}

template<class t_lcp, class t_cst>
void move_lcp(t_lcp& lcp, t_lcp& lcp_c, const t_cst&, lcp_plain_tag)
{
    lcp = std::move(lcp_c);
}

template<class t_lcp, class t_cst>
void move_lcp(t_lcp& lcp, t_lcp& lcp_c, const t_cst& cst, lcp_permuted_tag)
{
    lcp = std::move(lcp_c);
    lcp.set_csa(&(cst.csa));
}

template<class t_lcp, class t_cst>
void move_lcp(t_lcp& lcp, t_lcp& lcp_c, const t_cst& cst, lcp_tree_compressed_tag)
{
    lcp = std::move(lcp_c);
    lcp.set_cst(&cst);
}

template<class t_lcp, class t_cst>
void move_lcp(t_lcp& lcp, t_lcp& lcp_c, const t_cst& cst, lcp_tree_and_lf_compressed_tag)
{
    lcp = std::move(lcp_c);
    lcp.set_cst(&cst);
}


// swap lcp arrays
template<class t_lcp, class t_cst>
void swap_lcp(t_lcp& lcp1, t_lcp& lcp2, const t_cst& cst1, const t_cst& cst2)
{
    typename t_lcp::lcp_category tag;
    swap_lcp(lcp1, lcp2, cst1, cst2, tag);
}

template<class t_lcp, class t_cst>
void swap_lcp(t_lcp& lcp1, t_lcp& lcp2, const t_cst&, const t_cst&, lcp_plain_tag)
{
    lcp1.swap(lcp2);
}

template<class t_lcp, class t_cst>
void swap_lcp(t_lcp& lcp1, t_lcp& lcp2, const t_cst& cst1, const t_cst& cst2, lcp_permuted_tag)
{
    lcp1.swap(lcp2);
    lcp1.set_csa(&(cst1.csa));
    lcp2.set_csa(&(cst2.csa));
}

template<class t_lcp, class t_cst>
void swap_lcp(t_lcp& lcp1, t_lcp& lcp2, const t_cst& cst1, const t_cst& cst2, lcp_tree_compressed_tag)
{
    lcp1.swap(lcp2);
    lcp1.set_cst(&cst1);
    lcp2.set_cst(&cst2);
}

template<class t_lcp, class t_cst>
void swap_lcp(t_lcp& lcp1, t_lcp& lcp2, const t_cst& cst1, const t_cst& cst2, lcp_tree_and_lf_compressed_tag)
{
    lcp1.swap(lcp2);
    lcp1.set_cst(&cst1);
    lcp2.set_cst(&cst2);
}

// load lcp arrays
template<class t_lcp, class t_cst>
void load_lcp(t_lcp& lcp, std::istream& in, const t_cst& cst)
{
    typename t_lcp::lcp_category tag;
    load_lcp(lcp, in, cst, tag);
}

template<class t_lcp, class t_cst>
void load_lcp(t_lcp& lcp, std::istream& in, const t_cst&, lcp_plain_tag)
{
    lcp.load(in);
}

template<class t_lcp, class t_cst>
void load_lcp(t_lcp& lcp, std::istream& in, const t_cst& cst, lcp_permuted_tag)
{
    lcp.load(in, &(cst.csa));
}

template<class t_lcp, class t_cst>
void load_lcp(t_lcp& lcp, std::istream& in, const t_cst& cst, lcp_tree_compressed_tag)
{
    lcp.load(in, &cst);
}

template<class t_lcp, class t_cst>
void load_lcp(t_lcp& lcp, std::istream& in, const t_cst& cst, lcp_tree_and_lf_compressed_tag)
{
    lcp.load(in, &cst);
}

} // end namespace sdsl

#include "lcp_support_sada.hpp"     // type (b)
#include "lcp_byte.hpp"             // type (a)
#include "lcp_wt.hpp"               // type (a)
#include "lcp_vlc.hpp"              // type (a)
#include "lcp_dac.hpp"              // type (a)
#include "lcp_bitcompressed.hpp"    // type (a)
#include "lcp_support_tree.hpp"     // type (c)
#include "lcp_support_tree2.hpp"    // type (c)


#endif
