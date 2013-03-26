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
/*! \file test_index_performance.hpp
 *  \brief test_index_performance.hpp contains a set of benchmark methods
 *  \author Simon Gog
 */
#ifndef INCLUDE_SDSL_TEST_INDEX_PERFORMANCE
#define INCLUDE_SDSL_TEST_INDEX_PERFORMANCE

#include "int_vector.hpp"	// for bit_vector and int_vector
#include "util.hpp"			// for 
#include "algorithms.hpp"	// for backward_search
#include <cstdlib>			// for rand 
#include <algorithm>		// for swap
#include <vector>			// for std::vector	
#include <iostream>

namespace sdsl
{

//! Create 1^{log_s} random intergers mod m with seed x
/*
 */
int_vector<64> get_rnd_positions(uint8_t log_s, uint64_t& mask, uint64_t m=0, uint64_t x=17);

template<class Vector>
void test_int_vector_random_access(const Vector& v, bit_vector::size_type times=100000000)
{
    typedef bit_vector::size_type size_type;
    uint64_t mask;
    int_vector<64> rands = get_rnd_positions(20, mask, v.size());
    size_type cnt=0;
    write_R_output("int_vector","random access","begin",times,cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += v[rands[ i&mask ]];
    }
    write_R_output("int_vector","random access","end",times,cnt);
}

template<class Vector>
void test_int_vector_random_write(Vector& v, bit_vector::size_type times=100000000)
{
    typedef bit_vector::size_type size_type;
    uint64_t mask;
    int_vector<64> rands = get_rnd_positions(20, mask, v.size());
    size_type cnt=0;
    write_R_output("int_vector","random write","begin",times,cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += (v[rands[ i&mask ]] = i);
    }
    write_R_output("int_vector","random write","end",times,cnt);
}

template<class Vector>
void test_int_vector_sequential_write(Vector& v, bit_vector::size_type times=100000000)
{
    typedef bit_vector::size_type size_type;
    const uint64_t mask = (1ULL << bits::hi(v.size()))-1;
//	std::cout<<" mask="<<mask<<std::endl;
    size_type cnt=0;
    write_R_output("int_vector","seq write","begin",times,cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += (v[i&mask] = i);
    }
    write_R_output("int_vector","seq write","end",times,cnt);
}

//! Test random queries on rank data structure
/*
 */
template<class Rank>
void test_rank_random_access(const Rank& rank, bit_vector::size_type times=20000000)
{
    typedef bit_vector::size_type size_type;
    uint64_t mask;
    int_vector<64> rands = get_rnd_positions(20, mask, rank.size()+1);
    size_type cnt=0;
    write_R_output("rank","random access","begin",times,cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += rank.rank(rands[ i&mask ]);
    }
    write_R_output("rank","random access","end",times,cnt);
}

//! Test creation time for a rank data structure
/*
 * param size The size of the bit vector in bits for which to rank_support should be created
 */
template<class Rank>
void test_rank_construction(bit_vector::size_type size=838860800)
{
    typedef bit_vector::size_type size_type;
    bit_vector b(size);
    util::set_random_bits(b, 17);
    std::cout<<"# bit vector size of rank construct : "<<size<<std::endl;
    write_R_output("rank","construct","begin",1,0);
    Rank rs(&b); // construct rank_support data structure
    write_R_output("rank","construct","end",1,rs(size));
}



//! Test random queries on select data structure
/*
 */
template<class Select>
void test_select_random_access(const Select& select, bit_vector::size_type times=20000000)
{
    typedef bit_vector::size_type size_type;
    const int s = 20;
    const uint64_t mask = (1<<s)-1;
    int_vector<64> rands(1<<s ,0);
    util::set_random_bits(rands, 17);
    size_type args = util::get_one_bits(*(select.v));
    util::mod(rands, args);
    for (size_type i=0; i<rands.size(); ++i)
        rands[i] = rands[i]+1;
    size_type cnt=0;
    write_R_output("select","random access","begin",times,cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += select.select(rands[ i&mask ]);
    }
    write_R_output("select","random access","end",times,cnt);
}

//! Test random queries on select data structure
/*
 */
template<class Select>
void test_select_random_access(const Select& select, bit_vector::size_type args, bit_vector::size_type times)
{
    typedef bit_vector::size_type size_type;
    const int s = 20;
    const uint64_t mask = (1<<s)-1;
    int_vector<64> rands(1<<s ,0);
    util::set_random_bits(rands, 17);
    util::mod(rands, args);
    for (size_type i=0; i<rands.size(); ++i)
        rands[i] = rands[i]+1;
    size_type cnt=0;
    write_R_output("select","random access","begin",times,cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += select.select(rands[ i&mask ]);
    }
    write_R_output("select","random access","end",times,cnt);
}

//! Test creation time for a rank data structure
/*
 * param size The size of the bit vector in bits for which to rank_support should be created
 */
template<class Select>
void test_select_construction(bit_vector::size_type size=838860800)
{
    typedef bit_vector::size_type size_type;
    bit_vector b(size);
    util::set_random_bits(b, 17);
    std::cout<<"# bit vector size of select construct : "<<size<<std::endl;
    write_R_output("select","construct","begin",1,0);
    Select sls(&b); // construct rank_support data structure
    write_R_output("select","construct","end",1, sls(1));
    std::cout<<"# size of the select_support_mcl :  "<< ((double)util::get_size_in_bytes(sls))/util::get_size_in_bytes(b) <<std::endl;
}

template<class Csa>
void test_csa_access(const Csa& csa, typename Csa::size_type times=1000000)
{
    typedef typename Csa::size_type size_type;
    size_type cnt=0;
    uint64_t mask;
    int_vector<64> rands = get_rnd_positions(20, mask, csa.size());
    write_R_output("csa","csa[]","begin",times,cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += csa[ rands[ i&mask ] ];
    }
    write_R_output("csa","csa[]","end",times,cnt);
}

template<class Csa>
void test_icsa_access(const Csa& csa, typename Csa::size_type times=1000000)
{
    typedef typename Csa::size_type size_type;
    size_type cnt=0;
    uint64_t mask;
    int_vector<64> rands = get_rnd_positions(20, mask, csa.size());
    write_R_output("csa","csa()","begin",times,cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += csa(rands[ i&mask ]);
    }
    write_R_output("csa","csa()","end",times,cnt);
}

// Test random access on \f$\Psi\f$ function
template<class Csa>
void test_psi_access(const Csa& csa, typename Csa::size_type times=1000000)
{
    typedef typename Csa::size_type size_type;
    size_type cnt=0;
    uint64_t mask;
    int_vector<64> rands = get_rnd_positions(20, mask, csa.size());
    write_R_output("csa","psi[]","begin",times,cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += csa.psi[rands[ i&mask ]];
    }
    write_R_output("csa","psi[]","end",times,cnt);
}

//! Test random access on LF function
template<class Csa>
void test_lf_access(const Csa& csa, typename Csa::size_type times=1000000)
{
    typedef typename Csa::size_type size_type;
    size_type cnt=0;
    uint64_t mask;
    int_vector<64> rands = get_rnd_positions(20, mask, csa.size());
    write_R_output("csa","psi()","begin",times,cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += csa.psi(rands[ i&mask ]);
    }
    write_R_output("csa","psi()","end",times,cnt);
}

//! Test random access on bwt
template<class Csa>
void test_bwt_access(const Csa& csa, typename Csa::size_type times=1000000)
{
    typedef typename Csa::size_type size_type;
    uint64_t mask;
    int_vector<64> rands = get_rnd_positions(20, mask, csa.size());
    size_type cnt=0;
    write_R_output("csa","bwt[]","begin",times,cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += csa.bwt[ rands[ i&mask ] ];
    }
    write_R_output("csa","bwt[]","end",times,cnt);
}

//! Test speed for pattern matching
/*!
 * \par Selection of patterns: The pattern are selected from
        random position in the original
	Do we really need file_name? We can also extract text from the csa...
 */
template<class Csa>
void test_pattern_matching(const Csa& csa,
                           const char* file_name,
                           const typename Csa::size_type pattern_len=20,
                           typename Csa::size_type times=1000000)
{
    typedef typename Csa::size_type size_type;
    uint64_t mask;
    int_vector<64> rands = get_rnd_positions(15, mask, csa.size()-pattern_len);
    unsigned char* patterns = new unsigned char[rands.size()*pattern_len + 2];
    write_R_output("csa","extract patterns","begin",times,0);

    size_type file_size = 0;
    {
        int_vector<8> text;
        if (util::load_vector_from_file(text, file_name, 1)) {
            for (size_type i=0; i<rands.size(); ++i) {
                for (size_type j=rands[i]; j < rands[i] + pattern_len; ++j)
                    patterns[i*pattern_len + (j-rands[i])] = text[j];
            }
        } else {
            std::cerr << "ERROR: test_pattern_matching: could not open";
            std::cerr << "file `" << file_name << "`" << std::endl;
        }
    }
    write_R_output("csa","extract patterns","end",times,0);

    size_type cnt=0;
    write_R_output("csa","pattern_matching","begin",times,cnt);
    for (size_type i=0, l_res=0, r_res=0; i<times; ++i) {
        cnt += algorithm::backward_search(csa, 0, csa.size(),
                                          (patterns + ((i&mask)*pattern_len)),
                                          pattern_len, l_res, r_res);
    }
    write_R_output("csa","pattern_matching","end",times,cnt);
    delete [] patterns;
}


template<class Csa>
void test_rank_bwt_access(const Csa& csa, typename Csa::size_type times=1000000)
{
    typedef typename Csa::size_type size_type;
    size_type cnt=0;
    // precalc queries, distribution of chars is important!
    uint64_t mask, s=20;
    int_vector<64> rands = get_rnd_positions(s, mask, csa.size()+1);
    unsigned char* c_rands = new unsigned char[1<<s];
    for (size_type i=0; i<rands.size(); ++i)
        c_rands[i] = csa.bwt[ rands[i] ];

    write_R_output("csa","rank_bwt","begin",times,cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += csa.rank_bwt(rands[ i&mask ], c_rands[ i&mask ]);
    }
    write_R_output("csa","rank_bwt","end",times,cnt);
    delete [] c_rands;
}

template<class Csa>
void test_select_bwt_access(const Csa& csa, typename Csa::size_type times=500000)
{
    typedef typename Csa::size_type size_type;
    size_type cnt=0;

    // precalc queries, distribution of chars is important!
    uint64_t mask, s=20;
    int_vector<64> rands = get_rnd_positions(s, mask, csa.size());
    unsigned char* c_rands = new unsigned char[1<<s];
    for (size_type i=0; i<rands.size(); ++i) {
        c_rands[i] = csa.bwt[ rands[i] ];
        rands[i] = csa.rank_bwt(rands[i]+1, c_rands[i]);
    }
    write_R_output("csa","select_bwt","begin",times,cnt);
    for (size_type i=0; i<times; ++i) {
//		size_type c = i % csa.sigma;
//		cnt += csa.select_bwt( (i % (csa.C[c+1]-csa.C[c])) + 1, csa.comp2char[c]);
        cnt += csa.select_bwt(rands[ i&mask ], c_rands[ i&mask]);
    }
    write_R_output("csa","select_bwt","end",times,cnt);
}

//! Test performance of the depth first iterator of a CST
/*!
 * \param cst The CST that should be tested.
 * \param times The number of times the depth first iterator should be incremented.
 * \par Details
 *      This method increments the depth first search iterator of
 *      a CST n times.
 */
template<class Cst>
void test_cst_dfs_iterator(Cst& cst, typename Cst::size_type times=100000)
{
    if (times > cst.nodes())
        times = cst.nodes();
    typedef typename Cst::size_type size_type;
    size_type cnt=0;
    {
        // calc values for cnt
        typename Cst::const_iterator it = cst.begin();
        for (size_type i=0; i < std::min(times,(size_type)1000); ++i, ++it) {
            cnt += cst.depth(*it);
        }
    }
    write_R_output("cst","dfs","begin",times,cnt);
    typename Cst::const_iterator it = cst.begin();
    for (size_type i=0; i<times; ++i) {
        ++it;
    }
    write_R_output("cst", "dfs", "end", times, cnt + cst.depth(*it));
}

//! Test performance of the depth first iterator and the LCP array of a CST
/*!
 * \param cst The CST that should be tested.
 * \param times The number of times the depth first iterator should be incremented.
 * \par Details
 *      This method increments the depth first search iterator of
 *      a CST n times and calculates the depth for each visited node.
 */
template<class Cst>
void test_cst_dfs_iterator_and_depth(Cst& cst, typename Cst::size_type times=1000000, bool output=false)
{
    if (times > 2*cst.nodes()-cst.size())
        times = 2*cst.nodes()-cst.size();
    typedef typename Cst::size_type size_type;
    size_type cnt=0;
    write_R_output("cst","dfs and depth","begin",times,cnt);
    typename Cst::const_iterator it = cst.begin();
    if (!output) {
        for (size_type i=0; i<times; ++i, ++it) {
            if (!cst.is_leaf(*it))
                cnt += cst.depth(*it);
        }
    } else {
        for (size_type i=0; i<times; ++i, ++it) {
            if (!cst.is_leaf(*it)) {
                size_type d = cst.depth(*it);
                std::cerr << d << "-[" << cst.lb(*it) << "," << cst.rb(*it) << "] ";
                if (d < 60) {
                    for (int i=1; i<=d; ++i)
                        std::cerr<< cst.edge(*it, i);
                }
                std::cerr << std::endl;
                cnt += d;
            }
        }
    }
    write_R_output("cst","dfs and depth","end",times,cnt);
}

//! Test performance of the depth first iterator and the id method of the CST
/*!
 * \param cst The CST that should be tested.
 * \param times The number of times the depth first iterator should be incremented.
 * \par Details
 *      This method increments the depth first search iterator of
 *      a CST n times and calculates the function id for each visited node.
 */
template<class Cst>
void test_cst_dfs_iterator_and_id(Cst& cst, typename Cst::size_type times=1000000, bool output=false)
{
    if (times > 2*cst.nodes()-cst.size())
        times = 2*cst.nodes()-cst.size();
    typedef typename Cst::size_type size_type;
    size_type cnt=0;
    write_R_output("cst","dfs and id","begin",times,cnt);
    typename Cst::const_iterator it = cst.begin();
    if (!output) {
        for (size_type i=0; i<times; ++i, ++it) {
            cnt += cst.id(*it);
        }
    } else {
        for (size_type i=0; i<times; ++i, ++it) {
            size_type id = cst.id(*it);
            std::cerr << id << std::endl;
            cnt += id;
        }
    }
    write_R_output("cst","dfs and id","end",times,cnt);
}



// TODO: bottom up iterator

//! Make random accesse to an LCP array
template<class Lcp>
void test_lcp_random_access(Lcp& lcp, typename Lcp::size_type times=10000000)
{
    typedef typename Lcp::size_type size_type;
    size_type n = lcp.size();
    if (times > n)
        times = n;
    uint64_t mask;
    int_vector<64> rands = get_rnd_positions(20, mask, n);
    size_type cnt=0;
    write_R_output("lcp","random access","begin",times,cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += lcp[ rands[ i&mask ] ];
    }
    write_R_output("lcp","random access","end",times,cnt);
}

//! Make sequential accesses to an LCP array
/*!
 *  The access time for lcp sequential access is not important in practice,
 *  as sequential access should be realized by streaming the uncompressed
 *  lcp array from disk.
 */
template<class Lcp>
void test_lcp_sequential_access(Lcp& lcp, typename Lcp::size_type times=10000000)
{
    typedef typename Lcp::size_type size_type;
    size_type n = lcp.size();
    if (times > n)
        times = n;
    size_type cnt=0;
    write_R_output("lcp","sequential access","begin",times,cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += lcp[ i ];
    }
    write_R_output("lcp","sequential access","end",times,cnt);
}


//! Make random sequential accesse to an LCP array
/*!
 *  \param lcp   	The lcp data structure
 *  \param times 	The number of random accesses
 *  \param seq_len 	The number of sequential accesses after each random access
 *  This test was described in the master thesis of Rodrigo Canovas.
 */
template<class Lcp>
void test_lcp_random_sequential_access(Lcp& lcp, typename Lcp::size_type times=1000000, typename Lcp::size_type seq_len=64)
{
    typedef typename Lcp::size_type size_type;
    size_type n = lcp.size();
    if (times > n)
        times = n;
    if (seq_len >= n)
        seq_len = n-1;
    const int s = 20;
    const uint64_t mask = (1<<s)-1;
    int_vector<64> rands(1<<s ,0);
    util::set_random_bits(rands, 17);
    for (int_vector<64>::size_type i=0; i < rands.size(); ++i) {
        rands[i] = rands[i] % n;
        if (rands[i] + seq_len >= n)
            rands[i] = n - 1 - seq_len;
    }
    size_type cnt=0;
    write_R_output("lcp","random sequential access","begin",times*(seq_len+1),cnt);
    for (size_type i=0; i<times; ++i) {
        for (size_type j=rands[i&mask], k=0; k<=seq_len; ++k, ++j) {
            cnt += lcp[ j ];
        }
    }
    write_R_output("lcp","random sequential access","end",times*(seq_len+1),cnt);
}

//! Test the speed of the parent operation
/*!
 * \param Cst The compressed suffix tree
 * \param times Number of times a traversal from a random leaf to the root\
                is started to collect nodes for which the parent operation is performed.
 */
template<class Cst>
void test_cst_parent_operation(const Cst& cst, typename Cst::size_type times=100000, uint64_t x=17)
{
    typedef typename Cst::size_type size_type;
    typedef typename Cst::node_type node_type;

    srand(x);
    size_type n = cst.csa.size();
    // take \f$ time \f$ random leaves
    std::vector<node_type> rand_leaf(times);
    for (size_type i=0; i<rand_leaf.size(); ++i) {
        rand_leaf[i] = cst.select_leaf(1+ (rand() % n));
    }

    node_type p;
    size_type cnt=0;
    write_R_output("cst","parent","begin",times,cnt);
    for (size_type i=0; i<times; ++i, ++cnt) {
        p = cst.parent(rand_leaf[i]);
        while (p != cst.root()) {
            p = cst.parent(p);
            ++cnt;
        }
    }
    write_R_output("cst","parent","end",times,cnt);
}


//! Generate nodes of a cst by applying the child operation to each of \f$times\f$ random leaves until we get to the root
/*!
 * \param Cst The compressed suffix
 * \param times Numer of random leaves
 * \param nodes Reference to a vector which will contain the generated nodes
 * \param x     Seed for the random number generator for the generation of the leaves
 */
template<class Cst>
void generate_nodes_from_random_leaves(const Cst& cst, typename Cst::size_type times, std::vector<typename Cst::node_type>& nodes, uint64_t x=17)
{
    typedef typename Cst::size_type size_type;
    typedef typename Cst::node_type node_type;
    srand(x);
    size_type n = cst.csa.size();
    // generate nodes
    for (size_type i=0; i<times; ++i) {
        node_type p = cst.select_leaf(1+ (rand() % n));
        nodes.push_back(p);
        while (p != cst.root()) {
            p = cst.parent(p);
            nodes.push_back(p);
        }
    }
}

//! Test the speed of the child operation
/*!
 * \param Cst The compressed suffix tree
 * \param times Number of times a traversal from a random leaf to the root\
                is started to collect nodes for which the child operation is performed.
 */
template<class Cst>
void test_cst_child_operation(const Cst& cst, typename Cst::size_type times=5000, uint64_t x=17)
{
    typedef typename Cst::size_type size_type;
    typedef typename Cst::node_type node_type;

    std::vector<node_type> nodes;
    generate_nodes_from_random_leaves(cst, times, nodes, x);
//	for(size_type i=0; i<20; ++i){
//		std::cout<< cst.lb(nodes[i])<<" "<<cst.rb(nodes[i])<<std::endl;
//	}
    // choose some chars for the text
    unsigned char* letters = new unsigned char[nodes.size()+1];
    for (size_type i=0; i<nodes.size(); ++i) {
        letters[i] = cst.csa.bwt[i];
    }

    node_type c;  // for child node
    size_type char_pos=0;
    size_type cnt=0;
    write_R_output("cst","child","begin",nodes.size(),cnt);
    for (size_type i=0; i<nodes.size(); ++i) {
//		if(i<20){
//			std::cout<<"i="<<i<<" vl="<<cst.lb(nodes[i])<<" rb="<<cst.rb(nodes[i])<<std::endl;
//			std::cout<<cst.csa[cst.lb(nodes[i])]<<" "<<cst.depth(nodes[i])<<std::endl;
//		}
        c = cst.child(nodes[i], letters[i], char_pos);
        if (c==cst.root())
            ++cnt;
    }
    write_R_output("cst","child","end",nodes.size(),cnt);
    delete [] letters;
}

//! Test the speed of the 1th_child operation
/*!
 * \param Cst The compressed suffix tree
 * \param times Number of times a traversal from a random leaf to the root\
                is started to collect nodes for which the child operation is performed.
 */
template<class Cst>
void test_cst_1th_child_operation(const Cst& cst, typename Cst::size_type times=1000000, uint64_t x=17)
{
    typedef typename Cst::size_type size_type;
    typedef typename Cst::node_type node_type;

    std::vector<node_type> nodes;
    generate_nodes_from_random_leaves(cst, times, nodes, x);

    node_type c;  // for 1th_child node
    size_type cnt=0;
    write_R_output("cst","1th_child","begin",nodes.size(),cnt);
    for (size_type i=0; i<nodes.size(); ++i) {
        c = cst.select_child(nodes[i], 1);
        if (c==cst.root())
            ++cnt;
    }
    write_R_output("cst","1th_child","end",nodes.size(),cnt);
}

//! Test the speed of the sibling operation
/*!
 * \param Cst The compressed suffix tree
 * \param times Number of times a traversal from a random leaf to the root\
                is started to collect nodes for which the child operation is performed.
 */
template<class Cst>
void test_cst_sibling_operation(const Cst& cst, typename Cst::size_type times=100000, uint64_t x=17)
{
    typedef typename Cst::size_type size_type;
    typedef typename Cst::node_type node_type;

    std::vector<node_type> nodes;
    generate_nodes_from_random_leaves(cst, times, nodes, x);
    for (size_type i=0; i<nodes.size(); ++i) {
        nodes[i] = cst.sibling(nodes[i]);
    }

    node_type c;  // for sibling node
    size_type cnt=0;
    write_R_output("cst","sibling","begin",nodes.size(),cnt);
    for (size_type i=0; i<nodes.size(); ++i) {
        c = cst.sibling(nodes[i]);
        if (c==cst.root())
            ++cnt;
    }
    write_R_output("cst","sibling","end",nodes.size(),cnt);
}

//! Test id operation
template<class Cst>
void test_cst_id_operation(const Cst& cst, typename Cst::size_type times=100000, uint64_t x=17)
{
    typedef typename Cst::size_type size_type;
    typedef typename Cst::node_type node_type;
    std::vector<node_type> nodes;
    generate_nodes_from_random_leaves(cst, times, nodes, x);

    size_type cnt = 0;
    write_R_output("cst","id","begin",nodes.size(),cnt);
    for (size_type i=0; i < nodes.size(); ++i) {
        cnt += cst.id(nodes[i]);
    }
    write_R_output("cst","id","end",nodes.size(),cnt);
}

//! Test depth operations for leaves and inner nodes
template<class Cst>
void test_cst_depth_operation(const Cst& cst, typename Cst::size_type times=100000, uint64_t x=17)
{
    typedef typename Cst::size_type size_type;
    typedef typename Cst::node_type node_type;
    std::vector<node_type> nodes;
    generate_nodes_from_random_leaves(cst, times, nodes, x);

    size_type cnt = 0;
    write_R_output("cst","depth","begin",nodes.size(),cnt);
    for (size_type i=0; i < nodes.size(); ++i) {
        cnt += cst.depth(nodes[i]);
    }
    write_R_output("cst","depth","end",nodes.size(),cnt);
}


//! Test depth operations for inner nodes
template<class Cst>
void test_cst_depth_operation_for_inner_nodes(const Cst& cst, typename Cst::size_type times=100000, uint64_t x=17)
{
    typedef typename Cst::size_type size_type;
    typedef typename Cst::node_type node_type;
    std::vector<node_type> nodes;
    {
        std::vector<node_type> nodes2;
        generate_nodes_from_random_leaves(cst, times, nodes2, x);
        for (size_type i=0; i<nodes2.size(); ++i)
            if (!cst.is_leaf(nodes2[i])) {
                nodes.push_back(nodes2[i]);
            }
    }
    size_type cnt = 0;
    write_R_output("cst","depth of inner nodes","begin",nodes.size(),cnt);
    for (size_type i=0; i < nodes.size(); ++i) {
        cnt += cst.depth(nodes[i]);
    }
    write_R_output("cst","depth of inner nodes","end",nodes.size(),cnt);
}

//! Test lca operation
template<class Cst>
void test_cst_lca_operation(const Cst& cst, typename Cst::size_type times=1000000, uint64_t x=17)
{
    typedef typename Cst::size_type size_type;
    typedef typename Cst::node_type node_type;
    // 	generate \f$2^{19}\f$ random pairs of leafs
    size_type n = cst.csa.size();
    uint64_t mask = (1<<20)-1;
    std::vector<node_type> nodes(1<<20);
    srand(x);
    for (size_type i=0; i < nodes.size(); ++i) {
        nodes[i] = cst.select_leaf(rand()%n + 1);
    }

    size_type cnt=0;
    write_R_output("cst","lca","begin",times,cnt);
    for (size_type i=0; i<times; ++i) {
        node_type v = cst.lca(nodes[(2*i) & mask], nodes[(2*i+1) & mask]);
        if (v == cst.root())
            cnt++;
//		if(i<30)
//			std::cout<<"lca("<<cst.lb(nodes[(2*i)&mask])<<","<<cst.lb(nodes[(2*i+1)&mask])<<")=("<<cst.lb(v)<<","<<cst.rb(v)<<")"<<std::endl;
    }
    write_R_output("cst","lca","end",times,cnt);
}

//! Test suffix link operation
template<class Cst>
void test_cst_sl_operation(const Cst& cst, typename Cst::size_type times=500, uint64_t x=17)
{
    typedef typename Cst::size_type size_type;
    typedef typename Cst::node_type node_type;
    size_type n = cst.csa.size();
    if (times > n)
        times = n;

    std::vector<node_type> nodes(times);
    srand(x);
    // take \f$ times \f$ random leaves and calculate each parent
    for (size_type i=0; i<times; ++i) {
        nodes[i] = cst.parent(cst.select_leaf(rand()%n + 1));
    }

    size_type cnt=0;
    times = 0;
    write_R_output("cst","sl","begin",0,cnt);
    for (size_type i=0; i<nodes.size(); ++i) {
        node_type v = nodes[i];
//		std::cout<<"v="<<cst.lb(v)<<" "<<cst.rb(v)<<std::endl;
//		size_type d = cst.depth(v);
        while (v != cst.root()) { // while v is not the root
            ++cnt;
            v = cst.sl(v); // follow suffix link
//			if( cnt < 30 ){
//				std::cout<< cnt << " " << cst.lb(v) << " " << cst.rb(v) << " " << cst.depth(v) << std::endl;
//			}
//			size_type d2 = cst.depth(v);
//			if( d != d2+1 ){
//				std::cout<<"error at cnt "<<cnt<<" d="<<d<<" d2="<<d2<<std::endl;
//			}
//			d = d2;
        }
    }
    write_R_output("cst","sl","end",cnt,cnt);
}

//! Test matching statistics
/*! \param cst	Compressed suffix tree of sequence S1 of length n1
 *  \param S2	Pointer to the unsigned char array S2.
 *  \param n2	The length of S2.
 */
template<class Cst>
void test_cst_matching_statistics(const Cst& cst, unsigned char* S2, typename Cst::size_type n2)
{
    typedef typename Cst::size_type size_type;
    typedef typename Cst::node_type node_type;

    size_type cnt = 0;
    write_R_output("cst","mstats","begin",n2,cnt);
    size_type q  = 0;						// current match length
    size_type p2 = n2-1;              // position in S2
    size_type i  = 0, j = cst.csa.size()-1; // \f$ \epsilon \f$ matches all suffixes of S1
    while (p2+1 > 0) {
        size_type lb, rb;
        // perform backward search on interval \f$ [i,j] \f$
        size_type size = algorithm::backward_search(cst.csa, i, j, S2[p2], lb, rb);
        if (size > 0) {
            q = q + 1;
            i = lb; j = rb;
            p2 = p2 - 1;
        } else if (i==0 and j == cst.csa.size()) {
            p2 = p2 -1;
        } else {
            // map interval to a node of the cst and calculate parent
            node_type p = cst.parent(cst.node(i, j));
            q = cst.depth(p);	// update match length
            i = cst.lb(p); 		// update left bound
            j = cst.rb(p);		// update right bound
        }
        cnt += q;
    }
    write_R_output("cst","mstats","end",n2,cnt);
}

// test the speed of find_close at random opening parentheses
template<class Bps>
void test_bps_find_close_and_enclose(const Bps& bps, const bit_vector& b, uint64_t times=10000000, uint64_t x=17)
{
    typedef bit_vector::size_type size_type;
    uint64_t mask;
//	uint64_t n = bps.size();
    int_vector<64> rands = get_rnd_positions(20, mask, bps.size());
    for (size_type i=0; i<rands.size(); ++i) {
        if (!b[rands[i]]) { // if there is no opening parentheses at position rands[i]
            rands[i] = bps.find_open(rands[i]);
        }
    }
    uint64_t cnt = 0;
    write_R_output("bps","find_close","begin",times, cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += bps.find_close(rands[i&mask]);
    }
    write_R_output("bps","find_close","end",times, cnt);

    cnt = 0;
    write_R_output("bps","enclose rand","begin",times, cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += bps.enclose(rands[i&mask]);
    }
    write_R_output("bps","enclose rand","end",times, cnt);

    cnt = 0;
    size_type cnt2=0;
    write_R_output("bps","enclose tree","begin",times, cnt);
    for (size_type i=0,enc, size=bps.size(); cnt2 < times; ++i) {
        enc = rands[i&mask];
        while ((enc = bps.enclose(enc)) != size) {
            cnt += enc;
            ++cnt2;
        }
        ++cnt2;
    }
    write_R_output("bps","enclose tree","end",times, cnt);
}

// test the speed of find_close at random opening parentheses
template<class Bps>
void test_bps_find_open(const Bps& bps, const bit_vector& b, uint64_t times=10000000, uint64_t x=17)
{
    typedef bit_vector::size_type size_type;
    uint64_t mask;
    uint64_t n = bps.size();
    int_vector<64> rands = get_rnd_positions(20, mask, n);
    for (size_type i=0; i<rands.size(); ++i) {
        if (b[rands[i]]) {
            rands[i] = bps.find_close(rands[i]);
        }
    }
    uint64_t cnt = 0;
    write_R_output("bps","find_open","begin",times, cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += bps.find_open(rands[i&mask]);
    }
    write_R_output("bps","find_open","end",times, cnt);
}

// test the speed of the double enclose method for random opening parentheses i and i+1
template<class Bps>
void test_bps_double_enclose(const Bps& bps, const bit_vector& b, uint64_t times=10000000, uint64_t x=17)
{
    typedef bit_vector::size_type size_type;
    uint64_t mask;
    uint64_t n = bps.size();
    int_vector<64> rands = get_rnd_positions(20, mask, bps.size());
    for (size_type i=0; i<rands.size()/2; ++i) {
        if (!b[rands[2*i]]) { // if there is no opening parentheses at position rands[i]
            uint64_t pos = (rands[2*i]+1)%n;
            while (!b[pos] and pos != rands[2*i])  // go forward until we get an opening one
                pos = (pos+1) % n;
            if (pos + 1000 > rands.size())
                pos = 0;
            rands[2*i] = pos;
        }
        {
            uint64_t pos = (rands[2*i]+1)%n;
            while (!b[pos] and pos != rands[2*i])  // go forward until we get the next opening one
                pos = (pos+1) % n;
            rands[2*i+1] = pos;
        }
    }
    uint64_t cnt = 0;
    write_R_output("bps","double_enclose","begin",times, cnt);
    for (size_type i=0; i<times; ++i) {
        cnt += bps.double_enclose(rands[(i*2)&mask], rands[(i*2+1)&mask]);
    }
    write_R_output("bps","double_enclose","end",times, cnt);
}


}// end namespace sdsl

#endif
