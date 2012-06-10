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
/*! \file wt.hpp
    \brief wt.hpp contains a generic class for the wavelet tree proposed first by Grossi et al. 2003 and applied to the BWT in Foschini et al. 2004.
	\author Simon Gog

	\note Bug fix: m_first_symbol was not serialized 2010-07-14
*/
#ifndef INCLUDED_SDSL_WT
#define INCLUDED_SDSL_WT

#include "int_vector.hpp"
#include "bitmagic.hpp"
#include "util.hpp" // for util::ssign
#include <set> // for calculating the alphabet size
#include <map> // for mapping a symbol to its lexicographical index
#include <algorithm> // for std::swap
#include <stdexcept>
#include <vector>

#ifdef SDSL_DEBUG
#define SDSL_DEBUG_WT
#endif


//#define SDSL_DEBUG_WT

#ifdef SDSL_DEBUG_WT
#include "testutils.hpp"
#endif

//! Namespace for the succinct data structure library.
namespace sdsl
{

struct unsigned_char_map {
    unsigned char m_map[256];

    unsigned char& operator[](unsigned char i) {
        return *(m_map+i);
    }

    unsigned char operator[](unsigned char i)const {
        return m_map[i];
    }

    void clear() {
        for (uint16_t i=0; i<256; ++i)
            m_map[i] = 0;
    }

    uint16_t serialize(std::ostream& out)const {
        uint16_t written_bytes = 0;
        for (uint16_t i=0; i<256; ++i) {
            out.write((char*)&m_map[i], sizeof(m_map[i]));
            written_bytes += sizeof(m_map[256]);
        }
        return written_bytes;
    }

    void load(std::istream& in) {
        for (uint16_t i=0; i<256; ++i) {
            in.read((char*) &m_map[i], sizeof(m_map[i]));
        }
    }

    void swap(unsigned_char_map& map) {
        if (this != &map) {
            for (uint16_t i=0; i<256; ++i) {
                std::swap(m_map[i], map.m_map[i]);
            }
        }
    }
};


template<class RandomAccessContainer>
class wt_trait
{
    public:
        typedef typename RandomAccessContainer::size_type 		size_type;
        typedef typename RandomAccessContainer::value_type 		value_type;
        typedef RandomAccessContainer& 							reference_type;
        typedef std::map<value_type, size_type>					map_type;
        typedef std::map<value_type, size_type>					inv_map_type;
        enum { char_node_map_size=0 };

        static size_type alphabet_size_and_map(const reference_type rac, size_type n, map_type& map, inv_map_type& inv_map, value_type& first_symbol) {
            if (n > 0)
                first_symbol = rac[0];
            map.clear();
            inv_map.clear();
            size_type alphabet_size = 0;
            for (size_type i=0; i<n; ++i) {
                if (map.find(rac[i]) == map.end()) {
                    map[rac[i]]	= 1;
                }
            }
            for (typename map_type::iterator it = map.begin(); it != map.end(); ++it) { // this preserves the order
                it->second = alphabet_size;
                inv_map[alphabet_size] = it->first;
                ++alphabet_size;
            }
            return alphabet_size;
        }

        static bool symbol_available(const map_type& map, const value_type c, const value_type first_symbol, const size_type) {
            return map.find(c)!=map.end();
        }

        static size_type serialize_maps(std::ostream& out, const map_type& map, const inv_map_type& inv_map) {
            throw std::logic_error(util::demangle(typeid(wt_trait<RandomAccessContainer>).name())+": serialize not implemented");
            return 0;
        }

        static size_type load_maps(std::istream& in, map_type& map, inv_map_type& inv_map) {
            throw std::logic_error(util::demangle(typeid(wt_trait<RandomAccessContainer>).name())+": load not implemented");
            return 0;
        }
};

template<class character>
class wt_trait<character*>
{
    public:
        typedef bit_vector::size_type 			size_type;
        typedef character						value_type;
        typedef character*						reference_type;
        typedef std::map<value_type, size_type>	map_type;
        typedef std::map<value_type, size_type>	inv_map_type;
        enum { char_node_map_size=256 };

        static size_type alphabet_size_and_map(const reference_type rac, size_type n, map_type& map, inv_map_type& inv_map, value_type& first_symbol) {
            if (n > 0)
                first_symbol = *(rac+0);
            map.clear();
            inv_map.clear();
            size_type alphabet_size = 0;
            for (size_type i=0; i<n; ++i) {
                if (map.find(*(rac+i)) == map.end()) {
                    map[*(rac+i)]	= 1;
                }
            }

            for (typename map_type::iterator it = map.begin(); it != map.end(); ++it) { // this preserves the order
                it->second = alphabet_size;
                inv_map[alphabet_size] = it->first;
                ++alphabet_size;
            }
            return alphabet_size;
        }

        static bool symbol_available(const map_type& map, const value_type c, const value_type first_symbol, const size_type) {
            return map.find(c)!=map.end();
        }

        static size_type serialize_maps(std::ostream& out, const map_type& map, const inv_map_type& inv_map) {
            throw std::logic_error(util::demangle(typeid(wt_trait<character*>).name())+": serialize not implemented");
            return 0;
        }

        static void load_maps(std::istream& in, map_type& map, inv_map_type& inv_map) {
            throw std::logic_error(util::demangle(typeid(wt_trait<character*>).name())+": load not implemented");
            return 0;
        }
};

template<>
class wt_trait<unsigned char*>
{
    public:
        typedef bit_vector::size_type 		size_type;
        typedef unsigned char				value_type;
        typedef unsigned char*				reference_type;
        typedef unsigned_char_map			map_type;
        typedef unsigned_char_map			inv_map_type;
        enum { char_node_map_size=256 };

        static size_type alphabet_size_and_map(const reference_type rac, size_type n, map_type& map, inv_map_type& inv_map, value_type& first_symbol) {
            map.clear();
            inv_map.clear();
            if (n==0) {
                for (size_type i=0; i<256; ++i) {
                    map[i] = 255;    // mark each symbol as absent
                }
                return 0;
            }
            first_symbol	= *rac;
            map[*rac] = 0;
            inv_map[0] = *rac;
            size_type alphabet_size = 0;

            for (size_type i=0; i<256; ++i) {
                map[i] = 0;
            }

            for (size_type i=0; i<n; ++i) {
                value_type c = *(rac+i);
                map[c] = 1;
            }

            for (size_type i=0; i<256; ++i) {
                if (map[i]) {
                    map[i] = alphabet_size;
                    ++alphabet_size;
                } else {
                    map[i] = 255;
                }
                inv_map[map[i]] = i;
            }
            return alphabet_size;
        }

        static bool symbol_available(const map_type& map, const value_type c, const value_type first_symbol, const size_type sigma) {
            return sigma==256 or map[c] < 255;
        }

        static size_type serialize_maps(std::ostream& out, const map_type& map, const inv_map_type& inv_map) {
            size_type written_bytes = 0;
            written_bytes += map.serialize(out);
            written_bytes += inv_map.serialize(out);
            return written_bytes;
        }

        static void load_maps(std::istream& in, map_type& map, inv_map_type& inv_map) {
            map.load(in);
            inv_map.load(in);
        }
};

template<class size_type_class>
class wt_trait<int_vector_file_buffer<8, size_type_class> >
{
    public:
        typedef size_type_class		size_type;
        typedef unsigned char		value_type;
        typedef int_vector_file_buffer<8, size_type_class>&		reference_type;
        typedef unsigned_char_map 	map_type;
        typedef unsigned_char_map 	inv_map_type;
        enum { char_node_map_size=256 };

        static size_type alphabet_size_and_map(reference_type rac, size_type n, map_type& map, inv_map_type& inv_map, value_type& first_symbol) {
            map.clear();
            inv_map.clear();
            if (n==0) {
                for (size_type i=0; i<256; ++i) {
                    map[i] = 255;    // mark each symbol as absent
                }
                return 0;
            }
            rac.reset();
            if (rac.int_vector_size < n) {
                throw std::logic_error("wt<int_vector_file_buffer<8> >: n > rac.int_vector_size!");
                return 0;
            }

            for (size_type i=0; i<256; ++i) {
                map[i] = 0;
            }

            size_type alphabet_size = 0;
            size_type r = rac.load_next_block();
            first_symbol = rac[0];
            for (size_type i=0, r_sum=0; r_sum < n;) {
                if (r_sum +r > n) {  // make sure that not more than n characters are read
                    r = n-r_sum;
                }
                for (; i< r_sum+r; ++i) {
                    value_type c = rac[i-r_sum];
                    map[c] = 1;
                }
                r_sum += r; r = rac.load_next_block();
            }

            for (size_type i=0; i<256; ++i) {
                if (map[i]) {
                    map[i] = alphabet_size;
                    ++alphabet_size;
                } else {
                    map[i] = 255;
                }
                inv_map[map[i]] = i;
            }
            return alphabet_size;
        }

        static bool symbol_available(const map_type& map, const value_type c, const value_type first_symbol, const size_type sigma) {
            return sigma==256 or map[c] < 255;
        }

        static size_type serialize_maps(std::ostream& out, const map_type& map, const inv_map_type& inv_map) {
            size_type written_bytes = 0;
            written_bytes += map.serialize(out);
            written_bytes += inv_map.serialize(out);
            return written_bytes;
        }

        static void load_maps(std::istream& in, map_type& map, inv_map_type& inv_map) {
            map.load(in);
            inv_map.load(in);
        }
};

//! A wavelet tree class
/*!
 * A wavelet tree is build for a vector of characters over the alphabet \f$\Sigma\f$.
 * This class should be used only for small alphabets \f$\Sigma << n\f$ (see wt_int for a wavelet tree for big alphabets).
 * The wavelet tree \f$wt\f$ consists of a tree of  bitvectors and provides three efficient methods:
 *   - The "[]"-operator: \f$wt[i]\f$ returns the ith symbol of vector for which the wavelet tree was build for.
 *   - The rank method: \f$wt.rank(i,c)\f$ returns the number of occurences of symbol \f$c\f$ in the prefix [0..i-1] in the vector for which the wavelet tree was build for.
 *   - The select method: \f$wt.select(j,c)\f$ returns the index \f$i\in [0..size()-1]\f$ of the jth occurence of symbol \f$c\f$.
 *
 *
 *	\par Space complexity
 *		\f$\Order{n\log|\Sigma| + 2|\Sigma|\log n}\f$ bits, where \f$n\f$ is the size of the vector the wavelet tree was build for.
 *
 *  \tparam RandomAccessContainer 	Type of the input sequence.
 *  \tparam BitVector				Type of the bitvector used for representing the wavelet tree.
 *  \tparam RankSupport				Type of the support structure for rank on ones.
 *  \tparam SelectSupport			Type of the support structure for select on ones.
 *  \tparam SelectSupport			Type of the support structure for select on ones.
 *
 *   The wavelet tree was proposed first by Grossi et al. 2003 and applied to the BWT in Foschini et al. 2004.
 *   @ingroup wt
 */
template<class RandomAccessContainer=unsigned char*,
         class BitVector         = bit_vector,
         class RankSupport       = typename BitVector::rank_1_type,
         class SelectSupport     = typename BitVector::select_1_type,
         class SelectSupportZero = typename BitVector::select_0_type>
class wt
{
    public:
        typedef typename wt_trait<RandomAccessContainer>::size_type 		size_type;
        typedef typename wt_trait<RandomAccessContainer>::value_type 		value_type;
        typedef typename wt_trait<RandomAccessContainer>::map_type		     map_type;
        typedef typename wt_trait<RandomAccessContainer>::inv_map_type	inv_map_type;
    private:
        size_type 			m_size;
        size_type 			m_sigma; 		//<- \f$ |\Sigma| \f$
        BitVector 			m_tree;			// bit vector to store the wavelet tree
        RankSupport			m_tree_rank;	// rank support for the wavelet tree bit vector
        SelectSupport		m_tree_select1;	// select support for the wavelet tree bit vector
        SelectSupportZero	m_tree_select0;
        int_vector<64>		m_node_pointers;
        int_vector<64>		m_node_pointers_rank;
        size_type 			m_char_node_map[wt_trait<RandomAccessContainer>::char_node_map_size];

        value_type			m_first_symbol;

        mutable map_type			m_char_map;		// map the characters to integers; we use mutable as the std::map has no const version of the []-operator
        mutable inv_map_type		m_inv_char_map;	// map the integers back to characters; we use mutable " "  "  " "  "  "  "  "  "  "  "  "  "  "    "

        void copy(const wt& wt) {
            m_size 			= wt.m_size;
            m_sigma 		= wt.m_sigma;
            m_tree			= wt.m_tree;
            m_tree_rank 	= wt.m_tree_rank;
            m_tree_rank.set_vector(&m_tree);
            m_tree_select1	= wt.m_tree_select1;
            m_tree_select1.set_vector(&m_tree);
            m_tree_select0	= wt.m_tree_select0;
            m_tree_select0.set_vector(&m_tree);
            m_node_pointers	= wt.m_node_pointers;
            m_node_pointers_rank = wt.m_node_pointers_rank;
            m_first_symbol	= wt.m_first_symbol;
            m_char_map		= wt.m_char_map;
            m_inv_char_map	= wt.m_inv_char_map;
            if (wt_trait<RandomAccessContainer>::char_node_map_size == 256) {
                for (size_type i=0; i<256; ++i) m_char_node_map[i] = wt.m_char_node_map[i];
            }
        }

        void init_char_node_map() {
            if (wt_trait<RandomAccessContainer>::char_node_map_size == 256) {
                for (size_type i=0; i<256; ++i) m_char_node_map[i] = 0;
            }
        }

    public:

        const size_type& sigma;

        // Default constructor
        wt():m_size(0),m_sigma(0), sigma(m_sigma) {};



        //! Constructor
        /*!
         *	\param rac Reference to the vector (or unsigned char array) for which the wavelet tree should be build.
         *	\param size Size of the prefix of the vector (or unsigned char array) for which the wavelet tree should be build.
         *	\par Time complexity
         *		\f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         */
        wt(const typename wt_trait<RandomAccessContainer>::reference_type rac, size_type size):m_size(size), m_sigma(0), sigma(m_sigma) {
#ifdef SDSL_DEBUG_WT
            std::cerr<<"Wavelet Tree construct: size="<<size<<std::endl;
#endif
            // calculate alphabet size and the mappings for the symbols to the integers and back
            m_sigma = wt_trait<RandomAccessContainer>::alphabet_size_and_map(rac, m_size, m_char_map, m_inv_char_map, m_first_symbol);
            init_char_node_map();
            int_vector<> node_sizes = int_vector<>(2*m_sigma+1, 0, bit_magic::l1BP(m_size)+1);
            m_node_pointers = int_vector<64>(node_sizes.size()+1, 0);
            m_node_pointers_rank = int_vector<64>(node_sizes.size()+1, 0);

            if (m_sigma < 2) {
                if (m_sigma == 1) {  // ==> m_size > 0
                    if (wt_trait<RandomAccessContainer>::char_node_map_size == 256)
                        m_char_node_map[ m_first_symbol ] = 0; // first symbol corresponds to root node of the wavelet tree
                }
                return;
            } else {
                // TODO: slow compared to other constructor which uses
                // int_vector_file_buffer
                for (size_type i=0; i<m_size; ++i) {
                    size_type lex_idx = m_char_map[rac[(const size_type)i]]; // lex_idx in [0..m_sigma-1]
                    size_type sigma = m_sigma, node=0;
                    while (sigma >= 2) {
#ifdef SDSL_DEBUG_WT
                        if (node>=node_sizes.size()) {
                            std::cerr<<" node="<<node<<" node_sizes.size()="<<node_sizes.size()<<" rac["<<i<<"]="<<rac[i]<<" orig_lex_idx="<< (int)m_char_map[rac[i]] <<" sigma="<<sigma<<" lex_idx="<<lex_idx<<std::endl;
                        }
#endif
                        assert(node<node_sizes.size());
                        node_sizes[node] = node_sizes[node]+1;
                        if (lex_idx < (sigma+1)/2) {
                            sigma = (sigma+1)/2;
                            node = 2*node+1;
                        } else {
                            lex_idx -= (sigma+1)/2;
                            sigma -= (sigma+1)/2;
                            node = 2*node+2;
                        }
                    }

                    if (wt_trait<RandomAccessContainer>::char_node_map_size == 256)
                        m_char_node_map[rac[i]] = node;
                    if (lex_idx!=0) {
                        std::cout<<" _"<<rac[i]<<" ("<<lex_idx<<")";
                    }
                }
#ifdef SDSL_DEBUG_WT
                std::cerr<<"wavelet tree node sizes"<<std::endl;
                for (size_type i=0; i<node_sizes.size(); ++i) {
                    std::cerr<<node_sizes[i]<<" ";
                } std::cerr << std::endl;
#endif

                size_type max_node = 0;
                m_node_pointers[0] = 0;
                for (size_type i=1; i < m_node_pointers.size(); ++i) {
                    m_node_pointers[i] = m_node_pointers[i-1] + node_sizes[i-1];
                    if (node_sizes[i-1]!=0)
                        max_node = i-1;
                }
                node_sizes = int_vector<>(max_node+1, 0, 64);
#ifdef SDSL_DEBUG_WT
                std::cerr<<"all nodes size = "<<m_node_pointers[max_node+1]<<std::endl;
#endif
                bit_vector tree = bit_vector(m_node_pointers[max_node+1]);
                for (size_type i=0; i<m_size; ++i) {
                    size_type lex_idx = m_char_map[rac[i]]; // lex_idx in [0..m_sigma-1]
//std::cerr<<"lex_idx="<<lex_idx<<" "<<rac[i]<<std::endl;
                    size_type sigma = m_sigma, node=0;
                    while (sigma >= 2) {
                        if (lex_idx >= ((sigma+1)/2))
                            tree[m_node_pointers[node]+node_sizes[node]] = 1;
                        node_sizes[node] = node_sizes[node]+1;
                        if (lex_idx < (sigma+1)/2) {
                            sigma = (sigma+1)/2;
                            node = 2*node+1;
                        } else {
                            lex_idx -= (sigma+1)/2;
                            sigma -= (sigma+1)/2;
                            node = 2*node+2;
                        }
                    }
                }
                util::assign(m_tree, tree);
                util::init_support(m_tree_rank, &m_tree);
                for (size_type i=0; i < m_node_pointers.size(); ++i) {
                    m_node_pointers_rank[i] = m_tree_rank(m_node_pointers[i]);
                }
                util::init_support(m_tree_select0, &m_tree);
                util::init_support(m_tree_select1, &m_tree);
            }
        }

        template<class size_type_class>
        wt(int_vector_file_buffer<8, size_type_class>& rac, size_type size):m_size(size), m_sigma(0), sigma(m_sigma) {
            construct(rac, size);
        }

        //! Construct the wavelet tree from a random access container
        /*! \param rac A random access container
         *	\param size The length of the prefix of the random access container, for which the wavelet tree should be build
         */
        template<class size_type_class>
        void construct(int_vector_file_buffer<8, size_type_class>& rac, size_type size) {
            m_size = size;
            init_char_node_map();
            typedef int_vector_file_buffer<8, size_type_class> tIVFB;
//#define SDSL_DEBUG_WT
#ifdef SDSL_DEBUG_WT
            std::cerr<<"Wavelet Tree construct for int_vector_file_buffer: size="<<size<<std::endl;
            stop_watch sw; sw.start();
#endif

//		std::cerr << wt_trait<RandomAccessContainer>::char_node_map_size << std::endl;
            // calculate alphabet size and the mappings for the symbols to the integers and back
            m_sigma = wt_trait<tIVFB>::alphabet_size_and_map(rac, m_size, m_char_map, m_inv_char_map, m_first_symbol);

            int_vector<64> node_sizes = int_vector<64>(2*m_sigma+1, 0/*, bit_magic::l1BP(m_size)+1*/);
            m_node_pointers = int_vector<64>(node_sizes.size()+1, 0);
            m_node_pointers_rank = int_vector<64>(node_sizes.size()+1, 0);

            if (m_sigma < 2) {
                if (1 == m_sigma) {  // handle special case more efficient
                    m_char_node_map[m_first_symbol] = 0; // map the first symbol to the root node of the wavelet tree
                }
            } else {
                // O(n + |\Sigma|\log|\Sigma|) algorithm for calculating node sizes
                size_type C[256] = {0};
                rac.reset();
                //  1. Count occurrences of characters
                for (size_type i=0, r_sum=0, r = rac.load_next_block(); r_sum < m_size;) {
                    if (r_sum + r > m_size) {  // read not more than size chars in the next loop
                        r = m_size-r_sum;
                    }
                    for (; i < r_sum+r; ++i) {
                        ++C[rac[i-r_sum]];
                    }
                    r_sum += r; r = rac.load_next_block();
                }
                //  2. Sum up the node sizes for each character
                for (size_type i=0; i < 256; ++i) {
                    if (C[i]) {
                        size_type lex_idx = m_char_map[i]; // lex_idx in [0..m_sigma-1]
                        size_type sigma = m_sigma, node=0;
                        while (sigma >= 2) {
                            assert(node<node_sizes.size());
                            node_sizes[node] = node_sizes[node]+C[i];
                            if (lex_idx < ((sigma+1)>>1)) {
                                sigma = ((sigma+1)>>1);
                                node = (node<<1)+1;
                            } else {
                                lex_idx -= ((sigma+1)>>1);
                                sigma -= ((sigma+1)>>1);
                                node = (node<<1)+2;
                            }
                        }
                        m_char_node_map[i] = node;
                    } else {
                        m_char_node_map[i] = 0;
                    }
                }
#ifdef SDSL_DEBUG_WT
                sw.stop();
                std::cerr<<"Time for counting phase: "<< sw.get_real_time() << " ms real time , "<< sw.get_user_time()<<" ms user time"<< std::endl;
//			std::cerr<<"wavelet tree node sizes"<<std::endl;
//			for(size_type i=0; i<node_sizes.size(); ++i){ std::cerr<<node_sizes[i]<<" "; } std::cerr << std::endl;
#endif
                size_type max_node = 0;
                m_node_pointers[0] = 0;
                for (size_type i=1; i < m_node_pointers.size(); ++i) {
                    m_node_pointers[i] = m_node_pointers[i-1] + node_sizes[i-1];
                    if (node_sizes[i-1]!=0)
                        max_node = i-1;
                }
                node_sizes = int_vector<64>(max_node+1, 0);
#ifdef SDSL_DEBUG_WT
                std::cerr<<"all nodes size = "<<m_node_pointers[max_node+1]<<std::endl;
                sw.start();
#endif
                // initialize bit vector with 0's
                bit_vector tree = bit_vector(m_node_pointers[max_node+1], 0);
                // precalc paths in the tree for all symbols in the alphabet
                uint8_t path_len[256] = {0};
                uint8_t path[256] = {0};
                for (size_type i=0; i < 256; ++i) {
                    if (C[i]) {
                        size_type lex_idx = m_char_map[i]; // lex_idx in [0..m_sigma-1]
                        size_type sigma = m_sigma, depth=0;
                        while (sigma >= 2) {
                            if (lex_idx < ((sigma+1)>>1)) {
                                sigma = ((sigma+1)>>1);
                            } else {
                                lex_idx -= ((sigma+1)>>1);
                                sigma -= ((sigma+1)>>1);
                                path[i] |= (1<<depth);
                            }
                            ++path_len[i];
                            ++depth;
                        }
                    }
                }
///*
                rac.reset();
                for (size_type i=0, r_sum=0, r = rac.load_next_block(); r_sum < m_size;) {
                    if (r_sum + r > size) {  // read not more than size chars in the next loop
                        r = size-r_sum;
                    }
                    uint8_t old_chr = rac[i-r_sum];
                    uint8_t times = 0;
                    for (; i < r_sum+r; ++i) {
                        uint8_t chr = rac[i-r_sum];
                        if (chr	!= old_chr) {
                            uint8_t p = path[old_chr];
                            for (uint32_t l=0, node=0; l<path_len[old_chr]; ++l, p >>= 1) {
                                if (p&1) {
                                    tree.set_int(m_node_pointers[node]+node_sizes[node], 0xFFFFFFFFFFFFFFFFULL,times);
                                    node_sizes[node] += times; node = (node<<1)+2;
                                } else {
                                    node_sizes[node] += times; node = (node<<1)+1;
                                }
                            }
                            times = 1;
                            old_chr = chr;
                        } else { // chr == old_chr
                            ++times;
                            if (times == 64) {
                                uint8_t p = path[old_chr];
                                for (uint32_t l=0, node=0; l<path_len[old_chr]; ++l, p >>= 1) {
                                    if (p&1) {
                                        tree.set_int(m_node_pointers[node]+node_sizes[node], 0xFFFFFFFFFFFFFFFFULL,times);
                                        node_sizes[node] += times; node = (node<<1)+2;
                                    } else {
                                        node_sizes[node] += times; node = (node<<1)+1;
                                    }
                                }
                                times = 0;
                            }
                        }
                    }
                    if (times > 0) {
                        uint8_t p = path[old_chr];
                        for (uint32_t l=0, node=0; l<path_len[old_chr]; ++l, p >>= 1) {
                            if (p&1) {
                                tree.set_int(m_node_pointers[node]+node_sizes[node], 0xFFFFFFFFFFFFFFFFULL,times);
                                node_sizes[node] += times; node = (node<<1)+2;
                            } else {
                                node_sizes[node] += times; node = (node<<1)+1;
                            }
                        }
                    }
                    r_sum += r; r = rac.load_next_block();
                }
//*/
#ifdef SDSL_DEBUG_WT
                sw.stop();
                std::cerr<<"Time for building phase: "<< sw.get_real_time() << " ms real time , "<< sw.get_user_time()<<" ms user time"<< std::endl;
                sw.start();
#endif
                util::assign(m_tree, tree);
                util::init_support(m_tree_rank,&m_tree);
#ifdef SDSL_DEBUG_WT
                sw.stop();
                std::cerr<<"Time for rank init: "<< sw.get_real_time() << " ms real time , "<< sw.get_user_time()<<" ms user time"<< std::endl;
                sw.start();
#endif
                for (size_type i=0; i < m_node_pointers.size(); ++i) {
                    m_node_pointers_rank[i] = m_tree_rank(m_node_pointers[i]);
                }
                util::init_support(m_tree_select0,&m_tree);
#ifdef SDSL_DEBUG_WT
                std::cerr<<"select0 init ready!"<<std::endl;
#endif
                util::init_support(m_tree_select1,&m_tree);
#ifdef SDSL_DEBUG_WT
                sw.stop();
                std::cerr<<"Time for select init: "<< sw.get_real_time() << " ms real time , "<< sw.get_user_time()<<" ms user time"<< std::endl;
                sw.start();
#endif
            }
        }

        //! Copy constructor
        wt(const wt& wt):sigma(m_sigma) {
            copy(wt);
        }

        //! Assignment operator
        wt& operator=(const wt& wt) {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Swap operator
        void swap(wt& wt) {
            if (this != &wt) {
                std::swap(m_size, wt.m_size);
                std::swap(m_sigma,  wt.m_sigma);
                m_tree.swap(wt.m_tree);
                util::swap_support(m_tree_rank, wt.m_tree_rank, &m_tree, &(wt.m_tree));
                util::swap_support(m_tree_select1, wt.m_tree_select1, &m_tree, &(wt.m_tree));
                util::swap_support(m_tree_select0, wt.m_tree_select0, &m_tree, &(wt.m_tree));
                m_node_pointers.swap(wt.m_node_pointers);
                m_node_pointers_rank.swap(wt.m_node_pointers_rank);
                std::swap(m_first_symbol, wt.m_first_symbol);
                m_char_map.swap(wt.m_char_map);
                m_inv_char_map.swap(wt.m_inv_char_map);
                // swap char_node_map
                if (wt_trait<RandomAccessContainer>::char_node_map_size == 256) {
                    for (size_type i=0; i<256; ++i) {
                        std::swap(m_char_node_map[i],wt.m_char_node_map[i]);
                    }
                }
            }
        }

        //! Returns the size of the original vector.
        size_type size()const {
            return m_size;
        }

        //! Returns whether the wavelet tree contains no data.
        bool empty()const {
            return m_size == 0;
        }

        //! Recovers the ith symbol of the original vector.
        /*! \param i The index of the symbol in the original vector. \f$i \in [0..size()-1]\f$
         *	\returns The ith symbol of the original vector.
         */
        value_type operator[](size_type i)const {
            size_type lex_idx	= 0;
            size_type sigma		= m_sigma;
            size_type node		= 0;
            while (sigma >= 2) {
                if (m_tree[ m_node_pointers[node]+i ]) { // go to the right child
                    lex_idx += (sigma+1)/2;
                    i = m_tree_rank(m_node_pointers[node]+i) - m_node_pointers_rank[node];  //m_tree_rank( m_node_pointers[node] );
                    node = 2*node+2;
                    sigma -= (sigma+1)/2;
                } else { // go to the left child
                    i = i - (m_tree_rank(m_node_pointers[node]+i) - m_node_pointers_rank[node]);   //m_tree_rank( m_node_pointers[node] ) );
                    node = 2*node+1;
                    sigma = (sigma+1)/2;
                }
            }
            return m_inv_char_map[lex_idx];
        };

        //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
         *  \param c The symbol to count the occurences in the prefix.
         *	\returns The number of occurences of symbol c in the prefix [0..i-1] of the supported vector.
         *  \par Time complexity
         *		\f$ \Order{\log |\Sigma|} \f$
         */
        size_type rank(size_type i, value_type c)const {
            if (!wt_trait<RandomAccessContainer>::symbol_available(m_char_map, c, m_first_symbol, m_sigma)) {
                return 0;
            }
            size_type lex_idx 	= m_char_map[c]; // koennte man auch nur path, path_len ersetzen
            size_type sigma 	= m_sigma;
            size_type node		= 0;
            size_type result	= i;
            while (sigma >= 2 and result > 0) {
                if (lex_idx < (sigma+1)/2) {
                    result 	= result - (m_tree_rank(m_node_pointers[node]+result) -  m_node_pointers_rank[node]); //m_tree_rank(m_node_pointers[node]));
                    sigma 	= (sigma+1)/2;
                    node	= 2*node+1;
                } else {
                    result = m_tree_rank(m_node_pointers[node]+result) -  m_node_pointers_rank[node]; // m_tree_rank(m_node_pointers[node]);
                    lex_idx -= (sigma+1)/2;
                    sigma -= (sigma+1)/2;
                    node = 2*node+2;
                }
            }
            return result;
        };

        //! Calculates how many occurrences of symbol wt[i] are in the prefix [0..i-1] of the supported sequence.
        /*!
         *	\param i The index of the symbol.
         *  \return The number of occurrences of symbol wt[i] in the prefix [0..i-1]
         *	\par Time complexity
         *		\f$ \Order{\log |\Sigma|} \f$
         */

        size_type rank_ith_symbol(size_type i, value_type& c)const {
            assert(i>=0 and i < size());
            size_type lex_idx	= 0;
            size_type sigma		= m_sigma;
            size_type node		= 0;
            while (sigma >= 2) {
                if (m_tree[ m_node_pointers[node]+i ]) { // go to the right child
                    lex_idx += (sigma+1)/2;
                    i = m_tree_rank(m_node_pointers[node]+i) - m_node_pointers_rank[node];  //m_tree_rank( m_node_pointers[node] );
                    node = 2*node+2;
                    sigma -= (sigma+1)/2;
                } else { // go to the left child
                    i = i - (m_tree_rank(m_node_pointers[node]+i) - m_node_pointers_rank[node]);   //m_tree_rank( m_node_pointers[node] ) );
                    node = 2*node+1;
                    sigma = (sigma+1)/2;
                }
            }
            c = m_inv_char_map[lex_idx];
            return i;
        }

        // recursive internal version of the method interval_symbols
        void _interval_symbols(size_type i, size_type j, size_type& k,
                               std::vector<unsigned char>& cs,
                               std::vector<size_type>& rank_c_i,
                               std::vector<size_type>& rank_c_j,
                               size_type lex_idx, size_type sigma, size_type node) const {
            // not in a leaf
            if (sigma >= 2) {
                size_type i_new = (m_tree_rank(m_node_pointers[node] + i) - m_node_pointers_rank[node]);
                size_type j_new = (m_tree_rank(m_node_pointers[node] + j) - m_node_pointers_rank[node]);
                // goto left child
                i -= i_new; j -= j_new;
                if (i != j) {
                    _interval_symbols(i, j, k, cs, rank_c_i, rank_c_j, lex_idx, (sigma+1)>>1, (node<<1)|1 /*2*node+1*/);
                }
                // goto right child
                if (i_new != j_new) {
                    _interval_symbols(i_new, j_new, k, cs, rank_c_i, rank_c_j, lex_idx + ((sigma+1)>>1), sigma-((sigma+1)>>1), (node+1)<<1 /*2*node+2*/);
                }
            } else {
                rank_c_i[k] = i;
                rank_c_j[k] = j;
                cs[k++] = m_inv_char_map[lex_idx];
                return;
            }
        }

        //! Counts the characters in the range [0..i-1] which are smaller than character c.
        /* If the character c does not occur in the sequence 0 is returned.
         *
         */
        size_type count_lex_smaller(size_type i, value_type c)const {
            if (!wt_trait<RandomAccessContainer>::symbol_available(m_char_map, c, m_first_symbol, m_sigma)) {
                return 0;
            }
            size_type lex_idx 	= m_char_map[c];
            size_type sigma 	= m_sigma;  // start with the whole alphabet
            size_type node		= 0;
            size_type result	= 0;
            while (sigma >= 2) {
                if (lex_idx < (sigma+1)/2) { // symbols belongs to the left half of the alphabet
                    // calculate new i for the left child bit_vector
                    i 	= i - (m_tree_rank(m_node_pointers[node]+i) -  m_node_pointers_rank[node]);
                    sigma 	= (sigma+1)/2;
                    node	= 2*node+1;
                } else { // symbol belongs to the right half of the alphabet
                    size_type ones = m_tree_rank(m_node_pointers[node]+i) -  m_node_pointers_rank[node];
                    result += (i - ones); // all elements prefixed with 0 are lexicographic smaller than c
                    i = ones;
                    // calclate new i for the right child bit_vector
                    lex_idx -= (sigma+1)/2;
                    sigma -= (sigma+1)/2;
                    node = 2*node+2;
                }
            }
            return result;
        }

        //! Counts the characters in the range [i..j-1] which are smaller than character c.
        size_type count_lex_smaller(size_type i, size_type j, value_type c)const {
            if (i==j)
                return 0;
            if (i+1 == j) {
                return (*this)[i] < c;
            } else {
                return count_lex_smaller(j, c) - count_lex_smaller(i, c);
            }
        }

        //! Calculates for each symbol c in wt[i..j-1], how many times c occures in wt[0..i-1] and wt[0..j-1].
        /*!
         *	\param i The start index (inclusive) of the interval.
         *	\param j The end index (exclusive) of the interval.
         *	\param k Reference that will contain the number of different symbols in wt[i..j-1].
         *  \param cs Reference to a vector of size k that will contain all symbols that occur in wt[i..j-1] in ascending order.
         *  \param rank_c_i Reference to a vector which equals rank_c_i[p] = rank(cs[p],i), for \f$ 0 \leq p < k \f$
         *  \param rank_c_j Reference to a vector which equals rank_c_j[p] = rank(cs[p],j), for \f$ 0 \leq p < k \f$
         *	\par Time complexity
         *		\f$ \Order{\min{\sigma, k \log \sigma}} \f$
         *
         *  \par Precondition
         *       \f$ i\leq j \f$
         *       \f$ cs.size() \geq \sigma \f$
         *       \f$ rank_c_i.size() \geq \sigma \f$
         *       \f$ rank_c_j.size() \geq \sigma \f$
         */
        void interval_symbols(size_type i, size_type j, size_type& k,
                              std::vector<unsigned char>& cs,
                              std::vector<size_type>& rank_c_i,
                              std::vector<size_type>& rank_c_j) const {
            if (i==j) {
                k = 0;
                return;
            } else if ((j-i)==1) {
                k = 1;
                rank_c_i[0] = rank_ith_symbol(i, cs[0]);
                rank_c_j[0] = rank_c_i[0]+1;
                return;
            } else if ((j-i)==2) {
                rank_c_i[0] = rank_ith_symbol(i, cs[0]);
                rank_c_i[1] = rank_ith_symbol(i+1, cs[1]);
                if (cs[0]==cs[1]) {
                    k = 1;
                    rank_c_j[0] = rank_c_i[0]+2;
                    return;
                } else {
                    k = 2;
                    rank_c_j[0] = rank_c_i[0]+1;
                    rank_c_j[1] = rank_c_i[1]+1;
                    return;
                }
            } else {
                k = 0;
                _interval_symbols(i, j, k, cs, rank_c_i, rank_c_j, 0, m_sigma, 0);
            }
        }


        //! Calculates the ith occurence of the symbol c in the supported vector.
        /*!
         *  \param i The ith occurence. \f$i\in [1..rank(size(),c)]\f$.
         *  \param c The symbol c.
         *  \par Time complexity
         *		\f$ \Order{\log |\Sigma|} \f$
         */
        // TODO: was ist wenn c gar nicht vorkommt, oder es keine i Vorkommen gibt?
        size_type select(size_type i, value_type c)const {
            if (!wt_trait<RandomAccessContainer>::symbol_available(m_char_map, c, m_first_symbol, m_sigma)) {
                return size();
            }
            assert(i > 0);
#ifdef SDSL_DEBUG_WT
            std::cerr<<"wt: select("<<i<<", "<<c<<")"<<std::endl;
            std::cerr<<" i="<<i<< " <= "<< rank(size(),c) <<" ?" <<std::endl;
#endif
            assert(i <= rank(size(), c));

            size_type node		= 0;
            // first phase: go down to the node corresponding to the lex_idx of the character
            if (wt_trait<RandomAccessContainer>::char_node_map_size == 256) {
                node = m_char_node_map[c];
            } else {
                size_type lex_idx 	= m_char_map[c];
                size_type sigma 	= m_sigma;
                while (sigma >= 2) {
                    if (lex_idx < (sigma+1)/2) {
                        sigma = (sigma+1)/2;
                        node = 2*node+1;
                    } else {
                        lex_idx -= (sigma+1)/2;
                        sigma -= (sigma+1)/2;
                        node = 2*node+2;
                    }
                }
            }
            // second phase go up and select the right position
            size_type result = i-1;
            while (node != 0) {
                if (node&1) {// node is a left child
                    node = (node-1)/2;
                    result = m_tree_select0(m_node_pointers[node]-m_node_pointers_rank[node]+result+1)-m_node_pointers[node];
                } else { //node is a right child
                    node = (node-1)/2;
                    result = m_tree_select1(m_node_pointers_rank[node]+result+1)-m_node_pointers[node];
                }
            }
            return result;
        };

        //!
        void range_search_2d(size_type lb, size_type rb, value_type c1, value_type c2, std::vector<size_type>& result) const {
            size_type lex_idx_1 = m_char_map[c1];
            size_type lex_idx_2 = m_char_map[c2];
            _range_search_2d(0, lb, rb, lex_idx_1, lex_idx_2, m_sigma, result);
        }

        // add parameter path
        void _range_search_2d(size_type node, size_type lb, size_type rb, size_type lex_idx_1, size_type lex_idx_2, size_type sigma, std::vector<size_type>& result)const {
//		std::cerr<<"lb="<<lb<<" rb="<<rb<<" lex_idx_1="<<lex_idx_1<<" lex_idx_2="<<lex_idx_2<<" sigma="<<sigma<<" node="<<node<<std::endl;
            // [lex_idx_1..lex_idx_2] in [0..sigma]
            if (lb > rb)
                return;
            if (sigma <= 1) {
                // result[i]...
                size_type r = 0;
                while (node != 0) {
                    if (node&1) {// node is a left child
                        node = (node-1)/2;
                        size_type node_idx = m_node_pointers[node];
                        r = m_tree_select0(node_idx-m_tree_rank(node_idx)+r+1)-node_idx;
                    } else { //node is a right child
                        node = (node-1)/2;
                        r = m_tree_select1(m_tree_rank(m_node_pointers[node])+r+1)-m_node_pointers[node];
                    }
                }
                result.push_back(r);
                return;
            }
            size_type lex_mid = (sigma+1)/2;
            size_type lsigma = lex_mid;
            size_type rsigma = sigma - lex_mid;

            size_type ones_lb_1 = m_tree_rank(m_node_pointers[node]+lb) -  m_node_pointers_rank[node];  // ones in [0,lb)
            size_type ones_rb 	= m_tree_rank(m_node_pointers[node]+rb+1) -	m_node_pointers_rank[node]; // ones in [0,rb]
            size_type zeros_lb_1 = lb-ones_lb_1; // zeros in [0,lb)
            size_type zeros_rb	= rb+1-ones_rb; // zeros in [0,rb]

//		std::cerr<<"ones_lb_1="<<ones_lb_1<<" ones_rb="<<ones_rb<<std::endl;
//		std::cerr<<"zeros_lb_1="<<zeros_lb_1<<" zeros_rb="<<zeros_rb<<std::endl;

            if (lex_idx_1 < lex_mid and zeros_rb) {
                _range_search_2d(2*node+1, zeros_lb_1, zeros_rb-1, lex_idx_1, std::min(lex_idx_2,lex_mid-1), lsigma, result);
            }

            if (lex_idx_2 >= lex_mid  and  ones_rb) {
                size_type _lex_idx_1 = 0;
                if (lex_idx_1 > lex_mid)
                    _lex_idx_1 = lex_idx_1 - lex_mid;
                _range_search_2d(2*node+2, ones_lb_1, ones_rb-1, _lex_idx_1, lex_idx_2-lex_mid, rsigma, result);
            }
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out)const {
            size_type written_bytes = 0;
            written_bytes += util::write_member(m_size, out);
            written_bytes += util::write_member(m_sigma, out);
            written_bytes += m_tree.serialize(out);
            written_bytes += m_tree_rank.serialize(out);
            written_bytes += m_tree_select1.serialize(out);
            written_bytes += m_tree_select0.serialize(out);
            written_bytes += m_node_pointers.serialize(out);
            written_bytes += m_node_pointers_rank.serialize(out);
            written_bytes += wt_trait<RandomAccessContainer>::serialize_maps(out, m_char_map, m_inv_char_map);
            written_bytes += util::write_member(m_first_symbol, out);
            // serialize char_node_map
            if (wt_trait<RandomAccessContainer>::char_node_map_size == 256) {
                for (size_type i=0; i<256; ++i) {
                    written_bytes += util::write_member(m_char_node_map[i], out);
                }
            }
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            util::read_member(m_size, in);
            util::read_member(m_sigma, in);
            m_tree.load(in);
            m_tree_rank.load(in, &m_tree);
            m_tree_select1.load(in, &m_tree);
            m_tree_select0.load(in, &m_tree);
            m_node_pointers.load(in);
            m_node_pointers_rank.load(in);
            wt_trait<RandomAccessContainer>::load_maps(in, m_char_map, m_inv_char_map);
            util::read_member(m_first_symbol, in);
            // serialize char_node_map
            if (wt_trait<RandomAccessContainer>::char_node_map_size == 256) {
                for (size_type i=0; i<256; ++i) {
                    util::read_member(m_char_node_map[i], in);
                }
            }
        }

#ifdef MEM_INFO
        //! Print some infos about the size of the compressed suffix tree
        void mem_info(std::string label="")const {
            if (label=="")
                label="csa";
            size_type bytes = util::get_size_in_bytes(*this);
            std::cout << "list(label = \""<<label<<"\", size = "<< bytes/(1024.0*1024.0) <<"\n,";
            m_tree.mem_info("data");
            std::cout << ")\n";
        }
#endif
        /*
        	void print_info()const{
        		size_type rle_ed = 0;
        		for(size_type i=0; i < m_tree.size(); ++i){

        		}
        	}
        */

};

}// end namespace sds

#endif // end file 
