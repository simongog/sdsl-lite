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
/*! \file construct_lcp.hpp
    \brief construct_lcp.hpp contains a space and time efficient construction method for lcp arrays
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_CONSTRUCT_LCP
#define INCLUDED_SDSL_CONSTRUCT_LCP

#include "typedefs.hpp"  // includes definition of tMSS
#include "config.hpp"
#include "int_vector.hpp"
#include "rank_support_v.hpp"
#include "util.hpp"
#include "testutils.hpp"
#include "construct_isa.hpp"
#include "construct_bwt.hpp"
#include "wt_huff.hpp"
#include "construct_lcp_helper.hpp"

#include <iostream>
#include <stdexcept>
#include <list>
#include <algorithm>
#include <fstream>
#include <queue>
#include <list>

//#define STUDY_INFORMATIONS

namespace sdsl
{

//! Construct the LCP array for text over byte- or integer-alphabet.
/*!	The algorithm computes the lcp array and stores it to disk.
 *  \tparam t_width Width of the text. 0==integer alphabet, 8=byte alphabet.
 *  \param config	Reference to cache configuration
 *  \par Space complexity
 *		\f$ n (\log \sigma + \log n) \f$ bits
 *  \pre Text and Suffix array exist in the cache. Keys:
 *         * constants::KEY_TEXT for t_width=8  or constants::KEY_TEXT_INT for t_width=0
 *         * constants::KEY_SA
 *  \post LCP array exist in the cache. Key
 *         * constants::KEY_LCP
 *  \par Reference 
 *    Toru Kasai, Gunho Lee, Hiroki Arimura, Setsuo Arikawa, Kunsoo Park:
 *    Linear-Time Longest-Common-Prefix Computation in Suffix Arrays and Its Applications.
 *    CPM 2001: 181-192
 */
template<uint8_t t_width>
void construct_lcp_kasai(cache_config& config){
    int_vector<> lcp;
    typedef int_vector<>::size_type size_type;
    write_R_output("lcp", "construct LCP", "begin", 1, 0);
    construct_isa(config);  
    {
        write_R_output("lcp", "load text", "begin", 1, 0);
		int_vector<t_width> text;
        if (!util::load_from_cache(text, key_text_trait<t_width>::KEY_TEXT, config)) { return; }
        write_R_output("lcp", "load text", "end", 1, 0);
        int_vector_file_buffer<> isa_buf(config.file_map[constants::KEY_ISA], 1000000);   // init isa file_buffer
        int_vector<> sa;
        if (!util::load_from_cache(sa, constants::KEY_SA, config)) { return; }
        // use Kasai algorithm to compute the lcp values
        for (size_type i=0,j=0,sa_1=0,l=0, r_sum=0, r=isa_buf.load_next_block(), n=isa_buf.int_vector_size; r_sum < n;) {
            for (; i < r_sum+r; ++i) {
                sa_1 =  isa_buf[i-r_sum]; // = isa[i]
                if (sa_1) {
                    j = sa[sa_1-1];
                    if (l) --l;
                    assert(i!=j);
                    while (text[i+l]==text[j+l]) { // i+l < n and j+l < n are not necessary, since text[n]=0 and text[i]!=0 (i<n) and i!=j
                        ++l;
                    }
                    sa[ sa_1-1 ] = l; //overwrite sa array with lcp values
                } else {
                    l = 0;
                    sa[ n-1 ] = 0;
                }
            }
            r_sum += r;
            r = isa_buf.load_next_block();
        }

		for (size_type i=sa.size(); i>1; --i){
			sa[i-1] = sa[i-2];
		}
		sa[0] = 0;
		lcp.swap(sa);
    }
    write_R_output("lcp", "construct LCP", "end", 1, 0);
	util::store_to_cache(lcp, constants::KEY_LCP, config);
}

//! Construct the LCP array for text over byte- or integer-alphabet.
/*!	The algorithm computes the lcp array and stores it to disk.
 *  \par Space complexity
 *		\f$ n( \log \sigma + \log \n ) \f$ bits 
 *  \pre Text and Suffix array exist in the cache. Keys:
 *         * constants::KEY_TEXT_INT
 *         * constants::KEY_SA 
 *  \post LCP array exist in the cache. Key
 *         * constants::KEY_LCP
 *  \par Reference
 *     Juha Kärkkäinen, Giovanni Manzini, Simon J. Puglisi: 
 *     Permuted Longest-Common-Prefix Array. 
 *     CPM 2009: 181-192
 */
template<uint8_t t_width>
void construct_lcp_PHI(cache_config& config) {
    typedef int_vector<>::size_type size_type;
	typedef int_vector<t_width> text_type;
	const char * KEY_TEXT = key_text_trait<t_width>::KEY_TEXT;
    write_R_output("lcp", "construct LCP", "begin", 1, 0);
    int_vector_file_buffer<> sa_buf(config.file_map[constants::KEY_SA]);
    size_type n = sa_buf.int_vector_size; 

	assert( n > 0 );
	if ( 1 == n ){ // Handle special case: Input only the sentinel character.
		int_vector<> lcp(1, 0);
		util::store_to_cache( lcp, constants::KEY_LCP, config );
		return;
	}

//	(1) Calculate PHI (stored in array plcp)
    int_vector<> plcp(n, 0, sa_buf.width);
    for (size_type i=0, r_sum=0, r=sa_buf.load_next_block(), sai_1 = 0; r_sum < n;) {
        for (; i < r_sum+r; ++i) {
            size_type sai = sa_buf[i-r_sum];
            plcp[ sai ] = sai_1;
            sai_1 = sai;
        }
        r_sum += r; r = sa_buf.load_next_block();
    }

//  (2) Load text from disk
    write_R_output("lcp", "load text", "begin", 1, 0);
	text_type text;
	util::load_from_cache(text, KEY_TEXT, config);
    write_R_output("lcp", "load text", "end", 1, 0);

//  (3) Calculate permuted LCP array (text order), called PLCP
	size_type max_l = 0;
	for (size_type i=0, l=0; i < n-1; ++i) {
		size_type phii = plcp[i];
		while (text[i+l] == text[phii+l]) {
			++l;
		}
		plcp[i] = l;
		if ( l ){
			max_l = std::max(max_l, l);
			--l;
		}
	}
	util::clear(text);
	uint8_t lcp_width = bit_magic::l1BP(max_l)+1;

//	(4) Transform PLCP into LCP
	std::string lcp_file = util::cache_file_name(constants::KEY_LCP, config);
    std::ofstream lcp_out_buf(lcp_file.c_str(), std::ios::binary | std::ios::app | std::ios::out);   // open buffer for lcp

    size_type bit_size = n*lcp_width;
    lcp_out_buf.write((char*) &(bit_size), sizeof(bit_size));	// write size of vector
    lcp_out_buf.write((char*) &(lcp_width),sizeof(lcp_width));  // write int_width of vector
    size_type wb = 0;  // bytes written into lcp int_vector

    size_type buffer_size = 1000000; // buffer_size is a multiple of 8!

    int_vector<> lcp_buf(buffer_size, 0, lcp_width);
    lcp_buf[0] = 0;
    sa_buf.reset(buffer_size);
    size_type r = 0;// sa_buf.load_next_block();
    for (size_type i=1, r_sum=0; r_sum < n;) {
        for (; i < r_sum+r; ++i) {
            size_type sai = sa_buf[i-r_sum];
            lcp_buf[ i-r_sum ] = plcp[sai];
        }
        if (r > 0) {
            size_type cur_wb = (r*lcp_buf.width()+7)/8;
            lcp_out_buf.write((const char*)lcp_buf.data(), cur_wb);
            wb += cur_wb;
        }
        r_sum += r; r = sa_buf.load_next_block();
    }
    if (wb%8) {
        lcp_out_buf.write("\0\0\0\0\0\0\0\0", 8-wb%8);
    }
    lcp_out_buf.close();
	util::register_cache_file(constants::KEY_LCP, config);
    write_R_output("lcp", "construct LCP", "end", 1, 0);
}


//! Construct the LCP array out of the BWT (only for byte strings)
/*!	The algorithm computes the lcp array and stores it to disk. It needs only the Burrows and Wheeler transform.
 *  \param config	Reference to cache configuration
 *  \pre BWT exist in the cache. Keys:
 *         * constants::KEY_BWT
 *  \post LCP array exist in the cache. Key
 *         * constants::KEY_LCP
 *
 *  \par Time complexity
 *		  \f$ \Order{n \log{\sigma}} \f$
 *  \par Space complexity
 *		  Usually not more than \f$ 2.5n \f$ bytes
 *  \par Reference
 *        Timo Beller, Simon Gog, Enno Ohlebusch, Thomas Schnattinger: 
 *        Computing the Longest Common Prefix Array Based on the Burrows-Wheeler Transform.
 *        SPIRE 2011: 197-208
 */
void construct_lcp_bwt_based(cache_config& config);

//! Construct the LCP array out of the BWT (only for byte strings)
/*!	The algorithm computes the lcp array and stores it to disk. It needs only the Burrows and Wheeler transform.
 *  \param config	Reference to cache configuration
 *  \pre BWT exist in the cache. Keys:
 *         * constants::KEY_BWT
 *  \post LCP array exist in the cache. Key
 *         * constants::KEY_LCP
 *
 *  \par Time complexity
 *		  \f$ \Order{n \log{\sigma}} \f$
 *  \par Space complexity
 *		  Usually not more than \f$ 1.5n \f$ bytes
 *  \par Reference
 * 	      Timo Beller, Simon Gog, Enno Ohlebusch, Thomas Schnattinger: 
 *        Computing the longest common prefix array based on the Burrows-Wheeler transform. 
 *        J. Discrete Algorithms 18: 22-31 (2013)
 */
void construct_lcp_bwt_based2(cache_config& config);


// semi extern PHI algorithm of Karkainen, Manzini and Puglisi CPM 2009
//bool construct_lcp_semi_extern_PHI(tMSS& file_map, const std::string& dir, const std::string& id);

//class buffered_char_queue
//{
//        typedef bit_vector::size_type size_type;
//        typedef std::queue<uint8_t> tQ;
//    private:
//        static const uint32_t m_buffer_size =  10000;//409600;
//        uint8_t m_write_buf[m_buffer_size];
//        uint8_t m_read_buf[m_buffer_size];
//        size_type 	m_widx; // write index
//        size_type 	m_ridx; // read index
//        bool		m_sync; // are read and write buffer the same?
//        size_type 	m_disk_buffered_blocks; // number of blocks written to disk and not read again yet
//        char 		m_c;
//        size_type	m_rb; // read blocks
//        size_type	m_wb; // written blocks
//
//        std::string m_file_name;
//
//        std::fstream	m_stream;
//
//    public:
//
//        buffered_char_queue();
//        void init(const std::string& dir, char c);
//        ~buffered_char_queue();
//        void push_back(uint8_t x);
//        uint8_t pop_front();
//};
//
//typedef std::list<int_vector<>::size_type> tLI;
//typedef std::vector<int_vector<>::size_type> tVI;
//
//template<class size_type_class>
//void push_front_m_index(size_type_class i, uint8_t c, tLI(&m_list)[256], uint8_t (&m_chars)[256], size_type_class& m_char_count)
//{
//    if (m_list[c].empty()) {
//        m_chars[m_char_count++] = c;
//    }
//    m_list[c].push_front(i);
//}
//
//template<class size_type_class>
//void push_back_m_index(size_type_class i, uint8_t c, tLI(&m_list)[256], uint8_t (&m_chars)[256], size_type_class& m_char_count)
//{
//    if (m_list[c].empty()) {
//        m_chars[m_char_count++] = c;
//    }
//    m_list[c].push_back(i);
//}
//
//// only phase 1 of the new algorithm
//bool construct_lcp_simple_5n(tMSS& file_map, const std::string& dir, const std::string& id);
//
//// only phase 2 of the new algorithm
//// TODO: assert n > 0
//bool construct_lcp_simple2_9n(tMSS& file_map, const std::string& dir, const std::string& id);
//
////! Our new 2 phases lcp algorithm using usually 2n bytes.
///*!	The algorithm computes the lcp array and stores it to disk.
// *  \param file_map A map which contains the filenames of previous computed structures (like suffix array, Burrows and Wheeler transform, inverse suffix array, text,...)
// *  \param dir		Directory where the lcp array should be stored.
// *  \param id		Id for the file name of the lcp array.
// *  \par Space complexity
// *		Usually \f$ 2n \f$ bytes, worst case \f$5n bytes\f$
// */
//bool construct_lcp_go(tMSS& file_map, const std::string& dir, const std::string& id);
//
////! Our new 2 phases lcp algorithm using usually 2n bytes.
///*!	The algorithm computes the lcp array and stores it to disk.
// *  \param file_map A map which contains the filenames of previous computed structures (like suffix array, Burrows and Wheeler transform, inverse suffix array, text,...)
// *  \param dir		Directory where the lcp array should be stored.
// *  \param id		Id for the file name of the lcp array.
// *  \par Space complexity
// *		Usually \f$ 2n \f$ bytes, worst case \f$5n bytes\f$
// */
//bool construct_lcp_goPHI(tMSS& file_map, const std::string& dir, const std::string& id);
//
////! Our new 2 phases lcp algorithm using usually 1 n bytes.
///*!	The algorithm computes the lcp array and stores it to disk.
// *  \param file_map A map which contains the filenames of previous computed structures (like suffix array, Burrows and Wheeler transform, inverse suffix array, text,...)
// *  \param dir		Directory where the lcp array should be stored.
// *  \param id		Id for the file name of the lcp array.
// *  \par Space complexity:
// *		Usually \f$ n+\Order{1} \f$ bytes, worst case \f$ 5n \f$ bytes
// */
//bool construct_lcp_go2(tMSS& file_map, const std::string& dir, const std::string& id);
//



void lcp_info(tMSS& file_map);

}// end namespace

#endif
