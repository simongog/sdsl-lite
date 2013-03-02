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
/*! \file lcp_construct.hpp
    \brief lcp_construct.hpp contains a space and time efficient construction method for lcp arrays
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_LCP_CONSTRUCT
#define INCLUDED_SDSL_LCP_CONSTRUCT

#include "typedefs.hpp"  // includes definition of tMSS
#include "config.hpp"
#include "int_vector.hpp"
#include "rank_support_v.hpp"
#include "util.hpp"
#include "testutils.hpp"
#include "isa_construct.hpp"
#include "bwt_construct.hpp"
#include "wt_huff.hpp"

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

//! 5n byte variant of the algorithm of Kasai et al. (CPM 2001, "Linear-Time Longest-Common-Prefix Computation in Suffix Arrays and Its Applications")
/*!	The algorithm computes the lcp array and stores it to disk.
 *  \param file_map A map which contains the filenames of previous computed structures (like suffix array, Burrows and Wheeler transform, inverse suffix array, text,...)
 *  \param dir		Directory where the lcp array should be stored.
 *  \param id		Id for the file name of the lcp array.
 *  \par Space complexity
 *		\f$ 5n \f$ bytes
 */
void construct_lcp_kasai(int_vector<>& lcp, cache_config& config);
void construct_int_lcp_kasai(int_vector<>& lcp, const cache_config& config);


//
//
//// semi extern PHI algorithm of Karkainen, Manzini and Puglisi CPM 2009
//bool construct_lcp_semi_extern_PHI(tMSS& file_map, const std::string& dir, const std::string& id);
//
////! 5n byte variant of the algorithm of Kaerkkaeinen et al. (CPM 2009, "Permuted Longest Common Prefix Array")
///*!	The algorithm computes the lcp array and stores it to disk.
// *  \param file_map A map which contains the filenames of previous computed structures (like suffix array, Burrows and Wheeler transform, inverse suffix array, text,...)
// *  \param dir		Directory where the lcp array should be stored.
// *  \param id		Id for the file name of the lcp array.
// *  \param semi_external Boolean flag, which indicates if the algorithm should use only 4n bytes or 5n bytes
// *  \par Space complexity
// *		\f$ 5n \f$ bytes or \f$ 4n \f$ bytes if semi_external is set to true
// */
//bool construct_lcp_PHI(tMSS& file_map, const std::string& dir, const std::string& id, bool semi_external=false);
//
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
