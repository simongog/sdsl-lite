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
/*! \file construct_sa.hpp
    \brief construct_sa.hpp contains an interface to access suffix array construction algorithms
    \author Simon Gog
*/

#ifndef INCLUDED_SDSL_CONSTRUCT_SA
#define INCLUDED_SDSL_CONSTRUCT_SA

#include "config.hpp"
#include "int_vector.hpp"

#include "divsufsort.h"
#include "divsufsort64.h"

#include "qsufsort.hpp"

#include "construct_sa_se.hpp"
#include "construct_config.hpp"

namespace sdsl
{

//! Constructs the Suffix Array (SA) from text over byte-alphabet.
/*! The algorithm constructs the SA and stores it to disk.
 *  \param config Reference to cache configuration
 *  \par Space complexity
 *       Usually less than \f$1.5n \f$ bytes of main memory and
 *       \f$10n \f$ bytes of secondary memory
 *  \pre Text exist in the cache. Keys:
 *         * conf::KEY_TEXT
 *  \post SA exist in the cache. Key
 *         * conf::KEY_SA
 *
 *  This construction method uses less main memory, since data-structures
 *  are only kept in main memory, when random access to them is needed.
 *  Otherwise they are stored on disk. The disk-usage peak of this algorithm
 *  is about 10 times the input.
 *
 *  \par References
 *      [1] T. Beller, M. Zwerger, S. Gog and E. Ohlebusch:
 *          ,,Space-Efficient Construction of the Burrows-Wheeler Transform'',
 *          Proceedings of SPIRE 2013.
 *
 */
void construct_sa_se(cache_config& config);

namespace algorithm
{

//
// Forward declarations
//----------------------------------------------------------

//! Calculates the Suffix Array for a text.
/*!
 * \param c Text (c-string) to calculate the suffix array. The lex. order is given by the ascii-codes of the characters.
 * \param len Length of the text. *(c+len)=0 and for i<len *(c+len)!=0
 * \param sa Reference to a RandomAccessContainer which will contain the result of the calculation.
 * \pre sa.size() has to be equal to len.
 */
template<uint8_t fixedIntWidth>
void calculate_sa(const unsigned char* c, typename int_vector<fixedIntWidth>::size_type len, int_vector<fixedIntWidth>& sa)
{
    typedef typename int_vector<fixedIntWidth>::size_type size_type;
    if (len <= 1) { // handle special case
        sa = int_vector<fixedIntWidth>(len,0);
        return;
    }
    bool small_file = (sizeof(len) <= 4 or len < 0x7FFFFFFFULL);
    if (small_file) {
        uint8_t oldIntWidth = sa.width();
        if (32 == fixedIntWidth or(0==fixedIntWidth and 32 >= oldIntWidth)) {
            sa.width(32);
            sa.resize(len);
            divsufsort(c, (int32_t*)sa.data(), len);
            // copy integers back to the right positions
            if (oldIntWidth!=32) {
                for (size_type i=0; i<len; ++i) {
                    sa.set_int(i*oldIntWidth, sa.get_int(i<<5, 32), oldIntWidth);
                }
                sa.width(oldIntWidth);
                sa.resize(len);
            }
        } else {
            if (sa.width() < bits::hi(len)+1) {
                throw std::logic_error("width of int_vector is to small for the text!!!");
            }
            int_vector<> sufarray(len,0,32);
            divsufsort(c, (int32_t*)sufarray.data(), len);
            for (size_type i=0; i<len; ++i) {
                sa[i] = sufarray[i];
            }
        }
    } else {
        uint8_t oldIntWidth = sa.width();
        sa.width(64);
        sa.resize(len);
        divsufsort64(c, (int64_t*)sa.data(), len);
        // copy integers back to the right positions
        if (oldIntWidth!=64) {
            for (size_type i=0; i<len; ++i) {
                sa.set_int(i*oldIntWidth, sa.get_int(i<<6, 64), oldIntWidth);
            }
            sa.width(oldIntWidth);
            sa.resize(len);
        }
    }
}


} // end namespace algorithm

//! Constructs the Suffix Array (SA) from text over byte- or integer-alphabet.
/*!    The algorithm constructs the SA and stores it to disk.
 *  \tparam t_width Width of the text. 0==integer alphabet, 8=byte alphabet.
 *  \param config    Reference to cache configuration
 *  \par Space complexity
 *      \f$ 5n \f$ byte for t_width=8 and input < 2GB
 *      \f$ 9n \f$ byte for t_width=8 and input > 2GB
 *      \f$ n \log \sigma \f$ bits for t_width=0
 *  \pre Text exist in the cache. Keys:
 *         * conf::KEY_TEXT for t_width=8 or conf::KEY_TEXT_INT for t_width=0
 *  \post SA exist in the cache. Key
 *         * conf::KEY_SA
 *  \par Reference
 *    For t_width=8: DivSufSort (http://code.google.com/p/libdivsufsort/)
 *    For t_width=0: qsufsort (http://www.larsson.dogma.net/qsufsort.c)
 */
template<uint8_t t_width>
void construct_sa(cache_config& config)
{
    static_assert(t_width == 0 or t_width == 8 , "construct_sa: width must be `0` for integer alphabet and `8` for byte alphabet");
    const char* KEY_TEXT = key_text_trait<t_width>::KEY_TEXT;
    if (t_width == 8) {
        if (construct_config::byte_algo_sa == LIBDIVSUFSORT) {
            typedef int_vector<t_width> text_type;
            text_type text;
            load_from_cache(text, KEY_TEXT, config);
            // call divsufsort
            int_vector<> sa(text.size(), 0, bits::hi(text.size())+1);
            algorithm::calculate_sa((const unsigned char*)text.data(), text.size(), sa);
            store_to_cache(sa, conf::KEY_SA, config);
        } else if (construct_config::byte_algo_sa == SE_SAIS) {
            construct_sa_se(config);
        }
    } else if (t_width == 0) {
        // call qsufsort
        int_vector<> sa;
        sdsl::qsufsort::construct_sa(sa, cache_file_name(KEY_TEXT, config).c_str(), 0);
        store_to_cache(sa, conf::KEY_SA, config);
    } else {
        std::cerr << "Unknown alphabet type" << std::endl;
    }
}

} // end namespace sdsl

#endif
