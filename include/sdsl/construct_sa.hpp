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
#include "int_vector_mapper.hpp"

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
template<typename t_int_vec>
void calculate_sa(const unsigned char* c, typename t_int_vec::size_type len, t_int_vec& sa)
{
    typedef typename t_int_vec::size_type size_type;
    constexpr uint8_t t_width = t_int_vec::fixed_int_width;
    if (len <= 1) { // handle special case
        sa.width(1);
        sa.resize(len);
        if ( len > 0 )
            sa[0] = 0;
        return;
    }
    uint8_t sa_width = sa.width();
    bool small_file = (sizeof(len) <= 4 or len < 0x7FFFFFFFULL);
    if (small_file) {
        if (32 == t_width or (0==t_width and 32 >= sa_width)) {

auto sa1 = write_out_mapper<0>::create("@test", len, 32);
divsufsort(c, (int32_t*)(sa1.data()), len);
std::cout<<"hurray"<<std::endl;

            sa.width(32);
            sa.resize(len);
std::cout<<"sa.width()="<<(int)sa.width()<<std::endl;
std::cout<<"sa.size()="<<sa.size()<<" len="<<len<<std::endl;
{
//            int_vector<> tsa(sa.size(), 0, sa.width());
//            divsufsort(c, (int32_t*)tsa.data(), len);
    
              int32_t * p = (int32_t*)sa.data();
              for(int32_t i=0; i<(int32_t)sa.size(); ++i){
                p[i] = 0;  
                if ( p[i] != 0 ) {
                    std::cout<<"ERROR"<<i<<std::endl;
                }
                if ( i<1000)
                    std::cout<<c[i];
              }
              std::cout<<std::endl;
              
    
std::cout<<"start divsufsort"<<std::endl;
              divsufsort(c, p, len);
//              divsufsort(c, (int32_t*)(sa.data()), len);
//            std::copy(tsa.begin(), tsa.end(), sa.begin());
}
std::cout<<"finished divsufsort"<<std::endl;
            // copy integers back to the right positions
            if (sa_width!=32) {
                for (size_type i=0, j=0; i<len; ++i, j+=sa_width) {
                    sa.set_int(j, sa.get_int(i<<5, 32), sa_width);
                }
                sa.width(sa_width);
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
        sa.width(64);
        sa.resize(len);
        divsufsort64(c, (int64_t*)sa.data(), len);
        // copy integers back to the right positions
        if (sa_width!=64) {
            for (size_type i=0; i<len; ++i) {
                sa.set_int(i*sa_width, sa.get_int(i<<6, 64), sa_width);
            }
            sa.width(sa_width);
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
            if ( is_ram_file(cache_file_name(KEY_TEXT, config)) ) {
                read_only_mapper<t_width> text(KEY_TEXT, config);
                std::cout<<"text file size = "<<ram_fs::file_size(cache_file_name(KEY_TEXT, config))<<std::endl;
//                text_type text;
//                load_from_cache(text, KEY_TEXT, config);
                std::cout<<"text.size()="<<text.size()<<std::endl;
                std::cout<<"text.width()="<<(size_t)text.width()<<std::endl;
                auto sa = write_out_mapper<0>::create(cache_file_name(conf::KEY_SA, config),
                                                      text.size(), bits::hi(text.size())+1);
std::cout<<"created write_out_mapper: sa.width="<<(int)sa.width()<<std::endl;
std::cout<<"created write_out_mapper: sa.size()="<<(int)sa.size()<<std::endl;
                // call divsufsort
                algorithm::calculate_sa((const unsigned char*)text.data(), text.size(), sa);
            } else {
                text_type text;
                load_from_cache(text, KEY_TEXT, config);
                auto sa = int_vector<>(text.size(), 0, bits::hi(text.size())+1);
                // call divsufsort
                algorithm::calculate_sa((const unsigned char*)text.data(), text.size(), sa);
                store_to_cache(std::move(sa), conf::KEY_SA, config);
            }
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
