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
/*! \file construct_bwt.hpp
    \brief construct_bwt.hpp contains a space and time efficient construction method for the Burrows and Wheeler Transform (BWT).
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_CONSTRUCT_BWT
#define INCLUDED_SDSL_CONSTRUCT_BWT

#include "typedefs.hpp"
#include "int_vector.hpp"
#include "util.hpp"
#include "config.hpp" // for cache_config

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <list>

namespace sdsl
{

//! Constructs the Burrows and Wheeler Transform (BWT) from text over byte- or integer-alphabet and suffix array.
/*!	The algorithm constructs the BWT and stores it to disk.
 *  \tparam t_width Width of the text. 0==integer alphabet, 8=byte alphabet.
 *  \param config	Reference to cache configuration
 *  \par Space complexity
 *		\f$ n \log \sigma \f$ bits
 *  \pre Text and Suffix array exist in the cache. Keys:
 *         * constants::KEY_TEXT for t_width=8 or constants::KEY_TEXT_INT for t_width=0
 *         * constants::KEY_SA
 *  \post BWT exist in the cache. Key
 *         * constants::KEY_BWT for t_width=8 or constants::KEY_BWT_INT for t_width=0
 */
template<uint8_t t_width>
void construct_bwt(cache_config& config){
    typedef int_vector<>::size_type size_type;
    typedef int_vector<t_width> text_type;
    typedef int_vector<t_width> bwt_type;
    const char * KEY_TEXT = key_text_trait<t_width>::KEY_TEXT;
    const char * KEY_BWT = key_bwt_trait<t_width>::KEY_BWT;

    //  (1) Load text from disk
    write_R_output("bwt", "load text", "begin", 1, 0);
    text_type text;
    util::load_from_cache(text, KEY_TEXT, config);
    size_type n = text.size(); 
    uint8_t bwt_width = text.width();
    write_R_output("bwt", "load text", "end", 1, 0);

    //  (2) Prepare to stream SA from disc and BWT to disc
    write_R_output("bwt", "prepare io", "begin", 1, 0);
    size_type buffer_size = 1000000; // buffer_size is a multiple of 8!
    int_vector_file_buffer<> sa_buf(util::cache_file_name(constants::KEY_SA, config));
    sa_buf.reset(buffer_size);

    bwt_type bwt_buf(buffer_size, 0, bwt_width);

    std::string bwt_file = util::cache_file_name(KEY_BWT, config);
    std::ofstream bwt_out_buf(bwt_file.c_str(), std::ios::binary | std::ios::app | std::ios::out);   // open buffer for bwt
    size_type bit_size = n*bwt_width;
    bwt_out_buf.write((char*) &(bit_size), sizeof(bit_size));	// write size of vector
    if ( t_width != 8) {
        bwt_out_buf.write((char*) &(bwt_width),sizeof(bwt_width));  // write int_width of vector
    }
    write_R_output("bwt", "prepare io", "end", 1, 0);

    //  (3) Construct BWT sequentially by streaming SA and random access to text
    write_R_output("bwt", "construct BWT", "begin", 1, 0);
    size_type wb = 0;  // bytes written into bwt int_vector
    size_type to_add[2] = {(size_type)-1,n-1};
    for (size_type i=0, r_sum=0, r=0; r_sum < n;) {
        for (; i < r_sum+r; ++i) {
            bwt_buf[i-r_sum] = text[ sa_buf[i-r_sum]+to_add[sa_buf[i-r_sum]==0] ]; 
        }
        if (r > 0) {
            size_type cur_wb = (r*bwt_buf.width()+7)/8;
            bwt_out_buf.write((const char*)bwt_buf.data(), cur_wb);
            wb += cur_wb;
        }
        r_sum += r;
        r = sa_buf.load_next_block();
    }
    if (wb%8) {
        bwt_out_buf.write("\0\0\0\0\0\0\0\0", 8-wb%8);
    }
    bwt_out_buf.close();
    util::register_cache_file(KEY_BWT, config);
    write_R_output("bwt", "construct BWT", "end", 1, 0);
}

template <class T_charsize, class T_nsize, class int_vector_type >
void _construct_bwt_from_text(std::string text_filename, uint64_t sigma, int_vector_type &bwt, uint64_t rekursion);

//! Constructs the Burrows and Wheeler Transform (BWT) from text over byte-alphabet.
/*! The algorithm constructs the BWT and stores it to disk.
 *  \param config	Reference to cache configuration
 *  \par Space complexity
 *		\f Usually less than $1.5n \f$ bytes
 *  \pre Text exist in the cache. Keys:
 *         * constants::KEY_TEXT
 *  \post BWT exist in the cache. Key
 *         * constants::KEY_BWT
 */
/*
void construct_bwt_from_text(cache_config& config)
{
	std::string text_filename = util::cache_file_name(constants::KEY_TEXT, config);
	std::string bwt_filename = util::cache_file_name(constants::KEY_BWT, config);
	size_t sigma = 256;
	int_vector<8> bwt;
	// if(text.size() < (1ULL<<32))
	// {
		// _construct_bwt_from_text<uint8_t, uint32_t, int_vector<8> >(text, sigma, bwt);
	// }
	// else
	{
		_construct_bwt_from_text<uint8_t, uint64_t, int_vector<8> >(text_filename, sigma, bwt, 0);
	}
	util::store_to_cache(bwt, constants::KEY_BWT, config);
}

template <class T_charsize, class T_nsize, class int_vector_type >
void _construct_bwt_from_text(std::string text_filename, uint64_t sigma, int_vector_type &bwt, uint64_t rekursion)
{
	int_vector_type text;
	util::load_from_file(text, text_filename);
	uint64_t n = text.size();
	uint64_t buffersize = 1024*1024/8;
	uint8_t int_width = bit_magic::l1BP(n-1)+1;
	#if DEBUGLEVEL > 0
	std::cout << std::string(rekursion, '\t') << "recursion level=" << rekursion << " n=" << n << " text.int_width=" << (uint64_t)text.get_int_width() << " sigma=" << sigma << std::endl;
	uint64_t ts_begin = get_current_timestamp();
	uint64_t lasttime = ts_begin;
	#endif

	#if DEBUGLEVEL > 2
	std::cout << std::string(rekursion, '\t') << "Schritt 1 - Scan Text from right to left and count LMS, S and L characters" << std::endl;
	#endif

	// define variable
	size_t first_lms_pos=0;
	size_t number_of_lms_strings = 0;
	int_vector<> bkt_lms(sigma, 0, int_width);
	int_vector<> C(sigma, 0, int_width);
	size_t bkt_s_last = 0, bkt_s_sum=0, bkt_s_bound_char=0, bound_s=0;
	int_vector<> bkt_s(sigma, 0, int_width);
	size_t bkt_l_sum = 0, bkt_l_bound_char=0, bound_l=0;
	size_t parts = 0;
	if(rekursion == 0)
	{
		parts = 10;
		std::cerr << "parts=" << parts << "." << std::endl;
	}

	int_vector<> bkt_l(sigma, 0, int_width);
	bkt_s[text[n-1]] = 1;
	C[text[n-1]] = C[text[n-1]]+1;
	for(size_t i=n-2, was_s_typ = 1; i+1>0; --i)
	{
		C[text[i]] = C[text[i]]+1;
		if(text[i]>text[i+1])
		{
			if(was_s_typ)
			{
				bkt_lms[text[i+1]] = bkt_lms[text[i+1]] + 1;
				++number_of_lms_strings;
				first_lms_pos = i+1;
				was_s_typ = 0;
			}
		}
		else if(text[i]<text[i+1])
		{
			was_s_typ = 1;
		}

		if(was_s_typ)
		{
			bkt_s[text[i]] = bkt_s[text[i]]+1;
		}
		else
		{
			bkt_l[text[i]] = bkt_l[text[i]]+1;
		}
	}
	for(size_t i=1; i<C.size(); ++i)
	{
		C[i] = C[i]+C[i-1];
	}
	print_needed_time(1, lasttime, rekursion);

	#if DEBUGLEVEL > 2
	std::cout << std::string(rekursion, '\t') << "Schritt 2 - Scan Text from right to left and detect LMS-Positions. Sort and write them to disk" << std::endl;
	#endif
	cached_external_array<T_nsize> right("right_"+util::to_string(rekursion), buffersize, false, true);
	size_t right_pointer=0;
	cached_external_array<T_nsize> left("left_"+util::to_string(rekursion), buffersize, false, true);
	size_t left_pointer=0;
	{
		for(size_t i=0, tmp2=0, tmp=0; i<sigma; ++i)
		{
			tmp += bkt_lms[i];
			bkt_lms[i] = tmp2;
			tmp2 = tmp;
		}
		if(rekursion==0)
		{
			cached_external_array<T_nsize> lms_positions("lms_positions_"+util::to_string(rekursion), buffersize, true, true);
			size_t lms_positions_pointer=0;
			for(size_t i=n-2, was_s_typ = 1; i+1>0; --i)
			{
				if(text[i]>text[i+1])
				{
					if(was_s_typ)
					{
						lms_positions.write(lms_positions_pointer++, bkt_lms[text[i+1]]);
						lms_positions.write(lms_positions_pointer++, i+1);
						bkt_lms[text[i+1]] = bkt_lms[text[i+1]] + 1;
						was_s_typ = 0;
					}
				}
				else if(text[i]<text[i+1])
				{
					was_s_typ = 1;
				}
			}
			util::clear(text);
			{ 
				// Order lms_positions according to first character
				int_vector<> lms_strings(number_of_lms_strings, 0, int_width);
				for(size_t i=0; i<lms_positions_pointer; )
				{
					size_t idx = lms_positions.read(i++);
					size_t val = lms_positions.read(i++);
					lms_strings[idx] = val;
				}
				// Store it to file
				for(left_pointer=0; left_pointer<number_of_lms_strings; ++left_pointer)
				{
					left.write(left_pointer, lms_strings[number_of_lms_strings-left_pointer-1]);
				}
			}
			util::load_from_file(text, text_filename);
		}
		else
		{
			int_vector<> lms_strings(number_of_lms_strings, 0, int_width);
			for(size_t i=n-2, was_s_typ = 1; i+1>0; --i)
			{
				if(text[i]>text[i+1])
				{
					if(was_s_typ)
					{
						lms_strings[bkt_lms[text[i+1]]] = i+1;
						bkt_lms[text[i+1]] = bkt_lms[text[i+1]] + 1;
						was_s_typ = 0;
					}
				}
				else if(text[i]<text[i+1])
				{
					was_s_typ = 1;
				}
			}
			for(left_pointer=0; left_pointer<number_of_lms_strings; ++left_pointer)
			{
				left.write(left_pointer, lms_strings[number_of_lms_strings-left_pointer-1]);
			}
		}
	}
	left_pointer--;
	print_needed_time(2, lasttime, rekursion);

	#if DEBUGLEVEL > 2
	std::cout << std::string(rekursion, '\t') << "Schritt 3 - Scan virtual array from left to right" << std::endl;
	#endif
	print_R_line(util::to_string(rekursion), "blue");
	{
		// prepare bkt_lms and backup it into bkt_lms
		for(size_t i=0, tmp=0; i<sigma; ++i)
		{
			tmp = bkt_l[i];
			bkt_l[i] = bkt_l_sum;
			bkt_l_sum += tmp;
			bkt_lms[i] = bkt_l[i];
		}
		bound_l = bkt_l_sum;

		// determine splitting parameteres
		bkt_l_bound_char = sigma;
		for(size_t i=0; i<sigma; ++i)
		{
			if(bkt_l[i] > bkt_l_sum/2)
			{
				bkt_l_bound_char = i;
				bkt_l_sum = bkt_l[i];
				break;
			}
		}

		if(parts > 1)
		{
			size_t partsize = bound_l/parts+1;

			int_vector<> array(partsize, 0, int_width);
			vector< cached_external_array<T_nsize> > cached_array(parts-1);
			vector<uint64_t> cached_array_pointer(cached_array.size());
			for(size_t i=1; i<parts; ++i)
			{
				cached_array[i-1].init("rightbuffer_"+ util::to_string(i-1) +"_"+util::to_string(rekursion), buffersize, true, true);
				cached_array_pointer[i-1] = 0;
			}

			for(size_t c=0, pos=0, offset=0; c<sigma; ++c)
			{
				// begin with array
				for( ; pos<bkt_l[c]; ++pos)
				{
					// Load lazy values
					if(pos-offset >= partsize)
					{
						offset += partsize;
						for(size_t i=0, cur_part=pos/partsize-1; i<cached_array_pointer[cur_part]/2; ++i)
						{
							size_t src = cached_array[cur_part].read(2*i);
							size_t val = cached_array[cur_part].read(2*i+1);
							array[src-offset] = val;
						}
						cached_array[pos/partsize-1].reset();
					}

					size_t idx = array[pos-offset];
					if(idx == 0)
					{
						right.write(right_pointer++, idx);
					}
					else
					{
						size_t symbol = text[idx-1];
						if(symbol >= c)
						{
							size_t val = idx-1;
							size_t src = bkt_l[symbol];
							bkt_l[symbol] = bkt_l[symbol] + 1;
							if( (src-offset)/partsize == 0)
							{
								array[src-offset] = val;
							}
							else
							{
								size_t part = src/partsize-1;
								cached_array[part].write(cached_array_pointer[part]++, src);
								cached_array[part].write(cached_array_pointer[part]++, val);
							}
						}
						else
						{
							right.write(right_pointer++, idx);
						}
					}
				}
				// continue with stack
				while(left_pointer < number_of_lms_strings and text[left.read(left_pointer)] == c)
				{
					size_t idx = left.read(left_pointer--);
					if(idx == 0)
					{
						idx = n;
					}
					--idx;
					size_t symbol = text[idx];

					size_t val = idx;
					size_t src = bkt_l[symbol];
					bkt_l[symbol] = bkt_l[symbol] + 1;
					if( (src-offset)/partsize == 0)
					{
						array[src-offset] = val;
					}
					else
					{
						size_t part = src/partsize-1;
						cached_array[part].write(cached_array_pointer[part]++, src);
						cached_array[part].write(cached_array_pointer[part]++, val);
					}
				}
			}
		}
		else
		{
			int_vector<> array(bkt_l_sum, 0, int_width);
			cached_external_array<T_nsize> external_storage("rightbuffer_"+util::to_string(rekursion), buffersize, true, true);
			for(size_t c=0, pos=0, offset=0, written=0; c<sigma; ++c)
			{
				// Load lazy values
				if(c == bkt_l_bound_char)
				{
					offset = bkt_l_sum;
					for(size_t i=0; i<written; ++i)
					{
						size_t idx = external_storage.read(2*i);
						size_t val = external_storage.read(2*i+1);
						array[idx-offset] = val;
					}
				}

				// begin with array
				for( ; pos<bkt_l[c]; ++pos)
				{
					size_t idx = array[pos-offset];
					if(idx == 0)
					{
						right.write(right_pointer++, idx);
					}
					else
					{
						size_t symbol = text[idx-1];
						if(symbol >= c)
						{
							if(bkt_l[symbol]-offset < array.size())
							{
								array[bkt_l[symbol]-offset] = idx-1;
							}
							else
							{
								external_storage.write(2*written, bkt_l[symbol]);
								external_storage.write(2*written+1, idx-1);
								++written;
							}
							bkt_l[symbol] = bkt_l[symbol] + 1;
						}
						else
						{
							right.write(right_pointer++, idx);
						}
					}
				}
				// continue with stack
				while(left_pointer < number_of_lms_strings and text[left.read(left_pointer)] == c)
				{
					size_t idx = left.read(left_pointer--);
					if(idx == 0)
					{
						idx = n;
					}
					--idx;
					size_t symbol = text[idx];
					// array[bkt_l[symbol]] = idx;
					if(bkt_l[symbol]-offset < array.size())
					{
						array[bkt_l[symbol]-offset] = idx;
					}
					else
					{
						external_storage.write(2*written, bkt_l[symbol]);
						external_storage.write(2*written+1, idx);
						++written;
					}
					bkt_l[symbol] = bkt_l[symbol] + 1;
				}
			}
		}

		// Restore bkt_l from bkt_lms
		for(size_t i=0; i<sigma; ++i)
		{
			bkt_l[i] = bkt_lms[i];
		}
	}
	right_pointer--;
	print_needed_time(3, lasttime, rekursion);

	#if DEBUGLEVEL > 2
	std::cout << std::string(rekursion, '\t') << "Schritt 4 - Scan virtual array from right to left" << std::endl;
	#endif
	left_pointer = 0;
	left.reset();
	{
		// NEW

		// Prepare bkt_s and backup it into bkt_lms
		bkt_s_last=0, bkt_s_sum=0;
		for(size_t i=0; i<sigma; ++i)
		{
			bkt_s_sum += bkt_s[i];
			if(bkt_s[i])
			{
				bkt_s[i] = bkt_s_sum;
				bkt_s_last = bkt_s_sum;
			}
			else
			{
				bkt_s[i] = bkt_s_sum;
			}
			bkt_lms[i] = bkt_s[i];
		}
		bound_s = bkt_s_sum;

		// determine splitting parameters
		for(size_t i=0; i<sigma; ++i)
		{
			if(bkt_s[i] > bkt_s_sum/2)
			{
				bkt_s_bound_char = i;
				bkt_s_sum = bkt_s[i];
				break;
			}
		}

		if(parts > 1)
		{
			size_t partsize = bound_s/parts+1;
			int_vector<> array(partsize, 0, int_width);
			vector< cached_external_array<T_nsize> > cached_array(parts-1);
			vector<uint64_t> cached_array_pointer(cached_array.size());
			for(size_t i=0; i<parts-1; ++i)
			{
				cached_array[i].init("rightbuffer_"+ util::to_string(i) +"_"+util::to_string(rekursion), buffersize, true, true);
				cached_array_pointer[i] = 0;
			}

			for(size_t c=sigma-1, pos=bkt_s_last-1, offset=partsize*(parts-1); c<sigma; --c)
			{
				// begin with array
				for(; pos+1 > bkt_s[c]; --pos)
				{
					while(pos < offset)
					{
						// Load buffered values
						offset -= partsize;
						// offset = pos-partsize;
						for(size_t i=0, cur_part=offset/partsize; i<cached_array_pointer[cur_part]/2; ++i)
						{
							size_t src = cached_array[cur_part].read(2*i);
							size_t val = cached_array[cur_part].read(2*i+1);
							array[src-offset] = val;
						}
						cached_array[offset/partsize].reset();
					}

					size_t idx = array[pos-offset];
					if(idx==0)
					{
						idx = n;
					}
					--idx;
					size_t symbol = text[idx];
					if(symbol <= c)
					{
						bkt_s[symbol] = bkt_s[symbol] - 1;
						size_t val = idx;
						size_t src = bkt_s[symbol];
						if( src >= offset)
						{
							array[src-offset] = val;
						}
						else
						{
							size_t part = src/partsize;
							cached_array[part].write(cached_array_pointer[part]++, src);
							cached_array[part].write(cached_array_pointer[part]++, val);
						}
					}
					else
					{
						left.write(left_pointer++, array[pos-offset]);
					}
				}

				// continue with stack
				while(right_pointer < number_of_lms_strings and text[right.read(right_pointer)] == c)
				{
					size_t idx = right.read(right_pointer--);
					if(idx == 0)
					{
						idx = n;
					}
					--idx;
					size_t symbol = text[idx];
					bkt_s[symbol] = bkt_s[symbol] - 1;

					size_t val = idx;
					size_t src = bkt_s[symbol];
					if( src >= offset)
					{
						array[src-offset] = val;
					}
					else
					{
						size_t part = src/partsize;
						cached_array[part].write(cached_array_pointer[part]++, src);
						cached_array[part].write(cached_array_pointer[part]++, val);
					}
				}
			}
		}
		else
		{
			int_vector<> array(bkt_s_sum, 0, int_width);
			cached_external_array<T_nsize> external_storage("leftbuffer_"+util::to_string(rekursion), buffersize, true, true);
			for(size_t c=sigma-1, pos=bkt_s_last-1, offset=bkt_s_sum, written=0; c<sigma; --c)
			{
				// Load lazy values
				if(c == bkt_s_bound_char)
				{
					offset = 0;
					for(size_t i=0; i<written; ++i)
					{
						size_t idx = external_storage.read(2*i);
						size_t val = external_storage.read(2*i+1);
						array[idx-offset] = val;
					}
				}

				// begin with array
				for(; pos+1 > bkt_s[c]; --pos)
				{
					size_t idx = array[pos-offset];
					if(idx==0)
					{
						idx = n;
					}
					--idx;
					size_t symbol = text[idx];
					if(symbol <= c)
					{
						bkt_s[symbol] = bkt_s[symbol] - 1;
						if(bkt_s[symbol] >= offset)
						{
							array[bkt_s[symbol]-offset] = idx;
						}
						else
						{
							external_storage.write(2*written, bkt_s[symbol]);
							external_storage.write(2*written+1, idx);
							++written;
						}
					}
					else
					{
						left.write(left_pointer++, array[pos-offset]);
					}
				}
				while(right_pointer < number_of_lms_strings and text[right.read(right_pointer)] == c)
				{
					size_t idx = right.read(right_pointer--);
					if(idx == 0)
					{
						idx = n;
					}
					--idx;
					size_t symbol = text[idx];
					bkt_s[symbol] = bkt_s[symbol] - 1;
					if(bkt_s[symbol] >= offset)
					{
						array[bkt_s[symbol]-offset] = idx;
					}
					else
					{
						external_storage.write(2*written, bkt_s[symbol]);
						external_storage.write(2*written+1, idx);
						++written;
					}
				}
			}
		}

		// Restore bkt_s from bkt_lms
		for(size_t i=0; i<sigma; ++i)
		{
			bkt_s[i] = bkt_lms[i];
		}
	}
	print_R_line(util::to_string(rekursion), "blue");
	right.reset();
	right_pointer = 0;
	--left_pointer;
	print_needed_time(4, lasttime, rekursion);

	#if DEBUGLEVEL > 2
	std::cout << std::string(rekursion, '\t') << "Schritt 5 - Detect same lms-Strings, write text to file" << std::endl;
	#endif
	int_vector<1> same_lms(number_of_lms_strings, false);

	size_t last_end_pos = first_lms_pos, ordnung = number_of_lms_strings-1;
	same_lms[number_of_lms_strings-1] = true;
	for(size_t i=number_of_lms_strings-2, a=0, b=0, last_a=left.read(number_of_lms_strings-1); i<number_of_lms_strings; --i)
	{
		b = last_a;
		a = left.read(i);
		last_a = a;

		size_t end_pos = get_next_lms_position(text, a);
		if(end_pos-a == last_end_pos-b)
		{
			while(a < end_pos and text[a] == text[b])
			{
				++a;
				++b;
			}
			if(text[a] == text[b])
			{
				same_lms[i] = true;
				--ordnung;
			}
		}
		last_end_pos = end_pos;
	}
	util::clear(text);
	print_needed_time(5, lasttime, rekursion);

	#if DEBUGLEVEL > 2
	std::cout << std::string(rekursion, '\t') << "Schritt 6 - create renaming table and store to file" << std::endl;
	#endif
	std::string filename_translation = "translation_"+util::to_string(rekursion);
	{
		int_vector<> name_to_position(ordnung+1, 0, int_width);
		ordnung = 0;
		for(size_t i=number_of_lms_strings-1; i<number_of_lms_strings; --i)
		{
			if(!same_lms[i])
			{
				++ordnung;
			}
			name_to_position[ordnung] = left.read(i);
		}
		util::store_to_file(name_to_position, filename_translation);
	}
	print_needed_time(6, lasttime, rekursion);

	#if DEBUGLEVEL > 2
	std::cout << std::string(rekursion, '\t') << "Schritt 7 - create renamed string" << std::endl;
	#endif
	print_R_line(util::to_string(rekursion), "red");
	int_vector<> text_rek;
	if(rekursion==0)
	{
		text_rek.set_int_width((bit_magic::l1BP(ordnung+1)+1));
	}
	else
	{
		text_rek.set_int_width((bit_magic::l1BP(number_of_lms_strings+1)+1));
	}
	text_rek.resize(number_of_lms_strings);
std::cerr << "rekursion=" << rekursion << " text_rek.size()=" << text_rek.size() << " text_rek.get_int_width()=" << (uint64_t)text_rek.get_int_width() << std::endl;
	util::set_zero_bits(text_rek);
	{
		if(rekursion==0)
		{
			text_rek.resize(n/2+1);
			util::set_zero_bits(text_rek);
			ordnung = 0;
			for(size_t i=number_of_lms_strings-1; i<number_of_lms_strings; --i)
			{
				if(!same_lms[i])
				{
					++ordnung;
				}
				text_rek[left.read(left_pointer--)/2] = ordnung;
			}
			for(size_t i=0, pos=0; i<text_rek.size(); ++i)
			{
				if(text_rek[i]>0)
				{
					text_rek[pos++] = text_rek[i];
				}
			}
			text_rek[number_of_lms_strings-1] = 0;
			text_rek.resize(number_of_lms_strings);
		}
		else
		{
			int_vector<> names_in_correct_order(n/2+1, 0, (bit_magic::l1BP(ordnung)+1));
			ordnung = 0;
			for(size_t i=number_of_lms_strings-1; i<number_of_lms_strings; --i)
			{
				if(!same_lms[i])
				{
					++ordnung;
				}
				names_in_correct_order[left.read(left_pointer--)/2] = ordnung;
			}
			for(size_t i=0, pos=0; i<names_in_correct_order.size(); ++i)
			{
				if(names_in_correct_order[i]>0)
				{
					text_rek[pos++] = names_in_correct_order[i];
				}
			}
		}
	}
	util::clear(same_lms);
	left.reset();
	print_needed_time(7, lasttime, rekursion);

	#if DEBUGLEVEL > 2
	std::cout << std::string(rekursion, '\t') << "Schritt 8 - Determine complete LMS order (recursivly)" << std::endl;
	#endif
	int_vector<> bwt_rek;
	bwt_rek.set_int_width(text_rek.get_int_width());
	if(text_rek.size() > ordnung+1)
	{
		if(rekursion==0)
		{
			print_R_line(util::to_string(rekursion), "gray");
			std::string text_filename_rek = "text_rek_"+util::to_string(rekursion);
			util::store_to_file(text_rek, text_filename_rek);
			uint64_t n_rek = text_rek.size();
			util::clear(text_rek);
			if(ordnung+1 < (1ULL<<8))
			{
				if(n_rek < (1ULL<<32))
				{
					_construct_bwt_from_text<uint8_t, uint32_t, int_vector<> >(text_filename_rek, ordnung+1, bwt_rek, rekursion+1);
				}
				else
				{
					_construct_bwt_from_text<uint8_t, uint64_t, int_vector<> >(text_filename_rek, ordnung+1, bwt_rek, rekursion+1);
				}
			}
			else if(ordnung+1 < (1ULL<<16))
			{
				if(n_rek < (1ULL<<32))
				{
					_construct_bwt_from_text<uint16_t, uint32_t, int_vector<> >(text_filename_rek, ordnung+1, bwt_rek, rekursion+1);
				}
				else
				{
					_construct_bwt_from_text<uint16_t, uint64_t, int_vector<> >(text_filename_rek, ordnung+1, bwt_rek, rekursion+1);
				}
			}
			else if(ordnung+1 < (1ULL<<32))
			{
				if(n_rek < (1ULL<<32))
				{
					_construct_bwt_from_text<uint32_t, uint32_t, int_vector<> >(text_filename_rek, ordnung+1, bwt_rek, rekursion+1);
				}
				else
				{
					_construct_bwt_from_text<uint32_t, uint64_t, int_vector<> >(text_filename_rek, ordnung+1, bwt_rek, rekursion+1);
				}
			}
			else
			{
				_construct_bwt_from_text<uint64_t, uint64_t, int_vector<> >(text_filename_rek, ordnung+1, bwt_rek, rekursion+1);
			}
			remove(text_filename_rek.c_str());
			print_R_line(util::to_string(rekursion), "gray");
		}
		else
		{
			text_rek.resize(text_rek.size()*2);
			for(size_t i=0; i<number_of_lms_strings; ++i)
			{
				text_rek[number_of_lms_strings+i] = text_rek[i];
				text_rek[i] = 0;
			}

			print_R_line(util::to_string(rekursion), "green");
			construct_sa_simple(text_rek, text_rek, number_of_lms_strings, number_of_lms_strings, ordnung+1, rekursion+1);
			print_R_line(util::to_string(rekursion), "green");
			for(size_t i=0; i<number_of_lms_strings; ++i)
			{
				if(text_rek[i]>0)
				{
					text_rek[i] = text_rek[number_of_lms_strings+text_rek[i]-1];
				}
				else
				{
					text_rek[i] = 0;
				}
			}
			text_rek.resize(number_of_lms_strings);
			bwt_rek.resize(number_of_lms_strings);
			for(size_t i=0; i<number_of_lms_strings; ++i)
			{
				bwt_rek[i] = text_rek[i];
			}
		}
	}
	else
	{
		bwt_rek.resize(text_rek.size());
		bwt_rek[text_rek[0]] = text_rek[text_rek.size()-1];
		for(size_t i=1; i<text_rek.size(); ++i)
		{
			bwt_rek[text_rek[i]] = text_rek[i-1];
		}
	}
	util::clear(text_rek);
	print_needed_time(8, lasttime, rekursion);

	#if DEBUGLEVEL > 2
	std::cout << std::string(rekursion, '\t') << "Schritt 9 - retranslate bwt" << std::endl;
	#endif
	{
		int_vector<> name_to_position;
		util::load_from_file(name_to_position, filename_translation);
		remove(filename_translation.c_str());
		right_pointer = 0;
		for(size_t i=0; i<bwt_rek.size(); ++i)
		{
			right.write(right_pointer++, name_to_position[bwt_rek[i]]);
		}
	}
	print_needed_time(9, lasttime, rekursion);

	#if DEBUGLEVEL > 2
	std::cout << std::string(rekursion, '\t') << "Schritt 10 - Reload text" << std::endl;
	#endif
	util::clear(bwt_rek);
	util::load_from_file(text, text_filename);
	print_needed_time(10, lasttime, rekursion);

	// ToDo: Do this already in Step 6???
	#if DEBUGLEVEL > 2
	std::cout << std::string(rekursion, '\t') << "Schritt 11 - " << std::endl;
	#endif
	left_pointer = 0;
	while(right_pointer>0)
	{
		size_t idx = right.read(--right_pointer);
		if(idx == n-1)
		{
			idx = first_lms_pos;
		}
		else
		{
			idx = get_next_lms_position(text, idx);
		}
		left.write(left_pointer++, idx);
	}
	left_pointer--;
	print_needed_time(11, lasttime, rekursion);

	#if DEBUGLEVEL > 2
	std::cout << std::string(rekursion, '\t') << "Schritt 12 - Scan virtual array from left to right second time" << std::endl;
	#endif
	right.reset();
	right_pointer = 0;
	cached_external_array<T_charsize> charstack("charstack_"+util::to_string(rekursion), buffersize, false, true);
	size_t charstack_pointer = 0;
	cached_external_array<T_charsize> cached_bwt("bwt_"+util::to_string(rekursion), buffersize, false, true);
	size_t bwt_pointer = 0;
	print_R_line(util::to_string(rekursion), "blue");
	{
		if(parts>1)
		{
			size_t partsize = bound_l/parts+1;

			int_vector<> array(partsize, 0, int_width);
			vector< cached_external_array<T_nsize> > cached_array(parts-1);
			vector<uint64_t> cached_array_pointer(cached_array.size());
			for(size_t i=1; i<parts; ++i)
			{
				cached_array[i-1].init("rightbuffer_"+ util::to_string(i-1) +"_"+util::to_string(rekursion), buffersize, true, true);
				cached_array_pointer[i-1] = 0;
			}

			for(size_t c=0, pos=0, offset=0; c<sigma; ++c)
			{
				// begin with array
				for( ; pos<bkt_l[c]; ++pos)
				{
					// Load lazy values
					if(pos-offset >= partsize)
					{
						offset += partsize;
						for(size_t i=0, cur_part=pos/partsize-1; i<cached_array_pointer[cur_part]/2; ++i)
						{
							size_t src = cached_array[cur_part].read(2*i);
							size_t val = cached_array[cur_part].read(2*i+1);
							array[src-offset] = val;
						}
						cached_array[pos/partsize-1].reset();
					}

					size_t idx = array[pos-offset];
					if(idx == 0)
					{
						cached_bwt.write(bwt_pointer++, 0);
						right.write(right_pointer++, idx);
					}
					else
					{
						size_t symbol = text[idx-1];
						cached_bwt.write(bwt_pointer++, symbol);
						if(symbol >= c)
						{
							size_t val = idx-1;
							size_t src = bkt_l[symbol];
							bkt_l[symbol] = bkt_l[symbol] + 1;
							if( (src-offset)/partsize == 0)
							{
								array[src-offset] = val;
							}
							else
							{
								size_t part = src/partsize-1;
								cached_array[part].write(cached_array_pointer[part]++, src);
								cached_array[part].write(cached_array_pointer[part]++, val);
							}
						}
						else
						{
							right.write(right_pointer++, idx);
						}
					}
				}
				bwt_pointer = C[c];

				// continue with stack
				while(left_pointer < number_of_lms_strings and text[left.read(left_pointer)] == c)
				{
					size_t idx = left.read(left_pointer--);
					if(idx == 0)
					{
						idx = n;
					}
					--idx;
					size_t symbol = text[idx];
					charstack.write(charstack_pointer++, symbol);

					size_t val = idx;
					size_t src = bkt_l[symbol];
					bkt_l[symbol] = bkt_l[symbol] + 1;
					if( (src-offset)/partsize == 0)
					{
						array[src-offset] = val;
					}
					else
					{
						size_t part = src/partsize-1;
						cached_array[part].write(cached_array_pointer[part]++, src);
						cached_array[part].write(cached_array_pointer[part]++, val);
					}
				}
			}
		}
		else
		{
			int_vector<> array(bkt_l_sum, 0, int_width);
			cached_external_array<T_nsize> external_storage("rightbuffer_"+util::to_string(rekursion), buffersize, true, true);
			for(size_t c=0, pos=0, offset=0, written=0; c<sigma; ++c)
			{
				// Load lazy values
				if(c == bkt_l_bound_char)
				{
					offset = bkt_l_sum;
					for(size_t i=0; i<written; ++i)
					{
						size_t idx = external_storage.read(2*i);
						size_t val = external_storage.read(2*i+1);
						array[idx-offset] = val;
					}
				}

				// begin with array
				for( ; pos<bkt_l[c]; ++pos)
				{
					size_t idx = array[pos-offset];
					if(idx == 0)
					{
						cached_bwt.write(bwt_pointer++, 0);
						right.write(right_pointer++, idx);
					}
					else
					{
						size_t symbol = text[idx-1];
						cached_bwt.write(bwt_pointer++, symbol);
						if(symbol >= c)
						{
							if(bkt_l[symbol]-offset < array.size())
							{
								array[bkt_l[symbol]-offset] = idx-1;
							}
							else
							{
								external_storage.write(2*written, bkt_l[symbol]);
								external_storage.write(2*written+1, idx-1);
								++written;
							}
							bkt_l[symbol] = bkt_l[symbol] + 1;
						}
						else
						{
							right.write(right_pointer++, idx);
						}
					}
				}
				bwt_pointer = C[c];

				// continue with stack
				while(left_pointer < number_of_lms_strings and text[left.read(left_pointer)] == c)
				{
					size_t idx = left.read(left_pointer--);
					if(idx == 0)
					{
						idx = n;
					}
					--idx;
					size_t symbol = text[idx];
					charstack.write(charstack_pointer++, symbol);
					if(bkt_l[symbol]-offset < array.size())
					{
						array[bkt_l[symbol]-offset] = idx;
					}
					else
					{
						external_storage.write(2*written, bkt_l[symbol]);
						external_storage.write(2*written+1, idx);
						++written;
					}
					bkt_l[symbol] = bkt_l[symbol] + 1;
				}
			}
		}
	}
	left.reset();
	charstack_pointer--;
	right_pointer--;
	print_needed_time(12, lasttime, rekursion);

	#if DEBUGLEVEL > 2
	std::cout << std::string(rekursion, '\t') << "Schritt 13 - Scan virtual array from right to left second time " << std::endl;
	#endif
	cached_bwt.set_write_direction(false);
	{
		if(parts > 1)
		{
			size_t partsize = bound_s/parts+1;

			int_vector<> array(partsize, 0, int_width);
			vector< cached_external_array<T_nsize> > cached_array(parts-1);
			vector<uint64_t> cached_array_pointer(cached_array.size());
			for(size_t i=0; i<parts-1; ++i)
			{
				cached_array[i].init("leftbuffer_"+ util::to_string(i) +"_"+util::to_string(rekursion), buffersize, true, true);
				cached_array_pointer[i] = 0;
			}

			for(size_t c=sigma-1, pos=bkt_s_last-1, offset=partsize*(parts-1); c<sigma; --c)
			{
				// begin with array
				assert(c < C.size() );
				bwt_pointer = C[c]-1;
				for(; pos+1 > bkt_s[c]; --pos)
				{
					while(pos < offset)
					{
						// Load buffered values
						offset -= partsize;
						for(size_t i=0, cur_part=offset/partsize; i<cached_array_pointer[cur_part]/2; ++i)
						{
							size_t src = cached_array[cur_part].read(2*i);
							size_t val = cached_array[cur_part].read(2*i+1);
							assert( (src-offset) < array.size() );
							array[src-offset] = val;
						}
						assert( (offset/partsize) < cached_array.size() );
						cached_array[offset/partsize].reset();
					}

					assert( (pos-offset) < array.size() );
					size_t idx = array[pos-offset];
					if(idx==0)
					{
						idx = n;
					}
					--idx;
					assert( (idx) < text.size() );
					size_t symbol = text[idx];
					if(symbol <= c)
					{
						cached_bwt.write(bwt_pointer--, symbol);
						assert( (symbol) < bkt_s.size() );
						bkt_s[symbol] = bkt_s[symbol] - 1;
						size_t val = idx;
						size_t src = bkt_s[symbol];
						if( src >= offset)
						{
							assert( (src-offset) < array.size() );
							array[src-offset] = val;
						}
						else
						{
							size_t part = src/partsize;
							assert( part < cached_array.size() );
							cached_array[part].write(cached_array_pointer[part]++, src);
							cached_array[part].write(cached_array_pointer[part]++, val);
						}
					}
					else
					{
						cached_bwt.write(bwt_pointer--, charstack.read(charstack_pointer--));
					}
				}

				// continue with stack
				while(right_pointer < number_of_lms_strings and text[right.read(right_pointer)] == c)
				{
					size_t idx = right.read(right_pointer--);
					if(idx == 0)
					{
						idx = n;
					}
					--idx;
					size_t symbol = text[idx];
					assert( (symbol) < bkt_s.size() );
					bkt_s[symbol] = bkt_s[symbol] - 1;

					size_t val = idx;
					size_t src = bkt_s[symbol];
					if( src >= offset)
					{
						assert( (src-offset) < array.size() );
						array[src-offset] = val;
					}
					else
					{
						size_t part = src/partsize;
						assert( (part) < cached_array.size() );
						cached_array[part].write(cached_array_pointer[part]++, src);
						cached_array[part].write(cached_array_pointer[part]++, val);
					}
				}
			}
		}
		else
		{
			int_vector<> array(bkt_s_sum, 0, int_width);
			cached_external_array<T_nsize> external_storage("leftbuffer_"+util::to_string(rekursion), buffersize, true, true);
			for(size_t c=sigma-1, pos=bkt_s_last-1, offset=bkt_s_sum, written=0; c<sigma; --c)
			{
				// Load lazy values
				if(c == bkt_s_bound_char)
				{
					offset = 0;
					for(size_t i=0; i<written; ++i)
					{
						size_t idx = external_storage.read(2*i);
						size_t val = external_storage.read(2*i+1);
						array[idx-offset] = val;
					}
				}

				// begin with array
				bwt_pointer = C[c]-1;
				for(; pos+1 > bkt_s[c]; --pos)
				{
					size_t idx = array[pos-offset];
					if(idx==0)
					{
						idx = n;
					}
					--idx;
					size_t symbol = text[idx];
					if(symbol <= c)
					{
						cached_bwt.write(bwt_pointer--, symbol);
						bkt_s[symbol] = bkt_s[symbol] - 1;
						if(bkt_s[symbol] >= offset)
						{
							array[bkt_s[symbol]-offset] = idx;
						}
						else
						{
							external_storage.write(2*written, bkt_s[symbol]);
							external_storage.write(2*written+1, idx);
							++written;
						}
					}
					else
					{
						cached_bwt.write(bwt_pointer--, charstack.read(charstack_pointer--));
					}
				}

				// continue with stack
				while(right_pointer < number_of_lms_strings and text[right.read(right_pointer)] == c)
				{
					size_t idx = right.read(right_pointer--);
					if(idx == 0)
					{
						idx = n;
					}
					--idx;
					size_t symbol = text[idx];
					bkt_s[symbol] = bkt_s[symbol] - 1;
					if(bkt_s[symbol] >= offset)
					{
						array[bkt_s[symbol]-offset] = idx;
					}
					else
					{
						external_storage.write(2*written, bkt_s[symbol]);
						external_storage.write(2*written+1, idx);
						++written;
					}
				}
			}
		}
	}
	print_R_line(util::to_string(rekursion), "blue");
	print_needed_time(13, lasttime, rekursion);

	#if DEBUGLEVEL > 2
	std::cout << std::string(rekursion, '\t') << "Schritt 14 - Load bwt from file to memory" << std::endl;
	#endif
	cached_bwt.set_read_direction(true);
	text.swap(bwt);
	for(size_t i=0; i<n; ++i)
	{
		bwt[i] = cached_bwt.read(i);
	}
	print_needed_time(14, lasttime, rekursion);

	#if DEBUGLEVEL > 0
	std::cout << std::string(rekursion, '\t') << "recursion level=" << rekursion << " time=" << (get_current_timestamp()-ts_begin)/1000 << std::endl;
	#endif
	return;
}
*/

}// end namespace

#endif
