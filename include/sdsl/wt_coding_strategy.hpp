/* sdsl - succinct data structures library
    Copyright (C) 2012 Simon Gog

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
/*! \file wt_coding_strategy.hpp
    \brief wt_coding_strategy.hpp contains policy classes for the realization of prefix-free codes in wavelet trees.
	\author Simon Gog
*/

#ifndef INCLUDED_WT_CODING_STRATEGY
#define INCLUDED_WT_CODING_STRATEGY

#include "int_vector.hpp"
#include "int_vector_buffer.hpp"

namespace sdsl
{

//! Class for tree based Huffman solution on byte alphabets.
class fat_huff_byte
{
    public:
        typedef bit_vector::size_type size_type;
    public:
        // Default constructor
        fat_huff_byte();
        // Construct from byte stream
        template<class size_type_class>
        fat_huff_byte(int_vector_buffer<8>& text_buf);
        // The length of the code for symbol `c`
        uint8_t  code_len(uint8_t c);
        // The code for symbol `c`
        uint64_t code(uint8_t c);
        // The symbol for a codeword `code` of length `len`
        uint8_t  symbol(uint64_t code, uint8_t len);

        bool valid_code(uint64_t code, uint8_t len);

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr,
                            std::string name="")const;
        void load(std::istream& in);
};


//! Class for canonical Huffman solution on byte alphabets.
/*
class canonical_huff_tree_byte{
	private:
		int_vector<> freq_order; // decreasing order
		int_vector<> code_len;   //

	public:

		inline uint64_t code_len(uint8_t c){
			return code_len[c];
		}

}
*/

}

#endif
