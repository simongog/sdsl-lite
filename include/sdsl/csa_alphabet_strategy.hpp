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
/*! \file csa_alphabet_strategy.hpp
    \brief csa_alphabet_strategy.hpp includes different strategy classes for representing an alphabet of a CSA.
	\author Simon Gog
*/

#ifndef INCLUDED_CSA_ALPHABET_STRATEGY
#define INCLUDED_CSA_ALPHABET_STRATEGY

/* A alphabet strategy provides the following features:
 *   * Member `sigma` which contains the size (=number of unique symbols) of the alphabet.
 *   * Method `is_present(char_type c)` which indicates if character c occurs in the text.
 *   * Container `char2comp` which maps a symbol to a number [0..sigma-1]. The alphabetic
 *     order is preserved.
 *   * Container `comp2char` which is the inverse mapping of char2comp.
 *   * Container `C` contains the cumulative counts of occurrences. C[i] is the cumulative
 *     count of occurrences of symbols `comp2char[0]` to `comp2char[i-1]` in the text.
 *   * Typedefs for the four above members:
 *       * char2comp_type
 *       * comp2char_type
 *       * C_type
 *       * sigma_type
 *   * Constructor. Takes a int_vector_file_buffer<8, size_type_class> for byte-alphabets
 *     and int_vector_file_buffer<0, size_type_class> for integer-alphabets. 
 */

#include "int_vector.hpp"

namespace sdsl{

// forward declarations

class byte_alphabet_stategy;

template<class bit_vector_type     = bit_vector,
	     class rank_support_type   = typename bit_vector_type::rank_1_type,
	     class select_support_type = typename bit_vector_type::select_1_type,
		 class C_array_type		   = int_vector<>
	    >
class succinct_byte_alphabet_strategy;

//! A simple space greedy representation for byte alphabets.
/*!
 *  \par Space consumption:
 *       At least: 2.5 kB
 *       Details:  char2comp + comp2char take  2*256 + 2*8 bytes 
 *                 m_C                   takes       257*8 bytes
 *                 m_sigma               takes           2 bytes
 */
class byte_alphabet_stategy{
	public:
		typedef int_vector<>::size_type size_type;
		typedef int_vector<8>			char2comp_type;
		typedef int_vector<8>			comp2char_type;
		typedef int_vector<64>			C_type;
		typedef uint16_t				sigma_type;
	private:
	    char2comp_type	m_char2comp;
        comp2char_type	m_comp2char; 
        C_type  		m_C; 		 // cumulative counts for the compact alphabet [0..sigma-1]
        sigma_type		m_sigma;     
	
	public:

		const char2comp_type& char2comp;
		const comp2char_type& comp2char;
		const C_type&		  C;
		const sigma_type&     sigma;

		//! Default constructor
		byte_alphabet_stategy();

		//! Construct from a byte-stream
		/*!
		 *  \param text_buf	Byte stream.
		 *  \param len		Length of the byte stream. 
		 */
		template<class size_type_class>
		byte_alphabet_stategy(int_vector_file_buffer<8, size_type_class> &text_buf, size_type_class len): 
							   char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma)
		{
			m_sigma = 0;
    		text_buf.reset();
   			if (0 == text_buf.int_vector_size)
        		return;
			// initialize vectors 
			util::assign(m_C	    , int_vector<64>(257, 0));
			util::assign(m_char2comp, int_vector<8>(256,0));
    		util::assign(m_comp2char, int_vector<8>(256,0));
			// count occurrences of each symbol 
     		for (size_type i=0, r_sum=0, r = text_buf.load_next_block(); i < len;) {
        		for (; i < r_sum+r; ++i) {
            		++m_C[text_buf[i-r_sum]];
        		}
        		r_sum += r; r = text_buf.load_next_block();
    		}
    		assert(1 == m_C[0]); // null-byte should occur exactly once
    		m_sigma = 0;
			for (int i=0; i<256; ++i)
				if (m_C[i]) {
					m_char2comp[i] 	 	= m_sigma;
					m_comp2char[sigma]  = i; 
					m_C[m_sigma]		= m_C[i];
					++m_sigma;
				}
			comp2char.resize(m_sigma);
			m_C.resize(m_sigma+1);
			for (int i=(int)m_sigma; i > 0; --i) m_C[i] = m_C[i-1]; 
			m_C[0] = 0;
			for (int i=1; i <= (int)m_sigma; ++i) m_C[i] += m_C[i-1];
			assert(C[sigma]==len);
		}

		//! Swap operator
        void swap(byte_alphabet_stategy&);

		//! Serialize method
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const;

		//! Load method
        void load(std::istream& in);
};

} // end namespace sdsl

#endif
