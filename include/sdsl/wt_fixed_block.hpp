

/*
 *   template parameter block size 
 *
 * average access time \f$ \Order{H_k} \f$
 * References:
 *    Juha K"arkk"ainen and Simon J. Puglisi "Fixed Block Compression Boosting in FM-Indexes", SPIRE 2011
 *
 * Implementation for small alphabet size (<=256) 


 *
 */

#ifndef INCLUDED_SDSL_WT_FIXED_BLOCK
#define INCLUDED_SDSL_WT_FIXED_BLOCK

#include "int_vector.hpp"
#include "wt_helper.hpp"
#include "bp_support.hpp"

namespace sdsl{

template<uint32_t block_size = 1<<15, class bp_support_type = bp_support_sada<> >
class wt_fixed_block{
	public:
		typedef bit_vector::size_type size_type;
	private:
		size_type			m_size; // stores the size of the sequence
		size_type			m_blk_cnt; // stores the number of blocks
		uint16_t			m_sigma; // size of the effective alphabet
		uint16_t			m_c_to_nr; // map symbol c to a unqiue number in the range [0..m_sigma-1] 
		bit_vector			m_occ; // indicator bit vector for occurences of a character c
		                           // in a block i; if m_occ[m_c_to_nr[c]*(m_blk_cnt+1) + i +1] = 1
		bit_vector 			m_bp;  // the concatenated balanced parentheses sequences for the 
		                           // Huffman trees of the blocks
		bp_support_type 	m_bp_support; // supporting data structure for m_bp
		int_vector<>        m_bp_pointers; // for each block, we store a pointer to the start
		                                   // of the corrsponding start of the tree representation in m_bp
		bit_vector			m_tree; // store the concatenated bit vectors of the wavelet trees of
		                            // the blocks
		int_vector<>        m_tree_pointers; // for each block, we store a pointer to the start
		                                     // of the  
		
		// TODO: * for each node in a huffman shaped wt, we have to save the size of the right nodes
		//       * for each occurring character we have to store the number of its leaf
		//           or the position ... number of leaf is more space efficient but requires an additional select
		int_vector<>		m_leaf_nr;

		// if alphabet size is <= 4,text almost random, and block size not very larg, 
		// use scanning with broadword methods b1Cnt, b11Cnt,... to answer questions

		// if alphabet size equals 1, do not store m_tree, all answers can be 
	    // calculated easily

		// if the alphabet size is 2, we can decide with m_leaf_nr, which char is encoded as 1 or 0.
	public:

	//! Constructor
	template<class size_type_class>	
	wt_fixed_block(int_vector_file_buffer<8, size_type_class> &text, size_type size){
		m_size = size;
		if( m_size == 0 )
			return;
		// calculate number of blocks
		size_type m_blk_cnt = (m_size + block_size - 1) / block_size;
		// caculate the number of occurences of each character in the text 
		size_type C[256] = {0};
		calculate_character_occurences(text, m_size, C);
		// calculate the size of the effective alphabet
		calculate_effective_alphabet_size(C, m_sigma);	

		m_occ.resize( m_sigma * (m_blk_cnt+1) ); // +1 since we also store a dummy element for block[-1]

		// initialize m_c_to_nr
		uint16_t nr=0;
		for(size_type i=0; i<256; ++i){
			if( C[i]>0 ){
				m_c_to_nr[i] = nr++;
			}else{
				m_c_to_nr[i] = _undef_node;
			}
		}
		

	}
}

} // end namespace sdsl

#endif
