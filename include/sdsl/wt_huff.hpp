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
/*! \file wt_huff.hpp
    \brief wt_huff.hpp contains a class for the wavelet tree of byte sequences which is in Huffman shape.
	\author Simon Gog and Timo Beller
*/
#ifndef INCLUDED_SDSL_WT_HUFF
#define INCLUDED_SDSL_WT_HUFF

#include "int_vector.hpp"
#include "rank_support_v.hpp"
#include "rank_support_v5.hpp"
#include "select_support_mcl.hpp"
#include "rrr_vector.hpp"
#include "bitmagic.hpp"
#include "util.hpp"
#include "wt_helper.hpp"
#include "wt.hpp"
#include <algorithm> // for std::swap
#include <stdexcept>
#include <vector>
#include <utility> // for pair
#include <deque>
#include <queue>


#ifdef SDSL_DEBUG
	#define SDSL_DEBUG_WT_HUFF
#endif


//#define SDSL_DEBUG_WAVELET_TREE

#ifdef SDSL_DEBUG_WAVELET_TREE
#include "testutils.hpp"
#endif

//! Namespace for the succinct data structure library.
namespace sdsl{

const int_vector<>::size_type ZoO[2] = {0, (int_vector<>::size_type)-1};

template<class size_type>	
struct _node{
	size_type 	tree_pos; 		// pointer into the bit_vector, which represents the wavelet tree 
	size_type 	tree_pos_rank;	// precalculated rank for the prefix up to but not including tree_pos
	uint16_t	parent;			// pointer to the parent
	uint16_t	child[2];		// pointer to the children

	_node(	size_type tree_pos=0, size_type tree_pos_rank=0, uint16_t parent=_undef_node, 
			uint16_t child_left=_undef_node, uint16_t child_right=_undef_node):
		  	tree_pos(tree_pos), tree_pos_rank(tree_pos_rank), parent(parent)
	{
		child[0] = child_left;
		child[1] = child_right;	
	}

	_node& operator=(const _node& v){
		if( this != &v ){
			tree_pos 		= v.tree_pos;
			tree_pos_rank 	= v.tree_pos_rank;
			parent			= v.parent;
			child[0] 		= v.child[0];
			child[1] 		= v.child[1];
		}
		return *this;
	}

	size_type serialize(std::ostream &out)const{
		size_type written_bytes = 0;
		out.write((char*)&tree_pos, sizeof(tree_pos));
		written_bytes += sizeof(tree_pos);
		out.write((char*)&tree_pos_rank, sizeof(tree_pos_rank));
		written_bytes += sizeof(tree_pos_rank);
		out.write((char*)&parent, sizeof(parent));
		written_bytes += sizeof(parent);
		out.write((char*)child, 2*sizeof(child[0]));
		written_bytes += 2*sizeof(child[0]);
		return written_bytes;
	}

	void load(std::istream &in){
		in.read((char*) &tree_pos, sizeof(tree_pos));
		in.read((char*) &tree_pos_rank, sizeof(tree_pos_rank));
		in.read((char*) &parent, sizeof(parent));
		in.read((char*) child, 2*sizeof(child[0]));
	}
};

//! A Wavelet Tree class for byte sequences.
/*!
 * A wavelet tree is build for a vector of characters over the alphabet \f$\Sigma\f$. 
 * This class should be used only for small alphabets \f$\Sigma \ll n\f$ (see int_wavelet_tree for a wavelet tree for big alphabets). 
 * The wavelet tree \f$wt\f$ consists of a tree of bitvectors and provides three efficient methods:
 *   - The "[]"-operator: \f$wt[i]\f$ returns the ith symbol of vector for which the wavelet tree was build for.
 *   - The rank method: \f$wt.rank(i,c)\f$ returns the number of occurences of symbol \f$c\f$ in the prefix [0..i-1] in the vector for which the wavelet tree was build for.
 *   - The select method: \f$wt.select(j,c)\f$ returns the index \f$i\in [0..size()-1]\f$ of the jth occurence of symbol \f$c\f$. 
 *
 *	\par Space complexity
 *		\f$\Order{n H_0 + 2|\Sigma|\log n}\f$ bits, where \f$n\f$ is the size of the vector the wavelet tree was build for.
 *
 *   @ingroup wt
 */
template<class BitVector = bit_vector, 
	     class RankSupport = rank_support_v5<>, 
		 class SelectSupport=select_support_mcl<>, 
		 class SelectSupportZero=select_support_mcl<0>,
		 bool dfs_shape=0 >
class wt_huff{
	public:
	typedef int_vector<>::size_type	size_type;
	typedef unsigned char		 	value_type;
	typedef BitVector				bit_vector_type;
	typedef RankSupport				rank_1_type;
	typedef SelectSupport           select_1_type;
	typedef SelectSupportZero		select_0_type;

	private:
#ifdef WT_HUFF_CACHE
		mutable value_type m_last_access_answer;
		mutable size_type  m_last_access_i;
		mutable size_type  m_last_access_rl;
#endif		

		size_type 			m_size;
		size_type 			m_sigma; 		//<- \f$ |\Sigma| \f$
		bit_vector_type		m_tree;			// bit vector to store the wavelet tree
		RankSupport			m_tree_rank;	// rank support for the wavelet tree bit vector
		SelectSupport		m_tree_select1;	// select support for the wavelet tree bit vector
		SelectSupportZero	m_tree_select0; 

		_node<size_type> 	m_nodes[511];	 // nodes for the Huffman tree structure
		uint16_t			m_c_to_leaf[256];// map symbol c to a leaf in the tree structure
											 // if m_c_to_leaf[c] == _undef_node the char does
											 // not exists in the text
		uint64_t			m_path[256];	 // path information for each char; the bits at position 
											 // 0..55 hold path information; bits 56..63 the length
											 // of the path in binary representation
//		wavelet_tree<>		m_check;

		typedef std::pair<size_type, size_type> tPII;  // pair (frequency, node_number) for constructing the Huffman tree
		typedef std::priority_queue<tPII, std::vector<tPII>, std::greater<tPII> >  tMPQPII; // minimum priority queue 

		void copy(const wt_huff &wt){
			m_size 			= wt.m_size;
			m_sigma 		= wt.m_sigma;
			m_tree			= wt.m_tree;
			m_tree_rank 	= wt.m_tree_rank;
			m_tree_rank.set_vector(&m_tree);
			m_tree_select1	= wt.m_tree_select1;
			m_tree_select1.set_vector(&m_tree);
			m_tree_select0	= wt.m_tree_select0;
			m_tree_select0.set_vector(&m_tree);
			for(size_type i=0; i < 511; ++i)
				m_nodes[i] = wt.m_nodes[i];
			for(size_type i=0; i<256; ++i)
				m_c_to_leaf[i] = wt.m_c_to_leaf[i];
			for(size_type i=0; i<256; ++i){
				m_path[i] = wt.m_path[i];
			}
//			m_check = wt.m_check;
		}	

		// insert a character into the wavelet tree, see constuct method
		void insert_char(uint8_t old_chr, size_type *tree_pos, size_type times, bit_vector &f_tree){
			uint32_t path_len = (m_path[old_chr]>>56);
			uint64_t p = m_path[old_chr];
			for(uint32_t node=0, l=0; l<path_len; ++l, p >>= 1){
				if( p&1 ){
					f_tree.set_int(tree_pos[node], 0xFFFFFFFFFFFFFFFFULL,times);
				}
				tree_pos[node] += times;
				node = m_nodes[node].child[p&1];
			}
		}

		// calculates the Huffman tree and returns the size of the WT bit vector
		size_type construct_huffman_tree(size_type *C){
			tMPQPII pq; // priority queue
			std::vector<_node<size_type> > temp_nodes(2*m_sigma-1);  // vector for nodes of the Huffman tree
			size_type node_cnt=0;								// counter for the nodes
			for(size_type i=0; i < 256; ++i) // add leafs of Huffman tree
				if( C[i] > 0 ){
					pq.push( tPII(C[i], node_cnt) ); // push (frequency, pointer to node)
					temp_nodes[node_cnt++] = _node<size_type>( C[i], i ); // initial tree_pos with number of occurences 
																		  // and tree_pos_rank value with the code of the corresponding char
																		  // parent, child[0], and child[1] are set to _undef_node
				}
			while( pq.size() > 1 ){
				tPII v1, v2;
				v1 = pq.top(); pq.pop();
				v2 = pq.top(); pq.pop();
				temp_nodes[ v1.second ].parent = node_cnt; // parent is new node 
				temp_nodes[ v2.second ].parent = node_cnt; // parent is new node
				size_type frq_sum = v1.first + v2.first;
				pq.push( tPII(frq_sum, node_cnt) ); // push new node to the priority queue
				temp_nodes[ node_cnt++ ] = _node<size_type>(frq_sum, 0, _undef_node, v1.second, v2.second);
			}
			// Convert Huffman tree into breadth first search order in memory and
			// calculate tree_pos values
			m_nodes[0] = temp_nodes[node_cnt-1];  // insert root at index 0
			size_type tree_size = 0;
			node_cnt = 1;
			uint16_t last_parent = _undef_node;
				std::deque<size_type> q;
				q.push_back(0);
			while( !q.empty() ){
				size_type idx;
				if( !dfs_shape ){
					idx = q.front(); q.pop_front();
				}else{
					idx = q.back(); q.pop_back();
				}
				size_type frq = m_nodes[idx].tree_pos; // frq_sum was stored in tree_pos
				m_nodes[idx].tree_pos = tree_size;
				if( m_nodes[idx].child[0] != _undef_node ) // if node is not a leaf
					tree_size += frq;					   // add frequency, as leaves have size 0
				if( idx > 0 ){ // node is not the root
					if( last_parent != m_nodes[idx].parent )
						m_nodes[m_nodes[idx].parent].child[0] = idx;
					else
						m_nodes[m_nodes[idx].parent].child[1] = idx;
					last_parent = m_nodes[idx].parent;
				}
				if( m_nodes[idx].child[0] != _undef_node ){ // if node is not a leaf
					for(size_type k=0; k<2; ++k){			// add children to tree
						m_nodes[node_cnt] = temp_nodes[ m_nodes[idx].child[k] ];
						m_nodes[node_cnt].parent = idx;
						q.push_back(node_cnt);
						m_nodes[idx].child[k] = node_cnt++;
					}
				}
			}

			// initialize m_c_to_leaf
			for(size_type i=0; i<256; ++i)
				m_c_to_leaf[i] = _undef_node; // if c is not in the alphabet m_c_to_leaf[c] = _undef_node
			for(size_type i=0; i < 2*sigma-1; ++i){
				if( m_nodes[i].child[0] == _undef_node ) 				// if node is a leaf
					m_c_to_leaf[(uint8_t)m_nodes[i].tree_pos_rank] = i; // calculate value
			}
			// initialize path information
			// Note: In the case of a bfs search order, 
			// we can classify nodes as rigth child and left child with an easy criterion:
			//       node is a left child, if node%2==1
			//       node is a rigth child, if node%2==0
			for(size_type c=0; c<256; ++c){
				if( m_c_to_leaf[c] != _undef_node ){ // if char exists in the alphabet
					size_type node = m_c_to_leaf[c];
					uint64_t w = 0; // path
					uint64_t l = 0; // path len
					while( node != 0 ){ // while node is not the root
						w <<= 1;
//						if( (node & 1) == 0 )// if the node is a right child
//							w |= 1ULL;
						if( m_nodes[m_nodes[node].parent].child[1] == node )
							w |= 1ULL;
						++l;
						node = m_nodes[node].parent; // go up the tree
					}
					if( l > 56 ){
						std::cerr<<"Huffman tree has max depth > 56!!! ERROR"<<std::endl;
						throw std::logic_error("Huffman tree size is greater than 56!!!");
					}
					m_path[c] = w | (l << 56);
				}else{
					m_path[c] = 0; // i.e. len is also 0, good for special case in rank()
				}
			}
			return tree_size;
		}

		void construct_init_rank_select(){
			m_tree_rank.init(&m_tree);
			m_tree_select0.init(&m_tree);
			m_tree_select1.init(&m_tree);
		}

		void construct_precalc_node_ranks(){
			for(size_type i=0; i<2*m_sigma-1; ++i){
				if( m_nodes[i].child[0] != _undef_node ) // if node is not a leaf
					m_nodes[i].tree_pos_rank = m_tree_rank( m_nodes[i].tree_pos );
			}
		}


		// recursive internal version of the method interval_symbols
		void _interval_symbols(size_type i, size_type j, size_type &k, 
								std::vector<unsigned char> &cs, 
								std::vector<size_type> &rank_c_i, 
								std::vector<size_type> &rank_c_j, uint16_t node)
		{
			// invariant: j>i
			// goto right child
			size_type i_new = (m_tree_rank(m_nodes[node].tree_pos + i) - m_nodes[node].tree_pos_rank);
			size_type j_new = (m_tree_rank(m_nodes[node].tree_pos + j) - m_nodes[node].tree_pos_rank);
			if(i_new!=j_new){
				uint16_t node_new = m_nodes[node].child[1];
				// if node is not a leaf
				if(m_nodes[node_new].child[0] != _undef_node){
					_interval_symbols(i_new, j_new, k, cs, rank_c_i, rank_c_j, node_new);
				}
				else {
					rank_c_i[k] = i_new;
					rank_c_j[k] = j_new;
					cs[k++] = m_nodes[node_new].tree_pos_rank;
				}
			}
			// goto left child
			i -= i_new; j -= j_new; 
			if(i != j){
				uint16_t node_new = m_nodes[node].child[0];
				// if node is not a leaf
				if(m_nodes[node_new].child[0] != _undef_node) {
					_interval_symbols(i, j, k, cs, rank_c_i, rank_c_j, node_new);
				}
				else {
					rank_c_i[k] = i;
					rank_c_j[k] = j;
					cs[k++] = m_nodes[node_new].tree_pos_rank;
				}
			}	
		}


	public:	
	
	const size_type &sigma;	
	const bit_vector_type &tree;

	// Default constructor
	wt_huff():m_size(0),m_sigma(0), sigma(m_sigma),tree(m_tree){};



	//! Constructor
	/*!
	 *	\param rac Reference to the vector (or unsigned char array) for which the wavelet tree should be build.  
	 *	\param size Size of the prefix of the vector (or unsigned char array) for which the wavelet tree should be build.
	 *	\par Time complexity
	 *		\f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
	 */
	template<typename RandomAccessContainer>
	wt_huff(const RandomAccessContainer &rac, size_type size):m_size(size), m_sigma(0), sigma(m_sigma), tree(m_tree){
		construct( rac, size );
	}	

	template<uint8_t w>
	wt_huff(const int_vector<w> &rac):m_size(rac.size()), m_sigma(0), sigma(m_sigma), tree(m_tree){
		construct( rac, rac.size() );
	}
	

	template<typename RandomAccessContainer>
	void construct(const RandomAccessContainer &rac, size_type size){	
		m_size = size;
		if( m_size == 0 )
			return;
		// O(n + |\Sigma|\log|\Sigma|) algorithm for calculating node sizes
		size_type C[256] = {0};
		//  1. Count occurences of characters
		for(size_type i=0; i < size; ++i){
			++C[rac[i]];
		}
		// 2. Calculate effective alphabet size
		calculate_effective_alphabet_size(C, m_sigma);
		// 3. Generate Huffman tree 
		size_type tree_size = construct_huffman_tree(C);
		// 4. Generate wavelet tree bit sequence m_tree

		bit_vector tmp_tree(tree_size, 0);  // initialize bit_vector for the tree
	    //  Calculate starting position of wavelet tree nodes	
		size_type tree_pos[511];
		for(size_type i=0; i < 2*sigma-1; ++i){
			tree_pos[i] = m_nodes[i].tree_pos;
		}
		uint8_t old_chr = rac[0], times = 0; 
		for(size_type i=0; i < m_size; ++i){
			uint8_t chr = rac[i]; 
			if( chr	!= old_chr ){ 
				insert_char( old_chr, tree_pos, times, tmp_tree );
				times = 1;
				old_chr = chr;
			}
			else{ // chr == old_chr
				++times;
				if( times == 64 ){
					insert_char( old_chr, tree_pos, times, tmp_tree );
					times = 0;
				}
			}
		}
		if(times > 0 ){
			insert_char( old_chr, tree_pos, times, tmp_tree );
		}
		util::assign( m_tree, tmp_tree );
		// 5. Initialize rank and select data structures for m_tree
		construct_init_rank_select();
		// 6. Finish inner nodes by precalculating the tree_pos_rank values
		construct_precalc_node_ranks();
	}

	template<class size_type_class>
	wt_huff(int_vector_file_buffer<8, size_type_class> &rac, size_type size):m_size(size), m_sigma(0), sigma(m_sigma), tree(m_tree){
		construct(rac, size);
	}

	//! Construct the wavelet tree from a random access container
	/*! \param rac A random access container
	 *	\param size The length of the prefix of the random access container, for which the wavelet tree should be build 
	 */
	template<class size_type_class>
	void construct(int_vector_file_buffer<8, size_type_class> &rac, size_type size){
//		m_check.construct(rac, size);	
		m_size = size;
		if( m_size == 0 )
			return;
		// O(n + |\Sigma|\log|\Sigma|) algorithm for calculating node sizes
		size_type C[256] = {0};
//		rac.reset();
		//  1. Count occurences of characters
//		for(size_type i=0, r_sum=0, r = rac.load_next_block(); r_sum < m_size;){
//			for(; i < r_sum+r; ++i){
//				++C[rac[i-r_sum]];
//			}
//			r_sum += r; r = rac.load_next_block();
//		}
		// 1. Count occurences of characters
		calculate_character_occurences(rac, m_size, C);
		// 2. Calculate effective alphabet size
		calculate_effective_alphabet_size(C, m_sigma);
		// 3. Generate Huffman tree 
		size_type tree_size = construct_huffman_tree(C);
		// 4. Generate wavelet tree bit sequence m_tree
	
		bit_vector tmp_tree(tree_size, 0);  // initialize bit_vector for the tree
	    //  Calculate starting position of wavelet tree nodes	
		size_type tree_pos[511];
		for(size_type i=0; i < 2*sigma-1; ++i){
			tree_pos[i] = m_nodes[i].tree_pos;
		}
		rac.reset();
		for(size_type i=0, r_sum=0, r = rac.load_next_block(); r_sum < m_size; ){
			uint8_t old_chr = rac[i-r_sum], times = 0;
			for(; i < r_sum+r; ++i){
				uint8_t chr = rac[i-r_sum]; 
				if( chr	!= old_chr ){ 
				    insert_char( old_chr, tree_pos, times, tmp_tree );
					times = 1;
					old_chr = chr;
				}
				else{ // chr == old_chr
					++times;
					if( times == 64 ){
				    	insert_char( old_chr, tree_pos, times, tmp_tree );
						times = 0;
					}
				}
			}
			if(times > 0 ){
				insert_char( old_chr, tree_pos, times, tmp_tree );
			}
			r_sum += r; r = rac.load_next_block();
		}
		util::assign( m_tree, tmp_tree );
		// 5. Initialize rank and select data structures for m_tree
		construct_init_rank_select();
		// 6. Finish inner nodes by precalculating the tree_pos_rank values
		construct_precalc_node_ranks();
	}


	//! Copy constructor
	wt_huff(const wt_huff &wt):sigma(m_sigma), tree(m_tree){
		copy(wt);
	}	

	//! Assignment operator
	wt_huff& operator=(const wt_huff &wt){
		if( this != &wt ){
			copy(wt);
		}
		return *this;
	}

	//! Swap operator
	void swap(wt_huff &wt){
		if(this != &wt){
			std::swap(m_size, wt.m_size);
			std::swap(m_sigma,  wt.m_sigma);
			m_tree.swap(wt.m_tree);
			m_tree_rank.swap(wt.m_tree_rank); // rank swap after the swap of the bit vector m_tree
			m_tree_select1.swap(wt.m_tree_select1); // select1 swap after the swap of the bit vector m_tree
			m_tree_select0.swap(wt.m_tree_select0); // select0 swap after the swap of the bit vector m_tree

			for(size_type i=0; i < 511; ++i)
				std::swap(m_nodes[i], wt.m_nodes[i]);
			for(size_type i=0; i<256; ++i)
				std::swap(m_c_to_leaf[i], wt.m_c_to_leaf[i]);
			for(size_type i=0; i<256; ++i)
				std::swap(m_path[i], wt.m_path[i]);

//			m_check.swap( wt.m_check );
		}
	}

	//! Returns the size of the original vector.
	size_type size()const{
		return m_size;
	}

	//! Returns whether the wavelet tree contains no data.
	bool empty()const{
		return m_size == 0;
	}

	//! Recovers the ith symbol of the original vector.
	/*! \param i The index of the symbol in the original vector. \f$i \in [0..size()-1]\f$
	 *	\return The ith symbol of the original vector. 
	 *  \par Time complexity
	 *		\f$ \Order{H_0} \f$ on average, where \f$ H_0 \f$ is the zero order entropy
	 *      of the sequence.
	 */
	value_type operator[](size_type i)const{  // TODO: Maybe it is good to integrate a cache here
		                                      // which stores how many of the next symbols are equal
		                                      // with the current char
		size_type node = 0; // start at root node
		while( m_nodes[node].child[0] != _undef_node ){ // while node is not a leaf
			if( m_tree[ m_nodes[node].tree_pos + i]  ){ // goto the right child
				i = m_tree_rank( m_nodes[node].tree_pos + i ) - m_nodes[node].tree_pos_rank;
				node = m_nodes[node].child[1];
			}else{ // goto the left child
				i -=  (m_tree_rank( m_nodes[node].tree_pos + i ) - m_nodes[node].tree_pos_rank);
				node = m_nodes[node].child[0];
			}
		}	
		return m_nodes[node].tree_pos_rank;
	};

	//! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
	/*!
	 *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
   	 *  \param c The symbol to count the occurences in the prefix.
	 *	\return The number of occurences of symbol c in the prefix [0..i-1] of the supported vector. 
	 *  \par Time complexity
	 *		\f$ \Order{H_0} \f$
	 */
	size_type rank(size_type i, value_type c)const{
		assert(i>=0 and i <= size());
		uint64_t p = m_path[c]; 
		uint32_t path_len = (m_path[c]>>56); // equals zero if char was not present in the original text
		size_type result = i & ZoO[path_len>0]; // important: result has type size_type and ZoO has type size_type
		uint32_t node=0;
		for(uint32_t l=0; l<path_len and result; ++l, p >>= 1){	
//			result = (ZoO[1-(p&1)]&result) - (ZoO[1-(p&1)] ^ (m_tree_rank(m_nodes[node].tree_pos + result) - m_nodes[node].tree_pos_rank) );
			if( p&1 ){
				result 	= (m_tree_rank(m_nodes[node].tree_pos+result) -  m_nodes[node].tree_pos_rank); 
			}
			else{
				result -=  (m_tree_rank(m_nodes[node].tree_pos+result) -  m_nodes[node].tree_pos_rank); 
			}
			node = m_nodes[node].child[p&1]; // goto child
		}
//		size_type r2 = m_check.rank(i, c);
/*		if( r2 != result ){
			std::cerr<<"ERROR rank r="<<result<<" != "<<r2<<"=r2 for input ("<<i<<","<<c<<")"<<
				     " len="<<path_len<<" m_path[c]="<<m_path[c]<<" c="<<(char)m_nodes[node].tree_pos_rank<<std::endl;
			return r2;
		}
*/		return result;
	};

	//! Calculates how many occurences of symbol wt[i] are in the prefix [0..i-1] of the supported sequence.
	/*!
	 *	\param i The index of the symbol.
	 *  \param c Reference that will contain the symbol at position i after the execution of the method.
	 *  \return The number of occurences of symbol wt[i] in the prefix [0..i-1]
	 *	\par Time complexity
	 *		\f$ \Order{H_0} \f$ 
	 */
	size_type rank_ith_symbol(size_type i, value_type &c)const{
		assert(i>=0 and i < size());
		uint32_t node=0;
		while( m_nodes[node].child[0] != _undef_node ){ // while node is not a leaf
			if( m_tree[m_nodes[node].tree_pos + i] ){ // if bit is set at position goto rigth child
				i 	= (m_tree_rank(m_nodes[node].tree_pos + i) -  m_nodes[node].tree_pos_rank); 
				node = m_nodes[node].child[1]; 
			}else{ // goto left child
				i -= (m_tree_rank(m_nodes[node].tree_pos + i) -  m_nodes[node].tree_pos_rank); 
				node = m_nodes[node].child[0]; 
			}
		}
		c = m_nodes[node].tree_pos_rank;
		return i;
	}

	//! Calculates the ith occurence of the symbol c in the supported vector.
	/*!
	 *  \param i The ith occurence. \f$i\in [1..rank(size(),c)]\f$.
	 *  \param c The symbol c.
	 *  \par Time complexity
	 *		\f$ \Order{H_0} \f$
	 */
	size_type select(size_type i, value_type c)const{
		uint16_t node = m_c_to_leaf[c];
		if( node == _undef_node ){ // if c was not present in the original text 
			return m_size;		   // -> return a position right to the end 
		}
/*		if( i == 0 ){
			std::cerr<<"WARNING: i=0"<<std::endl;
			return m_size;
		}
		if( i > rank(size(), c) ){
			std::cerr<<"WARNING: i="<<i<<" > rank(size(),"<<c<<")="<<rank(size(),c)<<std::endl;
			return m_size;		
		}
*/		
		size_type result = i-1;		// otherwise 
		uint64_t p = m_path[c];
		uint32_t path_len = (p>>56);
		p <<=  (64-path_len);
		for(uint32_t l=0; l<path_len; ++l, p <<= 1){
//			if( node & 1 ){ // node was a left child, in the case of bfs order
			if( (p & 0x8000000000000000ULL)==0 ){ // node was a left child
				node = m_nodes[node].parent;
				result = m_tree_select0( m_nodes[node].tree_pos-m_nodes[node].tree_pos_rank + result + 1 )
						 - m_nodes[node].tree_pos;
			}else{ // node was a right child
				node = m_nodes[node].parent;
				result = m_tree_select1( m_nodes[node].tree_pos_rank + result + 1 )
				         - m_nodes[node].tree_pos;
			}
		}
/*		size_type r2 = m_check.select(i, c);
		if( r2 != result ){
			std::cerr<<"ERROR select r="<<result<<" != "<<r2<<"=r2 for input ("<<i<<","<<c<<")"<<std::endl;
			return r2;
		}
*/		
		return result;
	};


	//! Calculates for each symbol c in wt[i..j-1], how many times c occures in wt[0..i-1] and wt[0..j-1].
	/*!
	 *	\param i The start index (inclusive) of the interval.
	 *	\param j The end index (exclusive) of the interval.
	 *	\param k Reference that will contain the number of different symbols in wt[i..j-1].
	 *  \param cs Reference to a vector of size k that will contain all symbols that occur in wt[i..j-1] in arbitrary order.
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
	 void interval_symbols(size_type i, size_type j, size_type &k, 
		                    std::vector<unsigned char> &cs, 
						    std::vector<size_type> &rank_c_i, 
						    std::vector<size_type> &rank_c_j)
	{
		if( i==j ){
			k = 0;
			return;
		}
		else if((j-i)==1){
			k = 1;
			rank_c_i[0] = rank_ith_symbol(i, cs[0]);
			rank_c_j[0] = rank_c_i[0]+1;
			return;
		}
		else if((j-i)==2) {
			rank_c_i[0] = rank_ith_symbol(i, cs[0]);
			rank_c_i[1] = rank_ith_symbol(i+1, cs[1]);
			if(cs[0]==cs[1]){
				k = 1;
				rank_c_j[0] = rank_c_i[0]+2;
				return;
			}
			else{
				k = 2;
				rank_c_j[0] = rank_c_i[0]+1;
				rank_c_j[1] = rank_c_i[1]+1;
				return;
			}
		}
		else{
			k = 0;
			_interval_symbols(i, j, k, cs, rank_c_i, rank_c_j, 0);	
		}
	} 



	//! Serializes the data structure into the given ostream
	size_type serialize(std::ostream &out)const{
		size_type written_bytes = 0;
		out.write((char*)&m_size, sizeof(m_size));
		written_bytes += sizeof(m_size);
		out.write((char*)&m_sigma, sizeof(m_sigma));
		written_bytes += sizeof(m_sigma);
		written_bytes += m_tree.serialize(out);
		written_bytes += m_tree_rank.serialize(out);
		written_bytes += m_tree_select1.serialize(out);
		written_bytes += m_tree_select0.serialize(out);
		for(size_type i=0; i < 511; ++i){
			written_bytes += m_nodes[i].serialize(out);
		}
		out.write((char*) m_c_to_leaf, 256*sizeof(m_c_to_leaf[0]) );
		written_bytes += 256*sizeof(m_c_to_leaf[0]); // add written bytes from previous loop
		out.write((char*) m_path, 256*sizeof(m_path[0]) );
		written_bytes += 256*sizeof(m_path[0]); // add written bytes from previous loop
//		m_check.serialize(out);
		return written_bytes;
	}

	//! Loads the data structure from the given istream.
	void load(std::istream &in){
		in.read((char*) &m_size, sizeof(m_size));
		in.read((char*) &m_sigma, sizeof(m_sigma));
		m_tree.load(in);
		m_tree_rank.load(in, &m_tree);
		m_tree_select1.load(in, &m_tree);
		m_tree_select0.load(in, &m_tree);
		for(size_type i=0; i < 511; ++i){
			m_nodes[i].load(in);
		}
		in.read((char*) m_c_to_leaf, 256*sizeof(m_c_to_leaf[0]));
		in.read((char*) m_path, 256*sizeof(m_path[0]));
//		m_check.load(in);
	}

#ifdef MEM_INFO
		//! Print some infos about the size of the compressed suffix tree 
		void mem_info(std::string label="")const{
			if(label=="")
				label="wt_huff";
			size_type bytes = util::get_size_in_bytes(*this);
			std::cout << "list(label = \""<<label<<"\", size = "<< bytes/(1024.0*1024.0) <<"\n,";
			m_tree.mem_info("data"); std::cout<<",";
			m_tree_rank.mem_info("rank"); std::cout<<",";
			m_tree_select1.mem_info("select 1"); std::cout<<",";
			m_tree_select0.mem_info("select 0"); std::cout << ")\n";
			// TODO: add m_nodes, m_c_to_leaf, m_path?
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

typedef wt_huff<rrr_vector<>,
		        rrr_vector<>::rank_1_type,
				rrr_vector<>::select_1_type,
				rrr_vector<>::select_0_type, 0> wt_huff_rrr;

}// end namespace sdsl

#endif // end file 
