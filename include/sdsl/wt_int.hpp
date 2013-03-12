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
/*! \file wt_int.hpp
    \brief wt_int.hpp contains a specialized class for a wavelet tree of a permutation of the numbers from 0..n. This wavelet tree class takes
	       less memory than the generic class for wavelet trees.
	\author Simon Gog, Shanika Kuruppu

*/
#ifndef INCLUDED_SDSL_INT_WAVELET_TREE
#define INCLUDED_SDSL_INT_WAVELET_TREE

#include "sdsl_concepts.hpp"
#include "int_vector.hpp"
#include "rank_support_v.hpp"
#include "select_support_mcl.hpp"
#include "testutils.hpp"
#include "temp_write_read_buffer.hpp"
#include "util.hpp"
#include <set> // for calculating the alphabet size
#include <map> // for mapping a symbol to its lexicographical index
#include <algorithm> // for std::swap
#include <stdexcept>
#include <vector>

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A wavelet tree class for sequences of big alphabet size (like integer alphabet)
/*!
 * A wavelet tree is build for a vector of characters over the alphabet \f$\Sigma\f$.
 * The wavelet tree \f$wt\f$ consists of a tree of bit vector and provides three efficient methods:
 *   - The "[]"-operator: \f$wt[i]\f$ returns the ith symbol of vector for which the wavelet tree was build for.
 *   - The rank method: \f$wt.rank(i,c)\f$ returns the number of occurrences of symbol \f$c\f$ in the prefix [0..i-1] in the vector for which the wavelet tree was build for.
 *   - The select method: \f$wt.select(j,c)\f$ returns the index \f$i\in [0..size()-1]\f$ of the jth occurrence of symbol \f$c\f$.
 *
 *	\par Space complexity
 *		\f$\Order{n\log|\Sigma|}\f$ bits, where \f$n\f$ is the size of the vector the wavelet tree was build for.
 *
 *  \tparam BitVector				Type of the bitvector used for representing the wavelet tree.
 *  \tparam RankSupport				Type of the support structure for rank on ones.
 *  \tparam SelectSupport			Type of the support structure for select on ones.
 *  \tparam SelectSupport			Type of the support structure for select on ones.
 *
 *   @ingroup wt
 */
template<class BitVector   		 = bit_vector,
         class RankSupport 		 = typename BitVector::rank_1_type,
         class SelectSupport	 = typename BitVector::select_1_type,
         class SelectSupportZero = typename BitVector::select_0_type>
class wt_int {
    public:
        typedef int_vector<>::size_type 		size_type;
        typedef int_vector<>::value_type 		value_type;
        typedef BitVector						bit_vector_type;
        typedef RankSupport						rank_1_type;
        typedef SelectSupport					select_1_type;
        typedef SelectSupportZero				select_0_type;
		typedef wt_tag							index_category;
		typedef int_alphabet_tag				alphabet_category;
    protected:
        size_type 				m_size;
        size_type 				m_sigma; 		//<- \f$ |\Sigma| \f$
        bit_vector_type 		m_tree;			// bit vector to store the wavelet tree
        rank_1_type				m_tree_rank;	// rank support for the wavelet tree bit vector
        select_1_type			m_tree_select1;	// select support for the wavelet tree bit vector
        select_0_type			m_tree_select0;
        uint32_t				m_max_depth;
		mutable int_vector<64>	m_path_off;     // array keeps track of path offset in select-like methods
		mutable int_vector<64>  m_path_rank_off;// array keeps track of rank values for the offsets 

        void copy(const wt_int& wt) {
            m_size 			= wt.m_size;
            m_sigma 		= wt.m_sigma;
            m_tree			= wt.m_tree;
            m_tree_rank 	= wt.m_tree_rank;
            m_tree_rank.set_vector(&m_tree);
            m_tree_select1	= wt.m_tree_select1;
            m_tree_select1.set_vector(&m_tree);
            m_tree_select0	= wt.m_tree_select0;
            m_tree_select0.set_vector(&m_tree);
            m_max_depth			= wt.m_max_depth;
			m_path_off			= wt.m_path_off;
			m_path_rank_off		= wt.m_path_rank_off;
        }

	private:
		void init_buffers(uint32_t max_depth){
			util::assign(m_path_off, int_vector<64>(max_depth+1));
			util::assign(m_path_rank_off, int_vector<64>(max_depth+1));
		}	

    public:

        const size_type& sigma;	//!< Effective alphabet size of the wavelet tree.
        const bit_vector_type& tree; //!< A concatenation of all bit vectors of the wavelet tree.

        //! Default constructor
        wt_int():m_size(0),m_sigma(0), m_max_depth(0), sigma(m_sigma), tree(m_tree) { init_buffers(m_max_depth); };

        //! Semi-external constructor
        /*!	\param buf			File buffer of the int_vector for which the wt_int should be build.
		 *  \param size 		Size of the prefix of v, which should be indexed.  
         *	\param max_depth	Maximal depth of the wavelet tree. If set to 0, determined automatically.
         *	\param dir	Directory in which temporary files should be stored during the construction.
         *	\par Time complexity
         *		\f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         *		I.e. we need \Order{n\log n} if rac is a permutation of 0..n-1.
         *	\par Space complexity
         *		\f$ n\log|\Sigma| + O(1)\f$ bits, where \f$n=size\f$.
         */
        template<uint8_t int_width>
        wt_int(int_vector_file_buffer<int_width>& buf, size_type size, uint32_t max_depth=0, std::string dir="./")
            : m_size(size),m_sigma(0), m_max_depth(0), sigma(m_sigma), tree(m_tree) {
			init_buffers(m_max_depth);
			if( 0 == m_size )
				return; 	
            buf.reset();
            size_type n = buf.int_vector_size;  // set n
			if ( n < m_size ){
				throw std::logic_error("n="+util::to_string(n)+" < "+util::to_string(m_size)+"=m_size");
				return;
			}
            m_sigma = 0; // init sigma
            temp_write_read_buffer<> buf1(5000000, buf.width, dir);   // buffer for elements in the right node
            int_vector<int_width> rac(m_size, 0, buf.width);		  // initialize rac

            value_type x = 1;  // variable for the biggest value in rac
            for (size_type i=0,r=0,r_sum=0; i < m_size;) { // detect the largest value in rac
				if ( r_sum + r > m_size ){ // read not more than size chars in the next loop
					r = m_size - r_sum;
				}
                for (; i < r+r_sum; ++i) {
                    if (buf[i-r_sum] > x)
                        x = buf[i-r_sum];
                    rac[i] = buf[i-r_sum];
                }
                r_sum += r; r = buf.load_next_block();
            }

            if (max_depth == 0) {
                m_max_depth	= bit_magic::l1BP(x)+1; // we need max_depth bits to represent all values in the range [0..x]
            } else {
                m_max_depth = max_depth;
            }
			init_buffers(m_max_depth);

            std::string tree_out_buf_file_name = (dir+"m_tree"+util::to_string(util::pid())+"_"+util::to_string(util::id()));
            std::ofstream tree_out_buf(tree_out_buf_file_name.c_str(),
                                       std::ios::binary | std::ios::trunc | std::ios::out);   // open buffer for tree
            size_type bit_size = m_size*m_max_depth;
            tree_out_buf.write((char*) &bit_size, sizeof(bit_size));    // write size of bit_vector

            size_type tree_pos = 0;
            uint64_t tree_word = 0;

            uint64_t		mask_old = 1ULL<<(m_max_depth);
            for (uint32_t k=0; k<m_max_depth; ++k) {
                size_type	  	start	 = 0;
                const uint64_t	mask_new = 1ULL<<(m_max_depth-k-1);
                do {
                    buf1.reset();
                    size_type	i 		= start;
                    size_type	cnt0	=	0;
                    uint64_t	start_value = (rac[i]&mask_old);
                    uint64_t	x;
                    while (i < m_size and((x=rac[i])&mask_old)==start_value) {
                        if (x&mask_new) {
                            tree_word |= (1ULL << (tree_pos&0x3FULL));
                            buf1 << x;
                        } else {
                            rac[start + cnt0++ ] = x;
                        }
                        ++tree_pos;
                        if ((tree_pos & 0x3FULL) == 0) { // if tree_pos % 64 == 0 write old word
                            tree_out_buf.write((char*) &tree_word, sizeof(tree_word));
                            tree_word = 0;
                        }
                        ++i;
                    }
                    buf1.write_close();
                    size_type cnt1 = i-start-cnt0;
                    if (k+1 < m_max_depth) { // inner node
                        for (i=start + cnt0, start = start+cnt0+cnt1; i < start; ++i) {
                            buf1 >> x;
                            rac[ i ] = x;
                        }
                    } else { // leaf node
                        start += cnt0+cnt1;
                        ++m_sigma; // increase sigma for each leaf
                    }

                } while (start < m_size);
                mask_old += mask_new;
            }
            if ((tree_pos & 0x3FULL) != 0) { // if tree_pos % 64 > 0 => there are remaining entries we have to write
                tree_out_buf.write((char*) &tree_word, sizeof(tree_word));
            }
            tree_out_buf.close();
            rac.resize(0);
            bit_vector tree;
            util::load_from_file(tree, tree_out_buf_file_name);
            std::remove(tree_out_buf_file_name.c_str());
            util::assign(m_tree, tree);
            util::init_support(m_tree_rank, &m_tree);
            util::init_support(m_tree_select0, &m_tree);
            util::init_support(m_tree_select1, &m_tree);
        }

        //! Copy constructor
        wt_int(const wt_int& wt):sigma(m_sigma), tree(m_tree) {
            copy(wt);
        }

        //! Assignment operator
        wt_int& operator=(const wt_int& wt) {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Swap operator
        void swap(wt_int& wt) {
            if (this != &wt) {
                std::swap(m_size, wt.m_size);
                std::swap(m_sigma,  wt.m_sigma);
                m_tree.swap(wt.m_tree);
				util::swap_support(m_tree_rank, wt.m_tree_rank, &m_tree, &(wt.m_tree) );
				util::swap_support(m_tree_select1, wt.m_tree_select1, &m_tree, &(wt.m_tree) );
				util::swap_support(m_tree_select0, wt.m_tree_select0, &m_tree, &(wt.m_tree) );
                std::swap(m_max_depth,  wt.m_max_depth);
				m_path_off.swap(wt.m_path_off);
				m_path_rank_off.swap(wt.m_path_rank_off);
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

        //! Recovers the i-th symbol of the original vector.
        /*! \param i The index of the symbol in the original vector. \f$i \in [0..size()-1]\f$
         *	\returns The i-th symbol of the original vector.
         */
        value_type operator[](size_type i)const {
			assert( i < size() );
            size_type offset = 0;
            value_type res = 0;
            size_type node_size = m_size;
            for (uint32_t k=0; k < m_max_depth; ++k) {
                res <<= 1;
                size_type ones_before_o	  = m_tree_rank(offset);
                size_type ones_before_i   = m_tree_rank(offset + i) - ones_before_o;
                size_type ones_before_end = m_tree_rank(offset + node_size) - ones_before_o;
                if (m_tree[offset+i]) { // one at position i => follow right child
                    offset += (node_size - ones_before_end);
                    node_size = ones_before_end;
                    i = ones_before_i;
                    res |= 1;
                } else { // zero at position i => follow left child
                    node_size = (node_size - ones_before_end);
                    i = (i-ones_before_i);
                }
                offset += m_size;
            }
            return res;
        };

        //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
         *  \param c The symbol to count the occurrences in the prefix.
         *	\returns The number of occurrences of symbol c in the prefix [0..i-1] of the supported vector.
         *  \par Time complexity
         *		\f$ \Order{\log |\Sigma|} \f$
         */
        size_type rank(size_type i, value_type c)const {
			assert( i <= size() );
            size_type offset = 0;
            uint64_t mask	 = (1ULL) << (m_max_depth-1);
            size_type node_size = m_size;
            for (uint32_t k=0; k < m_max_depth and i; ++k) {
                size_type ones_before_o	  = m_tree_rank(offset);
                size_type ones_before_i   = m_tree_rank(offset + i) - ones_before_o;
                size_type ones_before_end = m_tree_rank(offset + node_size) - ones_before_o;
                if (c & mask) { // search for a one at this level
                    offset += (node_size - ones_before_end);
                    node_size = ones_before_end;
                    i = ones_before_i;
                } else { // search for a zero at this level
                    node_size = (node_size - ones_before_end);
                    i = (i-ones_before_i);
                }
                offset += m_size;
                mask >>= 1;
            }
            return i;
        };



        //! Calculates how many occurrences of symbol wt[i] are in the prefix [0..i-1] of the original sequence.
        /*!
         *	\param i The index of the symbol.
         *  \param c Reference that will contain symbol wt[i].
         *  \return The number of occurrences of symbol wt[i] in the prefix [0..i-1]
         */
		size_type inverse_select(size_type i, value_type&c )const{
			assert( i < size() );
			c = (*this)[i];
			return rank(i, c);
		}

        //! Calculates the i-th occurrence of the symbol c in the supported vector.
        /*!
         *  \param i The i-th occurrence. \f$i\in [1..rank(size(),c)]\f$.
         *  \param c The symbol c.
         *  \par Time complexity
         *		\f$ \Order{\log |\Sigma|} \f$
         */
        size_type select(size_type i, value_type c)const {
			assert( i > 0);
			assert( i <= rank(size(), c) );
            // possible optimization: if the array is a permutation we can start at the bottom of the tree
            size_type offset = 0;
            uint64_t mask	 = (1ULL) << (m_max_depth-1);
            size_type node_size = m_size;
            m_path_off[0] = m_path_rank_off[0] = 0;

            for (uint32_t k=0; k < m_max_depth and node_size; ++k) {
                size_type ones_before_o	  = m_tree_rank(offset);
                m_path_rank_off[k] = ones_before_o;
                size_type ones_before_end = m_tree_rank(offset + node_size) - ones_before_o;
                if (c & mask) { // search for a one at this level
                    offset += (node_size - ones_before_end);
                    node_size = ones_before_end;
                } else { // search for a zero at this level
                    node_size = (node_size - ones_before_end);
                }
                offset += m_size;
                m_path_off[k+1] = offset;
                mask >>= 1;
            }
            if (node_size < i) {
				throw std::logic_error("select("+util::to_string(i)+","+util::to_string(c)+"): c does not occur i times in the WT");
                return m_size;
            }
            mask = 1ULL;
            for (uint32_t k=m_max_depth; k>0; --k) {
                offset = m_path_off[k-1];
                size_type ones_before_o = m_path_rank_off[k-1];
                if (c & mask) { // right child => search i'th
                    i = m_tree_select1(ones_before_o + i) - offset + 1;
                } else { // left child => search i'th zero
                    i = m_tree_select0(offset - ones_before_o + i) - offset + 1;
                }
                mask <<= 1;
            }
            return i-1;
        };

        //! range_search_2d searches points in the index interval [lb..rb] and value interval [vlb..vrb].
        /*! \param lb         Left bound of index interval (inclusive)
         *  \param rb         Right bound of index interval (inclusive)
         *  \param vlb        Left bound of value interval (inclusive)
         *  \param vrb        Right bound of value interval (inclusive)
         *  \param idx_result Reference to a vector to which the resulting indices should be added
         *  \param val_result Reference to a vector to which the resulting values should be added
         */
        size_type range_search_2d(size_type lb, size_type rb, value_type vlb, value_type vrb,
                                  std::vector<size_type>* idx_result=NULL,
                                  std::vector<value_type>* val_result=NULL
                                 ) const {
            size_type offsets[m_max_depth+1];
            size_type ones_before_os[m_max_depth+1];
            offsets[0] = 0;
            if (vrb > (1ULL << m_max_depth))
                vrb = (1ULL << m_max_depth);
            if (vlb > vrb)
                return 0;
            size_type cnt_answers = 0;
            _range_search_2d(lb, rb, vlb, vrb, 0, 0, m_size, offsets, ones_before_os, 0, idx_result, val_result, cnt_answers);
            return cnt_answers;
        }

        // add parameter path
        // ilb interval left bound
        // irb interval right bound
        void _range_search_2d(size_type lb, size_type rb, value_type vlb, value_type vrb, size_type depth,
                              size_type ilb, size_type node_size, size_type offsets[], size_type ones_before_os[], size_type path,
                              std::vector<size_type>* idx_result, std::vector<size_type>* val_result, size_type& cnt_answers)
        const {
            if (lb > rb)
                return;
            if (depth == m_max_depth) {
                if (idx_result != NULL) {
                    for (size_type j=1; j <= node_size; ++j) {
                        size_type i = j;
                        size_type c = path;
                        for (uint32_t k=m_max_depth; k>0; --k) {
                            size_type offset = offsets[k-1];
                            size_type ones_before_o = ones_before_os[k-1];
                            if (c&1) {
                                i = m_tree_select1(ones_before_o + i) - offset + 1;
                            } else {
                                i = m_tree_select0(offset - ones_before_o + i) - offset + 1;
                            }
                            c >>= 1;
                        }
                        idx_result->push_back(i-1); // add resulting index; -1 cause of 0 based indexing
                    }
                }
                if (val_result != NULL) {
                    for (size_type j=1; j <= node_size; ++j) {
                        val_result->push_back(path);
                    }
                }
                cnt_answers += node_size;
                return;
            }
            size_type irb = ilb + (1ULL << (m_max_depth-depth));
            size_type mid = (irb + ilb)>>1;

            size_type offset = offsets[depth];

            size_type ones_before_o		= m_tree_rank(offset);
            ones_before_os[depth]		= ones_before_o;
            size_type ones_before_lb 	= m_tree_rank(offset + lb);
            size_type ones_before_rb 	= m_tree_rank(offset + rb + 1);
            size_type ones_before_end	= m_tree_rank(offset + node_size);
            size_type zeros_before_o	= offset - ones_before_o;
            size_type zeros_before_lb	= offset + lb - ones_before_lb;
            size_type zeros_before_rb	= offset + rb + 1 - ones_before_rb;
            size_type zeros_before_end	= offset + node_size - ones_before_end;
            if (vlb < mid and mid) {
                size_type nlb 			= zeros_before_lb - zeros_before_o;
                size_type nrb 			= zeros_before_rb - zeros_before_o;
                offsets[depth+1] 		= offset + m_size;
                if (nrb)
                    _range_search_2d(nlb, nrb-1, vlb, std::min(vrb,mid-1), depth+1, ilb, zeros_before_end - zeros_before_o, offsets, ones_before_os, path<<1, idx_result, val_result, cnt_answers);
            }
            if (vrb >= mid) {
                size_type nlb			= ones_before_lb - ones_before_o;
                size_type nrb			= ones_before_rb - ones_before_o;
                offsets[depth+1]		= offset + m_size + (zeros_before_end - zeros_before_o);
                if (nrb)
                    _range_search_2d(nlb, nrb-1, std::max(mid, vlb), vrb, depth+1, mid, ones_before_end - ones_before_o, offsets, ones_before_os, (path<<1)+1 ,idx_result, val_result, cnt_answers);
            }
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += util::write_member(m_size, out, child, "size");
            written_bytes += util::write_member(m_sigma, out, child, "sigma");
            written_bytes += m_tree.serialize(out, child, "tree");
            written_bytes += m_tree_rank.serialize(out, child, "tree_rank");
            written_bytes += m_tree_select1.serialize(out, child, "tree_select_1");
            written_bytes += m_tree_select0.serialize(out, child, "tree_select_0");
            written_bytes += util::write_member(m_max_depth, out, child, "max_depth");
            structure_tree::add_size(child, written_bytes);
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
            util::read_member(m_max_depth, in);
			init_buffers(m_max_depth);
        }
};

}// end namespace sdsl

#endif // end file 
