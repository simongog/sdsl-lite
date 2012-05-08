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

#include "int_vector.hpp"
#include "rank_support_v.hpp"
#include "select_support_mcl.hpp"
#include "bitmagic.hpp"
#include "testutils.hpp"
#include "temp_write_read_buffer.hpp"
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
 * The wavelet tree \f$wt\f$ consits of a tree of bit vector and provides three efficient methods:
 *   - The "[]"-operator: \f$wt[i]\f$ returns the ith symbol of vector for which the wavelet tree was build for.
 *   - The rank method: \f$wt.rank(i,c)\f$ returns the number of occurences of symbol \f$c\f$ in the prefix [0..i-1] in the vector for which the wavelet tree was build for.
 *   - The select method: \f$wt.select(j,c)\f$ returns the index \f$i\in [0..size()-1]\f$ of the jth occurence of symbol \f$c\f$.
 *
 *	\par Space complexity
 *		\f$\Order{n\log|\Sigma|}\f$ bits, where \f$n\f$ is the size of the vector the wavelet tree was build for.
 *
 *   @ingroup wt
 */
template<class RandomAccessContainer=int_vector<>, class RankSupport = rank_support_v<>, class SelectSupport=select_support_mcl<>, class SelectSupportZero=select_support_mcl<0> >
class wt_int
{
    public:
        typedef typename RandomAccessContainer::size_type 		size_type;
        typedef typename RandomAccessContainer::value_type 		value_type;
    protected:
        size_type 			m_size;
        size_type 			m_sigma; 		//<- \f$ |\Sigma| \f$
        bit_vector 			m_tree;			// bit vector to store the wavelet tree
        RankSupport			m_tree_rank;	// rank support for the wavelet tree bit vector
        SelectSupport		m_tree_select1;	// select support for the wavelet tree bit vector
        SelectSupportZero	m_tree_select0;
        uint32_t			m_logn;

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
            m_logn			= wt.m_logn;
        }

    public:

        const size_type& sigma;	//!< Effective alphabet size of the wavelet tree.
        const bit_vector& tree; //!< A concatenation of all bit vectors of the wavelet tree.

        //! Default constructor
        wt_int():m_size(0),m_sigma(0), m_logn(0), sigma(m_sigma), tree(m_tree) {};


        //! Constructor
        /*!	\param rac Reference a random access container of integer values for which the wavelet tree should be build.
         *	\param logn Let x > 0 be the biggest value in rac. logn should be bit_magic::l1BP(x-1)+1 to represent all values of rac.
         *	\par Time complexity
         *		\f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         *		I.e. we nee \Order{n\log n} if rac is a permutation of 0..n-1.
         *	\par Space complexity
         *		\f$ 3n\log|\Sigma| \f$ bits, where \f$n=size\f$.
         *		I.e. we need \f$3n\log n \f$ if rac is a permutation of 0..n-1.
         */
        wt_int(const RandomAccessContainer& rac, uint32_t logn=0):m_size(rac.size()), m_sigma(0), sigma(m_sigma), tree(m_tree) {
#ifdef SDSL_DEBUG_INT_WAVELET_TREE
            std::cerr<<"wt_int construct: size="<<m_size<<std::endl;
#endif
            size_type n	=	m_size;
            if (logn == 0) { // detect the biggest value in rac and set logn properly
                value_type x = 1;
                for (size_type i=0; i < n; ++i)
                    if (rac[i] > x)
                        x = rac[i];
                m_logn	= bit_magic::l1BP(x)+1; // we need logn bits to represent all values in the range [0..x]
            } else {
                m_logn = logn;
            }
            m_tree		=	bit_vector(n*m_logn, 0);  // initialize the tree

            int_vector<> 	perms[2];
            perms[0].set_int_width(m_logn); perms[1].set_int_width(m_logn);
            perms[0].resize(n);	perms[1].resize(n);

            for (size_type i=0; i < n; ++i) { // copy original array to perms[0]
                perms[0][i] = rac[i];
            }

            size_type tree_pos = 0;

            uint64_t		mask_old = 1ULL<<(m_logn);
            for (uint32_t k=0; k<m_logn; ++k) {
                int_vector<>& 	act_perm	 = perms[k%2];
                int_vector<>& 	next_perm	 = perms[1-(k%2)];
                size_type	  	start		 = 0;
                const uint64_t		mask_new = 1ULL<<(m_logn-k-1);
                do {
                    size_type	i 		= start;
                    size_type	cnt1	=	0;
                    uint64_t	start_value = (act_perm[i]&mask_old);
                    uint64_t	x;
                    while (i < n and((x=act_perm[i])&mask_old)==start_value) {
                        if (x&mask_new) {
                            ++cnt1; m_tree[tree_pos] = 1;
                        }
                        ++tree_pos;
                        ++i;
                    }
                    size_type cnt0 = i-start-cnt1;
                    if (k+1 < m_logn) { // inner node
                        size_type idxs[2] = {start, start+cnt0};
                        for (i=start, start = start+cnt0+cnt1; i < start; ++i) {
                            next_perm[idxs[(act_perm[i]&mask_new)>0 ]++] = act_perm[i];
                        }
                    } else { // leaf node
                        start += cnt0+cnt1;
                        ++m_sigma;  // increase sigma for each leaf
                    }
                } while (start < n);
                mask_old += mask_new;
            }
#ifdef SDSL_DEBUG_INT_WAVELET_TREE
            if (m_tree.size()<100) {
                std::cerr<<"tree="<<m_tree<<std::endl;
            }
#endif
            m_tree_rank.init(&m_tree);
            m_tree_select0.init(&m_tree);
            m_tree_select1.init(&m_tree);
        }

        //! Semi-external constructor
        /*!	\param buf	int_vector_file_buffer which contains the vector v for which a wt_int should be bould.
         *	\param logn Let x > 0 be the biggest value in v. logn should be bit_magic::l1BP(x-1)+1 to represent all values of v.
         *	\param dir	Derectory in which temporary files should be stored during the construction.
         *	\par Time complexity
         *		\f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         *		I.e. we nee \Order{n\log n} if rac is a permutation of 0..n-1.
         *	\par Space complexity
         *		\f$ n\log|\Sigma| + O(1)\f$ bits, where \f$n=size\f$.
         */
        template<uint8_t int_width>
        wt_int(int_vector_file_buffer<int_width> &buf, uint32_t logn=0, std::string dir="./")
            :m_size(0),m_sigma(0), m_logn(0), sigma(m_sigma), tree(m_tree) {
            construct(buf, logn, dir);
        }


        // TODO external construction method which occupies only additional constant memory
        template<uint8_t int_width>
        void construct(int_vector_file_buffer<int_width> &buf, uint32_t logn=0, std::string dir="./") {
            buf.reset();
            size_type n = buf.int_vector_size;  // set n
            m_size = n;				// set sigma and size
            m_sigma = 0;
            temp_write_read_buffer<> buf1(5000000, buf.int_width, dir);   // buffer for elements in the rigth node
            int_vector<int_width> rac(n, 0, buf.int_width);				  // initialze rac

            value_type x = 1;  // variable for the biggest value in rac
            for (size_type i=0,r=0,r_sum=0; i<n;) { // detect the biggest value in rac
                for (; i < r+r_sum; ++i) {
                    if (buf[i-r_sum] > x)
                        x = buf[i-r_sum];
                    rac[i] = buf[i-r_sum];
                }
                r_sum += r; r = buf.load_next_block();
            }

            if (logn == 0) {
                m_logn	= bit_magic::l1BP(x)+1; // we need logn bits to represent all values in the range [0..x]
            } else {
                m_logn = logn;
            }
            std::string tree_out_buf_file_name = (dir+"m_tree"+util::to_string(util::get_pid())+"_"+util::to_string(util::get_id()));
            std::ofstream tree_out_buf(tree_out_buf_file_name.c_str(),
                                       std::ios::binary | std::ios::trunc | std::ios::out);   // open buffer for tree
            size_type bit_size = n*m_logn;
            tree_out_buf.write((char*) &bit_size, sizeof(bit_size));    // write size of bit_vector

            size_type tree_pos = 0;
            uint64_t tree_word = 0;

            uint64_t		mask_old = 1ULL<<(m_logn);
            for (uint32_t k=0; k<m_logn; ++k) {
                size_type	  	start		 = 0;
                const uint64_t		mask_new = 1ULL<<(m_logn-k-1);
                do {
                    buf1.reset();
                    size_type	i 		= start;
                    size_type	cnt0	=	0;
                    uint64_t	start_value = (rac[i]&mask_old);
                    uint64_t	x;
                    while (i < n and((x=rac[i])&mask_old)==start_value) {
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
                    if (k+1 < m_logn) { // inner node
                        for (i=start + cnt0, start = start+cnt0+cnt1; i < start; ++i) {
                            buf1 >> x;
                            rac[ i ] = x;
                        }
                    } else { // leaf node
                        start += cnt0+cnt1;
                        ++m_sigma; // increase sigma for each leaf
                    }

                } while (start < n);
                mask_old += mask_new;
            }
            if ((tree_pos & 0x3FULL) != 0) { // if tree_pos % 64 > 0 => there are remaining entries we have to write
                tree_out_buf.write((char*) &tree_word, sizeof(tree_word));
            }
            tree_out_buf.close();
            rac.resize(0);
            util::load_from_file(m_tree, tree_out_buf_file_name.c_str());
            std::remove(tree_out_buf_file_name.c_str());
#ifdef SDSL_DEBUG_INT_WAVELET_TREE
            if (m_tree.size()<100) {
                std::cerr<<"tree="<<m_tree<<std::endl;
            }
#endif
            m_tree_rank.init(&m_tree);
            m_tree_select0.init(&m_tree);
            m_tree_select1.init(&m_tree);
        }

        //! Copy constructor
        wt_int(const wt_int& wt):sigma(m_sigma) {
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
                m_tree_rank.swap(wt.m_tree_rank); // rank swap after the swap of the bit vector m_tree
                m_tree_select1.swap(wt.m_tree_select1); // select1 swap after the swap of the bit vector m_tree
                m_tree_select0.swap(wt.m_tree_select0); // select0 swap after the swap of the bit vector m_tree
                std::swap(m_logn,  wt.m_logn);
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

        //! Recovers the ith symbol of the original vector.
        /*! \param i The index of the symbol in the original vector. \f$i \in [0..size()-1]\f$
         *	\returns The ith symbol of the original vector.
         */
        value_type operator[](size_type i)const {
//		size_type zeros_left_of_pos = 0;
            size_type offset = 0;
            value_type res = 0;
            size_type node_size = m_size;
            for (uint32_t k=0; k < m_logn; ++k) {
                res <<= 1;
                size_type ones_before_o	  = m_tree_rank(offset);
                size_type ones_before_i   = m_tree_rank(offset + i) - ones_before_o;
                size_type ones_before_end = m_tree_rank(offset + node_size) - ones_before_o;
//			std::cerr<<"k="<<k<<" i="<<i<<" offset="<<offset<<" node_size="<<node_size<<" obo="<<ones_before_o<<" obi="<<ones_before_i<<" obe="<<ones_before_end<<std::endl;
                if (m_tree[offset+i]) { // one at position i => follow right child
//				zeros_left_of_pos += (node_size - ones_before_end);
                    offset += (node_size - ones_before_end);
                    node_size = ones_before_end;
                    i = ones_before_i;
                    res |= 1;
                } else { // zero at position i => follow left child
                    node_size = (node_size - ones_before_end);
                    i = (i-ones_before_i);
                }
                //offset = (k+1)*m_size + zeros_left_of_pos;
                offset += m_size;
            }
//		std::cerr<<" offset="<<offset<<" node_size="<<node_size<<std::endl;
            return res;
        };

        //TODO: buchstaben, die nicht im orginaltext vorkommen behandeln!!!
        //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
         *  \param c The symbol to count the occurences in the prefix.
         *	\returns The number of occurences of symbol c in the prefix [0..i-1] of the supported vector.
         *  \par Time complexity
         *		\f$ \Order{\log |\Sigma|} \f$
         */
        size_type rank(size_type i, value_type c)const {
            size_type offset = 0;
            uint64_t mask	 = (1ULL) << (m_logn-1);
            size_type node_size = m_size;
            for (uint32_t k=0; k < m_logn and i; ++k) {
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

        //! Calculates the ith occurence of the symbol c in the supported vector.
        /*!
         *  \param i The ith occurence. \f$i\in [1..rank(size(),c)]\f$.
         *  \param c The symbol c.
         *  \par Time complexity
         *		\f$ \Order{\log |\Sigma|} \f$
         */
        // TODO: was ist wenn c gar nicht vorkommt, oder es keine i Vorkommen gibt?
        size_type select(size_type i, value_type c)const {
            // possible optimization: if the array is a permutation we can start at the bottom of the tree
            size_type offset = 0;
            uint64_t mask	 = (1ULL) << (m_logn-1);
            size_type node_size = m_size;
            size_type offsets[m_logn+1]; 			// offsets down the path
            size_type ones_before_os[m_logn+1];   // ones before the offsets
            offsets[0] = ones_before_os[0] = 0;

            for (uint32_t k=0; k < m_logn and node_size; ++k) {
                size_type ones_before_o	  = m_tree_rank(offset);
                ones_before_os[k] = ones_before_o;
                size_type ones_before_end = m_tree_rank(offset + node_size) - ones_before_o;
                if (c & mask) { // search for a one at this level
                    offset += (node_size - ones_before_end);
                    node_size = ones_before_end;
                } else { // search for a zero at this level
                    node_size = (node_size - ones_before_end);
                }
                offset += m_size;
                offsets[k+1] = offset;
                mask >>= 1;
            }
            if (node_size < i) {
                std::cerr<<"c="<<c<<" does not occure "<<i<<" times in the WT"<<std::endl;
                return m_size;
            }
            mask = 1ULL;
            for (uint32_t k=m_logn; k>0; --k) {
                offset = offsets[k-1];
                size_type ones_before_o = ones_before_os[k-1];
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
        /*! \param lb Left bound of index interval (inclusive)
         *  \param rb Right bound of index interval (inclusive)
         *  \param vlb Left bound of value interval (inclusive)
         *  \param vrb Right bound of value interval (inclusive)
         *  \param idx_result Reference to a vector to which the resulting indices should be added
         *  \param val_result Reference to a vector to which the resulting values should be added
         */
        size_type range_search_2d(size_type lb, size_type rb, value_type vlb, value_type vrb,
                                  std::vector<size_type> *idx_result=NULL,
                                  std::vector<value_type> *val_result=NULL
                                 ) const {
            size_type offsets[m_logn+1];
            size_type ones_before_os[m_logn+1];
            offsets[0] = 0;
            if (vrb > (1ULL << m_logn))
                vrb = (1ULL << m_logn);
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
                              std::vector<size_type> *idx_result, std::vector<size_type> *val_result, size_type& cnt_answers)
        const {
            if (lb > rb)
                return;
            if (depth == m_logn) {
                if (idx_result != NULL) {
                    for (size_type j=1; j <= node_size; ++j) {
                        size_type i = j;
                        size_type c = path;
                        for (uint32_t k=m_logn; k>0; --k) {
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
            size_type irb = ilb + (1ULL << (m_logn-depth));
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
        size_type serialize(std::ostream& out)const {
            size_type written_bytes = 0;
            out.write((char*)&m_size, sizeof(m_size));
            written_bytes += sizeof(m_size);
            out.write((char*)&m_sigma, sizeof(m_sigma));
            written_bytes += sizeof(m_sigma);
            written_bytes += m_tree.serialize(out);
            written_bytes += m_tree_rank.serialize(out);
            written_bytes += m_tree_select1.serialize(out);
            written_bytes += m_tree_select0.serialize(out);
            written_bytes += sizeof(m_logn);
            out.write((char*)&m_logn, sizeof(m_logn));
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            in.read((char*) &m_size, sizeof(m_size));
            in.read((char*) &m_sigma, sizeof(m_sigma));
            m_tree.load(in);
            m_tree_rank.load(in, &m_tree);
            m_tree_select1.load(in, &m_tree);
            m_tree_select0.load(in, &m_tree);
            in.read((char*) &m_logn, sizeof(m_logn));
        }
        /*
        	void print_info()const{
        		size_type rle_ed = 0;
        		for(size_type i=0; i < m_tree.size(); ++i){

        		}
        	}
        */

};

}// end namespace sds

#endif // end file 
