/* sdsl - succinct data structures library
    Copyright (C) 2011 Simon Gog

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
/*! \file wt_rlg8.hpp
    \brief wt_rlg8.hpp contains a class for the wavelet tree of byte sequences which is in Huffman shape and runs of character
	       are compressed.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_WT_RLG8
#define INCLUDED_SDSL_WT_RLG8

#include "int_vector.hpp"
#include "rank_support_v.hpp"
#include "rank_support_v5.hpp"
#include "select_support_mcl.hpp"
#include "select_support_bs.hpp"
#include "bitmagic.hpp"
#include "util.hpp"
#include "wt_huff.hpp"
#include <algorithm> // for std::swap
#include <stdexcept>
#include <vector>
#include <utility> // for pair
#include <queue>


#ifdef SDSL_DEBUG
#define SDSL_DEBUG_WAVELET_TREE_HUFFMAN_RLN4
#endif


//! Namespace for the succinct data structure library.
namespace sdsl
{

typedef wt_huff<bit_vector, rank_support_v5<>, select_support_bs<>, select_support_bs<> > wt_without_select;

//! A Wavelet Tree class for byte sequences.
/*!
 * A wavelet tree is build for a vector of characters over the alphabet \f$\Sigma\f$.
 * This class should be used only for small alphabets \f$\Sigma \ll n\f$ (see wt_int for a wavelet tree for big alphabets).
 * The wavelet tree \f$wt\f$ consists of a tree of bitvectors and provides three efficient methods:
 *   - The "[]"-operator: \f$wt[i]\f$ returns the ith symbol of vector for which the wavelet tree was build for.
 *   - The rank method: \f$wt.rank(i,c)\f$ returns the number of occurences of symbol \f$c\f$ in the prefix [0..i-1] in the vector for which the wavelet tree was build for.
 *   - The select method: \f$wt.select(j,c)\f$ returns the index \f$i\in [0..size()-1]\f$ of the jth occurence of symbol \f$c\f$.
 *
 *	\par Space complexity
 *		 \f$\Order{n\log|\Sigma| + 2|\Sigma|\log n}\f$ bits, where \f$n\f$ is the size of the vector the wavelet tree was build for.
 *
 *  \par Note
 *       We denote the length of the longest run in the sequence with \f$ L \f$
 */
template<class RankSupport = rank_support_v5<>, class WaveletTree = wt_without_select >
class wt_rlg8
{
    public:
        typedef int_vector<>::size_type	size_type;
        typedef unsigned char		 	value_type;
        typedef RankSupport				rank_support_type;
        typedef WaveletTree             wt_type;

    private:
        size_type 				m_size;         // size of the original input sequence
        wt_type					m_wt;	        // wavelet tree for all levels
        bit_vector				m_b;	        // bit vector which indicates if a pair consists of
        // two equal chars
        rank_support_type		m_b_rank;       // rank support for vector b
        int_vector<64>          m_b_border_rank;// Vector in which we store the rank values of m_b at the
        // border positions.
        int_vector<64>          m_b_border;       // Vector in which we store the borders of the different levels
        // Takes \f$\Order{\max(1, \log L)\log n}\f$ bits.
        int_vector<64>          m_wt_rank;  // Vector in which we store the rank value for each character
        // and each border.
        // Takes \f$\Order{\sigma\max(1, \log L)\log n}\f bits
        int_vector<8>			m_char2comp;    //
        int_vector<64>          m_char_occ;     //
        uint16_t				m_sigma;

        void copy(const wt_rlg8& wt) {
            m_size          = wt.m_size;
            m_wt            = wt.m_wt;
            m_b             = wt.m_b;
            m_b_rank        = wt.m_b_rank;
            m_b_rank.set_vector(&m_b);
            m_b_border_rank = wt.m_b_border_rank;
            m_b_border      = wt.m_b_border;
            m_wt_rank       = wt.m_wt_rank;
            m_char2comp     = wt.m_char2comp;
            m_char_occ      = wt.m_char_occ;
            m_size			= wt.m_size;
        }

    public:

        const uint16_t& sigma;

        // Default constructor
        wt_rlg8():m_size(0), m_sigma(0), sigma(m_sigma) {};

        // Constructor
        /*
         *	\param rac Reference to the vector (or unsigned char array) for which the wavelet tree should be build.
         *	\param size Size of the prefix of the vector (or unsigned char array) for which the wavelet tree should be build.
         *	\par Time complexity
         *		\f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         */
        wt_rlg8(const unsigned char* rac, size_type size):m_size(size), m_sigma(0), sigma(m_sigma) {
            std::cerr << "ERROR: Constructor of wt_rlg8 not implemented yet!!!" << std::endl;
            throw std::logic_error("This constructor of wt_rlg8 is not yet implemented!");
        }

        //! Construct the wavelet tree from a file_buffer
        /*! \param text_buf	A int_vector_file_buffer to the original text.
         *	\param size The length of the prefix of the text, for which the wavelet tree should be build.
         */
        wt_rlg8(int_vector_file_buffer<8>& rac, size_type size):m_size(size), m_sigma(0), sigma(m_sigma) {
            // TODO: remove absolute file name
            std::string temp_file = "wt_rlg8_" + util::to_string(util::get_pid()) + "_" + util::to_string(util::get_id());
            std::ofstream wt_out(temp_file.c_str(), std::ios::binary | std::ios::trunc);
            size_type bit_cnt=0;
            wt_out.write((char*)&bit_cnt, sizeof(bit_cnt)); // initial dummy write

            m_b = bit_vector(size/4+1,0);
            uint8_t* next_bwt = new uint8_t[size/4+4];
//		bit_vector same_prev_char(size/2+1,0); //

            m_b_border.resize(bit_magic::l1BP(size) + 1);
            m_b_border[0] = 0;

            int m=0;

            rac.reset();
            bit_vector b_sigma(256, 0);
            uint8_t last_c[9] = {0,1,0,1,0,1,0,1};
            uint8_t c = '\0';
            size_type b_cnt = 0, pair1cnt=0, pair0cnt=0;
            for (size_type i=0, r=0, r_sum=0; r_sum < size;) {
                if (r_sum + r > size) {  // read not more than size chars in the next loop
                    r = size-r_sum;
                }
                for (; i < r+r_sum; ++i) {
                    c = rac[i-r_sum];
                    b_sigma[c] = 1;
                    last_c[i&7ULL] = c;
                    if ((i & 7ULL)==7) {
                        if (last_c[0] == last_c[1] and last_c[2] == last_c[3] and last_c[0] == last_c[2] and
                            last_c[4] == last_c[5] and last_c[6] == last_c[7] and last_c[4] == last_c[6]  and
                            last_c[0] == last_c[4]
                           ) { // join octuple
                            m_b[b_cnt] = 1;
                            next_bwt[pair1cnt] = c;
                            ++pair1cnt;
                        } else { // write pair to stream
                            //m_b[b_cnt] = 0; // since m_b is initialized to zero, this is not necessary
                            wt_out.write((char*)last_c, 8*sizeof(c));
                            ++pair0cnt;
                        }
                        ++b_cnt;
                    }
                }
                r_sum += r;
                r = rac.load_next_block();
            }
            if (size & 7ULL) {
                wt_out.write((char*)&last_c, 8*sizeof(c));
                ++pair0cnt;
                ++b_cnt;
            }
            m_sigma = 0;
            for (size_type i=0; i<b_sigma.size(); ++i) {
                m_sigma += b_sigma[i];
            }

            uint32_t level = 0;
            //  handle remaining levels
            while (pair1cnt > 0) {
                ++m;
                m_b_border[++level] = b_cnt;
                size_type level_size = pair1cnt;
                pair1cnt = 0;
                for (size_type i=7; i < level_size; i+=8) {
                    if (next_bwt[i] == next_bwt[i-1] and next_bwt[i-2] == next_bwt[i-3] and next_bwt[i] == next_bwt[i-2] and
                        next_bwt[i-4] == next_bwt[i-5] and next_bwt[i-6] == next_bwt[i-7] and next_bwt[i-4] == next_bwt[i-6] and
                        next_bwt[i] == next_bwt[i-4]) {
                        m_b[b_cnt] = 1;
                        next_bwt[pair1cnt] = next_bwt[i];
                        ++pair1cnt;
                    } else {
                        //m_b[b_cnt] = 0; // since m_b is initialized to zero, this is not necessary
                        wt_out.write((char*)next_bwt+i-7, 8*sizeof(c));
                        ++pair0cnt;
                    }
                    ++b_cnt;
                }
                if (level_size & 7ULL) { // handle last element
                    wt_out.write((char*)next_bwt + (level_size/8)*8, (level_size&7ULL)*sizeof(c));
                    wt_out.write((char*)"\0\0\0\0\0\0\0\0", 8-(level_size&7ULL));
                    ++pair0cnt;
                    ++b_cnt;
                }
            }
            delete [] next_bwt;
            m_b.resize(b_cnt);
            m_b_border.resize(level+1);

            wt_out.seekp(0, std::ios::beg);
            bit_cnt = (8*pair0cnt)*8;
            wt_out.write((char*)&bit_cnt, sizeof(bit_cnt));
            wt_out.close();

            {
                int_vector_file_buffer<8> temp_bwt_buf(temp_file.c_str());
				util::assign(m_wt, wt_type(temp_bwt_buf, temp_bwt_buf.int_vector_size));
            }

            util::init_support(m_b_rank, &m_b);
            m_b_border_rank.resize(m_b_border.size());

            for (size_type i=0; i<m_b_border.size(); ++i) {
                m_b_border_rank[i] = m_b_rank.rank(m_b_border[i]);
            }

            m_char2comp = int_vector<8>(256,255);
            for (uint16_t c=0, cnt=0; c<256; ++c) {
                if (b_sigma[c]) {
                    m_char2comp[c] = cnt++;
                }
            }


            m_wt_rank.resize(m_sigma * m_b_border.size());
            m_char_occ.resize(m_sigma);
            for (size_type c=0; c < 256; ++c) {
                uint16_t cc = m_char2comp[c];
                if (cc < m_sigma) {
                    for (size_type i=0; i < m_b_border.size(); ++i) {
                        size_type zeros  = m_b_border[i] - m_b_border_rank[i];
                        m_wt_rank[cc * m_b_border.size() + i] = m_wt.rank(8*zeros, c);
                    }
                    m_char_occ[cc] = m_wt.rank(m_wt.size(), c);
                }
            }
            std::remove(temp_file.c_str());
        }

        //! Copy constructor
        wt_rlg8(const wt_rlg8& wt):sigma(m_sigma) {
            copy(wt);
        }

        //! Assignment operator
        wt_rlg8& operator=(const wt_rlg8& wt) {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Swap operator
        void swap(wt_rlg8& wt) {
            if (this != &wt) {
                std::swap(m_size, wt.m_size);
                m_wt.swap(wt.m_wt);
                m_b.swap(wt.m_b);
                util::swap_support(m_b_rank, wt.m_b_rank, &m_b, &(wt.m_b));
                m_b_border.swap(wt.m_b_border);
                m_b_border_rank.swap(wt.m_b_border_rank);
                m_wt_rank.swap(wt.m_wt_rank);
                m_char2comp.swap(wt.m_char2comp);
                m_char_occ.swap(wt.m_char_occ);
                std::swap(m_sigma, wt.m_sigma);
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
         *	\return The ith symbol of the original vector.
         *  \par Time complexity
         *		\f$ \Order{H_0 + \log L} \f$ on average, where \f$ H_0 \f$ is the zero order entropy of
         *      the sequence and \f$L\f$ the maximal length of a run in the sequence.
         */
        value_type operator[](size_type i)const {
            size_type level = 0;
            while (m_b[(i>>3) + m_b_border[level]]) {
                i = m_b_rank((i>>3) + m_b_border[level]) - m_b_border_rank[level];
                ++level;
            }
            size_type zeros = (i>>3) + m_b_border[level] - m_b_rank((i>>3) + m_b_border[level]);
            return m_wt[(zeros<<3) + (i&7ULL)];
        };

        //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
         *  \param c The symbol to count the occurrences in the prefix.
         *	\return The number of occurrences of symbol c in the prefix [0..i-1] of the supported vector.
         *  \par Time complexity
         *		\f$ \Order{H_0 \log L} \f$ on average, where \f$ H_0 \f$ is the zero order entropy of
         *      the sequence and \f$L\f$ the maximal length of a run of \f$c\f$s in the sequence.
         */
        size_type rank(size_type i, value_type c)const {
            value_type cc    = m_char2comp[c];
            if (((size_type)cc) >= m_char_occ.size()) { // char does not occur
                return 0;
            }
            size_type  res   = 0;
            size_type  level = 0;
            size_type  added = 0;
            size_type cs     = 0;
            while (i>0 and cs != m_char_occ[cc]) {
                size_type ones  = m_b_rank((i>>3) + m_b_border[level]); // # of ones till this position
                size_type zeros = m_b_border[level] + (i>>3) - ones;  // # of zeros till this position
                res += ((cs=m_wt.rank((zeros<<3), c)) - (m_wt_rank[cc*m_b_border.size() + level])) << ((level<<1)+level);
                if (i & 7ULL) {//  i not a multiple of 8
                    if (m_b[(i>>3) + m_b_border[level]]) {
                        added += (1<<((level<<1)+level)) * (8-(i & 7ULL));
                        i += (8-(i & 7ULL));
                        ++ones;
                    } else {
                        if (m_wt[(zeros<<3) + (i&7ULL) - 1] == c) {
                            res += (1<<((level<<1)+level)) - added;
                        }
                        if ((i&7ULL) > 1) {
                            size_type cnt = m_wt.rank((zeros<<3) + (i&7ULL) - 1, c) - cs;
                            res += (cnt<<((level<<1)+level));
                        }
                        added = 0;
                    }
                } else { // i is a multiple of 8
                    if (added > 0 and m_b[(i>>3) + m_b_border[level] - 1] == 0) {
                        if (m_wt[(zeros<<3)-1] == c) {
                            res -= added;
                        }
                        added = 0;
                    }
                }
                i = ones - m_b_border_rank[level];
                ++level;
            }
            return res;
        };

        //! Calculates how many occurrences of symbol wt[i] are in the prefix [0..i-1] of the supported sequence.
        /*!
         *	\param i The index of the symbol.
         *  \param c Reference that will contain the symbol at position i after the execution of the method.
         *  \return The number of occurrences of symbol wt[i] in the prefix [0..i-1]
         *	\par Time complexity
         *		\f$ \Order{H_0 \log L} \f$
         */
        size_type rank_ith_symbol(size_type i, value_type& c)const {
            return rank(i, c=(*this)[i]);
        }

        //! Calculates the ith occurrence of the symbol c in the supported vector.
        /*!
         *  \param i The ith occurrence. \f$i\in [1..rank(size(),c)]\f$.
         *  \param c The symbol c.
         *  \par Time complexity
         *		\f$ \Order{\log n H_0 \log L} \f$ on average, where \f$ H_0 \f$ is the zero order
         *      entropy of the sequence and \f$L\f$ the maximal length of a run
         *      of \f$c\f$s in the sequence.
         */
        size_type select(size_type i, value_type c)const {
            if (((size_type)m_char2comp[c]) >= m_char_occ.size()) { // char does not occur
                return size();
            }
            size_type lb = 0, rb = m_size;  // lb inclusive, rb exclusive
            while (rb > lb) {
                size_type m = (lb+rb)>>1;
                if (rank(m+1, c) < i) {
                    lb = m+1;
                } else { // rank(m+1,c) >= i
                    rb = m;
                }
            }
            return lb;
        };

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += util::write_member(m_size, out, child, "size");
            written_bytes += m_wt.serialize(out, child, "wt");
            written_bytes += m_b.serialize(out, child, "b");
            written_bytes += m_b_rank.serialize(out, child, "b_rank");
            written_bytes += m_b_border.serialize(out, child, "b_border");
            written_bytes += m_b_border_rank.serialize(out, child, "b_border_rank");
            written_bytes += m_wt_rank.serialize(out, child, "wt_rank");
            written_bytes += m_char2comp.serialize(out, child, "char2comp");
            written_bytes += m_char_occ.serialize(out, child, "char_occ");
            written_bytes += util::write_member(m_sigma, out, child, "sigma");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            util::read_member(m_size, in);
            m_wt.load(in);
            m_b.load(in);
            m_b_rank.load(in, &m_b);
            m_b_border.load(in);
            m_b_border_rank.load(in);
            m_wt_rank.load(in);
            m_char2comp.load(in);
            m_char_occ.load(in);
            util::read_member(m_sigma, in);
        }
};


}// end namespace sdsl

#endif // end file 
