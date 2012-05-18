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
/*! \file wt_rlmn.hpp
    \brief wt_rlmn.hpp contains a class for the wavelet tree of byte sequences which is in Huffman shape and runs of character are compressed.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_WT_RLMN
#define INCLUDED_SDSL_WT_RLMN

#include "int_vector.hpp"
#include "rank_support_v.hpp"
#include "rank_support_v5.hpp"
#include "select_support_mcl.hpp"
#include "bitmagic.hpp"
#include "util.hpp"
#include "wt_huff.hpp"
#include <algorithm> // for std::swap
#include <stdexcept>
#include <vector>
#include <utility> // for pair
#include <queue>
#include <iostream>


#ifdef SDSL_DEBUG
#define SDSL_DEBUG_WAVELET_TREE_HUFFMAN_RL
#endif


//! Namespace for the succinct data structure library.
namespace sdsl
{



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
 *		 \f$ nH_0 + 2|\Sigma|\log n + 2n + o(n) \f$ bits, where \f$n\f$ is the size of the vector the wavelet tree was build for.
 *
 *  \par Note
 *       This implementation is based on the idea of Veli Makinen and Gonzalo Navarro presented in the paper
 *       "Succint Suffix Arrays Based on Run-Length Encoding" (CPM 2005)
 *
 *   @ingroup wt
 *
 * TODO: make it possible to replace the bit_vector class for m_bl and m_bf by the rrr_vector class
 */
template<class BitVector = bit_vector, class RankSupport = rank_support_v5<>, class SelectSupport = select_support_mcl<>, class WaveletTree = wt_huff<> >
class wt_rlmn
{
    public:
        typedef int_vector<>::size_type	size_type;
        typedef unsigned char		 	value_type;
        typedef BitVector				bit_vector_type;
        typedef RankSupport				rank_support_type;
        typedef SelectSupport           select_support_type;
        typedef WaveletTree             wt_type;

    private:
        size_type 				m_size;         // size of the original input sequence
        bit_vector_type			m_bl;	        // bit vector which indicates the starts of runs in
        // the BWT (or last column), i.e. _b_ _l_ast
        bit_vector_type         m_bf;           // bit vector which indicates the starts of runs in
        // the first column of the sorted suffixes, i.e _b_ _f_irst
        wt_type					m_wt;	        // wavelet tree for all levels
        // two equal chars
        rank_support_type		m_bl_rank;      // rank support for bit vector bl
        rank_support_type		m_bf_rank;	    // rank support for bit vector bf
        select_support_type     m_bl_select;    // select support for bit vector bl
        select_support_type     m_bf_select;    // select support for bit vector bf
        int_vector<64>          m_C;     		//
        int_vector<64>          m_C_bf_rank;    // stores the number of 1s in m_bf for the prefixes
        // m_bf[0..m_C[0]],m_bf[0..m_C[1]],....,m_bf[0..m_C[255]];
        // named C_s in the original paper

        void copy(const wt_rlmn& wt) {
            m_size          = wt.m_size;
            m_bl            = wt.m_bl;
            m_bf            = wt.m_bf;
            m_wt            = wt.m_wt;
            m_bl_rank       = wt.m_bl_rank;
            m_bl_rank.set_vector(&m_bl);
            m_bf_rank       = wt.m_bf_rank;
            m_bf_rank.set_vector(&m_bf);
            m_bl_select     = wt.m_bl_select;
            m_bl_select.set_vector(&m_bl);
            m_bf_select     = wt.m_bf_select;
            m_bf_select.set_vector(&m_bf);
            m_C             = wt.m_C;
            m_C_bf_rank     = wt.m_C_bf_rank;
        }

    public:

        const size_type& sigma;

        // Default constructor
        wt_rlmn():m_size(0), sigma(m_wt.sigma) {};



        //! Constructor
        /*!
         *	\param rac Reference to the vector (or unsigned char array) for which the wavelet tree should be build.
         *	\param size Size of the prefix of the vector (or unsigned char array) for which the wavelet tree should be build.
         *	\par Time complexity
         *		\f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         */
        wt_rlmn(const unsigned char* rac, size_type size):m_size(size), sigma(m_wt.sigma) {
            // TODO
            std::cerr << "ERROR: Constructor of wt_rlmn not implemented yet!!!" << std::endl;
            throw std::logic_error("This constructor of wt_rlmn is not yet implemented!");
        }

        template<class size_type_class>
        wt_rlmn(int_vector_file_buffer<8, size_type_class>& rac, size_type size):m_size(size), sigma(m_wt.sigma) {
            construct(rac, size);
        }

        //! Construct the wavelet tree from a random access container
        /*! \param rac A random access container
         *	\param size The length of the prefix of the random access container, for which the wavelet tree should be build
         */
        template<class size_type_class>
        void construct(int_vector_file_buffer<8, size_type_class>& rac, size_type size) {
            m_size = size;
            typedef size_type_class size_type;
            // TODO: remove absolute file name
            std::string temp_file = "/tmp/wt_rlmn_" + util::to_string(util::get_pid()) + "_" + util::to_string(util::get_id());
            std::ofstream wt_out(temp_file.c_str(), std::ios::binary | std::ios::trunc);
            size_type bit_cnt=0;
            wt_out.write((char*)&bit_cnt, sizeof(bit_cnt)); // initial dummy write
            {
                // scope for bl and bf
                bit_vector bl = bit_vector(size, 0);
                m_C  = int_vector<64>(256, 0);

                rac.reset();
                uint8_t last_c = '\0';
                for (size_type i=0, r=0, r_sum=0; r_sum < size;) {
                    for (; i < r+r_sum; ++i) {
                        uint8_t c = rac[i-r_sum];
                        if (last_c != c or i==0) {
                            bl[i] = 1;
                            wt_out.write((char*)&c, sizeof(c));
                            bit_cnt += 8;
                        }
                        ++m_C[c];
                        last_c = c;
                    }
                    r_sum += r;
                    r = rac.load_next_block();
                }

                wt_out.seekp(0, std::ios::beg);
                wt_out.write((char*)&bit_cnt, sizeof(bit_cnt));
                wt_out.close();

                for (size_type i=0, prefix_sum=0, t=0; i<256; ++i) {
                    t = m_C[i];
                    m_C[i] = prefix_sum;
                    prefix_sum += t;
                }

                int_vector<64> lf_map = m_C;
                bit_vector bf = bit_vector(size+1, 0);
                bf[size] = 1; // initialize last element
                rac.reset();
                for (size_type i=0, r=0, r_sum=0; r_sum < size;) {
                    for (; i < r+r_sum; ++i) {
                        uint8_t c = rac[i-r_sum];
                        if (bl[i]) {
                            bf[lf_map[c]] = 1;
                        }
                        ++lf_map[c];
                    }
                    r_sum += r;
                    r = rac.load_next_block();
                }
                {
                    int_vector_file_buffer<8, size_type> temp_bwt_buf(temp_file.c_str());
                    m_wt.construct(temp_bwt_buf, temp_bwt_buf.int_vector_size);

                }
                util::assign(m_bl, bl);
                util::assign(m_bf, bf);
//			m_bl = bit_vector_type(bl);
//			m_bf = bit_vector_type(bf);
            }

            m_bl_rank.init(&m_bl);
            m_bf_rank.init(&m_bf);
            m_bf_select.init(&m_bf);
            m_bl_select.init(&m_bl);
            m_C_bf_rank = int_vector<64>(256,0);
            for (size_type i=0; i<256; ++i) {
                m_C_bf_rank[i] = m_bf_rank(m_C[i]);
            }
            std::remove(temp_file.c_str());
        }

        //! Copy constructor
        wt_rlmn(const wt_rlmn& wt):sigma(wt.sigma) {
            copy(wt);
        }

        //! Assignment operator
        wt_rlmn& operator=(const wt_rlmn& wt) {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Swap operator
        void swap(wt_rlmn& wt) {
            if (this != &wt) {
                std::swap(m_size, wt.m_size);
                m_bl.swap(wt.m_bl);
                m_bf.swap(wt.m_bf);
                m_wt.swap(wt.m_wt);

                m_bl_rank.swap(wt.m_bl_rank);
                m_bl_rank.set_vector(&m_bl);
                wt.m_bl_rank.set(&(wt.m_bl));
                m_bf_rank.swap(wt.m_bf_rank);
                m_bf_rank.set_vector(&m_bf);
                wt.m_bf_rank.set(&(wt.m_bf));

                m_bl_select.swap(wt.m_bl_select);
                m_bl_select.set_vector(&m_bl);
                wt.m_bl_select.set(&(wt.m_bl));
                m_bf_select.swap(wt.m_bf_select);
                m_bf_select.set_vector(&m_bf);
                wt.m_bf_select.set(&(wt.m_bf));

                m_C.swap(wt.m_C);
                m_C_bf_rank.swap(wt.m_C_bf_rank);
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
         *		\f$ \Order{H_0} \f$ on average, where \f$ H_0 \f$ is the zero order entropy of
         *      the sequence
         */
        value_type operator[](size_type i)const {
            return m_wt[m_bl_rank(i+1)-1];
        };

        //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
         *  \param c The symbol to count the occurences in the prefix.
         *	\return The number of occurences of symbol c in the prefix [0..i-1] of the supported vector.
         *  \par Time complexity
         *		\f$ \Order{H_0} \f$ on average, where \f$ H_0 \f$ is the zero order entropy of
         *      the sequence
         */
        size_type rank(size_type i, value_type c)const {
            if (i == 0)
                return 0;
            size_type wt_ex_pos = m_bl_rank(i);
            size_type c_runs = m_wt.rank(wt_ex_pos, c);
            if (c_runs == 0)
                return 0;
            if (m_wt[wt_ex_pos-1] == c) {
                size_type c_run_begin = m_bl_select(wt_ex_pos);
                return m_bf_select(m_C_bf_rank[c] + c_runs) - m_C[c] + i - c_run_begin;
            } else {
                return m_bf_select(m_C_bf_rank[c] + c_runs + 1) - m_C[c];
            }
        };

        //! Calculates how many occurences of symbol wt[i] are in the prefix [0..i-1] of the supported sequence.
        /*!
         *	\param i The index of the symbol.
         *  \param c Reference that will contain the symbol at position i after the execution of the method.
         *  \return The number of occurences of symbol wt[i] in the prefix [0..i-1]
         *	\par Time complexity
         *		\f$ \Order{H_0} \f$
         */
        size_type rank_ith_symbol(size_type i, value_type& c)const {
            if (i == 0) {
                c = m_wt[0];
                return 0;
            }
            size_type wt_ex_pos = m_bl_rank(i+1);
            size_type c_runs = m_wt.rank_ith_symbol(wt_ex_pos-1, c)+1;
            if (c_runs == 0)
                return 0;
            if (m_wt[wt_ex_pos-1] == c) {
                size_type c_run_begin = m_bl_select(wt_ex_pos);
                return m_bf_select(m_C_bf_rank[c] + c_runs) - m_C[c] + i - c_run_begin;
            } else {
                return m_bf_select(m_C_bf_rank[c] + c_runs + 1) - m_C[c];
            }
        }

        //! Calculates the ith occurence of the symbol c in the supported vector.
        /*!
         *  \param i The ith occurence. \f$i\in [1..rank(size(),c)]\f$.
         *  \param c The symbol c.
         *  \par Time complexity
         *		\f$ \Order{H_0} \f$ on average, where \f$ H_0 \f$ is the zero order
         *      entropy of the sequence
         */
        size_type select(size_type i, value_type c)const {
            size_type c_runs = m_bf_rank(m_C[c]+i) - m_C_bf_rank[c];
            size_type offset = m_C[c] + i - 1 - m_bf_select(c_runs + m_C_bf_rank[c]);
            return m_bl_select(m_wt.select(c_runs, c)+1) + offset;
        };

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out)const {
            size_type written_bytes = 0;
            out.write((char*)&m_size, sizeof(m_size));
            written_bytes += sizeof(m_size);
            written_bytes += m_bl.serialize(out);
            written_bytes += m_bf.serialize(out);
            written_bytes += m_wt.serialize(out);
            written_bytes += m_bl_rank.serialize(out);
            written_bytes += m_bf_rank.serialize(out);
            written_bytes += m_bl_select.serialize(out);
            written_bytes += m_bf_select.serialize(out);
            written_bytes += m_C.serialize(out);
            written_bytes += m_C_bf_rank.serialize(out);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            in.read((char*) &m_size, sizeof(m_size));
            m_bl.load(in);
            m_bf.load(in);
            m_wt.load(in);
            m_bl_rank.load(in, &m_bl);
            m_bf_rank.load(in, &m_bf);
            m_bl_select.load(in, &m_bl);
            m_bf_select.load(in, &m_bf);
            m_C.load(in);
            m_C_bf_rank.load(in);
        }

        void print_info()const {
            std::cout<<"# m_wt.size in MB="<<util::get_size_in_bytes(m_wt)/(1024.0*1024.0)<<std::endl;
            std::cout<<"# m_bf.size in MB="<<util::get_size_in_bytes(m_bf)/(1024.0*1024.0)<<std::endl;
            std::cout<<"# m_bl.size in MB="<<util::get_size_in_bytes(m_bl)/(1024.0*1024.0)<<std::endl;
        }

#ifdef MEM_INFO
        void mem_info(std::string label="")const {
            if (label=="")
                label = "wt_rlmn";
            size_type bytes = util::get_size_in_bytes(*this);
            std::cout << "list(label=\""<<label<<"\", size = "<< bytes/(1024.0*1024.0) << "\n,";
            m_bl.mem_info("bl"); std::cout << ",\n";
            m_bf.mem_info("bf"); std::cout << ",\n";
            m_wt.mem_info("wt"); std::cout << ",\n";
            m_bl_rank.mem_info("rank bl"); std::cout << ",\n";
            m_bf_rank.mem_info("rank bf"); std::cout << ",\n";
            m_bl_select.mem_info("select bl"); std::cout << ",\n";
            m_bf_select.mem_info("select bf"); std::cout << ",\n";
            m_C.mem_info("C"); std::cout << ",\n";
            m_C_bf_rank.mem_info("rank C bf"); std::cout << ")\n";
        }
#endif

};


}// end namespace sdsl

#endif // end file 
