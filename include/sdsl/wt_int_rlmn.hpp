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
/*! \file wt_int_rlmn.hpp
    \brief wt_int_rlmn.hpp contains a class for a compressed wavelet tree. Compression is achieved by exploiting runs in the input sequence.
	\author Simon Gog
TODO: merge with wt_rlmn
*/
#ifndef INCLUDED_SDSL_WT_INT_RLMN
#define INCLUDED_SDSL_WT_INT_RLMN

#include "sdsl_concepts.hpp"
#include "int_vector.hpp"
#include "sd_vector.hpp"  // for standard initialisation of template parameters 
#include "util.hpp"
#include "wt_int.hpp"
#include "rrr_vector.hpp"
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
 *   - The rank method: \f$wt.rank(i,c)\f$ returns the number of occurrences of symbol \f$c\f$ in the prefix [0..i-1] in the vector for which the wavelet tree was build for.
 *   - The select method: \f$wt.select(j,c)\f$ returns the index \f$i\in [0..size()-1]\f$ of the jth occurrence of symbol \f$c\f$.
 *
 *	\par Space complexity
 *		 \f$ nH_0 + 2|\Sigma|\log n + 2n + o(n) \f$ bits, where \f$n\f$ is the size of the vector the wavelet tree was build for.
 *
 *  \par Note
 *       This implementation is based on the idea of Veli Mäkinen and Gonzalo Navarro presented in the paper
 *       "Succint Suffix Arrays Based on Run-Length Encoding" (CPM 2005)
 *
 *   @ingroup wt
 *
 *  \tparam BitVector 		Type of the bitvector which is used to represent bf and bl which mark the head of each run in the original sequence.
 *  \tparam RankSupport		Type of the rank support for bitvectors bf and bl.
 *  \tparam SelectSupport   Type of the select support for bitvectors bf and lb.
 *  \tparam WaveletTree		Type of the wavelet tree for the string consisting of the heads of the runs of the original sequence.
 */
template<class BitVector = sd_vector<>,
         class RankSupport = typename BitVector::rank_1_type,
         class SelectSupport = typename BitVector::select_1_type,
         class WaveletTree = wt_int<rrr_vector<63> > >
class wt_int_rlmn
{
    public:
        typedef int_vector<>::size_type		size_type;
        typedef int_vector<>::value_type	value_type;
        typedef BitVector					bit_vector_type;
        typedef RankSupport					rank_support_type;
        typedef SelectSupport           	select_support_type;
        typedef WaveletTree             	wt_type;
        typedef wt_tag						index_category;
        typedef int_alphabet_tag			alphabet_category;
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
        int_vector<>            m_C;     		//
        int_vector<>            m_C_bf_rank;    // stores the number of 1s in m_bf for the prefixes
        // m_bf[0..m_C[0]],m_bf[0..m_C[1]],....,m_bf[0..m_C[sigma-1]];
        // named C_s in the original paper

        void copy(const wt_int_rlmn& wt) {
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
        wt_int_rlmn():m_size(0), sigma(m_wt.sigma) {};

        // Construct the wavelet tree from a random access container
        /*
         *	\param rac Reference to the vector (or unsigned char array) for which the wavelet tree should be build.
         *	\param size Size of the prefix of the vector (or unsigned char array) for which the wavelet tree should be build.
         *	\par Time complexity
         *		\f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         */
        wt_int_rlmn(const unsigned char* rac, size_type size):m_size(size), sigma(m_wt.sigma) {
            // TODO: Delegate this to the file_buffer constructor using a wrapper for the file_buffer
            std::cerr << "ERROR: Constructor of wt_int_rlmn not implemented yet!!!" << std::endl;
            throw std::logic_error("This constructor of wt_int_rlmn is not yet implemented!");
        }

        //! Construct the wavelet tree from a file_buffer
        /*! \param text_buf	A int_vector_file_buffer to the original text.
         *	\param size The length of the prefix of the text, for which the wavelet tree should be build.
         */
        wt_int_rlmn(int_vector_file_buffer<>& text_buf, size_type size):m_size(size), sigma(m_wt.sigma) {
            if (0 == text_buf.int_vector_size or 0 == size)
                return;
            int_vector<> condensed_bwt;
            {
                // scope for bl and bf
                bit_vector bl = bit_vector(size, 0);
                std::map<uint64_t, uint64_t> C;
                text_buf.reset();
                uint64_t last_c = 0;
                size_type runs = 0;
                for (size_type i=0, r=0, r_sum=0; r_sum < size;) {
                    if (r_sum + r > size) {  // read not more than size chars in the next loop
                        r = size-r_sum;
                    }
                    for (; i < r+r_sum; ++i) {
                        uint64_t c = text_buf[i-r_sum];
                        if (last_c != c or i==0) {
                            bl[i] = 1;
                            ++runs;
                        }
                        ++C[c];
                        last_c = c;
                    }
                    r_sum += r;
                    r = text_buf.load_next_block();
                }
                uint64_t max_symbol = (--C.end())->first;
                util::assign(m_C, int_vector<>(max_symbol+1, 0, bit_magic::l1BP(size)+1));
                for (size_type i=0, prefix_sum=0; i<=max_symbol; ++i) {
                    m_C[i] = prefix_sum;
                    prefix_sum += C[i];
                }

                int_vector<> lf_map = m_C;
                bit_vector bf = bit_vector(size+1, 0);
                bf[size] = 1; // initialize last element
                text_buf.reset();
                util::assign(condensed_bwt, int_vector<>(runs, 0, bit_magic::l1BP(max_symbol)+1));
                runs = 0;
                for (size_type i=0, r=0, r_sum=0; r_sum < size;) {
                    if (r_sum + r > size) {  // read not more than size chars in the next loop
                        r = size-r_sum;
                    }
                    for (; i < r+r_sum; ++i) {
                        uint64_t c = text_buf[i-r_sum];
                        if (bl[i]) {
                            bf[lf_map[c]] = 1;
                            condensed_bwt[runs++] = c;
                        }
                        ++lf_map[c];
                    }
                    r_sum += r;
                    r = text_buf.load_next_block();
                }
                {
                    // TODO: remove absolute file name
                    std::string temp_file = "tmp_wt_int_rlmn_" + util::to_string(util::pid()) + "_" + util::to_string(util::id());
                    util::store_to_file(condensed_bwt, temp_file);
                    util::clear(condensed_bwt);
                    int_vector_file_buffer<> temp_bwt_buf(temp_file);
                    util::assign(m_wt, wt_type(temp_bwt_buf, temp_bwt_buf.int_vector_size));
                    std::remove(temp_file.c_str());
                }
                util::assign(m_bl, bl);
                util::assign(m_bf, bf);
            }

            util::init_support(m_bl_rank, &m_bl);
            util::init_support(m_bf_rank, &m_bf);
            util::init_support(m_bf_select, &m_bf);
            util::init_support(m_bl_select, &m_bl);
            util::assign(m_C_bf_rank, int_vector<>(m_C.size(), 0, bit_magic::l1BP(size)+1));
            for (size_type i=0; i<m_C.size(); ++i) {
                m_C_bf_rank[i] = m_bf_rank(m_C[i]);
            }
        }

        //! Copy constructor
        wt_int_rlmn(const wt_int_rlmn& wt):sigma(wt.sigma) {
            copy(wt);
        }

        //! Assignment operator
        wt_int_rlmn& operator=(const wt_int_rlmn& wt) {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Swap operator
        void swap(wt_int_rlmn& wt) {
            if (this != &wt) {
                std::swap(m_size, wt.m_size);
                m_bl.swap(wt.m_bl);
                m_bf.swap(wt.m_bf);
                m_wt.swap(wt.m_wt);

                m_bl_rank.swap(wt.m_bl_rank);
                m_bl_rank.set_vector(&m_bl);
                wt.m_bl_rank.set_vector(&(wt.m_bl));
                m_bf_rank.swap(wt.m_bf_rank);
                m_bf_rank.set_vector(&m_bf);
                wt.m_bf_rank.set_vector(&(wt.m_bf));

                m_bl_select.swap(wt.m_bl_select);
                m_bl_select.set_vector(&m_bl);
                wt.m_bl_select.set_vector(&(wt.m_bl));
                m_bf_select.swap(wt.m_bf_select);
                m_bf_select.set_vector(&m_bf);
                wt.m_bf_select.set_vector(&(wt.m_bf));

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
            assert(i < size());
            return m_wt[m_bl_rank(i+1)-1];
        };

        //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
         *  \param c The symbol to count the occurrences in the prefix.
         *	\return The number of occurrences of symbol c in the prefix [0..i-1] of the supported vector.
         *  \par Time complexity
         *		\f$ \Order{H_0} \f$ on average, where \f$ H_0 \f$ is the zero order entropy of
         *      the sequence
         */
        size_type rank(size_type i, value_type c)const {
            assert(i <= size());
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

        //! Calculates how many occurrences of symbol wt[i] are in the prefix [0..i-1] of the supported sequence.
        /*!
         *	\param i The index of the symbol.
         *  \param c Reference that will contain the symbol at position i after the execution of the method.
         *  \return The number of occurrences of symbol wt[i] in the prefix [0..i-1]
         *	\par Time complexity
         *		\f$ \Order{H_0} \f$
         */
        size_type inverse_select(size_type i, value_type& c)const {
            assert(i < size());
            if (i == 0) {
                c = m_wt[0];
                return 0;
            }
            size_type wt_ex_pos = m_bl_rank(i+1);
            size_type c_runs = m_wt.inverse_select(wt_ex_pos-1, c)+1;
            if (c_runs == 0)
                return 0;
            if (m_wt[wt_ex_pos-1] == c) {
                size_type c_run_begin = m_bl_select(wt_ex_pos);
                return m_bf_select(m_C_bf_rank[c] + c_runs) - m_C[c] + i - c_run_begin;
            } else {
                return m_bf_select(m_C_bf_rank[c] + c_runs + 1) - m_C[c];
            }
        }

        //! Calculates the ith occurrence of the symbol c in the supported vector.
        /*!
         *  \param i The ith occurrence. \f$i\in [1..rank(size(),c)]\f$.
         *  \param c The symbol c.
         *  \par Time complexity
         *		\f$ \Order{H_0} \f$ on average, where \f$ H_0 \f$ is the zero order
         *      entropy of the sequence
         */
        size_type select(size_type i, value_type c)const {
            assert(i > 0);
            assert(i <= rank(size(), c));
            size_type c_runs = m_bf_rank(m_C[c]+i) - m_C_bf_rank[c];
            size_type offset = m_C[c] + i - 1 - m_bf_select(c_runs + m_C_bf_rank[c]);
            return m_bl_select(m_wt.select(c_runs, c)+1) + offset;
        };

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += util::write_member(m_size, out, child, "size");
            written_bytes += m_bl.serialize(out, child, "bl");
            written_bytes += m_bf.serialize(out, child, "bf");
            written_bytes += m_wt.serialize(out, child, "wt");
            written_bytes += m_bl_rank.serialize(out, child, "bl_rank");
            written_bytes += m_bf_rank.serialize(out, child, "bf_rank");
            written_bytes += m_bl_select.serialize(out, child, "bl_select");
            written_bytes += m_bf_select.serialize(out, child, "bf_select");
            written_bytes += m_C.serialize(out, child, "C");
            written_bytes += m_C_bf_rank.serialize(out, child, "C_bf_rank");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            util::read_member(m_size, in);
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
};


}// end namespace sdsl

#endif // end file 
