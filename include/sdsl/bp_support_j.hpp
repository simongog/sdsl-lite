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
/*! \file bp_support_j.hpp
    \brief bp_support_j.hpp contains an implementation of a balanced parentheses support data structure.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_BP_SUPPORT_J
#define INCLUDED_SDSL_BP_SUPPORT_J

#include "int_vector.hpp"
#include "rank_support.hpp"
#include "select_support.hpp"
#include "algorithms.hpp"
#include <stack>
#include <map>
#include <set>
#include <utility>
#include <iostream>
#include <sstream> // for get_info method
#include <fstream>

namespace sdsl
{

//! A class that provides support for bit_vectors that represent a balanced parentheses sequence. Implementation was proposed by Jacobson (1989) and Geary et al. (CPM 2004).
/*! An opening parenthesis is represented by a 1 in the bit_vector and a closing parenthesis by a 0.
 *  This class could be parametrized by a rank_support and select_support.
 *  @ingroup bps
 */
template<class RankSupport = rank_support_v<>, class SelectSupport = select_support_mcl<1> >
class bp_support_j
{
    public:
        typedef bit_vector::size_type size_type;
    private:
        const bit_vector* m_bp; // the supported balanced parentheses sequence as bit_vector
        RankSupport 	m_rank;     // a rank dictionary for calculation of excess
        SelectSupport   m_select;   // additional select dictionary

        bit_vector m_pioneer_bitmap;       // bitmap for pioneer positions
        RankSupport m_rank_pioneer_bitmap; //
        int_vector<> m_pioneer;            //

        bit_vector m_enclose_pioneer_bitmap;       // bitmap for enclose pioneer positions
        RankSupport m_rank_enclose_pioneer_bitmap; //
        int_vector<> m_enclose_pioneer;            //

        uint32_t m_block_size;
        uint32_t m_blocks; // number of blocks
        size_type m_size;


        // TODO: implement this space efficient!!! with the sorted_stack_support!!! 2009-12-03
        //TODO: replace this by a call of algorithm::calculate_pioneer...
        void calculate_pioneers_bitmap() {
            std::map<size_type, size_type> pioneer_matches;
            m_pioneer_bitmap.resize(m_size);      // resize pioneers bitmap
            util::set_zero_bits(m_pioneer_bitmap);  // initialize bitmap with zeros

            std::map<size_type, size_type> pioneer_matches_for_enclose;
            m_enclose_pioneer_bitmap.resize(m_size);	  // resize pioneers bitmap for enclose
            util::set_zero_bits(m_enclose_pioneer_bitmap);  // initialize bitmap with zeros

//			algorithm::calculate_pioneers_bitmap(*m_bp, m_block_size, m_pioneer_bitmap);
//			algorithm::calculate_matches_for_pioneers(*m_bp, m_pioneer_bitmap, m_pioneer);

            std::stack<size_type> opening_parenthesis;

            // calculate positions of findclose and findopen pioneers
            for (size_type block_nr = 0; block_nr < m_blocks; ++block_nr) {
                std::map<size_type, size_type> block_and_position; // for find_open and find_close
                std::map<size_type, size_type> matching_position;  // for find_open and find_close
                std::map<size_type, size_type> block_and_position_for_enclose; // for enclose
                std::map<size_type, size_type> matching_position_for_enclose; // for enclose
                for (size_type i=0, j=block_nr*m_block_size; i < m_block_size and j < m_size; ++i, ++j) {
                    if ((*m_bp)[j]) {//opening parenthesis
                        if (!opening_parenthesis.empty()) {
                            size_type position = opening_parenthesis.top();
                            size_type blockpos = position/m_block_size;
//								if( block_and_position_for_enclose.find(blockpos) == block_and_position_for_enclose.end() ){ // smallest j is pioneer
                            block_and_position_for_enclose[blockpos] = position;
                            matching_position_for_enclose[blockpos] = j;
//								}
                        }
                        opening_parenthesis.push(j);
                    } else { // closing parenthesis
                        size_type position = opening_parenthesis.top();
                        size_type blockpos = position/m_block_size;
                        opening_parenthesis.pop();
                        block_and_position[blockpos] = position;
                        matching_position[blockpos]  = j; // greatest j is pioneer
                    }
                }
                for (std::map<size_type, size_type>::const_iterator it = block_and_position.begin(),
                     end = block_and_position.end(),
                     mit = matching_position.begin(); it != end and it->first != block_nr; ++it, ++mit) {
                    // opening and closing pioneers are symmetric
                    m_pioneer_bitmap[it->second] = 1;
                    pioneer_matches[it->second] = mit->second;
                    m_pioneer_bitmap[mit->second] = 1;
                    pioneer_matches[mit->second] = it->second;
                }
                for (std::map<size_type, size_type>::const_iterator it = block_and_position_for_enclose.begin(),
                     end = block_and_position_for_enclose.end(),
                     mit = matching_position_for_enclose.begin(); it != end and it->first != block_nr; ++it, ++mit) {
                    m_enclose_pioneer_bitmap[mit->second] = 1;
                    pioneer_matches_for_enclose[mit->second] = it->second ;
                }

            }
            // assert that the sequence is balanced
            assert(opening_parenthesis.empty());
            // store matching positions of pioneers
            m_pioneer.set_int_width(bit_magic::l1BP(m_size)+1);
            m_pioneer.resize(pioneer_matches.size());
            size_type cnt=0;
            for (std::map<size_type, size_type>::const_iterator mit = pioneer_matches.begin(); mit!= pioneer_matches.end();  ++mit) {
                m_pioneer[cnt++] = mit->second;
            }

            // initialize the rank dictionary for the pioneer bitmap
            util::init_support(m_rank_pioneer_bitmap, &m_pioneer_bitmap);

            // store matching positions of enclose pioneers
            m_enclose_pioneer.set_int_width(bit_magic::l1BP(m_size)+1);
            m_enclose_pioneer.resize(pioneer_matches_for_enclose.size());
            cnt = 0;
            for (std::map<size_type, size_type>::const_iterator mit = pioneer_matches_for_enclose.begin(); mit != pioneer_matches_for_enclose.end(); ++mit) {
                m_enclose_pioneer[cnt++] = mit->second;
            }
            // initialize the rank dictionary for the enclose pioneer bitmap
            util::init_support(m_rank_enclose_pioneer_bitmap, &m_enclose_pioneer_bitmap);

            /*				if(m_size<120 and m_size>0){
            				std::cerr<<"bp"<<std::endl;
            				for(size_type i=0; i<m_size; ++i){ std::cerr<<(*m_bp)[i]; }
            				std::cerr<<std::endl<<"pioneer bitmap"<<std::endl;
            				for(size_type i=0; i<m_size; ++i){ std::cerr<<m_pioneer_bitmap[i]; }
            				std::cerr<<std::endl;
            				for(std::map<size_type, size_type>::const_iterator mit = pioneer_matches.begin(); mit!= pioneer_matches.end();  ++mit){
            					std::cerr<<"_"<<mit->first<<" "<<mit->second<<std::endl;
            				}
            				for(size_type i=0;i<m_pioneer.size();++i){
            					std::cerr<<"__"<<i<<" "<<m_pioneer[i]<<std::endl;
            				}
            			}
            */
        }

        void copy(const bp_support_j& bp_support) {
            m_bp = bp_support.m_bp;
            m_block_size = bp_support.m_block_size;
            m_size = bp_support.m_size;
            m_blocks = bp_support.m_blocks;

            m_rank = bp_support.m_rank;
            m_rank.set_vector(m_bp);

            m_select = bp_support.m_select;
            m_select.set_vector(m_bp);

            m_pioneer_bitmap = bp_support.m_pioneer_bitmap;
            m_rank_pioneer_bitmap = bp_support.m_rank_pioneer_bitmap;
            m_rank_pioneer_bitmap.set_vector(&m_pioneer_bitmap);
            m_pioneer = bp_support.m_pioneer;

            m_enclose_pioneer_bitmap = bp_support.m_enclose_pioneer_bitmap;
            m_rank_enclose_pioneer_bitmap = bp_support.m_rank_enclose_pioneer_bitmap;
            m_rank_enclose_pioneer_bitmap.set_vector(&m_enclose_pioneer_bitmap);
            m_enclose_pioneer = bp_support.m_enclose_pioneer;
        }

    public:

        //! Constructor
        bp_support_j(const bit_vector* bp = NULL, uint32_t ignore=0/*, uint32_t used_block_size = 64*/):m_bp(bp), m_block_size(64/*block_size==0?used_block_size:block_size*/), m_size(bp==NULL?0:bp->size()) {
//				assert(m_block_size > 0 and m_block_size <= 8192 and m_block_size%64 == 0);
            m_blocks = (m_size+m_block_size-1)/m_block_size;
            util::init_support(m_rank, m_bp);
            util::init_support(m_select, m_bp);
            calculate_pioneers_bitmap();
        }

        //! Copy constructor
        bp_support_j(const bp_support_j& bp_support) {
            copy(bp_support);
        }

        //! Assignment operator
        bp_support_j& operator=(const bp_support_j& bp_support) {
            if (this != &bp_support) {
                copy(bp_support);
            }
            return *this;
        }

        void set_vector(const bit_vector* bp) {
            m_bp = bp;
            m_rank.set_vector(bp);
            m_select.set_vector(bp);
        }

        /*! Calculates the excess value at index i.
         * \param i The index of which the excess value should be calculated.
         */
        size_type excess(size_type i)const {
            return (m_rank(i+1)<<1)-i-1;
        }

        /*! Returns the number of opening parentheses up to and including index i.
         * \pre{ \f$ 0\leq i < size() \f$ }
         */
        size_type rank(size_type i)const {
            return m_rank(i+1);
        }

        /*! Returns the index of the i-th opening parenthesis.
         * \param i Number of the parenthesis to select.
         * \pre{ \f$1\leq i < rank(size())\f$ }
         * \post{ \f$ 0\leq select(i) < size() \f$ }
         */
        size_type select(size_type i)const {
            return m_select(i);
        }

        /*! Calculate the index of the matching closing parenthesis to the parenthesis at index i.
         * \param i Index of an parenthesis. 0 <= i < size().
         * \return * i, if the parenthesis at index i is closing,
         *         * the position j of the matching closing parenthesis, if a matching parenthesis exists,
         *         * size() if no matching closing parenthesis exists.
         */
        size_type find_close(size_type i)const {
#ifdef SDSL_DEBUG_BP
            if (i >= m_size) {
                throw std::out_of_range("OUT_OF_RANGE: bp_support_j::find_close");
            }
#endif
            if (!(*m_bp)[i]) {// if there is a closing parenthesis at index i return i
                return i;
            }
            if (m_pioneer_bitmap[i]) {
                return m_pioneer[m_rank_pioneer_bitmap.rank(i)];
            } else {
                size_type ip1 = i+1; // move one position to the left
                uint64_t byte_prefix_sums_x_2;
                const uint64_t* data = m_bp->data() + (ip1>>6);
                uint8_t offset = ip1&0x3F;
                uint64_t w = (~(*data))>>offset; // add $offset$ opening parenthesis at the beginning of the word

                uint8_t pos = bit_magic::first_excess_position(w, 1, byte_prefix_sums_x_2);
                if (pos != 64) {// found position of matching parenthesis
                    return ip1+pos;
                } else { // the excess value at the end of w is lower or equal zero
                    assert((i%64) != 0);
                    offset = (offset + 63)&0x3F;
                    // search the pioneer bitmap for the closesd preceding pioneer
                    pos = bit_magic::l1BP(*(m_pioneer_bitmap.data() + ((i-1)>>6)) & *(m_bp->data() + ((i-1)>>6))  & bit_magic::Li1Mask[offset]);
                    size_type pioneer_index = ((i-1)&0xFFFFFFFFFFFFFFC0ULL)+pos;
                    assert(m_pioneer_bitmap[pioneer_index]  == 1);
                    size_type match_index = m_pioneer[m_rank_pioneer_bitmap.rank(pioneer_index)];
                    assert((match_index&0xFFFFFFFFFFFFFFC0ULL) > i);
                    uint8_t  excess_difference = excess((match_index&0xFFFFFFFFFFFFFFC0ULL)-1)-(excess(i)-1);
                    assert(excess_difference <=63);
                    return (match_index&0xFFFFFFFFFFFFFFC0ULL)
                           + bit_magic::first_excess_position(~*(m_bp->data()+(match_index>>6)) , excess_difference, byte_prefix_sums_x_2);
                }
            }
        }

        //! Calculate the matching opening parenthesis to the closing parenthesis at position i
        /*! \param i Index of a closing parenthesis.
          * \return * i, if the parenthesis at index i is closing,
          *         * the position j of the matching opening parenthesis, if a matching parenthesis exists,
          *         * size() if no matching closing parenthesis exists.
          */
        size_type find_open(size_type i)const {
#ifdef SDSL_DEBUG_BP
            if (i >= m_size) {
                throw std::out_of_range("OUT_OF_RANGE: bp_support_j::find_open");
            }
#endif
            if ((*m_bp)[i]) {// if there is a opening parenthesis at index i return i
                return i;
            }
            if (m_pioneer_bitmap[i]) {
                return m_pioneer[m_rank_pioneer_bitmap.rank(i)];
            } else {
                size_type im1 = i-1; // move one position to the right
                const uint64_t* data = m_bp->data() + (im1>>6);
                uint8_t close_parenthesis_index = (i&0x3F) + ((i==0)<<6);
                uint8_t pos = bit_magic::find_open(*data, close_parenthesis_index);
                if (pos!=64) { // found position of the matching parenthesis
                    return (im1&0xFFFFFFFFFFFFFFC0ULL)+pos;
                } else {
                    assert((i%64)!=63);
                    // search the pioneer bitmap for the closest succeeding pioneer
                    pos = bit_magic::r1BP(*(m_pioneer_bitmap.data() + ((i+1)>>6)) & ~*(m_bp->data() + ((i+1)>>6))  & bit_magic::Li0Mask[i&0x3F]);
                    size_type pioneer_index = ((i+1)&0xFFFFFFFFFFFFFFC0ULL)+pos;
                    assert(m_pioneer_bitmap[pioneer_index] == 1);
                    size_type match_index = m_pioneer[m_rank_pioneer_bitmap.rank(pioneer_index)];
                    assert(match_index < i);
                    int8_t excess_difference = excess(i);
                    if (match_index >= 64) {
                        excess_difference -= excess((match_index&0xFFFFFFFFFFFFFFC0ULL)-1);
                    }
                    assert(excess_difference >=-64 and excess_difference <= 64);
                    uint64_t dummy;
                    if (excess_difference >= 0) {
                        return (match_index&0xFFFFFFFFFFFFFFC0ULL)
                               + bit_magic::last_excess_position(*(m_bp->data()+(match_index>>6)), excess_difference, dummy) + 1;
                    } else {
                        return (match_index&0xFFFFFFFFFFFFFFC0ULL)
                               + bit_magic::last_excess_position(~*(m_bp->data()+(match_index>>6)), -excess_difference, dummy) + 1;
                    }
                }
            }
        }

        //! Calculate the index of the opening parenthesis corresponding to the closest matching parenthesis pair enclosing i.
        /*! \param i Index of an opening parenthesis.
         *  \return The index of the opening parenthesis corresponding to the closest matching parenthesis pair enclosing i,
         *          or size() if no such pair exists.
         */
        size_type enclose(size_type i)const {
#ifdef SDSL_DEBUG_BP
            if (i >= m_size) {
                throw std::out_of_range("OUT_OF_RANGE: bp_support_j::enclose.");
            }
#endif
            if (!(*m_bp)[i]) { // if there is closing parenthesis at position i
                throw std::logic_error("LOGIC_ERROR: bp_support_j::enclose. A opening parenthesis is expected as argument.");
            }
            if (m_enclose_pioneer_bitmap[i]) {
                return m_enclose_pioneer[m_rank_enclose_pioneer_bitmap.rank(i)];
            } else {
                if (i==0)
                    return size();
                size_type im1 = i-1; // move one position to the right
                const uint64_t* data = m_bp->data() + (im1>>6);
                uint8_t open_parenthesis_index = (i&0x3F) + ((i==0)<<6);
                uint8_t pos = bit_magic::find_open(*data, open_parenthesis_index);
                if (pos!=64) {
                    return (im1&0xFFFFFFFFFFFFFFC0ULL)+pos;
                } else {
                    assert((i%64)!= 63);
                    // search the pioneer bitmap for the closest succeeding pioneer
                    uint64_t w = *(m_enclose_pioneer_bitmap.data() + ((i+1)>>6)) & bit_magic::Li0Mask[i&0x3F];
                    if (w) {
                        pos = bit_magic::r1BP(*(m_enclose_pioneer_bitmap.data() + ((i+1)>>6)) & bit_magic::Li0Mask[i&0x3F]);
                        size_type pioneer_index = ((i+1)&0xFFFFFFFFFFFFFFC0ULL)+pos;
                        assert(m_enclose_pioneer_bitmap[pioneer_index] == 1);
                        size_type match_index = m_enclose_pioneer[m_rank_enclose_pioneer_bitmap.rank(pioneer_index)];
                        assert(match_index < i);
                        int8_t excess_difference = excess(i)-2;
                        if (match_index >= 64) {
                            excess_difference -= excess((match_index&0xFFFFFFFFFFFFFFC0ULL)-1);
                        }
                        assert(excess_difference >=-64 and excess_difference <= 64);
                        uint64_t dummy;
                        if (excess_difference >= 0) {
                            pos = bit_magic::last_excess_position(*(m_bp->data()+(match_index>>6)), excess_difference, dummy);
                            if (pos==64)
                                return match_index&0xFFFFFFFFFFFFFFC0ULL;
                            else
                                return (match_index&0xFFFFFFFFFFFFFFC0ULL) + pos + 1;
//								return (match_index&0xFFFFFFFFFFFFFFC0ULL)
//									  + bit_magic::last_excess_position( *(m_bp.data()+(match_index>>6)), excess_difference, dummy) + 1;
                        } else {
                            pos = bit_magic::last_excess_position(~*(m_bp->data()+(match_index>>6)), -excess_difference, dummy);
                            if (pos==64)
                                return match_index&0xFFFFFFFFFFFFFFC0ULL;
                            else
                                return (match_index&0xFFFFFFFFFFFFFFC0ULL) + pos + 1;
//								return (match_index&0xFFFFFFFFFFFFFFC0ULL)
//									  + bit_magic::last_excess_position( ~*(m_bp.data()+(match_index>>6)), -excess_difference, dummy) + 1;
                        }
                    } else { // there exists no pioneer => there exists no enclosing parentheses pair
                        return m_size;
                    }
                }
            }
        }

        //! The range restricted enclose operation.
        /*! \param i Index of an opening parenthesis.
        	\param j Index of an opening parenthesis/ \f$ i<j \wedge findclose(i) < j \f$.
        	\return The smallest index, say k, of an opening parenthesis such that findclose(i) < k < j and
        	findclose(j) < findclose(k). If such a k does not exists, restricted_enclose(i,j) returns size().
        */
        size_type rr_enclose(size_type i, size_type j)const {
            assert(j > i and j < m_size);
            size_type mi = find_close(i); // matching parenthesis to i
            assert(mi > i and mi < j);
            assert(find_close(j) > j);
            size_type k = enclose(j);
            if (k == m_size or k < i) // there exists no opening parenthesis at position mi<k<j.
                return m_size;
            size_type kk;
            do {
                kk = k;
                k = enclose(k);
            } while (k != m_size and k > mi);
            return kk;
        }

        size_type rr_enclose_naive(size_type i, size_type j)const {
            return rr_enclose(i, j);
        }

        //! The double enclose operation
        /*! \param i Index of an opening parenthesis.
         *  \param j Index of an opening parenthesis \f$ i<j \wedge findclose(i) < j \f$.
         *  \return The maximal opening parenthesis, say k, such that \f$ k<j \wedge k>findclose(j) \f$.
         *          If such a k does not exists, double_enclose(i,j) returns size().
         */
        size_type double_enclose(size_type i, size_type j)const {
            assert(j > i);
            assert((*m_bp)[i]==1 and (*m_bp)[j]==1);
            size_type k = rr_enclose(i, j);
            if (k == size())
                return enclose(j);
            else
                return enclose(k);
        }

        //! Return the number of zeros which procede position i in the balanced parentheses sequence.
        /*! \param i Index of an parenthesis.
         */
        size_type preceding_closing_parentheses(size_type i)const {
            assert(i < m_size);
            if (!i) return 0;
            size_type ones = m_rank(i);
            if (ones) { // ones > 0
                assert(m_select(ones) < i);
                return i - m_select(ones) - 1;
            } else {
                return i;
            }
            /*
            			size_type im1 = i-1;
            			const uint64_t *data = m_bp->data()+(im1>>6);
            			uint8_t offset = (im1&0x3F)+1;
            			uint64_t w = *data << (64-offset);//& bit_magic::Li1Mask[offset];
            			if( offset != 64 ){
            				w |= *(--data)>>(offset);
            			}
            			size_type result = 0;
            */
        }
        /*
        		size_type restricted_enclose2(size_type i, size_type j)const{
        			assert( j > i );
        			size_type mi = find_close(i); // matching parenthesis to i
        			assert(  find_close(i) > i and mi < j  );
        			//
        //				m_rank_pioneer_bitmap.rank(mi)

        			assert( find_close(j) > j );
        			size_type k = enclose(j);
        			if( k == m_size or k < i )// there exists no opening parenthesis at position mi<k<j.
        				return m_size;
        			size_type kk;
        			do{
        				kk = k;
        				k = enclose(k);
        			}while( k != m_size and k > mi );
        			return kk;
        		}
        */

        /*! The size of the supported balanced parentheses sequence.
         * \return the size of the supported balanced parentheses sequence.
         */
        size_type size() const {
            return m_size;
        }

        //! Serializes the bp_support_j to a stream.
        /*!
         * \param out The outstream to which the data structure is written.
         * \return The number of bytes written to out.
         */
        size_type serialize(std::ostream& out) const {
            size_type written_bytes = 0;
            written_bytes += m_rank.serialize(out);
            written_bytes += m_select.serialize(out);

            written_bytes += m_pioneer_bitmap.serialize(out);
            written_bytes += m_rank_pioneer_bitmap.serialize(out);
            written_bytes += m_pioneer.serialize(out);

            written_bytes += m_enclose_pioneer_bitmap.serialize(out);
            written_bytes += m_rank_enclose_pioneer_bitmap.serialize(out);
            written_bytes += m_enclose_pioneer.serialize(out);

            return written_bytes;
        }

        //! Load the bp_support_j for a bit_vector v.
        /*!
         * \param in The instream from which the data strucutre is read.
         * \param bp Bit vector representing a balanced parentheses sequence that is supported by this data structure.
         */
        void load(std::istream& in, const bit_vector* bp) {
            m_bp = bp;
            if (m_bp == NULL)
                return;
            m_size = m_bp->size();
            m_block_size = 64;
            m_blocks = (m_size+m_block_size-1)/m_block_size;

            m_rank.load(in, m_bp);
            m_select.load(in, m_bp);

            m_pioneer_bitmap.load(in);
            m_rank_pioneer_bitmap.load(in, &m_pioneer_bitmap);
            m_pioneer.load(in);

            m_enclose_pioneer_bitmap.load(in);
            m_rank_enclose_pioneer_bitmap.load(in, &m_enclose_pioneer_bitmap);
            m_enclose_pioneer.load(in);
        }

        std::string get_info()const {
            std::stringstream ss;
            if (m_bp == NULL) {
                ss<<"ERROR: bp_support_j is unitialized!"<<std::endl;
                return ss.str();
            }
            ss<<"number of parentheses: "<<m_bp->size()<<std::endl;
            ss<<"number of pioneers: "<< m_pioneer.size()<<std::endl;
            std::ofstream out("/dev/null");
            ss<<"size of parentheses sequence: "<< m_bp->serialize(out) << std::endl;
            ss<<"size of pioneer bitmap: "<< m_pioneer_bitmap.serialize(out) << std::endl;
            ss<<"size of pioneer data  : "<< m_pioneer.serialize(out) << std::endl;
            ss<<"size of rank for pioneers: "<< m_rank_pioneer_bitmap.serialize(out) << std::endl;

            ss<<"size of enclose  pioneer bitmap: "<< m_enclose_pioneer_bitmap.serialize(out) << std::endl;
            ss<<"size of pioneer data  : "<< m_pioneer.serialize(out) << std::endl;
            ss<<"size of rank for pioneers: "<< m_rank_enclose_pioneer_bitmap.serialize(out) << std::endl;

            return ss.str();
        }
};
}// end namespace sdsl




#endif
