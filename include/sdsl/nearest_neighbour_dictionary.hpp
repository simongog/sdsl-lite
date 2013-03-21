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
/*! \file nearest_neighbour_dictionary.hpp
    \brief nearest_neighbour_dictionary.hpp contains a class which supports rank/select for sparse populated sdsl::bit_vectors.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_NEAREST_NEIGHBOUR_DICTIONARY
#define INCLUDED_SDSL_NEAREST_NEIGHBOUR_DICTIONARY

#include "int_vector.hpp"
#include "rank_support.hpp"
#include "util.hpp"
#include <stdexcept>
#include <string>

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! Nearest neighbour dictionary for sparse uniform sets (described in Geary et al., A Simple Optimal Representation for Balanced Parentheses, CPM 2004).
/*!
 * Template parameter sample_dens corresponds to parameter t in the paper.
 * The data structure the following methods:
 *  - rank
 *  - select
 *  - prev
 *  - next
 * @ingroup rank_support_group
 * @ingroup select_support_group
 *
*/
// TODO: implement an iterator for the ones in the nearest neighbour dictionary!!! used in the construction of the balanced parentheses support
template<uint8_t sample_dens>
class nearest_neighbour_dictionary
{
    public:
        typedef bit_vector::size_type size_type;
    private:
        int_vector<> m_abs_samples; // absolute samples array corresponds to array \f$ A_1 \f$ in the paper
        int_vector<> m_differences; // vector for the differences in between the samples; corresponds to array \f$ A_2 \f$ in the paper
        size_type    m_ones; // corresponds to N in the paper
        size_type    m_size; // corresponds to M in the paper
        bit_vector   m_contains_abs_sample; // vector which stores for every block of length sample_dens of the original bit_vector if an absolute sample lies in this block.
        // Corresponds to array \f$ A_3 \f$ in the paper.
        rank_support_v<>  m_rank_contains_abs_sample; // rank support for m_contains_abs_sample. Corresponds to array \f$ A_4 \f$ in the paper.
        // NOTE: A faster version should store the absolute samples and the differences interleaved

        void copy(const nearest_neighbour_dictionary& nnd) {
            // copy all members of the data structure
            m_abs_samples 	= 		nnd.m_abs_samples;
            m_differences 	= 		nnd.m_differences;
            m_ones 			= 		nnd.m_ones;
            m_size			= 		nnd.m_size;
            m_contains_abs_sample =	nnd.m_contains_abs_sample;
            m_rank_contains_abs_sample = nnd.m_rank_contains_abs_sample;
            m_rank_contains_abs_sample.set_vector(&m_contains_abs_sample);
        }

    public:

        //! Default constructor
        nearest_neighbour_dictionary():m_ones(0),m_size(0) { }

        //! Constructor
        /*! \param v The supported bit_vector.
         */
        nearest_neighbour_dictionary(const bit_vector& v):m_ones(0), m_size(0) {
            if (sample_dens==0) { // first logical error check
                throw std::logic_error(util::demangle(typeid(this).name())+": sample_dens should not be equal 0!");
            }
            size_type max_distance_between_two_ones = 0;
            size_type ones = 0; // counter for the ones in v

            // get maximal distance between to ones in the bit vector
            // speed this up by broadword computing
            for (size_type i=0, last_one_pos_plus_1=0; i < v.size(); ++i) {
                if (v[i]) {
                    if (i+1-last_one_pos_plus_1 > max_distance_between_two_ones)
                        max_distance_between_two_ones = i+1-last_one_pos_plus_1;
                    last_one_pos_plus_1 = i+1;
                    ++ones;

                }
            }
            m_ones = ones;
            m_size = v.size();
//			std::cerr<<ones<<std::endl;
            // initialize absolute samples m_abs_samples[0]=0
            m_abs_samples = int_vector<>(m_ones/sample_dens + 1, 0,  bits::l1BP(v.size())+1);
            // initialize different values
            m_differences = int_vector<>(m_ones - m_ones/sample_dens, 0, bits::l1BP(max_distance_between_two_ones)+1);
            // initialize m_contains_abs_sample
            m_contains_abs_sample = bit_vector((v.size()+sample_dens-1)/sample_dens, 0);
            ones = 0;
            for (size_type i=0, last_one_pos=0; i < v.size(); ++i) {
                if (v[i]) {
                    ++ones;
                    if ((ones % sample_dens) == 0) {  // insert absolute samples
                        m_abs_samples[ones/sample_dens] = i;
                        m_contains_abs_sample[i/sample_dens] = 1;
                    } else {
                        m_differences[ones - ones/sample_dens - 1] = i - last_one_pos;
                    }
                    last_one_pos = i;
                }
            }
            util::init_support(m_rank_contains_abs_sample, &m_contains_abs_sample);
        }

        //! Copy constructor
        nearest_neighbour_dictionary(const nearest_neighbour_dictionary& nnd) {
            // copy all members of the data structure
            copy(nnd);
        }

        //! Destructor
        ~nearest_neighbour_dictionary() {}

        nearest_neighbour_dictionary& operator=(const nearest_neighbour_dictionary& nnd) {
            if (this != &nnd) {
                copy(nnd);
            }
            return *this;
        }

        void swap(nearest_neighbour_dictionary& nnd) {
            // copy all members of the data structure
            m_abs_samples.swap(nnd.m_abs_samples);
            m_differences.swap(nnd.m_differences);
            std::swap(m_ones, nnd.m_ones);
            std::swap(m_size, nnd.m_size);
            m_contains_abs_sample.swap(nnd.m_contains_abs_sample);
            util::swap_support(m_rank_contains_abs_sample, nnd.m_rank_contains_abs_sample,
                               &m_contains_abs_sample, &(nnd.m_contains_abs_sample));
        }

        //! Answers rank queries for the supported bit_vector
        /*! \param idx Argument for the length of the prefix v[0..idx-1].
         *  \return Number of 1-bits in the prefix [0..idx-1] of the supported bit_vector.
         *  \par Time complexity \f$ \Order{1} \f$
         */
        size_type rank(size_type idx)const {
            assert(idx <= m_size);
            size_type r = m_rank_contains_abs_sample.rank(idx/sample_dens); //
            size_type result = r*sample_dens;
            size_type i = m_abs_samples[r];
            while (++result <= m_ones) {
                if ((result % sample_dens) == 0) {
                    i = m_abs_samples[result/sample_dens];
                } else {
                    i = i+m_differences[result - result/sample_dens-1];
                }
                if (i >= idx)
                    return result-1;
            }
            return result-1;
        };

        //! Answers select queries for the supported bit_vector
        /*! \param i Select the \f$i\f$th 1 in the supported bit_vector. \f$i\in [1..ones()]\f$
         *  \return The position of the \f$i\f$th 1 in the supported bit_vector.
         *  \par Time complexity \f$ \Order{1} \f$
         */
        size_type select(size_type i)const {
            assert(i > 0 and i <= m_ones);
            size_type j = i/sample_dens;
            size_type result = m_abs_samples[j];
            j = j*sample_dens - j;
            for (size_type end = j + (i%sample_dens); j < end; ++j) {
                result += m_differences[j];
            }
            return result;
        }

        //! Answers "previous occurence of one" queries for the supported bit_vector.
        /*! \param i Position \f$ i \in [0..size()-1] \f$.
         *  \return The maximal position \f$j \leq i\f$ where the supported bit_vector v equals 1.
         *  \pre rank(i+1)>0
         *  \par Time complexity \f$ \Order{1} \f$
         */
        size_type prev(size_type i)const {
            size_type r = rank(i+1);
            assert(r>0);
            return select(r);
        }
        /*! Answers "next occurence of one" queries for the supported bit_vector.
         * \param i Position \f$ i \in [0..size()-1] \f$.
         * \return The minimal position \f$ j \geq i \f$ where the supported bit_vector v equals 1.
         * \pre rank(i) < ones()
         *  \par Time complexity \f$ \Order{1} \f$
         */
        size_type next(size_type i)const {
            size_type r = rank(i);
            assert(r < m_ones);
            return select(r+1);
        }

        size_type size()const {
            return m_size;
        }

        size_type ones()const {
            return m_ones;
        }

        //! Serializes the nearest_neighbour_dictionary.
        /*! \param out Out-Stream to serialize the data to.
        */
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            size_type written_bytes = 0;
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            written_bytes += m_abs_samples.serialize(out, child, "absolute_samples");
            written_bytes += m_differences.serialize(out, child, "differences");
            written_bytes += util::write_member(m_ones, out, child, "ones");
            written_bytes += util::write_member(m_size,out,  child, "size");
            written_bytes += m_contains_abs_sample.serialize(out, child, "contains_abs_sample");
            written_bytes += m_rank_contains_abs_sample.serialize(out, child, "rank_contains_abs_sample");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the nearest_neighbour_dictionary.
        /*! \param in In-Stream to load the rank_support data from.
         */
        void load(std::istream& in) {
            m_abs_samples.load(in);
            m_differences.load(in);
            util::read_member(m_ones, in);
            util::read_member(m_size, in);
            m_contains_abs_sample.load(in);
            m_rank_contains_abs_sample.load(in, &m_contains_abs_sample);
        }
};


}// end namespace sdsl


#endif // end file 
