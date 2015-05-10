/* sdsl - succinct data structures library
    Copyright (C) 2014 Simon Gog

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
/*! \file wt_ap.hpp
    \brief wt_ap.hpp contains a space-efficient class to support select, 
            rank and access on inputs with potentially large alphabets.
    \author Johannes Bader, Simon Gog
*/
#ifndef INCLUDED_SDSL_WT_AP
#define INCLUDED_SDSL_WT_AP

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/vectors.hpp>

//! Namespace for the succinct data structure library.
namespace sdsl
{
    
//! A wavelet tree class for integer sequences.
/*!
 *    \par Space complexity
 *        \f$\Order{n\log|\Sigma|}\f$ bits, where \f$n\f$ is the size of the vector the wavelet tree was build for.
 *
 *  \tparam t_bitvector   Type of the bitvector used for representing the wavelet tree.
 *  \tparam t_rank        Type of the support structure for rank on pattern `1`.
 *  \tparam t_select      Type of the support structure for select on pattern `1`.
 *  \tparam t_select_zero Type of the support structure for select on pattern `0`.
 *
 *   @ingroup wt
 */
template<class t_wt = wt_int<>>
class wt_ap
{
    public:

        typedef int_vector<>::size_type              size_type;
        typedef int_vector<>::value_type             value_type;
        typedef random_access_const_iterator<wt_ap>  const_iterator;
        typedef const_iterator                       iterator;
        typedef t_wt                                 wt_type;
        typedef wt_tag                               index_category;
        typedef int_alphabet_tag                     alphabet_category;
        enum 	{lex_ordered=0};

    protected:

        size_type              m_size  = 0;
        value_type             m_sigma = 0;     //<- \f$ |\Sigma| \f$
        value_type             m_singleton_classes = 0;
        value_type             m_classes = 0;
        wt_type                m_char2class;
        wt_type                m_class;
        std::vector<wt_type>   m_offset;

        void copy(const wt_ap& wt) {
            m_size              = wt.m_size;
            m_sigma             = wt.m_sigma;
            m_singleton_classes = wt.m_singleton_classes;
            m_classes           = wt.m_classes;
            m_char2class        = wt.m_char2class;
            m_class             = wt.m_class;
            m_offset            = wt.m_offset;
        }

    private:

    public:

        const size_type&       sigma = m_sigma;         //!< Effective alphabet size of the wavelet tree.

        //! Default constructor
        wt_ap() {}

        //! Semi-external constructor
        /*! \param buf         File buffer of the int_vector for which the wt_ap should be build.
         *  \param size        Size of the prefix of v, which should be indexed.
         *  \param max_level   Maximal level of the wavelet tree. If set to 0, determined automatically.
         *    \par Time complexity
         *        \f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         *        I.e. we need \Order{n\log n} if rac is a permutation of 0..n-1.
         *    \par Space complexity
         *        \f$ n\log|\Sigma| + O(1)\f$ bits, where \f$n=size\f$.
         */
        template<uint8_t int_width>
        wt_ap(int_vector_buffer<int_width>& buf, size_type size) : m_size(size) {
            size_type n = buf.size();  // set n
            if (n < m_size) {
                throw std::logic_error("n="+util::to_string(n)+" < "+util::to_string(m_size)+"=m_size");
                return;
            }

            // calculate effective sigma and character frequencies
            m_sigma = 0;
            int_vector<int_width> rac(m_size, 0, buf.width());
            std::vector<std::pair<size_type, value_type>> char_freq;
            value_type pseudo_entries = 0;
            for (size_type i=0; i < m_size; ++i) {
                auto element = rac[i] = buf[i];
                while (element >= m_sigma)
                {
                    char_freq.push_back(std::make_pair(0, m_sigma));
                    m_sigma++;
                    pseudo_entries++;
                }
                if (char_freq[element].first == 0) {
                    pseudo_entries--;
                }
                char_freq[element].first++;
            }
            std::sort(char_freq.rbegin(), char_freq.rend());
            value_type m_compact_sigma = m_sigma - pseudo_entries;

            m_singleton_classes = 1; // + bits::hi(m_sigma); // OR m_compact_sigma?
            m_classes = bits::hi(m_compact_sigma - m_singleton_classes + 1) + m_singleton_classes;
            
            std::vector<std::pair<size_type, int_vector<>>> m_offset_buffer;
            
            // assign character classes
            int_vector<> m_char2class_buffer(m_sigma, m_classes, bits::hi(m_classes+1)+1);
            for (value_type i=0; i < m_singleton_classes; ++i) {
                m_char2class_buffer[char_freq[i].second] = i;
            }
            value_type current_symbol = m_singleton_classes;
            value_type class_size = 1;
            for (value_type i=m_singleton_classes; i < m_classes; ++i) {
                class_size <<= 1;
                size_type class_frequency = 0;
                value_type offset=0;
                for (; offset < class_size && current_symbol < m_compact_sigma; ++offset, ++current_symbol) {
                    m_char2class_buffer[char_freq[current_symbol].second] = i;
                    class_frequency += char_freq[current_symbol].first;
                }
                m_offset_buffer.push_back(std::make_pair(0, int_vector<>(class_frequency, 0, bits::hi(offset)+1)));
            }
            
            construct_im(m_char2class, m_char2class_buffer);
            
            // calculate text-order classes and offsets
            int_vector<> m_class_buffer(m_size, 0, bits::hi(m_classes)+1);
            for (size_type i=0; i < m_size; ++i) {
                value_type ch = rac[i];
                value_type cl = m_class_buffer[i] = m_char2class_buffer[ch];
                if (cl >= m_singleton_classes) {
                    value_type offset = m_char2class.rank(ch, cl);
                    cl -= m_singleton_classes;
                    m_offset_buffer[cl].second[m_offset_buffer[cl].first++] = offset;
                }
            }
            
            construct_im(m_class, m_class_buffer);
            m_offset.resize(m_classes-m_singleton_classes);
            for (value_type i=0; i < m_classes-m_singleton_classes; ++i) {
                construct_im(m_offset[i], m_offset_buffer[i].second);
            }
        }

        //! Copy constructor
        wt_ap(const wt_ap& wt) {
            copy(wt);
        }

        //! Copy constructor
        wt_ap(wt_ap&& wt) {
            *this = std::move(wt);
        }

        //! Assignment operator
        wt_ap& operator=(const wt_ap& wt) {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Assignment move operator
        wt_ap& operator=(wt_ap&& wt) {
            if (this != &wt) {
                m_size              = wt.m_size;
                m_sigma             = wt.m_sigma;
                m_singleton_classes = wt.m_singleton_classes;
                m_classes           = wt.m_classes;
                m_char2class        = std::move(wt.m_char2class);
                m_class             = std::move(wt.m_class);
                m_offset            = std::move(wt.m_offset);
            }
            return *this;
        }

        //! Swap operator
        void swap(wt_ap& wt) {
            if (this != &wt) {
                std::swap(m_size, wt.m_size);
                std::swap(m_sigma, wt.m_sigma);
                std::swap(m_singleton_classes, wt.m_singleton_classes);
                std::swap(m_classes, wt.m_classes);
                m_char2class.swap(wt.m_char2class);
                m_class.swap(wt.m_class);
                std::swap(m_offset,  wt.m_offset);
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
        /*! \param i The index of the symbol in the original vector.
         *  \returns The i-th symbol of the original vector.
         *  \par Precondition
         *       \f$ i < size() \f$
         */
        value_type operator[](size_type i)const {
            assert(i < size());
            auto textoffset_class = m_class.inverse_select(i);
            auto cl = textoffset_class.second;
            value_type offset = cl < m_singleton_classes 
                ? 0 
                : m_offset[cl-m_singleton_classes][textoffset_class.first];
            return m_char2class.select(offset+1, cl);
        };

        //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
         *  \param c The symbol to count the occurrences in the prefix.
         *    \returns The number of occurrences of symbol c in the prefix [0..i-1] of the supported vector.
         *  \par Time complexity
         *       \f$ \Order{\log |\Sigma|} \f$
         *  \par Precondition
         *       \f$ i \leq size() \f$
         */
        size_type rank(size_type i, value_type c)const {
            assert(i <= size());
            if (c >= m_sigma) { // c is greater than any symbol in text
                return 0;
            }
            auto offset_class = m_char2class.inverse_select(c);
            auto cl = offset_class.second;
            if (cl == m_classes) { // c never occurs in text
                return 0;
            }
            size_type count = m_class.rank(i, cl);
            return cl < m_singleton_classes
                ? count
                : m_offset[cl-m_singleton_classes].rank(count, offset_class.first);
        };



        //! Calculates how many occurrences of symbol wt[i] are in the prefix [0..i-1] of the original sequence.
        /*!
         *  \param i The index of the symbol.
         *  \return  Pair (rank(wt[i],i),wt[i])
         *  \par Precondition
         *       \f$ i < size() \f$
         */
        std::pair<size_type, value_type>
        inverse_select(size_type i)const {
            assert(i < size());

            value_type val = operator [](i);
            size_type cnt = rank(i, val);
            return std::make_pair(val,cnt);
        }

        //! Calculates the i-th occurrence of the symbol c in the supported vector.
        /*!
         *  \param i The i-th occurrence.
         *  \param c The symbol c.
         *  \par Time complexity
         *       \f$ \Order{\log |\Sigma|} \f$
         *  \par Precondition
         *       \f$ 1 \leq i \leq rank(size(), c) \f$
         */
        size_type select(size_type i, value_type c)const {
            assert(1 <= i and i <= rank(size(), c));
            if (c >= m_sigma) { // c is greater than any symbol in text
                return m_size;
            }
            auto offset_class = m_char2class.inverse_select(c);
            auto cl = offset_class.second;
            if (cl == m_classes) { // c never occurs in text
                return m_size;
            }
            size_type offset = cl < m_singleton_classes
                ? i
                : 1 + m_offset[cl-m_singleton_classes].select(i, offset_class.first);
            return m_class.select(offset, cl);
        };

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += write_member(m_sigma, out, child, "sigma");
            written_bytes += write_member(m_singleton_classes, out, child, "singleton_classes");
            written_bytes += write_member(m_classes, out, child, "classes");
            written_bytes += m_char2class.serialize(out, child, "char2class");
            written_bytes += m_class.serialize(out, child, "class");
            written_bytes += write_member(m_size, out, child, "offset_size");
            for (int i=0; i<m_offset.size(); ++i) {
                written_bytes += m_offset[i].serialize(out, child, "offset");
            }
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            read_member(m_size, in);
            read_member(m_sigma, in);
            read_member(m_singleton_classes, in);
            read_member(m_classes, in);
            m_char2class.load(in);
            m_class.load(in);
            typename std::vector<wt_type>::size_type offset_size;
            read_member(offset_size, in);
            m_offset.resize(offset_size);
            for (int i=0; i<offset_size; ++i) {
                m_offset[i].load(in);
            }
        }
};

}// end namespace sdsl
#endif
