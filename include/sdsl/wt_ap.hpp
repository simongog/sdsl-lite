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
 *        \f$\order{n} (H_0 + 1)\f$ bits, where \f$n\f$ is the size of the vector the wavelet tree was build for.1
 *
 *  \tparam t_wt     Type of the wavelet trees used for internal representation.
 *
 *    \par References
 *    [1] J. Barbay, F. Claude, T. Gagie, G. Navarro and Y. Nekrich:
 *        ,,Efficient Fully-Compressed Sequence Representations''
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

        // retrieves a character's class and offset - if the character exists in the text
        bool try_get_char_class_offset(value_type c, value_type& cl, value_type& offset)const {
            assert(i <= size());
            if (c >= m_char2class.size()) { // c is greater than any symbol in text
                return false;
            }
            auto offset_class = m_char2class.inverse_select(c);
            cl = offset_class.second;
            if (cl == m_classes) { // c never occurs in text
                return false;
            }
            offset = offset_class.first;
            return true;
        }

    public:

        const size_type&       sigma = m_sigma;

        //! Default constructor
        wt_ap() {}

        //! Semi-external constructor
        /*! \param buf         File buffer of the int_vector for which the wt_ap should be build.
         *  \param size        Size of the prefix of v, which should be indexed.
         */
        template<uint8_t int_width>
        wt_ap(int_vector_buffer<int_width>& buf, size_type size) : m_size(size) {
            size_type n = buf.size();  // set n
            if (n < m_size) {
                throw std::logic_error("n="+util::to_string(n)+" < "+util::to_string(m_size)+"=m_size");
                return;
            }

            // calculate effective sigma and character frequencies
            value_type max_symbol = 0;
            int_vector<int_width> rac(m_size, 0, buf.width());
            std::vector<std::pair<size_type, value_type>> char_freq;
            value_type pseudo_entries = 0;
            for (size_type i=0; i < m_size; ++i) {
                auto element = rac[i] = buf[i];
                while (element >= max_symbol)
                {
                    char_freq.push_back(std::make_pair(0, max_symbol));
                    max_symbol++;
                    pseudo_entries++;
                }
                if (char_freq[element].first == 0) {
                    pseudo_entries--;
                }
                char_freq[element].first++;
            }
            std::sort(char_freq.rbegin(), char_freq.rend());
            m_sigma = max_symbol - pseudo_entries;

            m_singleton_classes = std::min(max_symbol, (value_type)bits::hi(m_sigma)); 
            //                                                           OR max_symbol?
            m_classes = bits::hi(m_sigma - m_singleton_classes + 1) + m_singleton_classes;            

            std::vector<std::pair<size_type, int_vector<>>> m_offset_buffer;

            // assign character classes
            int_vector<> m_char2class_buffer(max_symbol, m_classes, bits::hi(m_classes+1)+1);
            for (value_type i=0; i < m_singleton_classes; ++i) {
                m_char2class_buffer[char_freq[i].second] = i;
            }
            value_type current_symbol = m_singleton_classes;
            value_type class_size = 1;
            for (value_type i=m_singleton_classes; i < m_classes; ++i) {
                class_size <<= 1;
                size_type class_frequency = 0;
                value_type offset=0;
                for (; offset < class_size && current_symbol < m_sigma; ++offset, ++current_symbol) {
                    m_char2class_buffer[char_freq[current_symbol].second] = i;
                    class_frequency += char_freq[current_symbol].first;
                }
                m_offset_buffer.push_back(std::make_pair(0, int_vector<>(class_frequency, 0, bits::hi(offset)+1)));
            }

            construct_im(m_char2class, m_char2class_buffer);

            // calculate text-order classes and offsets
            for (size_type i=0; i < m_size; ++i) {
                value_type ch = rac[i];
                value_type cl = rac[i] = m_char2class_buffer[ch];
                if (cl >= m_singleton_classes) {
                    value_type offset = m_char2class.rank(ch, cl);
                    cl -= m_singleton_classes;
                    m_offset_buffer[cl].second[m_offset_buffer[cl].first++] = offset;
                }
            }
            
            construct_im(m_class, rac);
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
         *  \par Worst case time complexity
         *       \f$ \Order{\log \log |\Sigma|} \f$
         *  \par Average case time complexity
         *       \f$ \Order{\log H_0} \f$
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
         *  \par Worst case time complexity
         *       \f$ \Order{\log \log |\Sigma|} \f$
         *  \par Average case time complexity
         *       \f$ \Order{\log H_0} \f$
         *  \par Precondition
         *       \f$ i \leq size() \f$
         */
        size_type rank(size_type i, value_type c)const {
            assert(i <= size());
            value_type offset;
            value_type cl;
            if (!try_get_char_class_offset(c, cl, offset)) {
                return 0;
            }
            size_type count = m_class.rank(i, cl);
            return cl < m_singleton_classes
                ? count
                : m_offset[cl-m_singleton_classes].rank(count, offset);
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

            auto textoffset_class = m_class.inverse_select(i);
            auto textoffset = textoffset_class.first;
            auto cl = textoffset_class.second;
            if (cl < m_singleton_classes) {
                return std::make_pair(textoffset, m_char2class.select(1, cl));
            }
            auto class_result = m_offset[cl-m_singleton_classes].inverse_select(textoffset);
            return std::make_pair(class_result.first, m_char2class.select(class_result.second+1, cl));
        }

        //! Calculates the i-th occurrence of the symbol c in the supported vector.
        /*!
         *  \param i The i-th occurrence.
         *  \param c The symbol c.
         *  \par Worst case time complexity
         *       \f$ \Order{\log \log |\Sigma|} \f$
         *  \par Average case time complexity
         *       \f$ \Order{\log H_0} \f$
         *  \par Precondition
         *       \f$ 1 \leq i \leq rank(size(), c) \f$
         */
        size_type select(size_type i, value_type c)const {
            assert(1 <= i and i <= rank(size(), c));
            value_type offset;
            value_type cl;
            if (!try_get_char_class_offset(c, cl, offset)) {
                return m_size;
            }
            size_type text_offset = cl < m_singleton_classes
                ? i
                : 1 + m_offset[cl-m_singleton_classes].select(i, offset);
            return m_class.select(text_offset, cl);
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
            value_type offset_size = m_classes - m_singleton_classes;
            m_offset.resize(offset_size);
            for (int i=0; i<offset_size; ++i) {
                m_offset[i].load(in);
            }
        }
};

}// end namespace sdsl
#endif
