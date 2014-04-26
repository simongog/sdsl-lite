/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

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
/*! \file enc_vector.hpp
   \brief enc_vector.hpp contains the sdsl::enc_vector class.
   \author Simon Gog
*/
#ifndef SDSL_ENC_VECTOR
#define SDSL_ENC_VECTOR

#include "int_vector.hpp"
#include "coder.hpp"
#include "iterators.hpp"


//! Namespace for the succinct data structure library.
namespace sdsl
{

template<uint8_t t_width>
struct enc_vector_trait {
    typedef int_vector<0> int_vector_type;
};

template<>
struct enc_vector_trait<32> {
    typedef int_vector<32> int_vector_type;
};

template<>
struct enc_vector_trait<64> {
    typedef int_vector<64> int_vector_type;
};

//! A generic immutable space-saving vector class for unsigned integers.
/*! A vector v is stored more space-efficiently by self-delimiting coding
 *  the deltas v[i+1]-v[i] (v[-1]:=0). Space of the structure and random
 *  access time to it can be controlled by a sampling parameter t_dens.
 *
 *  \tparam t_coder  Self-delimiting coder.
 *  \tparam t_dens   Every t_dens-th element of v is sampled.
 *  \tparam t_width  Width of the int_vector used to store the samples and pointers.
 *  This class is a parameter of csa_sada.
 * @ingroup int_vector
 */
template<class t_coder=coder::elias_delta,
         uint32_t t_dens = 128, uint8_t t_width=0>
class enc_vector
{
    private:
        static_assert(t_dens > 1 , "enc_vector: sample density must be larger than `1`");
    public:
        typedef uint64_t                                 value_type;
        typedef random_access_const_iterator<enc_vector> iterator;
        typedef iterator                                 const_iterator;
        typedef const value_type                         reference;
        typedef const value_type                         const_reference;
        typedef const value_type*                        const_pointer;
        typedef ptrdiff_t                                difference_type;
        typedef int_vector<>::size_type                  size_type;
        typedef t_coder                                  coder;
        typedef typename enc_vector_trait<t_width>::int_vector_type int_vector_type;
        typedef iv_tag                                   index_category;
        static  const uint32_t                           sample_dens    = t_dens;

        int_vector<0>     m_z;                       // storage for encoded deltas
    private:
        int_vector_type   m_sample_vals_and_pointer; // samples and pointers
        size_type         m_size = 0;                // number of vector elements

        void clear() {
            m_z.resize(0);
            m_size = 0;
            m_sample_vals_and_pointer.resize(0);
        }

    public:
        enc_vector() = default;
        enc_vector(const enc_vector&) = default;
        enc_vector(enc_vector&&) = default;
        enc_vector& operator=(const enc_vector&) = default;
        enc_vector& operator=(enc_vector&&) = default;

        //! Constructor for a Container of unsigned integers.
        /*! \param c A container of unsigned integers.
          */
        template<class Container>
        enc_vector(const Container& c);

        //! Constructor for an int_vector_buffer of unsigned integers.
        /*
            \param v_buf A int_vector_buf.
        */
        template<uint8_t int_width>
        enc_vector(int_vector_buffer<int_width>& v_buf);

        //! Default Destructor
        ~enc_vector() { }

        //! The number of elements in the enc_vector.
        size_type size()const {
            return m_size;
        }

        //! Return the largest size that this container can ever have.
        static size_type max_size() {
            return int_vector<>::max_size()/2;
        }

        //!    Returns if the enc_vector is empty.
        bool empty() const {
            return 0==m_size;
        }

        //! Swap method for enc_vector
        void swap(enc_vector& v);

        //! Iterator that points to the first element of the enc_vector.
        const const_iterator begin()const {
            return const_iterator(this, 0);
        }

        //! Iterator that points to the position after the last element of the enc_vector.
        const const_iterator end()const {
            return const_iterator(this, this->m_size);
        }

        //! operator[]
        /*! \param i Index. \f$ i \in [0..size()-1]\f$.
         */
        value_type operator[](size_type i)const;

        //! Serialize the enc_vector to a stream.
        /*! \param out Out stream to write the data structure.
            \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const;

        //! Load the enc_vector from a stream.
        void load(std::istream& in);

        //! Returns the i-th sample of enc_vector
        /*! \param i The index of the sample. 0 <= i < size()/get_sample_dens()
         *  \return The value of the i-th sample.
         */
        value_type sample(const size_type i) const;

        uint32_t get_sample_dens() const {
            return t_dens;
        }

        /*!
         * \param i The index of the sample for which all values till the next sample should be decoded. 0 <= i < size()/get_sample_dens()
         * \param it A pointer to a uint64_t vector, whereto the values should be written
         */
        void get_inter_sampled_values(const size_type i, uint64_t* it)const {
            *(it++) = 0;
            if (i*t_dens + t_dens - 1 < size()) {
                t_coder::template decode<true, true>(m_z.data(), m_sample_vals_and_pointer[(i<<1)+1], t_dens - 1, it);
            } else {
                assert(i*t_dens < size());
                t_coder::template decode<true, true>(m_z.data(), m_sample_vals_and_pointer[(i<<1)+1], size()-i*t_dens - 1, it);
            }
        };
};

template<class t_coder, uint32_t t_dens, uint8_t t_width>
inline typename enc_vector<t_coder, t_dens,t_width>::value_type enc_vector<t_coder, t_dens,t_width>::operator[](const size_type i)const
{
    assert(i+1 != 0);
    assert(i < m_size);
    size_type idx = i/get_sample_dens();
    return m_sample_vals_and_pointer[idx<<1] + t_coder::decode_prefix_sum(m_z.data(), m_sample_vals_and_pointer[(idx<<1)+1], i-t_dens*idx);
}

template<class t_coder, uint32_t t_dens, uint8_t t_width>
inline typename enc_vector<t_coder, t_dens,t_width>::value_type enc_vector<t_coder, t_dens,t_width>::sample(const size_type i)const
{
    assert(i*get_sample_dens()+1 != 0);
    assert(i*get_sample_dens() < m_size);
    return m_sample_vals_and_pointer[i<<1];
}

template<class t_coder, uint32_t t_dens, uint8_t t_width>
void enc_vector<t_coder, t_dens,t_width>::swap(enc_vector<t_coder, t_dens,t_width>& v)
{
    if (this != &v) { // if v and _this_ are not the same object
        m_z.swap(v.m_z);
        m_sample_vals_and_pointer.swap(v.m_sample_vals_and_pointer);
        std::swap(m_size, v.m_size);
    }
}

template<class t_coder, uint32_t t_dens, uint8_t t_width>
template<class Container>
enc_vector<t_coder, t_dens,t_width>::enc_vector(const Container& c)
{
    // clear bit_vectors
    clear();

    if (c.empty())  // if c is empty there is nothing to do...
        return;
    typename Container::const_iterator    it             = c.begin(), end = c.end();
    typename Container::value_type         v1            = *it, v2, max_sample_value=0, x;
    size_type samples=0;
    size_type z_size = 0;
//  (1) Calculate maximal value of samples and of deltas
    for (size_type i=0, no_sample=0; it != end; ++it,++i, --no_sample) {
        v2 = *it;
        if (!no_sample) { // add a sample
            no_sample = get_sample_dens();
            if (max_sample_value < v2) max_sample_value = v2;
            ++samples;
        } else {
            z_size += t_coder::encoding_length(v2-v1);
        }
        v1=v2;
    }
//    (2) Write sample values and deltas
    {
        if (max_sample_value > z_size+1)
            m_sample_vals_and_pointer.width(bits::hi(max_sample_value) + 1);
        else
            m_sample_vals_and_pointer.width(bits::hi(z_size+1) + 1);
        m_sample_vals_and_pointer.resize(2*samples+2); // add 2 for last entry
        util::set_to_value(m_sample_vals_and_pointer, 0);
        typename int_vector_type::iterator sv_it = m_sample_vals_and_pointer.begin();
        z_size = 0;
        size_type no_sample=0;
        for (it = c.begin(); it != end; ++it, --no_sample) {
            v2 = *it;
            if (!no_sample) { // add a sample
                no_sample = get_sample_dens();
                *sv_it = v2; ++sv_it;
                *sv_it = z_size; ++sv_it;
            } else {
                x = v2-v1;
                z_size += t_coder::encoding_length(x);
            }
            v1=v2;
        }
        *sv_it = 0; ++sv_it;        // initialize
        *sv_it = z_size+1; ++sv_it; // last entry

        m_z = int_vector<>(z_size, 0, 1);
        uint64_t* z_data = t_coder::raw_data(m_z);
        uint8_t offset = 0;
        no_sample = 0;
        for (it = c.begin(); it != end; ++it, --no_sample) {
            v2 = *it;
            if (!no_sample) { // add a sample
                no_sample = get_sample_dens();
            } else {
                t_coder::encode(v2-v1, z_data, offset);
            }
            v1=v2;
        }
    }
    m_size = c.size();
}

template<class t_coder, uint32_t t_dens, uint8_t t_width>
template<uint8_t int_width>
enc_vector<t_coder, t_dens,t_width>::enc_vector(int_vector_buffer<int_width>& v_buf)
{
    // clear bit_vectors
    clear();
    size_type n = v_buf.size();
    if (n == 0)  // if c is empty there is nothing to do...
        return;
    value_type     v1=0, v2=0, max_sample_value=0;
    size_type samples=0, z_size=0;
    const size_type sd = get_sample_dens();
//  (1) Calculate maximal value of samples and of deltas
    for (size_type i=0, no_sample = 0; i < n; ++i, --no_sample) {
        v2 = v_buf[i];
        if (!no_sample) { // is sample
            no_sample = sd;
            if (max_sample_value < v2) max_sample_value = v2;
            ++samples;
        } else {
            z_size += t_coder::encoding_length(v2-v1);
        }
        v1 = v2;
    }

//    (2) Write sample values and deltas
//    (a) Initialize array for sample values and pointers
    if (max_sample_value > z_size+1)
        m_sample_vals_and_pointer.width(bits::hi(max_sample_value) + 1);
    else
        m_sample_vals_and_pointer.width(bits::hi(z_size+1) + 1);
    m_sample_vals_and_pointer.resize(2*samples+2); // add 2 for last entry
    util::set_to_value(m_sample_vals_and_pointer, 0);

//    (b) Initilize bit_vector for encoded data
    m_z = int_vector<>(z_size, 0, 1);
    uint64_t* z_data = t_coder::raw_data(m_z);
    uint8_t offset = 0;

//    (c) Write sample values and deltas
    z_size = 0;
    for (size_type i=0, j=0, no_sample = 0; i < n; ++i, --no_sample) {
        v2 = v_buf[i];
        if (!no_sample) { // is sample
            no_sample = sd;
            m_sample_vals_and_pointer[j++] = v2;    // write samples
            m_sample_vals_and_pointer[j++] = z_size;// write pointers
        } else {
            z_size += t_coder::encoding_length(v2-v1);
            t_coder::encode(v2-v1, z_data, offset);   // write encoded values
        }
        v1 = v2;
    }
    m_size = n;
}

template<class t_coder, uint32_t t_dens, uint8_t t_width>
enc_vector<>::size_type enc_vector<t_coder, t_dens,t_width>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += write_member(m_size, out, child, "size");
    written_bytes += m_z.serialize(out, child, "encoded deltas");
    written_bytes += m_sample_vals_and_pointer.serialize(out, child, "samples_and_pointers");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<class t_coder, uint32_t t_dens, uint8_t t_width>
void enc_vector<t_coder, t_dens,t_width>::load(std::istream& in)
{
    read_member(m_size, in);
    m_z.load(in);
    m_sample_vals_and_pointer.load(in);
}

} // end namespace sdsl
#endif
