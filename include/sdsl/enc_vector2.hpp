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
/*! \file enc_vector2.hpp
   \brief enc_vector2.hpp contains the sdsl::enc_vector2 class.
   \author Simon Gog
*/
#ifndef SDSL_ENC_VECTORII
#define SDSL_ENC_VECTORII

#include "int_vector.hpp"
#include "coder.hpp"
#include "iterators.hpp"


//! Namespace for the succinct data structure library.
namespace sdsl
{

template<uint8_t t_width>
struct enc_vector2_trait {
    typedef int_vector<0> int_vector_type;
};

template<>
struct enc_vector2_trait<32> {
    typedef int_vector<32> int_vector_type;
};

template<>
struct enc_vector2_trait<64> {
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
class enc_vector2
{
    private:
        static_assert(t_dens > 1 , "enc_vector2: sample density must be larger than `1`");
    public:
        typedef uint64_t                                 value_type;
        typedef random_access_const_iterator<enc_vector2> iterator;
        typedef iterator                                 const_iterator;
        typedef const value_type                         reference;
        typedef const value_type                         const_reference;
        typedef const value_type*                        const_pointer;
        typedef ptrdiff_t                                difference_type;
        typedef int_vector<>::size_type                  size_type;
        typedef t_coder                                  coder;
        typedef typename enc_vector2_trait<t_width>::int_vector_type int_vector_type;
        typedef iv_tag                                   index_category;
        static  constexpr uint32_t                       sample_dens    = t_dens;
        typedef enc_vector2                              enc_vec_type;

        int_vector<0>     m_z;                       // storage for encoded deltas
    private:
        int_vector_type            m_samples;        // samples
        sd_vector<>                m_pointers;
        sd_vector<>::select_1_type m_pointers_sel;
        size_type                  m_size = 0;       // number of vector elements

        void clear()
        {
            m_z.resize(0);
            m_size = 0;
            m_samples.resize(0);
            m_pointers = sd_vector<>();
        }

    public:
        enc_vector2() = default;
        enc_vector2(const enc_vector2&) = default;
        enc_vector2(enc_vector2&&) = default;
        enc_vector2& operator=(const enc_vector2&) = default;
        enc_vector2& operator=(enc_vector2&&) = default;

        //! Constructor for a Container of unsigned integers.
        /*! \param c A container of unsigned integers.
          */
        template<class Container>
        enc_vector2(const Container& c);

        //! Constructor for an int_vector_buffer of unsigned integers.
        /*
            \param v_buf A int_vector_buf.
        */
        template<uint8_t int_width>
        enc_vector2(int_vector_buffer<int_width>& v_buf);

        //! Default Destructor
        ~enc_vector2() { }

        //! The number of elements in the enc_vector2.
        size_type size()const
        {
            return m_size;
        }

        //! Return the largest size that this container can ever have.
        static size_type max_size()
        {
            return int_vector<>::max_size()/2;
        }

        //!    Returns if the enc_vector2 is empty.
        bool empty() const
        {
            return 0==m_size;
        }

        //! Swap method for enc_vector2
        void swap(enc_vector2& v);

        //! Iterator that points to the first element of the enc_vector2.
        const const_iterator begin()const
        {
            return const_iterator(this, 0);
        }

        //! Iterator that points to the position after the last element of the enc_vector2.
        const const_iterator end()const
        {
            return const_iterator(this, this->m_size);
        }

        //! operator[]
        /*! \param i Index. \f$ i \in [0..size()-1]\f$.
         */
        value_type operator[](size_type i)const;

        //! Serialize the enc_vector2 to a stream.
        /*! \param out Out stream to write the data structure.
            \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const;

        //! Load the enc_vector2 from a stream.
        void load(std::istream& in);

        //! Returns the i-th sample of enc_vector2
        /*! \param i The index of the sample. 0 <= i < size()/get_sample_dens()
         *  \return The value of the i-th sample.
         */
        value_type sample(const size_type i) const;

        uint32_t get_sample_dens() const
        {
            return t_dens;
        }

        /*!
         * \param i The index of the sample for which all values till the next sample should be decoded. 0 <= i < size()/get_sample_dens()
         * \param it A pointer to a uint64_t vector, whereto the values should be written
         */
        void get_inter_sampled_values(const size_type i, uint64_t* it)const
        {
            // TODO: this will not work for blocks with m_pointers_sel(i+1)+t_dens==m_pointers_sel(i+2)
            *(it++) = 0;
            if (i*t_dens + t_dens - 1 < size()) {
                if (i+1 < m_samples.size() and m_samples[i] + t_dens == m_samples[i+1]) {
                    if (m_pointers_sel(i+1) != m_pointers_sel(i+2)) {
                        throw std::logic_error("Should not be here");
                    }
                    uint64_t x = 1;
                    while (x < t_dens) {
                        *(it++) = x;
                        ++x;
                    }
//                    throw std::logic_error("Should not be here");
                } else {
                    t_coder::template decode<true, true>(m_z.data(), m_pointers_sel(i+1), t_dens - 1, it);
                }
            } else {
                assert(i*t_dens < size());
                t_coder::template decode<true, true>(m_z.data(), m_pointers_sel(i+1), size()-i*t_dens - 1, it);
            }
        };
};

template<class t_coder, uint32_t t_dens, uint8_t t_width>
inline typename enc_vector2<t_coder, t_dens,t_width>::value_type enc_vector2<t_coder, t_dens,t_width>::operator[](const size_type i)const
{
    assert(i+1 != 0);
    assert(i < m_size);
    size_type idx = i/get_sample_dens();
    if (idx+1 < m_samples.size() and m_samples[idx]+t_dens == m_samples[idx+1]) {
        return m_samples[idx] + i-t_dens*idx;
    }
    return m_samples[idx] + t_coder::decode_prefix_sum(m_z.data(), m_pointers_sel(idx+1), i-t_dens*idx);
}

template<class t_coder, uint32_t t_dens, uint8_t t_width>
inline typename enc_vector2<t_coder, t_dens,t_width>::value_type enc_vector2<t_coder, t_dens,t_width>::sample(const size_type i)const
{
    assert(i*get_sample_dens()+1 != 0);
    assert(i*get_sample_dens() < m_size);
    return m_samples[i];
}

template<class t_coder, uint32_t t_dens, uint8_t t_width>
void enc_vector2<t_coder, t_dens,t_width>::swap(enc_vector2<t_coder, t_dens,t_width>& v)
{
    if (this != &v) { // if v and _this_ are not the same object
        m_z.swap(v.m_z);
        m_samples.swap(v.m_samples);
        m_pointers.swap(v.m_pointers);
        util::swap_support(m_pointers_sel, v.m_pointers_sel, &m_pointers, &v.m_pointers);
        std::swap(m_size, v.m_size);
    }
}

template<class t_coder, uint32_t t_dens, uint8_t t_width>
template<uint8_t int_width>
enc_vector2<t_coder, t_dens,t_width>::enc_vector2(int_vector_buffer<int_width>& v_buf)
{
    // clear bit_vectors
    clear();
    size_type n = v_buf.size();
    if (n == 0)  // if c is empty there is nothing to do...
        return;
    value_type     v1=0, v2=0, max_sample_value=0;
    size_type samples=0, z_size=0;
    const size_type sd = get_sample_dens();
    size_type tmp_z = 0;
    bool uniform = true;
//  (1) Calculate maximal value of samples and of deltas
    for (size_type i=0, no_sample = 0; i < n; ++i, --no_sample) {
        v2 = v_buf[i];
        if (!no_sample) { // is sample
            uniform &= (v2==v1+1);
            if (!uniform) {
                z_size += tmp_z;
            }
            uniform = true;
            tmp_z = 0;
            no_sample = sd;
            if (max_sample_value < v2) max_sample_value = v2;
            ++samples;
        } else {
            uniform &= (v2==v1+1);
            tmp_z += t_coder::encoding_length(v2-v1);
        }
        v1 = v2;
    }
    z_size += tmp_z;

//  (2) Write sample values and deltas
//  (a) Initialize array for sample values and pointers
    m_samples = int_vector<>(samples+1, 0, bits::hi(max_sample_value)+1);

    sd_vector_builder builder(z_size, samples);

//  (b) Initilize bit_vector for encoded data
    m_z = int_vector<>(z_size, 0, 1);
    uint64_t* z_data = t_coder::raw_data(m_z);
    uint8_t offset = 0;

//  (c) Write sample values and deltas
    z_size = 0;
    tmp_z = 0;
    uniform = true;
    std::vector<uint64_t> delta;
    for (size_type i=0, j=0, no_sample = 0; i < n; ++i, --no_sample) {
        v2 = v_buf[i];
        if (!no_sample) { // is sample
            uniform &= (v2==v1+1);
            if (!uniform) {
                for (size_t k=0; k<delta.size(); ++k) {
                    z_size += t_coder::encoding_length(delta[k]);
                    t_coder::encode(delta[k], z_data, offset); // write encoded data
                }
            }
            delta.clear();
            uniform = true;
            no_sample = sd;
            m_samples[j++] = v2;    // write samples
            builder.set(z_size);
        } else {
            uniform &= (v2==v1+1);
            delta.push_back(v2-v1);
        }
        v1 = v2;
    }
    for (size_t k=0; k<delta.size(); ++k) {
        t_coder::encode(delta[k], z_data, offset); // write encoded data
        z_size += t_coder::encoding_length(delta[k]);
    }
    m_size = n;
    m_pointers = sd_vector<>(builder);
    m_pointers_sel.set_vector(&m_pointers);
}

template<class t_coder, uint32_t t_dens, uint8_t t_width>
enc_vector2<>::size_type enc_vector2<t_coder, t_dens,t_width>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += write_member(m_size, out, child, "size");
    written_bytes += m_z.serialize(out, child, "encoded deltas");
    written_bytes += m_samples.serialize(out, child, "samples");
    written_bytes += m_pointers.serialize(out, child, "pointers");
    written_bytes += m_pointers_sel.serialize(out, child, "pointers_sel");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<class t_coder, uint32_t t_dens, uint8_t t_width>
void enc_vector2<t_coder, t_dens,t_width>::load(std::istream& in)
{
    read_member(m_size, in);
    m_z.load(in);
    m_samples.load(in);
    m_pointers.load(in);
    m_pointers_sel.load(in);
    m_pointers_sel.set_vector(&m_pointers);
}

} // end namespace sdsl
#endif
