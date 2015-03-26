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
/*! \file vlc_vector.hpp
   \brief vlc_vector.hpp contains a vector which stores the values with variable length codes.
   \author Simon Gog
*/
#ifndef SDSL_VLC_VECTOR
#define SDSL_VLC_VECTOR

#include "int_vector.hpp"
#include "coder_elias_delta.hpp"
#include "iterators.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{

template<uint8_t t_width>
struct vlc_vector_trait {
    typedef int_vector<0> int_vector_type;
};

template<>
struct vlc_vector_trait<32> {
    typedef int_vector<32> int_vector_type;
};

//! A generic immutable space-saving vector class for unsigned integers.
/*! The values of a vlc_vector are immutable after the constructor call. The class
 *   could be parametrized with a self-delimiting code t_coder and the sample density.
 *  \tparam t_coder Type of self-delimiting coder.
 *  \tparam t_dens  Sampling density of pointers into the stream of self-delimiting coded numbers.
 *  \tparam t_width Width of the underlying int_vector for the pointers.
 */
template<class    t_coder = coder::elias_delta,
         uint32_t t_dens  = 128,
         uint8_t  t_width = 0>
class vlc_vector
{
    private:
        static_assert(t_dens > 1 , "vlc_vector: Sampling density must be larger than 1");
    public:
        typedef uint64_t                                 value_type;
        typedef random_access_const_iterator<vlc_vector> iterator;
        typedef iterator                                 const_iterator;
        typedef const value_type                         reference;
        typedef const value_type                         const_reference;
        typedef const value_type*                        const_pointer;
        typedef ptrdiff_t                                difference_type;
        typedef int_vector<>::size_type                  size_type;
        typedef t_coder                                  coder;
        typedef iv_tag                                   index_category;
        typedef typename
        vlc_vector_trait<t_width>::int_vector_type      int_vector_type;

        static  const uint32_t                           sample_dens = t_dens;
        bit_vector      m_z;     // compressed bit stream
    private:
        int_vector_type m_sample_pointer;
        size_type       m_size        = 0;    // number of elements
        uint32_t        m_sample_dens = t_dens;

        void copy(const vlc_vector& v);

        void clear()
        {
            m_z.resize(0);
            m_size = 0;
            m_sample_pointer.resize(0);
        }

    public:
        vlc_vector() = default;
        vlc_vector(const vlc_vector&) = default;
        vlc_vector(vlc_vector&&) = default;
        vlc_vector& operator=(const vlc_vector&) = default;
        vlc_vector& operator=(vlc_vector&&) = default;

        //! Constructor for a Container of unsigned integers.
        /*! \param c A container of unsigned integers.
            \pre No two adjacent values should be equal.
          */
        template<class Container>
        vlc_vector(const Container& c);

        //! Constructor for an int_vector_buffer of unsigned integers.
        template<uint8_t int_width>
        vlc_vector(int_vector_buffer<int_width>& v_buf);

        //! The number of elements in the vlc_vector.
        size_type size()const
        {
            return m_size;
        }
        //! Return the largest size that this container can ever have.
        static size_type max_size()
        {
            return int_vector<>::max_size()/2;
        }

        //!    Returns if the vlc_vector is empty.
        bool empty() const
        {
            return 0 == m_size;
        }

        //! Swap method for vlc_vector
        void swap(vlc_vector& v);

        //! Iterator that points to the first element of the vlc_vector.
        const const_iterator begin()const
        {
            return const_iterator(this, 0);
        }


        //! Iterator that points to the position after the last element of the vlc_vector.
        const const_iterator end()const
        {
            return const_iterator(this, this->m_size);
        }

        //! []-operator
        value_type operator[](size_type i)const;

        //! Serializes the vlc_vector to a stream.
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const;

        //! Load the vlc_vector from a stream.
        void load(std::istream& in);

        //! Returns the ith sample of vlc_vector
        value_type sample(const size_type i) const;

        uint32_t get_sample_dens() const;
        void set_sample_dens(const uint32_t sdens);
};

template<class t_coder, uint32_t t_dens, uint8_t t_width>
inline uint32_t vlc_vector<t_coder, t_dens, t_width>::get_sample_dens() const
{
    if (t_dens == 0)
        return m_sample_dens;
    else
        return t_dens;
}

template<class t_coder, uint32_t t_dens, uint8_t t_width>
inline void vlc_vector<t_coder, t_dens, t_width>::set_sample_dens(const uint32_t sdens)
{
    m_sample_dens = sdens;
}

template<class t_coder, uint32_t t_dens, uint8_t t_width>
inline typename vlc_vector<t_coder, t_dens,t_width>::value_type vlc_vector<t_coder, t_dens,t_width>::operator[](const size_type i)const
{
    assert(i+1 != 0);
    assert(i < m_size);
    size_type idx = i/get_sample_dens();
    return (t_coder::template decode<false, false, int*>(m_z.data(), m_sample_pointer[idx], i-t_dens*idx+1)) - 1;
}

template<class t_coder, uint32_t t_dens, uint8_t t_width>
void vlc_vector<t_coder, t_dens,t_width>::swap(vlc_vector<t_coder, t_dens,t_width>& v)
{
    if (this != &v) { // if v and _this_ are not the same object
        m_z.swap(v.m_z);                    // swap compressed bit streams
        m_sample_pointer.swap(v.m_sample_pointer);
        std::swap(m_size, v.m_size);// swap the number of elements
    }
}


template<class t_coder, uint32_t t_dens, uint8_t t_width>
template<class Container>
vlc_vector<t_coder, t_dens, t_width>::vlc_vector(const Container& c)
{
    clear(); // clear bit_vectors

    if (c.empty())  // if c is empty there is nothing to do...
        return;
    size_type samples = 0, z_size = 0;
//  (1) Calculate size of z
    for (size_type i=0; i < c.size(); ++i) {
        if (c[i]+1<1) {
            throw std::logic_error("vlc_vector cannot decode values smaller than 1!");
        }
        z_size += t_coder::encoding_length(c[i]+1);
    }
    samples = (c.size()+get_sample_dens()-1)/get_sample_dens();
//    (2) Write z
    m_sample_pointer = int_vector<>(samples+1, 0, bits::hi(z_size+1)+1);

    m_z.bit_resize(z_size);
    z_size = 0;
    uint64_t* z_data = t_coder::raw_data(m_z);
    uint8_t offset = 0;
    size_type no_sample = 0;
    for (size_type i=0, sample_cnt=0; i < c.size(); ++i, --no_sample) {
        if (!no_sample) { // add a sample pointer
            no_sample = get_sample_dens();
            m_sample_pointer[sample_cnt++] = z_size;
        }
        t_coder::encode(c[i]+1, z_data, offset);
        z_size += t_coder::encoding_length(c[i]+1);
    }
    m_size = c.size();
}

template<class t_coder, uint32_t t_dens, uint8_t t_width>
template<uint8_t int_width>
vlc_vector<t_coder, t_dens, t_width>::vlc_vector(int_vector_buffer<int_width>& v_buf)
{
    clear(); // clear bit_vectors
    size_type n = v_buf.size();
    if (n == 0)  // if c is empty there is nothing to do...
        return;
    size_type samples=0, z_size=0;
//  (1) Calculate size of z
    for (size_type i=0; i < n; ++i) {
        size_type x = v_buf[i]+1;
        if (x < 1) {
            throw std::logic_error("vlc_vector cannot decode values smaller than 1!");
        }
        z_size += t_coder::encoding_length(x);
    }
    samples = (n+get_sample_dens()-1)/get_sample_dens();
//    (2) Write z

    m_sample_pointer = int_vector<>(samples+1, 0, bits::hi(z_size+1)+1);  // add 1 for last entry

//     (b) Initilize bit_vector for encoded data
    m_z.bit_resize(z_size);
    z_size = 0;
    uint64_t* z_data = t_coder::raw_data(m_z);
    uint8_t offset = 0;

//     (c) Write sample values and deltas
    size_type no_sample = 0;
    for (size_type i=0, sample_cnt = 0; i < n; ++i, --no_sample) {
        if (!no_sample) { // add a sample pointer
            no_sample = get_sample_dens();
            m_sample_pointer[sample_cnt++] = z_size;
        }
        size_type x = v_buf[i]+1;
        t_coder::encode(x, z_data, offset);   // write encoded values
        z_size += t_coder::encoding_length(x);
    }
    m_size = n;
}

template<class t_coder, uint32_t t_dens, uint8_t t_width>
vlc_vector<>::size_type vlc_vector<t_coder, t_dens,t_width>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += write_member(m_size, out, child, "m_size");
    written_bytes += m_z.serialize(out, child, "m_z");
    written_bytes += m_sample_pointer.serialize(out, child, "m_sample_pointer");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<class t_coder, uint32_t t_dens, uint8_t t_width>
void vlc_vector<t_coder, t_dens,t_width>::load(std::istream& in)
{
    read_member(m_size, in);
    m_z.load(in);
    m_sample_pointer.load(in);
}

} // end namespace sdsl
#endif
