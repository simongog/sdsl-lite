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
#include "int_vector.hpp"
#include "elias_delta_coder.hpp"
#include "iterators.hpp"

#ifndef SDSL_VLC_VECTOR
#define SDSL_VLC_VECTOR

//! Namespace for the succinct data structure library.
namespace sdsl
{

template<uint8_t fixedIntWidth>
struct vlc_vector_trait {
    typedef int_vector<0> int_vector_type;
};

template<>
struct vlc_vector_trait<32> {
    typedef int_vector<32> int_vector_type;
};

//! A generic immutable space-saving vector class for unsigned positive integers.
/*! The values of a vlc_vector are immutable after the constructor call. The class
 *   could be parametrized with a self-delimiting code Coder and the sample density.
 *
 */
template<class Coder=coder::elias_delta,
         uint32_t SampleDens = 16, 
		 uint8_t fixedIntWidth=0>
class vlc_vector
{
    public:
        typedef uint64_t 							value_type;  	// STL Container requirement
        typedef random_access_const_iterator<vlc_vector> iterator;// STL Container requirement
        typedef iterator							const_iterator; // STL Container requirement
        typedef const value_type		 			reference;
        typedef const value_type 					const_reference;
        typedef const value_type*					const_pointer;
        typedef ptrdiff_t 							difference_type;// STL Container requirement
        typedef int_vector<>::size_type				size_type;		// STL Container requirement
        typedef Coder								coder;
        typedef typename vlc_vector_trait<fixedIntWidth>::int_vector_type int_vector_type;
        static  const uint32_t 						sample_dens	= SampleDens;


        bit_vector 		m_z; 		// compressed bit stream
    private:
        int_vector_type m_sample_pointer;
        size_type		m_elements;    // number of elements
        uint32_t		m_sample_dens;

        void copy(const vlc_vector& v);

        void clear() {
            m_z.resize(0);
            m_elements = 0;
            m_sample_pointer.resize(0);
        }

    public:
        //! Default Constuctor
        vlc_vector() : m_elements(0), m_sample_dens(16) {}
        //! Copy constructor
        vlc_vector(const vlc_vector& v);

        //! Constructor for a Container of positive integers.
        /*! \param c A container of positive integers.
        	\pre No two adjacent values should be equal.
          */
        template<class Container>
        vlc_vector(const Container& c);

        //! Constructor for an int_vector_file_buffer of positive integers.
        template<uint8_t int_width>
        vlc_vector(int_vector_file_buffer<int_width>& v_buf);

        //! The number of elements in the vlc_vector.
        size_type size()const;

        //! Return the largest size that this container can ever have.
        static size_type max_size();

        //!	Returns if the vlc_vector is empty.
        bool empty() const;

        //! Swap method for vlc_vector
        void swap(vlc_vector& v);

        //! Iterator that points to the first element of the vlc_vector.
        const const_iterator begin()const;

        //! Iterator that points to the position after the last element of the vlc_vector.
        const const_iterator end()const;

        //! []-operator
        value_type operator[](size_type i)const;

        //! Assignment Operator
        /*!
         *	Required for the Assignable Concept of the STL.
         */
        vlc_vector& operator=(const vlc_vector& v);

		//! Equality Operator
        bool operator==(const vlc_vector& v)const;

        //! Inequality Operator
        bool operator!=(const vlc_vector& v)const;

        //! Serialzes the vlc_vector to a stream.
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const;

        //! Load the vlc_vector from a stream.
        void load(std::istream& in);

        //! Returns the ith sample of vlc_vector
        value_type sample(const size_type i) const;

        uint32_t get_sample_dens() const;
        void set_sample_dens(const uint32_t sample_dens);
};


template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
inline uint32_t vlc_vector<Coder, SampleDens, fixedIntWidth>::get_sample_dens() const {
    if (SampleDens == 0)
        return m_sample_dens;
    else
        return SampleDens;
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
inline void vlc_vector<Coder, SampleDens, fixedIntWidth>::set_sample_dens(const uint32_t sample_dens) {
    m_sample_dens = sample_dens;
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
inline typename vlc_vector<Coder, SampleDens,fixedIntWidth>::value_type vlc_vector<Coder, SampleDens,fixedIntWidth>::operator[](const size_type i)const {
    assert(i+1 != 0);
#ifdef SDSL_DEBUG
    if (i >= m_elements) {
        throw std::out_of_range("OUT_OF_RANGE_ERROR: vlc_vector::operator[](size_type); i >= size()!");
        return 0;
    }
#endif
    size_type idx = i/get_sample_dens();
    return (Coder::template decode<false, false, int*>(m_z.data(), m_sample_pointer[idx], i-SampleDens*idx+1)) - 1;
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
inline vlc_vector<>::size_type vlc_vector<Coder, SampleDens,fixedIntWidth>::size()const {
    return m_elements;
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
inline vlc_vector<>::size_type vlc_vector<Coder, SampleDens,fixedIntWidth>::max_size() {
    return int_vector<>::max_size()/2; // each element could possible occupy double space with selfdelimiting codes
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
inline bool vlc_vector<Coder, SampleDens,fixedIntWidth>::empty()const {
    return 0==m_elements;
}


template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
void vlc_vector<Coder, SampleDens,fixedIntWidth>::copy(const vlc_vector<Coder, SampleDens,fixedIntWidth>& v) {
    m_z					= v.m_z;				// copy compressed bit stream
    m_sample_pointer	= v.m_sample_pointer;   // copy sample values
    m_elements			= v.m_elements;			// copy number of stored elements
	m_sample_dens       = v.m_sample_dens;
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
vlc_vector<Coder, SampleDens,fixedIntWidth>::vlc_vector(const vlc_vector& v) : m_elements(0), m_sample_dens(16) {
    copy(v);
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
vlc_vector<Coder, SampleDens,fixedIntWidth>& vlc_vector<Coder, SampleDens,fixedIntWidth>::operator=(const vlc_vector<Coder, SampleDens,fixedIntWidth>& v) {
    if (this != &v) {// if v and _this_ are not the same object
        copy(v);
    }
    return *this;
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
bool vlc_vector<Coder, SampleDens,fixedIntWidth>::operator==(const vlc_vector<Coder, SampleDens,fixedIntWidth>& v)const {
    if (this == &v)
        return true;
    return	 	m_elements == v.m_elements
                and	m_z == v.m_z
                and	m_sample_pointer == v.m_sample_pointer;
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
bool vlc_vector<Coder, SampleDens,fixedIntWidth>::operator!=(const vlc_vector<Coder, SampleDens,fixedIntWidth>& v)const {
    return !(*this == v);
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
void vlc_vector<Coder, SampleDens,fixedIntWidth>::swap(vlc_vector<Coder, SampleDens,fixedIntWidth>& v) {
    if (this != &v) { // if v and _this_ are not the same object
        m_z.swap(v.m_z);					// swap compressed bit streams
        m_sample_pointer.swap(v.m_sample_pointer);
        std::swap(m_elements, v.m_elements);// swap the number of elements
    }
}


template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
template<class Container>
vlc_vector<Coder, SampleDens, fixedIntWidth>::vlc_vector(const Container& c) : m_elements(0), m_sample_dens(16) {
    clear(); // clear bit_vectors

    if (c.empty())  // if c is empty there is nothing to do...
        return;
    size_type samples = 0, z_size = 0;
//  (1) Calculate size of z
    for (size_type i=0; i < c.size(); ++i) {
        if (c[i]+1<1) {
            throw std::logic_error("vlc_vector cannot decode values smaller than 1!");
        }
        z_size += Coder::encoding_length(c[i]+1);
    }
    samples = (c.size()+get_sample_dens()-1)/get_sample_dens();
//	(2) Write z

    m_sample_pointer.set_int_width(bit_magic::l1BP(z_size+1) + 1);
    m_sample_pointer.resize(samples+1); // add 1 for last entry

    std::cout<<"z_size="<<z_size<<std::endl;
    std::cout<<"int_width = "<<bit_magic::l1BP(z_size+1)+1<<std::endl;
    std::cout<<samples<<std::endl;

    m_z.bit_resize(z_size);
    z_size = 0;
    uint64_t* z_data = Coder::raw_data(m_z);
    uint8_t offset = 0;
    size_type no_sample = 0;
    for (size_type i=0, sample_cnt=0; i < c.size(); ++i, --no_sample) {
        if (!no_sample) { // add a sample pointer
            no_sample = get_sample_dens();
            m_sample_pointer[sample_cnt++] = z_size;
        }
        Coder::encode(c[i]+1, z_data, offset);
        z_size += Coder::encoding_length(c[i]+1);
    }
    m_elements = c.size();
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
template<uint8_t int_width>
vlc_vector<Coder, SampleDens, fixedIntWidth>::vlc_vector(int_vector_file_buffer<int_width>& v_buf) 
	: m_elements(0), m_sample_dens(16) {
    clear(); // clear bit_vectors
    size_type n = v_buf.int_vector_size;
    if (n == 0)  // if c is empty there is nothing to do...
        return;
    v_buf.reset();
    size_type samples=0, z_size=0;
//  (1) Calculate size of z
    for (size_type i=0, r_sum=0, r = v_buf.load_next_block(); r_sum < n;) {
        for (; i < r_sum+r; ++i) {
            size_type x = v_buf[i-r_sum]+1;
            if (x < 1) {
                throw std::logic_error("vlc_vector cannot decode values smaller than 1!");
            }
            z_size += Coder::encoding_length(x);
        }
        r_sum += r; r = v_buf.load_next_block();
    }
    samples = (n+get_sample_dens()-1)/get_sample_dens();
//	(2) Write z

    m_sample_pointer.set_int_width(bit_magic::l1BP(z_size+1) + 1);
    m_sample_pointer.resize(samples+1); // add 1 for last entry

//     (b) Initilize bit_vector for encoded data
    m_z.bit_resize(z_size);
    z_size = 0;
    uint64_t* z_data = Coder::raw_data(m_z);
    uint8_t offset = 0;

//     (c) Write sample values and deltas
    v_buf.reset();
    size_type no_sample = 0;
    for (size_type i=0, sample_cnt = 0, r_sum=0, r = v_buf.load_next_block(); r_sum < n;) {
        for (; i < r_sum+r; ++i, --no_sample) {
            if (!no_sample) { // add a sample pointer
                no_sample = get_sample_dens();
                m_sample_pointer[sample_cnt++] = z_size;
            }
            size_type x = v_buf[i-r_sum]+1;
            Coder::encode(x, z_data, offset);   // write encoded values
            z_size += Coder::encoding_length(x);
        }
        r_sum += r; r = v_buf.load_next_block();
    }
    m_elements = n;
}


template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
vlc_vector<>::size_type vlc_vector<Coder, SampleDens,fixedIntWidth>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const {
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
	written_bytes += util::write_member(m_elements, out, child, "m_elements");
    written_bytes += m_z.serialize(out, child, "m_z");
    written_bytes += m_sample_pointer.serialize(out, child, "m_sample_pointer");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
void vlc_vector<Coder, SampleDens,fixedIntWidth>::load(std::istream& in) {
    in.read((char*) &m_elements, sizeof(m_elements));
    m_z.load(in);
    m_sample_pointer.load(in);
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
const typename vlc_vector<Coder,SampleDens,fixedIntWidth>::const_iterator vlc_vector<Coder, SampleDens,fixedIntWidth>::begin()const {
    return const_iterator(this, 0);
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
const typename vlc_vector<Coder,SampleDens,fixedIntWidth>::const_iterator vlc_vector<Coder, SampleDens,fixedIntWidth>::end()const {
    return const_iterator(this, this->m_elements);
}


} // end namespace sdsl

#endif
