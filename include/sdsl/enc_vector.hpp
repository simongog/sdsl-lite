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
#include "int_vector.hpp"
#include "elias_delta_coder.hpp"
#include "iterators.hpp"

#ifndef SDSL_ENC_VECTOR
#define SDSL_ENC_VECTOR

//! Namespace for the succinct data structure library.
namespace sdsl
{

template<uint8_t fixedIntWidth>
struct enc_vector_trait {
    typedef int_vector<0> int_vector_type;
};

template<>
struct enc_vector_trait<32> {
    typedef int_vector<32> int_vector_type;
};

//! A generic immutable space-saving vector class for unsigned positiv integers. It encodes each integer with its self-delimiting code and still provides constant time access.
/*! The values of a enc_vector are immutable after the constructor call. The class
 *   can be parametrized with a self-delimiting codes (parameter Coder)
 *   and the sample density.
 *
 *  This class is a parameter of csa_sada.
 * @ingroup int_vector
 */
template<class Coder=coder::elias_delta,
         uint32_t SampleDens = 8, uint8_t fixedIntWidth=0>
class enc_vector
{
    public:
        typedef uint64_t 							value_type;  	// STL Container requirement
        typedef random_access_const_iterator<enc_vector> iterator;// STL Container requirement
        typedef iterator							const_iterator; // STL Container requirement
        typedef const value_type		 			reference;
        typedef const value_type 					const_reference;
        typedef const value_type*					const_pointer;
        typedef ptrdiff_t 							difference_type;// STL Container requirement
        typedef int_vector<>::size_type				size_type;		// STL Container requirement
        typedef Coder								coder;
        typedef typename enc_vector_trait<fixedIntWidth>::int_vector_type int_vector_type;
        static  const uint32_t 						sample_dens	= SampleDens;

        int_vector<0> 	m_z; 		// compressed bit stream
    private:
        int_vector_type   m_sample_vals_and_pointer;
        size_type		m_elements;    // number of elements
        uint32_t		m_sample_dens;

        // workaround function for the constructor
        void construct() {
            m_elements = 0;
            m_sample_dens = 8;
        }
        void copy(const enc_vector& v);

        void clear() {
            m_z.resize(0);
            m_elements = 0;
            m_sample_vals_and_pointer.resize(0);
        }

    public:
        //! Default Constuctor
        enc_vector() {
            construct();
        };
        //! Copy constructor
        /*! \param v The enc_vector to copy.
          	Required for the Assignable Concept of the STL
         */
        enc_vector(const enc_vector& v);

        //! Constructor for a Container of positive integers.
        /*! \param c A container of positive integers.
        	\pre No two adjacent values should be equal.
          */
        template<class Container>
        enc_vector(const Container& c) {
            construct();
            init(c);
        }

        //! Constructor for an int_vector_file_buffer of positive integers.
        /*
            \param v_buf A int_vector_file_buf.
        	\pre No two adjacent values should be equal.
        */
        template<uint8_t int_width, class size_type_class>
        enc_vector(int_vector_file_buffer<int_width, size_type_class>& v_buf) {
            construct();
            init(v_buf);
        }

        template<class Container>
        void init(const Container& c);

        template<uint8_t int_width, class size_type_class>
        void init(int_vector_file_buffer<int_width, size_type_class>& v_buf);

        //! Default Destructor
        ~enc_vector() {
        }

        //! The number of elements in the enc_vector.
        /*!
         	Required for the Container Concept of the STL.
        	\sa max_size
         */
        size_type size()const;

        //! Return the largest size that this container can ever have.
        /*! Required for the Container Concept of the STL.
         */
        static size_type max_size();

        //!	Returns if the enc_vector is empty.
        /*! Equivalent to size() == 0.
         *
         * 	Required for the STL Container Concept.
         *  \sa size()
         */
        bool empty() const;

        //! Swap method for enc_vector
        /*! The swap method can be defined in terms of assignment.
        	This requires three assignments, each of which, for a container type, is linear
        	in the container's size. In a sense, then, a.swap(b) is redundant.
        	This implementation guaranties a run-time complexity that is constant rather than linear.
        	\param v enc_vector to swap.

        	Required for the Assignable Conecpt of the STL.
          */
        void swap(enc_vector& v);

        //! Iterator that points to the first element of the enc_vector.
        /*!
         * 	Required for the Container Concept of the STL.
         *  \sa end()
         */
        const const_iterator begin()const;

        //! Iterator that points to the position after the last element of the enc_vector.
        /*!
         *	Required for the Container Concept of the STL
         *  \sa begin()
         */
        const const_iterator end()const;

        // Iterator that points to the last element of the enc_vector.
        /*
         * 	Required for the Container Concept of the STL.
         *  \sa rend()
         */
//		reverse_iterator rbegin()const;

        // Iterator that points to the position before the first element of the enc_vector.
        /*
         *	Required for the Container Concept of the STL
         *  \sa rbegin()
         */
//		reverse_iterator rend()const;

        //! []-operator
        /*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
         *
         *  Required for the STL Random Access Container Concept.
         */
        value_type operator[](size_type i)const;

        //! Assignment Operator
        /*!
         *	Required for the Assignable Concept of the STL.
         */
        enc_vector& operator=(const enc_vector& v);

        //! Equality Operator
        /*! Two enc_vectors are equal if all member variables are equal
         *  (including the sample density of the enc_vectors).
         *  \note If the sample density is not equal you should use
         *  SDSAlgorithm::equal_container_values to compare two enc_vectors.
         *
         * 	Required for the Equality Comparable Concept of the STL.
         *  \sa operator!=
         */
        bool operator==(const enc_vector& v)const;

        //! Unequality Operator
        /*! Two enc_vectors are unuequal if not all member variables are equal
         *  (including the sample density of the enc_vectors).
         *  \note If the sample density is not equal you should use
         *  SDSAlgorithm::equal_container_values to compare two enc_vectors.
         *
         * 	Required for the Equality Comparable Concept of the STL.
         *  \sa operator==
         */
        bool operator!=(const enc_vector& v)const;

        //! Serialzes the enc_vector to a stream.
        /*! \param out Outstream to write the data structure.
            \return The number of written bytes.
         */
        size_type serialize(std::ostream& out) const;

        //! Load the enc_vector from a stream.
        void load(std::istream& in);

        //! Returns the ith sample of enc_vector
        /*! \param i The index of the sample. 0 <= i < size()/get_sample_dens()
         *  \return The value of the ith sample.
         */
        value_type sample(const size_type i) const;

        uint32_t get_sample_dens() const;
        void set_sample_dens(const uint32_t sample_dens);

        /*!
         * \param i The index of the sample for which all values till the next sample should be decoded. 0 <= i < size()/get_sample_dens()
         * \param it A pointer to a uint64_t vector, whereto the values should be written
         */
//		template<class Iterator>
        void get_inter_sampled_values(const size_type i, uint64_t*& it)const {
            *(it++) = 0;
            if (i*SampleDens + SampleDens - 1 < size()) {
                Coder::template decode<true, true>(m_z.data(), m_sample_vals_and_pointer[(i<<1)+1], SampleDens - 1, it);
            } else {
                assert(i*SampleDens < size());
                Coder::template decode<true, true>(m_z.data(), m_sample_vals_and_pointer[(i<<1)+1], size()-i*SampleDens - 1, it);
            }
        };

#ifdef MEM_INFO
        void mem_info(std::string label="")const {
            if (label=="")
                label = "enc vector";
            size_type bytes = util::get_size_in_bytes(*this);
            std::cout << "list(label = \""<<label<<"\", size = "<< bytes/(1024.0*1024.0) <<"\n,";
            m_z.mem_info("variable-length code");
            std::cout<<",";
            m_sample_vals_and_pointer.mem_info("samples");
            std::cout << ")\n";
        }
#endif
};


template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
inline uint32_t enc_vector<Coder, SampleDens, fixedIntWidth>::get_sample_dens() const
{
    if (SampleDens == 0)
        return m_sample_dens;
    else
        return SampleDens;
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
inline void enc_vector<Coder, SampleDens, fixedIntWidth>::set_sample_dens(const uint32_t sample_dens)
{
    m_sample_dens = sample_dens;
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
inline typename enc_vector<Coder, SampleDens,fixedIntWidth>::value_type enc_vector<Coder, SampleDens,fixedIntWidth>::operator[](const size_type i)const
{
    assert(i+1 != 0);
#ifdef SDSL_DEBUG
    if (i >= m_elements) {
        throw std::out_of_range("OUT_OF_RANGE_ERROR: enc_vector::operator[](size_type); i >= size()!");
        return 0;
    }
#endif
    size_type idx = i/get_sample_dens();
    return m_sample_vals_and_pointer[idx<<1] + Coder::decode_prefix_sum(m_z.data(), m_sample_vals_and_pointer[(idx<<1)+1], i-SampleDens*idx);
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
inline typename enc_vector<Coder, SampleDens,fixedIntWidth>::value_type enc_vector<Coder, SampleDens,fixedIntWidth>::sample(const size_type i)const
{
    assert(i*get_sample_dens()+1 != 0);
#ifdef SDSL_DEBUG
    if (i*get_sample_dens() >= m_elements) {
        throw std::out_of_range("OUT_OF_RANGE_ERROR: enc_vector::sample(size_type); i*get_sample_dens() >= size()!");
        return 0;
    }
#endif
    return m_sample_vals_and_pointer[i<<1];
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
inline enc_vector<>::size_type enc_vector<Coder, SampleDens,fixedIntWidth>::size()const
{
    return m_elements;
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
inline enc_vector<>::size_type enc_vector<Coder, SampleDens,fixedIntWidth>::max_size()
{
    return int_vector<>::max_size()/2; // each element could possible occupy double space with selfdelimiting codes
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
inline bool enc_vector<Coder, SampleDens,fixedIntWidth>::empty()const
{
    return 0==m_elements;
}


template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
void enc_vector<Coder, SampleDens,fixedIntWidth>::copy(const enc_vector<Coder, SampleDens,fixedIntWidth>& v)
{
    m_z					= v.m_z;				// copy compressed bit stream
    m_sample_vals_and_pointer		= v.m_sample_vals_and_pointer;      // copy sample values
    m_elements			= v.m_elements;			// copy number of stored elements
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
enc_vector<Coder, SampleDens,fixedIntWidth>::enc_vector(const enc_vector& v)
{
    copy(v);
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
enc_vector<Coder, SampleDens,fixedIntWidth>& enc_vector<Coder, SampleDens,fixedIntWidth>::operator=(const enc_vector<Coder, SampleDens,fixedIntWidth>& v)
{
    if (this != &v) {// if v and _this_ are not the same object
        copy(v);
    }
    return *this;
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
bool enc_vector<Coder, SampleDens,fixedIntWidth>::operator==(const enc_vector<Coder, SampleDens,fixedIntWidth>& v)const
{
    if (this == &v)
        return true;
    return	 	m_elements == v.m_elements
                and	m_z == v.m_z
                and	m_sample_vals_and_pointer == v.m_sample_vals_and_pointer;
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
bool enc_vector<Coder, SampleDens,fixedIntWidth>::operator!=(const enc_vector<Coder, SampleDens,fixedIntWidth>& v)const
{
    return !(*this == v);
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
void enc_vector<Coder, SampleDens,fixedIntWidth>::swap(enc_vector<Coder, SampleDens,fixedIntWidth>& v)
{
    if (this != &v) { // if v and _this_ are not the same object
        m_z.swap(v.m_z);					// swap compressed bit streams
        m_sample_vals_and_pointer.swap(v.m_sample_vals_and_pointer);
        std::swap(m_elements, v.m_elements);// swap the number of elements
    }
}


template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
template<class Container>
void enc_vector<Coder, SampleDens,fixedIntWidth>::init(const Container& c)
{
    // clear bit_vectors
    clear();

    if (c.empty())  // if c is empty there is nothing to do...
        return;
    typename Container::const_iterator	it		 	= c.begin(), end = c.end();
    typename Container::value_type 		v1			= *it, v2, max_value=0, max_sample_value=0, x;
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
            if (max_value < v2-v1) max_value = v2 - v1;
            z_size += Coder::encoding_length(v2-v1);
        }
        v1=v2;
    }
//	(2) Write sample values and deltas
    {
        if (max_sample_value > z_size+1)
            m_sample_vals_and_pointer.set_int_width(bit_magic::l1BP(max_sample_value) + 1);
        else
            m_sample_vals_and_pointer.set_int_width(bit_magic::l1BP(z_size+1) + 1);
        m_sample_vals_and_pointer.resize(2*samples+2); // add 2 for last entry
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
                if (v2 == v1) {
                    throw std::logic_error("enc_vector cannot decode adjacent equal values!");
                }
                z_size += Coder::encoding_length(x);
            }
            v1=v2;
        }
        *sv_it = 0; ++sv_it;        // initialize
        *sv_it = z_size+1; ++sv_it; // last entry

        m_z.bit_resize(z_size);
        uint64_t* z_data = Coder::raw_data(m_z);
        uint8_t offset = 0;
        no_sample = 0;
        for (it = c.begin(); it != end; ++it, --no_sample) {
            v2 = *it;
            if (!no_sample) { // add a sample
                no_sample = get_sample_dens();
            } else {
                Coder::encode(v2-v1, z_data, offset);
            }
            v1=v2;
        }
//		Coder::encode(delta_c, m_z); // encode delta_c to m_z
    }
//	delta_c.resize(0);
//std::cerr<<"Finished "<<std::endl;,
    m_elements = c.size();
}


template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
template<uint8_t int_width, class size_type_class>
void enc_vector<Coder, SampleDens,fixedIntWidth>::init(int_vector_file_buffer<int_width, size_type_class>& v_buf)
{
    // clear bit_vectors
    clear();
    size_type n = v_buf.int_vector_size;
    if (n == 0)  // if c is empty there is nothing to do...
        return;
    v_buf.reset();
    value_type 	v1=0, v2=0, max_value=0, max_sample_value=0;
    size_type samples=0, z_size=0;
    const size_type sd = get_sample_dens();
//  (1) Calculate maximal value of samples and of deltas
    for (size_type i=0, r_sum=0, r = v_buf.load_next_block(), no_sample = 0; r_sum < n;) {
        for (; i < r_sum+r; ++i, --no_sample) {
            v2 = v_buf[i-r_sum];
            if (!no_sample) { // is sample
                no_sample = sd;
                if (max_sample_value < v2) max_sample_value = v2;
                ++samples;
            } else {
                if (max_value < v2-v1) max_value = v2-v1;
                if (v2 == v1) {
                    throw std::logic_error("enc_vector cannot decode adjacent equal values!");
                }
                z_size += Coder::encoding_length(v2-v1);
            }
            v1 = v2;
        }
        r_sum += r; r = v_buf.load_next_block();
    }

//	(2) Write sample values and deltas
//     (a) Initialize array for sample values and pointers
    if (max_sample_value > z_size+1)
        m_sample_vals_and_pointer.set_int_width(bit_magic::l1BP(max_sample_value) + 1);
    else
        m_sample_vals_and_pointer.set_int_width(bit_magic::l1BP(z_size+1) + 1);
    m_sample_vals_and_pointer.resize(2*samples+2); // add 2 for last entry

//     (b) Initilize bit_vector for encoded data
    m_z.bit_resize(z_size);
    uint64_t* z_data = Coder::raw_data(m_z);
    uint8_t offset = 0;

//     (c) Write sample values and deltas
    v_buf.reset();
    z_size = 0;
    for (size_type i=0, j=0, r_sum=0, r = v_buf.load_next_block(), no_sample = 0; r_sum < n;) {
        for (; i < r_sum+r; ++i, --no_sample) {
            v2 = v_buf[i-r_sum];
            if (!no_sample) { // is sample
                no_sample = sd;
                m_sample_vals_and_pointer[j++] = v2;	// write samples
                m_sample_vals_and_pointer[j++] = z_size;// write pointers
            } else {
                z_size += Coder::encoding_length(v2-v1);
                Coder::encode(v2-v1, z_data, offset);   // write encoded values
            }
            v1 = v2;
        }
        r_sum += r; r = v_buf.load_next_block();
    }
    m_elements = n;
}


template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
enc_vector<>::size_type enc_vector<Coder, SampleDens,fixedIntWidth>::serialize(std::ostream& out) const
{
    size_type written_bytes = 0;
    out.write((char*) &m_elements, sizeof(m_elements));
    written_bytes += sizeof(m_elements);
    written_bytes += m_z.serialize(out);
    written_bytes += m_sample_vals_and_pointer.serialize(out);
    return written_bytes;
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
void enc_vector<Coder, SampleDens,fixedIntWidth>::load(std::istream& in)
{
    in.read((char*) &m_elements, sizeof(m_elements));
    m_z.load(in);
    m_sample_vals_and_pointer.load(in);
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
const typename enc_vector<Coder,SampleDens,fixedIntWidth>::const_iterator enc_vector<Coder, SampleDens,fixedIntWidth>::begin()const
{
    return const_iterator(this, 0);
}

template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
const typename enc_vector<Coder,SampleDens,fixedIntWidth>::const_iterator enc_vector<Coder, SampleDens,fixedIntWidth>::end()const
{
    return const_iterator(this, this->m_elements);
}


} // end namespace sdsl

#endif
