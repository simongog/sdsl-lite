/* sdsl - succinct data structures library
    Copyright (C) 2008-2013 Simon Gog

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
/*! \file int_vector.hpp
    \brief int_vector.hpp contains the sdsl::int_vector class.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_INT_VECTOR
#define INCLUDED_SDSL_INT_VECTOR

#include <sys/mman.h>

#define HUGE_LEN 1073741824
#define HUGE_PROTECTION (PROT_READ | PROT_WRITE)
#define HUGE_FLAGS (MAP_HUGETLB | MAP_ANONYMOUS | MAP_PRIVATE)

#include "compatibility.hpp"
#include "bits.hpp"
#include "structure_tree.hpp"
#include "util.hpp"
#include "io.hpp"
#include "config.hpp"
#include "uintx_t.hpp"

#include "memory_management.hpp"
#include "ram_fs.hpp"
#include "sfstream.hpp"

#include <iosfwd>    // forward declaration of ostream
#include <stdexcept> // for exceptions
#include <iostream>  // for cerr
#include <typeinfo>
#include <cassert>
#include <iterator>
#include <cstdlib>
#include <cstddef>
#include <ctime>    // for rand initialization
#include <cstring>  // for memcpy
#include <ostream>
#include <istream>
#include <string>
#include <initializer_list>


//! Namespace for the succinct data structure library.
namespace sdsl
{

typedef uint64_t std_size_type_for_int_vector;

template<uint8_t t_width=0>
class int_vector;

template<class int_vector_type>
class mm_item;

namespace algorithm
{
template<uint8_t t_width>
static void calculate_sa(const unsigned char* c,
                         typename int_vector<t_width>::size_type len,
                         int_vector<t_width>& sa);
}


//! bit_vector is a specialization of the int_vector.
typedef int_vector<1> bit_vector;

template<class t_int_vector>
class int_vector_reference;

template<class t_int_vector>
class int_vector_iterator_base;

template<class t_int_vector>
class int_vector_iterator;

template<class t_int_vector>
class int_vector_const_iterator;

template<uint8_t b, uint8_t t_patter_len>  // forward declaration
class rank_support_v;

class rank_support;

class select_support;

template<uint8_t t_bit_pattern, uint8_t t_pattern_len>
class select_support_mcl;

namespace coder
{
class fibonacci;
class elias_delta;
}

template<uint8_t=0>
class int_vector_file_buffer;

class char_array_serialize_wrapper;

template<uint8_t t_width>
struct int_vector_trait {
    typedef uint64_t                                    value_type;
    typedef int_vector<t_width>                         int_vector_type;
    typedef int_vector_reference<int_vector_type>       reference;
    typedef const uint64_t                              const_reference;
    typedef uint8_t                                     int_width_type;
    typedef int_vector_iterator<int_vector_type>        iterator;
    typedef int_vector_const_iterator<int_vector_type>  const_iterator;
    // Sets int_width to new_int_width
    static void width(int_width_type& int_width, const uint8_t new_int_width) {
        if (t_width==1)
            return;
        if (new_int_width>0 and new_int_width<=64)
            int_width = new_int_width;
        else
            int_width = 64;
    }

    // read the size and int_width of a bit_vector
    static void read_header(int_vector_size_type& size, int_width_type& int_width, std::istream& in) {
        read_member(size, in);
        if (0 == t_width) {
            read_member(int_width, in);
            width(int_width, int_width);
        }
    }

    static iterator begin(int_vector_type* v, uint64_t*) {
        return iterator(v, 0);
    }
    static iterator end(int_vector_type* v, uint64_t*, int_vector_size_type) {
        return iterator(v, v->size()*v->width()) ;
    }
    static const_iterator begin(const int_vector_type* v, const uint64_t*) {
        return const_iterator(v, 0);
    }
    static const_iterator end(const int_vector_type* v, const uint64_t*, int_vector_size_type) {
        return const_iterator(v, v->size()*v->width());
    }
};

template<>
struct int_vector_trait<64> {
    typedef uint64_t        value_type;
    typedef int_vector<64>  int_vector_type;
    typedef uint64_t&       reference;
    typedef const uint64_t  const_reference;
    typedef const uint8_t   int_width_type;
    typedef uint64_t*       iterator;
    typedef const uint64_t* const_iterator;


    static void width(int_width_type&, const uint8_t) {}

    // read the size and int_width of a bit_vector
    static void read_header(int_vector_size_type& size, int_width_type&, std::istream& in) {
        read_member(size, in);
    }

    static iterator begin(int_vector_type*, uint64_t* begin) {
        return begin;
    }
    static iterator end(int_vector_type*, uint64_t* begin, int_vector_size_type size) {
        return begin+size;
    }
    static const_iterator begin(const int_vector_type*, const uint64_t* begin) {
        return begin;
    }
    static const_iterator end(const int_vector_type*, const uint64_t* begin, int_vector_size_type size) {
        return begin+size;
    }
};

template<>
struct int_vector_trait<32> {
    typedef uint32_t        value_type;
    typedef int_vector<32>  int_vector_type;
    typedef uint32_t&       reference;
    typedef const uint32_t  const_reference;
    typedef const uint8_t   int_width_type;
    typedef uint32_t*       iterator;
    typedef const uint32_t* const_iterator;
    static void width(int_width_type&, const uint8_t) {}

    // read the size and int_width of a bit_vector
    static void read_header(int_vector_size_type& size, int_width_type&, std::istream& in) {
        read_member(size, in);
    }

    static iterator begin(int_vector_type*, uint64_t* begin) {
        return (uint32_t*)begin;
    }
    static iterator end(int_vector_type*, uint64_t* begin, int_vector_size_type size) {
        return ((uint32_t*)begin)+size;
    }
    static const_iterator begin(const int_vector_type*, const uint64_t* begin) {
        return (uint32_t*)begin;
    }
    static const_iterator end(const int_vector_type*, const uint64_t* begin, int_vector_size_type size) {
        return ((uint32_t*)begin)+size;
    }
};

template<>
struct int_vector_trait<16> {
    typedef uint16_t        value_type;
    typedef int_vector<16>  int_vector_type;
    typedef uint16_t&       reference;
    typedef const uint16_t  const_reference;
    typedef const uint8_t   int_width_type;
    typedef uint16_t*       iterator;
    typedef const uint16_t* const_iterator;
    static void width(int_width_type&, const uint8_t) {}

    // read the size and int_width of a bit_vector
    static void read_header(int_vector_size_type& size, int_width_type&, std::istream& in) {
        read_member(size, in);
    }

    static iterator begin(int_vector_type*, uint64_t* begin) {
        return (uint16_t*)begin;
    }
    static iterator end(int_vector_type*, uint64_t* begin, int_vector_size_type size) {
        return ((uint16_t*)begin)+size;
    }
    static const_iterator begin(const int_vector_type*, const uint64_t* begin) {
        return (uint16_t*)begin;
    }
    static const_iterator end(const int_vector_type*, const uint64_t* begin, int_vector_size_type size) {
        return ((uint16_t*)begin)+size;
    }
};

template<>
struct int_vector_trait<8> {
    typedef uint8_t         value_type;
    typedef int_vector<8>   int_vector_type;
    typedef uint8_t&        reference;
    typedef const uint8_t   const_reference;
    typedef const uint8_t   int_width_type;
    typedef uint8_t*        iterator;
    typedef const uint8_t*  const_iterator;
    static void width(int_width_type&, const uint8_t) {}

    // read the size and int_width of a bit_vector
    static void read_header(int_vector_size_type& size, int_width_type&, std::istream& in) {
        read_member(size, in);
    }

    static iterator begin(int_vector_type*, uint64_t* begin) {
        return (uint8_t*)begin;
    }
    static iterator end(int_vector_type*, uint64_t* begin, int_vector_size_type size) {
        return ((uint8_t*)begin)+size;
    }
    static const_iterator begin(const int_vector_type*, const uint64_t* begin) {
        return (uint8_t*)begin;
    }
    static const_iterator end(const int_vector_type*, const uint64_t* begin, int_vector_size_type size) {
        return ((uint8_t*)begin)+size;
    }
};

//! A generic vector class for integers of width \f$w\in [1..64]\f$.
/*! \author Simon Gog
 *
 *    This generic vector class could be used to generate a vector
 *    that contains integers of fixed width \f$w\in [1..64]\f$.
 *
 *  \tparam t_width Width of the integer. If set to `0` it is variable
 *          during runtime, otherwise fixed at compile time.
 *  @ingroup int_vector
 */
template<uint8_t t_width>
class int_vector
{
    private:
        static_assert(t_width <= 64 , "int_vector: width of must be at most 64bits.");
    public:
        typedef typename int_vector_trait<t_width>::value_type      value_type;
        typedef typename int_vector_trait<t_width>::iterator        iterator;
        typedef typename int_vector_trait<t_width>::const_iterator  const_iterator;
        typedef typename int_vector_trait<t_width>::reference       reference;
        typedef typename int_vector_trait<t_width>::const_reference const_reference;
        typedef int_vector_reference<int_vector>*                   pointer;
        typedef const value_type*                                   const_pointer;
        typedef ptrdiff_t                                           difference_type;
        typedef int_vector_size_type                                size_type;
        typedef typename int_vector_trait<t_width>::int_width_type  int_width_type;
        typedef rank_support_v<1,1>                                 rank_1_type;
        typedef rank_support_v<0,1>                                 rank_0_type;
        typedef select_support_mcl<1,1>                             select_1_type;
        typedef select_support_mcl<0,1>                             select_0_type;
        typedef int_vector_file_buffer<t_width>                     file_buffer;

        friend struct int_vector_trait<t_width>;
        friend class  int_vector_iterator_base<int_vector>;
        friend class  int_vector_iterator<int_vector>;
        friend class  int_vector_const_iterator<int_vector>;
        friend class  coder::elias_delta;
        friend class  coder::fibonacci;
        friend class  int_vector_file_buffer<t_width>;
        friend class  mm_item<int_vector>;

        friend void mm::realloc<int_vector>(int_vector&, const size_type);
        friend void util::set_random_bits<int_vector>(int_vector& v, int);
        friend void util::_set_zero_bits<int_vector>(int_vector&);
        friend void util::_set_one_bits<int_vector>(int_vector&);
        friend void util::bit_compress<int_vector>(int_vector&);
        friend void util::set_to_value<int_vector>(int_vector&, uint64_t);
        friend bool load_vector_from_file<int_vector>(int_vector&, const std::string&,uint8_t,uint8_t);
        friend void algorithm::calculate_sa<t_width>(const unsigned char* c, typename int_vector<t_width>::size_type len, int_vector<t_width>& sa);

        enum { fixed_int_width = t_width }; // make template parameter accessible

    private:

        size_type      m_size;  //!< Number of bits needed to store int_vector.
        uint64_t*      m_data;  //!< Pointer to the memory for the bits.
        int_width_type m_width; //!< Width of the integers.

    public:

        //! Constructor for int_vector.
        /*! \param size          Number of elements. Default value is 0.
            \param default_value Initialize all value to `default value`.
            \param int_width     The width of each integer.
            \sa resize, width
         */
        int_vector(size_type size = 0, value_type default_value = 0,
                   uint8_t int_width = t_width);

        //! Constructor for initializer_list.
        template<class t_T>
        int_vector(std::initializer_list<t_T> il) : int_vector() {
            resize(il.size());
            size_type idx = 0;
            for (auto x : il) {
                (*this)[idx++] = x;
            }
        }

        //! Move constructor.
        int_vector(int_vector&& v);

        //! Copy constructor.
        int_vector(const int_vector& v);

        //! Destructor.
        ~int_vector();

        //! Equivalent to size() == 0.
        bool empty() const {
            return 0==m_size;
        }

        //! Swap method for int_vector.
        void swap(int_vector& v);

        //! Resize the int_vector in terms of elements.
        /*! \param size The size to resize the int_vector in terms of elements.
         */
        void resize(const size_type size) {
            bit_resize(size * width());
        }

        //! Resize the int_vector in terms of bits.
        /*! \param size The size to resize the int_vector in terms of bits.
         */
        void bit_resize(const size_type size);

        //! The number of elements in the int_vector.
        /*! \sa max_size, bit_size, capacity
         */
        size_type size() const {
            return m_size/m_width;
        }

        //! Maximum size of the int_vector.
        /*! \sa size, bit_size, capacity
        */
        static size_type max_size() {
            return ((size_type)1)<<(sizeof(size_type)*8-6);
        }

        //! The number of bits in the int_vector.
        /*!  \sa size, max_size, bit_size, capacity
         */
        size_type bit_size() const {
            return m_size;
        }

        //! Returns the size of the occupied bits of the int_vector.
        /*! The capacity of a int_vector is greater or equal to the
            bit_size of the vector: capacity() >= bit_size().
            \sa size, bit_size, max_size, capacity
         */
        size_type capacity() const {
            return ((m_size+63)>>6)<<6;
        }

        //! Pointer to the raw data of the int_vector
        /*! \returns Const pointer to the raw data of the int_vector
         */
        const uint64_t* data() const {
            return m_data;
        }

        //! Get the integer value of the binary string of length len starting at position idx in the int_vector.
        /*! \param idx Starting index of the binary representation of the integer.
            \param len Length of the binary representation of the integer. Default value is 64.
            \returns The integer value of the binary string of length len starting at position idx.
            \sa setInt, getBit, setBit
        */
        const value_type get_int(size_type idx, const uint8_t len=64) const;

        //! Set the bits from position idx to idx+len-1 to the binary representation of integer x.
        /*! The bit at position idx represents the least significant bit(lsb), and the bit at
            position idx+len-1 the most significant bit (msb) of x.
            \param idx Starting index of the binary representation of x.
            \param x   The integer to store in the int_vector.
            \param len The length used to store x in the int_vector. Default value is 64.
            \sa getInt, getBit, setBit
        */
        void set_int(size_type idx, value_type x, const uint8_t len=64);

        //! Returns the width of the integers which are accessed via the [] operator.
        /*! \returns The width of the integers which are accessed via the [] operator.
            \sa width
        */
        const uint8_t width() const {
            return m_width;
        }

        //! Sets the width of the integers which are accessed via the [] operator, if t_width equals 0.
        /*! \param intWidth New width of the integers accessed via the [] operator.
            \note This method has no effect if t_width is in the range [1..64].
              \sa width
        */
        void width(uint8_t f_width) {
            int_vector_trait<t_width>::width(m_width, f_width); // delegate to trait function
        }

        //! Serializes the int_vector to a stream.
        /*! \return The number of bytes written to out.
         *  \sa load
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr,
                            std::string name = "", bool write_fixed_as_variable=false) const;

        //! Load the int_vector for a stream.
        void load(std::istream& in);

        //! non const version of [] operator
        /*! \param i Index the i-th integer of length width().
         *  \return A reference to the i-th integer of length width().
         */
        inline reference operator[](const size_type& i);

        //! const version of [] operator
        /*! \param i Index the i-th integer of length width().
         *  \return The value of the i-th integer of length width().
         */
        inline const_reference operator[](const size_type& i) const;

        //! Assignment operator.
        /*! \param v The vector v which should be assigned
         */
        int_vector& operator=(const int_vector& v);

        //! Move assignment operator.
        int_vector& operator=(int_vector&& v);

        //! Equality operator for two int_vectors.
        /*! Two int_vectors are equal if
         *    - capacities and sizes are equal and
         *    - width are equal and
         *    - the bits in the range [0..bit_size()-1] are equal.
         */
        bool operator==(const int_vector& v) const;

        //! Inequality operator for two int_vectors.
        /*! Two int_vectors are not equal if
         *    - capacities and sizes are not equal or
         *    - int widths are not equal or
         *    - the bits in the range [0..bit_size()-1] are not equal.
         */
        bool operator!=(const int_vector& v) const;

        //! Less operator for two int_vectors
        /*! int_vector w is less than v if
         *    - w[i]==v[i] for i<j and w[j]<v[j] with j in [0, min(w.size(), v.size()) )
         *    - or w[i]==v[i] for all i < min(w.size(), v.size()) and w.size()<v.size().
         *  \sa operator>
        */
        bool operator<(const int_vector& v) const;

        //! Greater operator for two int_vectors
        /*! int_vector w is greater than v if
         *    - w[i]==v[i] for i<j and w[j]>v[j] with j in [0, min(w.size(), v.size()) )
         *    - or w[i]==v[i] for all i < min(w.size(), v.size()) and w.size()>v.size().
        */
        bool operator>(const int_vector& v) const;

        //! Less or equal operator
        bool operator<=(const int_vector& v) const;

        //! Greater of equal operator
        bool operator>=(const int_vector& v) const;

        //! Iterator that points to the first element of the int_vector.
        /*!  Time complexity guaranty is O(1).
         */
        const iterator begin() {
            return int_vector_trait<t_width>::begin(this, m_data);
        }

        //! Iterator that points to the element after the last element of int_vector.
        /*! Time complexity guaranty is O(1).
         */
        const iterator end() {
            return int_vector_trait<t_width>::end(this, m_data, (m_size/m_width));
        }

        //! Const iterator that points to the first element of the int_vector.
        const const_iterator begin() const {
            return int_vector_trait<t_width>::begin(this, m_data);
        }

        //! Const iterator that points to the element after the last element of int_vector.
        const const_iterator end() const {
            return int_vector_trait<t_width>::end(this, m_data, (m_size/m_width));
        }
};

//! A proxy class that acts as a reference to an integer of length \p len bits in a int_vector.
/*! \tparam t_int_vector The specific int_vector class.
 */
template<class t_int_vector>
class int_vector_reference
{
    public:
        typedef typename t_int_vector::value_type value_type;
    private:
        typename t_int_vector::value_type* const m_word;
        const uint8_t m_offset;
        const uint8_t m_len; //!< Length of the integer referred to in bits.
    public:
        //! Constructor for the reference class
        /*! \param word Pointer to the corresponding 64bit word in the int_vector.
            \param offset Offset to the starting bit (offset in [0..63])
            \param len length of the integer, should be v->width()!!!
        */
        int_vector_reference(value_type* word, uint8_t offset, uint8_t len):
            m_word(word),m_offset(offset),m_len(len) {};

        //! Assignment operator for the proxy class
        /*!
            The integer x is assign to the referenced
            position in the t_int_vector with the specified width
            of the int_vector
            \param x 64bit integer to assign
            \return A const_reference to the assigned reference
         */
        int_vector_reference& operator=(value_type x) {
            bits::write_int(m_word, x, m_offset, m_len);
            return *this;
        };

        int_vector_reference& operator=(const int_vector_reference& x) {
            return *this = value_type(x);
        };

        //! Cast the reference to a int_vector<>::value_type
        operator value_type()const {
            return bits::read_int(m_word, m_offset, m_len);
        }

        //! Prefix increment of the proxy object
        int_vector_reference& operator++() {
            value_type x = bits::read_int(m_word, m_offset, m_len);
            bits::write_int(m_word, x+1, m_offset, m_len);
            return *this;
        }

        //! Postfix increment of the proxy object
        value_type operator++(int) {
            value_type val = (typename t_int_vector::value_type)*this;
            ++(*this);
            return val;
        }

        //! Prefix decrement of the proxy object
        int_vector_reference& operator--() {
            value_type x = bits::read_int(m_word, m_offset, m_len);
            bits::write_int(m_word, x-1, m_offset, m_len);
            return *this;
        }

        //! Postfix decrement of the proxy object
        value_type operator--(int) {
            value_type val = (value_type)*this;
            --(*this);
            return val;
        }

        //! Add assign from the proxy object
        int_vector_reference& operator+=(const value_type x) {
            value_type w = bits::read_int(m_word, m_offset, m_len);
            bits::write_int(m_word, w+x, m_offset, m_len);
            return *this;
        }

        //! Subtract assign from the proxy object
        int_vector_reference& operator-=(const value_type x) {
            value_type w = bits::read_int(m_word, m_offset, m_len);
            bits::write_int(m_word, w-x, m_offset, m_len);
            return *this;
        }

        bool operator==(const int_vector_reference& x)const {
            return value_type(*this) == value_type(x);
        }

        bool operator<(const int_vector_reference& x)const {
            return value_type(*this) < value_type(x);
        }
};

// For C++11
template<class t_int_vector>
inline void swap(int_vector_reference<t_int_vector> x,
                 int_vector_reference<t_int_vector> y)
{
    // TODO: more efficient solution?
    typename int_vector_reference<t_int_vector>::value_type tmp = x;
    x = y;
    y = tmp;
}

// For C++11
template<class t_int_vector>
inline void swap(typename int_vector_reference<t_int_vector>::value_type& x,
                 int_vector_reference<t_int_vector> y)
{
    // TODO: more efficient solution?
    typename int_vector_reference<t_int_vector>::value_type tmp = x;
    x = y;
    y = tmp;
}

// For C++11
template<class t_int_vector>
inline void swap(int_vector_reference<t_int_vector> x,
                 typename int_vector_reference<t_int_vector>::value_type& y)
{
    // TODO: more efficient solution?
    typename int_vector_reference<t_int_vector>::value_type tmp = x;
    x = y;
    y = tmp;
}

// specialization for int_vector_reference for int_vector == bit_vector
// special thanks to Timo Beller, who pointed out that the specialization is missing
// Same implementation as in stl_bvector.h.
template<>
class int_vector_reference<bit_vector>
{
    public:
        typedef bool value_type;
    private:
        uint64_t* const m_word;
        uint64_t m_mask;
    public:
        //! Constructor for the reference class
        /*! \param word Pointer to the corresponding 64bit word in the int_vector.
            \param offset Offset to the starting bit (offset in [0..63])
        */
        int_vector_reference(uint64_t* word, uint8_t offset, uint8_t):
            m_word(word),m_mask(1ULL<<offset) {};

        //! Assignment operator for the proxy class
        int_vector_reference& operator=(bool x) {
            if (x)
                *m_word |= m_mask;
            else
                *m_word &= ~m_mask;
            return *this;
        };

        int_vector_reference& operator=(const int_vector_reference& x) {
            return *this = bool(x);
        };

        //! Cast the reference to a bool
        operator bool()const {
            return !!(*m_word & m_mask);
        }

        bool operator==(const int_vector_reference& x)const {
            return bool(*this) == bool(x);
        }

        bool operator<(const int_vector_reference& x)const {
            return !bool(*this) && bool(x);
        }
};

// For C++11
template<>
inline void swap(int_vector_reference<bit_vector> x,
                 int_vector_reference<bit_vector> y)
{
    // TODO: more efficient solution?
    bool tmp = x;
    x = y;
    y = tmp;
}

// For C++11
template<>
inline void swap(bool& x,
                 int_vector_reference<bit_vector> y)
{
    // TODO: more efficient solution?
    bool tmp = x;
    x = y;
    y = tmp;
}

// For C++11
template<>
inline void swap(int_vector_reference<bit_vector> x,
                 bool& y)
{
    // TODO: more efficient solution?
    bool tmp = x;
    x = y;
    y = tmp;
}



template<class t_int_vector>
class int_vector_iterator_base: public std::iterator<std::random_access_iterator_tag, typename t_int_vector::value_type, typename t_int_vector::difference_type>
{
    public:
        typedef uint64_t  size_type;
    protected:
        uint8_t           m_offset;
        uint8_t           m_len;

    public:
        int_vector_iterator_base(uint8_t offset, uint8_t len):
            m_offset(offset),m_len(len) {}

        int_vector_iterator_base(const t_int_vector* v=nullptr, size_type idx=0):
            m_offset(idx&0x3F), m_len(v==nullptr ? 0 : v->m_width) {}
};

template<class t_int_vector>
class int_vector_iterator : public int_vector_iterator_base<t_int_vector>
{
    public:

        typedef int_vector_reference<t_int_vector>     reference;
        typedef uint64_t                               value_type;
        typedef int_vector_iterator                    iterator;
        typedef reference*                             pointer;
        typedef typename t_int_vector::size_type       size_type;
        typedef typename t_int_vector::difference_type difference_type;

    private:

        using int_vector_iterator_base<t_int_vector>::m_offset; // make m_offset easy usable
        using int_vector_iterator_base<t_int_vector>::m_len;    // make m_len easy usable

        typename t_int_vector::value_type* m_word;

    public:

        int_vector_iterator(t_int_vector* v=nullptr, size_type idx=0):
            int_vector_iterator_base<t_int_vector>(v, idx),
            m_word((v != nullptr) ? v->m_data + (idx>>6) : nullptr) {}


        int_vector_iterator(const int_vector_iterator<t_int_vector>& it) :
            m_word(it.m_word) {
            m_offset = it.m_offset;
            m_len = it.m_len;
        }

        reference operator*() const {
            return reference(m_word, m_offset, m_len);
        }

        //! Prefix increment of the Iterator
        iterator& operator++() {
            m_offset+=m_len;
            if (m_offset >= 64) {
                m_offset &= 0x3F;
                ++m_word;
            }
            return *this;
        }

        //! Postfix increment of the Iterator
        iterator operator++(int) {
            int_vector_iterator it = *this;
            ++(*this);
            return it;
        }

        //! Prefix decrement of the Iterator
        iterator& operator--() {
            m_offset-=m_len;
            if (m_offset >= 64) {
                m_offset &= 0x3F;
                --m_word;
            }
            return *this;
        }

        //! Postfix decrement of the Iterator
        iterator operator--(int) {
            int_vector_iterator it = *this;
            --(*this);
            return it;
        }

        iterator& operator+=(difference_type i) {
            if (i<0)
                return *this -= (-i);
            difference_type t = i*m_len;
            m_word += (t>>6);
            if ((m_offset+=(t&0x3F))&~0x3F) {  // if new offset is >= 64
                ++m_word;       // add one to the word
                m_offset&=0x3F; // offset = offset mod 64
            }
            return *this;
        }

        iterator& operator-=(difference_type i) {
            if (i<0)
                return *this += (-i);
            difference_type t = i*m_len;
            m_word -= (t>>6);
            if ((m_offset-=(t&0x3F))&~0x3F) {  // if new offset is < 0
                --m_word;
                m_offset&=0x3F;
            }
            return *this;
        }

        iterator& operator=(const int_vector_iterator<t_int_vector>& it) {
            if (this != &it) {
                m_word   = it.m_word;
                m_offset = it.m_offset;
                m_len    = it.m_len;
            }
            return *this;
        }

        iterator operator+(difference_type i) const {
            iterator it = *this;
            return it += i;
        }

        iterator operator-(difference_type i) const {
            iterator it = *this;
            return it -= i;
        }

        reference operator[](difference_type i) const {
            return *(*this + i);
        }

        bool operator==(const int_vector_iterator& it)const {
            return it.m_word == m_word && it.m_offset == m_offset;
        }

        bool operator!=(const int_vector_iterator& it)const {
            return !(*this==it);
        }

        bool operator<(const int_vector_iterator& it)const {
            if (m_word == it.m_word)
                return m_offset < it.m_offset;
            return m_word < it.m_word;
        }

        bool operator>(const int_vector_iterator& it)const {
            if (m_word == it.m_word)
                return m_offset > it.m_offset;
            return m_word > it.m_word;
        }

        bool operator>=(const int_vector_iterator& it)const {
            return !(*this < it);
        }

        bool operator<=(const int_vector_iterator& it)const {
            return !(*this > it);
        }
        inline difference_type operator-(const int_vector_iterator& it) {
            return (((m_word - it.m_word)<<6) + m_offset - it.m_offset) / m_len;
        }
};

//template<class t_int_vector>
//void swap(const int_vector_iterator<t_int_vector> &x, const int_vector_iterator<t_int_vector> &y){
//  x->swap(*y);
//}

template<class t_int_vector>
inline int_vector_iterator<t_int_vector> operator+(typename int_vector_iterator<t_int_vector>::difference_type n, const int_vector_iterator<t_int_vector>& it)
{
    return it+n;
}

template<class t_int_vector>
class int_vector_const_iterator : public int_vector_iterator_base<t_int_vector>
{
    public:

        typedef typename t_int_vector::value_type        const_reference;
        typedef const typename t_int_vector::value_type* pointer;
        typedef int_vector_const_iterator                const_iterator;
        typedef typename t_int_vector::size_type         size_type;
        typedef typename t_int_vector::difference_type   difference_type;

        template<class X>
        friend typename int_vector_const_iterator<X>::difference_type
        operator-(const int_vector_const_iterator<X>& x, const int_vector_const_iterator<X>& y);
        friend class int_vector_iterator<t_int_vector>;
        friend class int_vector_iterator_base<t_int_vector>;

    private:

        using int_vector_iterator_base<t_int_vector>::m_offset; // make m_offset easy usable
        using int_vector_iterator_base<t_int_vector>::m_len;    // make m_len easy usable

        const typename t_int_vector::value_type* m_word;

    public:

        int_vector_const_iterator(const t_int_vector* v=nullptr, size_type idx=0):
            int_vector_iterator_base<t_int_vector>(v, idx),
            m_word((v != nullptr) ? v->m_data + (idx>>6) : nullptr) {}

        int_vector_const_iterator(const int_vector_const_iterator& it):
            m_word(it.m_word) {
            m_offset = it.m_offset;
            m_len = it.m_len;
        }

        int_vector_const_iterator(const int_vector_iterator<t_int_vector>& it):
            m_word(it.m_word) {
            m_offset = it.m_offset;
            m_len = it.m_len;
        }

        const_reference operator*() const {
            if (m_offset+m_len <= 64) {
                return ((*m_word)>>m_offset)&bits::lo_set[m_len];
            }
            return ((*m_word)>>m_offset) |
                   ((*(m_word+1) & bits::lo_set[(m_offset+m_len)&0x3F])<<(64-m_offset));
        }

        //! Prefix increment of the Iterator
        const_iterator& operator++() {
            m_offset+=m_len;
            if (m_offset >= 64) {
                m_offset &= 0x3F;
                ++m_word;
            }
            return *this;
        }

        //! Postfix increment of the Iterator
        const_iterator operator++(int) {
            int_vector_const_iterator it = *this;
            ++(*this);
            return it;
        }

        //! Prefix decrement of the Iterator
        const_iterator& operator--() {
            m_offset-=m_len;
            if (m_offset >= 64) {
                m_offset &= 0x3F;
                --m_word;
            }
            return *this;
        }

        //! Postfix decrement of the Iterator
        const_iterator operator--(int) {
            int_vector_const_iterator it = *this;
            --(*this);
            return it;
        }

        const_iterator& operator+=(difference_type i) {
            if (i<0)
                return *this -= (-i);
            difference_type t = i*m_len;
            m_word += (t>>6);
            if ((m_offset+=(t&0x3F))&~0x3F) {// if new offset >= 64
                ++m_word;       // add one to the word
                m_offset&=0x3F; // offset = offset mod 64
            }
            return *this;
        }

        const_iterator& operator-=(difference_type i) {
            if (i<0)
                return *this += (-i);
            difference_type t = i*m_len;
            m_word -= (t>>6);
            if ((m_offset-=(t&0x3F))&~0x3F) {// if new offset is < 0
                --m_word;
                m_offset&=0x3F;
            }
            return *this;
        }

        const_iterator operator+(difference_type i) const {
            const_iterator it = *this;
            return it += i;
        }

        const_iterator operator-(difference_type i) const {
            const_iterator it = *this;
            return it -= i;
        }

        const_reference operator[](difference_type i) const {
            return *(*this + i);
        }

        bool operator==(const int_vector_const_iterator& it)const {
            return it.m_word == m_word && it.m_offset == m_offset;
        }

        bool operator!=(const int_vector_const_iterator& it)const {
            return !(*this==it);
        }

        bool operator<(const int_vector_const_iterator& it)const {
            if (m_word == it.m_word)
                return m_offset < it.m_offset;
            return m_word < it.m_word;
        }

        bool operator>(const int_vector_const_iterator& it)const {
            if (m_word == it.m_word)
                return m_offset > it.m_offset;
            return m_word > it.m_word;
        }

        bool operator>=(const int_vector_const_iterator& it)const {
            return !(*this < it);
        }

        bool operator<=(const int_vector_const_iterator& it)const {
            return !(*this > it);
        }

};

template<class t_int_vector>
inline typename int_vector_const_iterator<t_int_vector>::difference_type
operator-(const int_vector_const_iterator<t_int_vector>& x,
          const int_vector_const_iterator<t_int_vector>& y)
{
    return (((x.m_word - y.m_word)<<6) + x.m_offset - y.m_offset) / x.m_len;
}

template<class t_int_vector>
inline int_vector_const_iterator<t_int_vector>
operator+(typename int_vector_const_iterator<t_int_vector>::difference_type n,
          const int_vector_const_iterator<t_int_vector>& it)
{
    return it + n;
}

inline std::ostream& operator<<(std::ostream& os, const int_vector<1>& v)
{
    for (auto it=v.begin(), end = v.end(); it != end; ++it) {
        os << *it;
    }
    return os;
}

template<uint8_t t_width>
inline std::ostream& operator<<(std::ostream& os, const int_vector<t_width>& v)
{
    for (auto it=v.begin(), end = v.end(); it != end; ++it) {
        os << *it;
        if (it+1 != end) os << " ";
    }
    return os;
}


// ==== int_vector implementation  ====

template<uint8_t t_width>
inline int_vector<t_width>::int_vector(size_type size, value_type default_value, uint8_t intWidth):
    m_size(0), m_data(nullptr), m_width(t_width)
{
    mm::add(this);
    int_vector_trait<t_width>::width(m_width, intWidth);
    resize(size);
    util::set_to_value(*this, default_value); // new initialization
}

template<uint8_t t_width>
inline int_vector<t_width>::int_vector(int_vector&& v) :
    m_size(v.m_size), m_data(v.m_data), m_width(v.m_width)
{
    v.m_size = 0;       // has to be set for mm::remove
    v.m_data = nullptr; // ownership of v.m_data now transfered
    mm::add(this, true);
}

template<uint8_t t_width>
inline int_vector<t_width>::int_vector(const int_vector& v):
    m_size(0), m_data(nullptr), m_width(v.m_width)
{
    mm::add(this);
    bit_resize(v.bit_size());
    if (v.capacity() > 0) {
        if (memcpy(m_data, v.data() ,v.capacity()/8)==nullptr) {
            throw std::bad_alloc();
        }
    }
    int_vector_trait<t_width>::width(m_width, v.m_width);
}

template<uint8_t t_width>
int_vector<t_width>& int_vector<t_width>::operator=(const int_vector& v)
{
    if (this != &v) {// if v is not the same object
        bit_resize(v.bit_size());
        if (v.bit_size()>0) {
            if (memcpy(m_data, v.data() ,v.capacity()/8)==nullptr) {
                throw std::bad_alloc();
            }
        }
        width(v.width());
    }
    return *this;
}

template<uint8_t t_width>
int_vector<t_width>& int_vector<t_width>::operator=(int_vector&& v)
{
    swap(v);
    return *this;
}

// Destructor
template<uint8_t t_width>
int_vector<t_width>::~int_vector()
{
    mm::remove(this);
    free(m_data);
}

template<uint8_t t_width>
void int_vector<t_width>::swap(int_vector& v)
{
    if (this != &v) { // if v and _this_ are not the same object
        size_type size     = m_size;
        uint64_t* data     = m_data;
        uint8_t  int_width = m_width;
        m_size   = v.m_size;
        m_data   = v.m_data;
        int_vector_trait<t_width>::width(m_width, v.m_width);
        v.m_size = size;
        v.m_data = data;
        int_vector_trait<t_width>::width(v.m_width, int_width);
    }
}

template<uint8_t t_width>
void int_vector<t_width>::bit_resize(const size_type size)
{
    mm::realloc(*this, size);
}

template<uint8_t t_width>
auto int_vector<t_width>::get_int(size_type idx, const uint8_t len)const -> const value_type
{
#ifdef SDSL_DEBUG
    if (idx+len > m_size) {
        throw std::out_of_range("OUT_OF_RANGE_ERROR: int_vector::get_int(size_type, uint8_t); idx+len > size()!");
    }
    if (len > 64) {
        throw std::out_of_range("OUT_OF_RANGE_ERROR: int_vector::get_int(size_type, uint8_t); len>64!");
    }
#endif
    return bits::read_int(m_data+(idx>>6), idx&0x3F, len);
}

template<uint8_t t_width>
inline void int_vector<t_width>::set_int(size_type idx, value_type x, const uint8_t len)
{
#ifdef SDSL_DEBUG
    if (idx+len > m_size) {
        throw std::out_of_range("OUT_OF_RANGE_ERROR: int_vector::set_int(size_type, uint8_t); idx+len > size()!");
    }
    if (len > 64) {
        throw std::out_of_range("OUT_OF_RANGE_ERROR: int_vector::set_int(size_type, uint8_t); len>64!");
    }
#endif
    bits::write_int(m_data+(idx>>6), x, idx&0x3F, len);
}

template<uint8_t t_width>
inline auto int_vector<t_width>::operator[](const size_type& idx) -> reference {
    assert(idx < this->size());
    size_type i = idx * m_width;
    return reference(this->m_data + (i>>6), i&0x3F, m_width);
}

// specialized [] operator for 64 bit access.
template<>
inline auto int_vector<64>::operator[](const size_type& idx) -> reference {
    assert(idx < this->size());
    return *(this->m_data+idx);
}

// specialized [] operator for 32 bit access.
template<>
inline auto int_vector<32>::operator[](const size_type& idx) -> reference {
    assert(idx < this->size());
    return *(((uint32_t*)(this->m_data))+idx);
}

// specialized [] operator for 16 bit access.
template<>
inline auto int_vector<16>::operator[](const size_type& idx) -> reference {
    assert(idx < this->size());
    return *(((uint16_t*)(this->m_data))+idx);
}

// specialized [] operator for 8 bit access.
template<>
inline auto int_vector<8>::operator[](const size_type& idx) -> reference {
    assert(idx < this->size());
    return *(((uint8_t*)(this->m_data))+idx);
}

template<uint8_t t_width>
inline auto
int_vector<t_width>::operator[](const size_type& idx)const -> const_reference
{
    assert(idx < this->size());
    return get_int(idx * t_width, t_width);
}

template<>
inline auto
int_vector<0>::operator[](const size_type& idx)const -> const_reference
{
    assert(idx < this->size());
    return get_int(idx * m_width, m_width);
}

template<>
inline auto
int_vector<64>::operator[](const size_type& idx)const -> const_reference
{
    assert(idx < this->size());
    return *(this->m_data+idx);
}

template<>
inline auto
int_vector<32>::operator[](const size_type& idx)const -> const_reference
{
    assert(idx < this->size());
    return *(((uint32_t*)this->m_data)+idx);
}

template<>
inline auto
int_vector<16>::operator[](const size_type& idx)const -> const_reference
{
    assert(idx < this->size());
    return *(((uint16_t*)this->m_data)+idx);
}

template<>
inline auto
int_vector<8>::operator[](const size_type& idx)const -> const_reference
{
    assert(idx < this->size());
    return *(((uint8_t*)this->m_data)+idx);
}

template<>
inline auto
int_vector<1>::operator[](const size_type& idx)const -> const_reference
{
    assert(idx < this->size());
    return ((*(m_data+(idx>>6)))>>(idx&0x3F))&1;
}
template<uint8_t t_width>
bool int_vector<t_width>::operator==(const int_vector& v)const
{
    if (capacity() != v.capacity())
        return false;
    if (bit_size() != v.bit_size())
        return false;
    if (empty())
        return true;
    const uint64_t* data1 = v.data();
    const uint64_t* data2 = data();
    for (size_type i=0; i < (capacity()>>6)-1; ++i) {
        if (*(data1++) != *(data2++))
            return false;
    }
    int8_t l = 64-(capacity()-bit_size());
    return ((*data1)&bits::lo_set[l])==((*data2)&bits::lo_set[l]);
}

template<uint8_t t_width>
bool int_vector<t_width>::operator<(const int_vector& v)const
{
    size_type min_size = size();
    if (min_size > v.size())
        min_size = v.size();
    for (auto it = begin(), end = begin()+min_size, it_v = v.begin(); it!=end; ++it, ++it_v) {
        if (*it == *it_v)
            continue;
        else
            return *it < *it_v;
    }
    return  size() < v.size();
}

template<uint8_t t_width>
bool int_vector<t_width>::operator>(const int_vector& v)const
{
    size_type min_size = size();
    if (min_size > v.size())
        min_size = v.size();
    for (auto it = begin(), end = begin()+min_size, it_v = v.begin(); it!=end; ++it, ++it_v) {
        if (*it == *it_v)
            continue;
        else
            return *it > *it_v;
    }
    return  size() > v.size();
}

template<uint8_t t_width>
bool int_vector<t_width>::operator<=(const int_vector& v)const
{
    return *this==v or *this<v;
}

template<uint8_t t_width>
bool int_vector<t_width>::operator>=(const int_vector& v)const
{
    return *this==v or *this>v;
}

template<uint8_t t_width>
bool int_vector<t_width>::operator!=(const int_vector& v)const
{
    return !(*this==v);
}

template<class size_type_class>
size_type_class _sdsl_serialize_size_and_int_width(std::ostream& out, uint8_t fixed_int_width, uint8_t int_width, size_type_class size)
{
    size_type_class written_bytes = 0;
    out.write((char*) &size, sizeof(size));
    written_bytes += sizeof(size);
    if (fixed_int_width == 0) {
        out.write((char*) &int_width, sizeof(int_width));
        written_bytes += sizeof(int_width);
    }
    return written_bytes;
}

template<uint8_t t_width>
typename int_vector<t_width>::size_type int_vector<t_width>::serialize(std::ostream& out,
        structure_tree_node* v,
        std::string name,
        bool write_fixed_as_variable) const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    if (t_width > 0 and write_fixed_as_variable) {
        written_bytes += _sdsl_serialize_size_and_int_width(out, 0, t_width, m_size);
    } else {
        written_bytes += _sdsl_serialize_size_and_int_width(out, t_width, m_width, m_size);
    }
    uint64_t* p = m_data;
    size_type idx = 0;
    while (idx+constants::SDSL_BLOCK_SIZE < (capacity()>>6)) {
        out.write((char*) p, constants::SDSL_BLOCK_SIZE*sizeof(uint64_t));
        written_bytes += constants::SDSL_BLOCK_SIZE*sizeof(uint64_t);
        p     += constants::SDSL_BLOCK_SIZE;
        idx    += constants::SDSL_BLOCK_SIZE;
    }
    out.write((char*) p, ((capacity()>>6)-idx)*sizeof(uint64_t));
    written_bytes += ((capacity()>>6)-idx)*sizeof(uint64_t);
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<uint8_t t_width>
void int_vector<t_width>::load(std::istream& in)
{
    size_type size;
    int_vector_trait<t_width>::read_header(size, m_width, in);

    bit_resize(size);
    uint64_t* p = m_data;
    size_type idx = 0;
    while (idx+constants::SDSL_BLOCK_SIZE < (capacity()>>6)) {
        in.read((char*) p, constants::SDSL_BLOCK_SIZE*sizeof(uint64_t));
        p     += constants::SDSL_BLOCK_SIZE;
        idx += constants::SDSL_BLOCK_SIZE;
    }
    in.read((char*) p, ((capacity()>>6)-idx)*sizeof(uint64_t));
}

//! A wrapper class which allows us to serialize an char array as an int_vector.
class char_array_serialize_wrapper
{
    public:
        typedef int_vector_size_type    size_type;
    private:
        size_type m_n;  // number of char
        const unsigned char* m_cp;
    public:
        inline char_array_serialize_wrapper(const unsigned char* c=nullptr, size_type n=0):m_n(n), m_cp(c) {}
        size_type serialize(std::ostream& out)const;
};

//! A class for reading an int_vector buffered from a file.
template<uint8_t t_width>
class int_vector_file_buffer
{
    public:

        typedef typename int_vector<t_width>::size_type        size_type;
        typedef typename int_vector<t_width>::value_type       value_type;
        typedef typename int_vector<t_width>::const_reference  const_reference;
        typedef typename int_vector<t_width>::int_width_type   int_width_type;

    private:

        isfstream      m_in;
        uint64_t*      m_buf = nullptr;
        size_type      m_off = 0; // offset in the first 64bit word of the buffer
        size_type      m_read_values = 0; // number of values read in the last buffer operation
        size_type      m_len = 0;
        size_type      m_int_vector_size = 0;
        size_type      m_read_values_sum = 0;
        int_width_type m_width = 0;
        std::string    m_file = "";
        bool           m_load_from_plain = false;

        void load_size_and_width() {
            int_vector_trait<t_width>::read_header(m_int_vector_size, m_width, m_in);
            m_int_vector_size/=m_width;
        }

        size_type words_to_read(size_type len) {
            if (len == 0)
                return 0;
            size_type last_off = (m_read_values_sum*m_width)%64;
            assert(((64-last_off)%64) < m_width);
            return ((len*m_width - ((64-last_off)%64))+63)/64;
        }

        void init() {
            m_int_vector_size = 0;
            m_off             = 0;
            m_read_values     = 0;
            m_read_values_sum = 0;
        }

    public:

        const size_type& int_vector_size  = m_int_vector_size;
        const uint8_t& width              = m_width;
        const std::string& file_name      = m_file;

        //! Constructor
        /*
         * \param f_file     File which contains the int_vector.
         * \param len             Length of the buffer in elements.
         */
        int_vector_file_buffer(std::string f_file="", size_type len=1000000,
                               uint8_t int_width=0): m_in(), m_width(t_width) {
            m_load_from_plain = false;
            int_vector_trait<t_width>::width(m_width, int_width);
            m_len             = len;
            init();
            if (f_file == "") {
                return;
            }
            m_file = f_file;
            m_in.open(m_file);
            if (m_in.is_open()) {
                load_size_and_width();
                m_buf = new uint64_t[(m_len*m_width+63)/64 + 2];
            } else {
                m_buf = nullptr;
                std::cerr<< "WARNING: int_vector_file_buffer:\
                     could not open `"<<m_file<<std::endl;
            }
        }

        // initialize int_vector_file_buffer from a plain file
        // works only for t_width = 8 // TODO extent to 1,16,32,64
        bool load_from_plain(const std::string& f_file, size_type len=1000000,
                             uint8_t int_width=0) {
            if (t_width != 8) {
                std::logic_error("ERROR: int_vector_file_buffer: load_from_plain\
                          is only implemented for t_width=8.");
                return false;
            }
            m_load_from_plain = true;
            int_vector_trait<t_width>::width(m_width, int_width);
            m_len = len;
            m_file = f_file;
            m_int_vector_size = util::file_size(m_file);
            m_in.open(m_file.c_str());
            if (m_in.is_open()) {
                m_buf = new uint64_t[(m_len*m_width+63)/64 + 2];
            } else {
                m_buf = nullptr;
                return false;
            }
            return true;
        }

        bool reset(size_type new_buf_len=0) {
            m_in.clear();
            if (!m_in.seekg(0, std::ios::beg)) {
                std::cout<<"m_file="<<m_file<<std::endl;
                throw std::ios_base::failure("int_vector_file_buffer: reset()");
                return false;
            };
            init();
            if (!m_load_from_plain) {
                load_size_and_width();
            } else {
                m_int_vector_size = util::file_size(m_file);
            }
            if (new_buf_len > 0 and new_buf_len != m_len) {
                if (m_buf != nullptr)
                    delete [] m_buf;
                m_len = new_buf_len;
                m_buf = new uint64_t[(m_len*m_width+63)/64 + 2];
            }
            return true;
        }

        bool valid() {
            return m_buf!=nullptr;
        }

        size_type load_next_block() {
            if (m_read_values_sum == m_int_vector_size) {
                return 0;
            }
            assert(m_read_values_sum < m_int_vector_size);
            size_type values_to_read     = m_len;
            if (values_to_read + m_read_values_sum > m_int_vector_size) {
                values_to_read = m_int_vector_size - m_read_values_sum;
            }
            if (((m_read_values_sum*m_width)%64) == 0) { //if the new offset == 0
                m_in.read((char*) m_buf, words_to_read(values_to_read) * sizeof(uint64_t));
                m_off = 0;
            } else { // the new offset != 0
                m_buf[0] = m_buf[(m_read_values*m_width+m_off)/64 ];
                m_in.read((char*)(m_buf+1), words_to_read(values_to_read) * sizeof(uint64_t));
                m_off = (m_read_values_sum*m_width)%64;
            }
            m_read_values_sum     += values_to_read;
            m_read_values         =  values_to_read;
            return m_read_values;
        }

        //! Returns the i-th element in the buffer
        /* \pre \f$ i < len \f$
         */
        value_type operator[](const size_type i)const {
            assert(i<m_len);
            size_type idx = i*m_width+m_off;
            return bits::read_int(m_buf + (idx>>6), idx&0x3F, m_width);
        }

        void set_int(const size_type i, uint64_t x) {
            size_type idx = i*m_width+m_off;
            bits::write_int(m_buf + (idx>>6), x, idx&0x3F, m_width);
        }

        const uint64_t* data()const {
            return m_buf;
        }

        ~int_vector_file_buffer() {
            if (m_in.is_open()) {
                m_in.close(); // close stream
                delete [] m_buf;
            }
        }
};

}// end namespace sdsl
#endif
