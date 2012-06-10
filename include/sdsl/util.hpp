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
/*! \file util.hpp
    \brief util.hpp contains some helper methods for int_vector and other stuff like demangle class names.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_UTIL
#define INCLUDED_SDSL_UTIL

#include "bitmagic.hpp"
#include "typedefs.hpp"
#include <iosfwd> // forward declaration of ostream
#include <stdint.h> // for uint64_t uint32_t declaration
#include <iostream>// for cerr
#include <cassert>
#include <fstream> // file stream for storeToFile and loadFromFile
#include <ctime>  // for rand initialization
#include <string>
#include <string.h> // for strlen and strdup
#include <libgen.h> // for basename
#include <cstdlib>
#include <unistd.h> // for getpid 
#include <sstream> // for to_string method
#include <stdexcept>   // for std::logic_error


// macros to transform a defined name to a string
#define SDSL_STR(x) #x
#define SDSL_XSTR(s) SDSL_STR(s)

//! Namespace for the succinct data structure library.
namespace sdsl
{

template<uint8_t, class size_type_class>
class int_vector_file_buffer; // forward declaration

template<uint8_t, class size_type_class>
class int_vector;	 // forward declaration

//! A namespace for helper functions
namespace util
{

static bool verbose = false;

//! Returns the basename of a file_name
std::string basename(const std::string& file_name);

//! Returns the directory of a file_name. Trailing / are removed.
std::string dirname(const std::string& file_name);

//! Sets all bits of the int_vector to pseudo-random bits.
/*! \param v The int_vector whose bits should be set to random bits
 *  \param seed If seed = 0, the time is used to initialize the
 *              pseudo random number generator, otherwise the seed
 *              parameter is used.
 */
template<class int_vector_type>
void set_random_bits(int_vector_type& v, int seed=0);
//! Sets all bits of the int_vector to 0-bits.
template<class int_vector_type>
void set_zero_bits(int_vector_type& v);
//! Sets all bits of the int_vector to 1-bits.
template<class int_vector_type>
void set_one_bits(int_vector_type& v);

//! Bit compress the int_vector
/*! Determine the biggest value X and then set the
 *  int_width to the smallest possible so that we
 *  still can represent X
 */
template<class int_vector_type>
void bit_compress(int_vector_type& v);

//! All elements of v modulo m
template<class int_vector_type, class size_type_class>
void all_elements_mod(int_vector_type& v, size_type_class m);


//! Set all entries of int_vector to value k
/*! \param  v The int_vector which should be set
 *	\param  k The value which should be inserted into v.
 *  \par Details
 *   This method precalculates the content of at most 64
 *   words and then repeatedly inserts these words into v.
 */
template<class int_vector_type>
void set_all_values_to_k(int_vector_type& v, uint64_t k);

//! Sets each entry of the numerical vector v at position \$fi\f$ to value \$fi\$f
template<class int_vector_type>
void set_to_id(int_vector_type& v);

//! Counts and returns the 1-bits an int_vector contains.
/*! \param v The int_vector to count the 1-bits.
  	\return The number of 1-bits in v.
 */
template<class int_vector_type>
typename int_vector_type::size_type get_one_bits(const int_vector_type& v);

//! Counts 10 bit pair occurencies.
/*! \sa getOneBits, getOneZeroBits
 */
template<class int_vector_type>
typename int_vector_type::size_type get_onezero_bits(const int_vector_type& v);

//! Counts 01 bit pair occurencies.
/*! \sa getOneBits, getZeroOneBits
 */
template <class int_vector_type>
typename int_vector_type::size_type get_zeroone_bits(const int_vector_type& v);

//! Load a data structure from a file.
/*! The data structure has to provide a load function.
 * \param v Data structure to load.
   \param file_name Name of the serialized file.
 */
template<class T>
bool load_from_file(T& v, const char* file_name);

template<>
bool load_from_file(void*&, const char* file_name);

template<class size_type_class>
bool load_from_int_vector_buffer(unsigned char*& text, int_vector_file_buffer<8, size_type_class>& text_buf);

//! Specialization of load_from_file for a char array
/*  \pre v=NULL
 *
 */
bool load_from_file(char*& v, const char* file_name);

//! Store a data structure to a file.
/*! The data structure has to provide a serialize function.
	\param v Data structure to store.
	\param file_name Name of the file where to store the data structure.
	\param Return if the data structure was stored successfully
 */
template<class T>
bool store_to_file(const T& v, const char* file_name);

//! Specialization of store_to_file for a char array
bool store_to_file(const char* v, const char* file_name);

//! Specialization of store_to_file for int_vector
template<uint8_t fixed_int_width, class size_type_class>
bool store_to_file(const int_vector<fixed_int_width, size_type_class>& v, const char* file_name, bool write_fixed_as_variable=false);

//! Demangle the class name of typeid(...).name()
/*!
 *	\param name A pointer to the the result of typeid(...).name()
 */
std::string demangle(const char* name);

//! Demangle the class name of typeid(...).name() and remove the "sdsl::"-prefix, "unsigned int",...
std::string demangle2(const char* name);

template<class T>
std::string get_class_name(const T& t)
{
    std::string result = demangle2(typeid(t).name());
    size_t template_pos = result.find("<");
    if (template_pos != std::string::npos) {
        result = result.erase(template_pos);
    }
    return result;
}

//! Get the size of a data structure in bytes.
/*!
 *	\param v A reference to the data structure for which the size in bytes should be calculated.
 */
template<class T>
typename T::size_type get_size_in_bytes(const T& t);

//! Get the size of a data structure in mega bytes (MB).
/*!
 *	\param t A reference to the data structure for which the size in bytes should be calculated.
 */
template<class T>
double get_size_in_mega_bytes(const T& t);

struct nullstream : std::ostream {
    struct nullbuf: std::streambuf {
        int overflow(int c) {
            return traits_type::not_eof(c);
        }
    } m_sbuf;
    nullstream(): std::ios(&m_sbuf), std::ostream(&m_sbuf) {}
};

// Writes primitive-typed variable t to stream out
template<class T>
size_t write_member(const T& t, std::ostream& out)
{
    out.write((char*)&t, sizeof(t));
    return sizeof(t);
}

// Specialization for std::string
template<>
size_t write_member<std::string>(const std::string& t, std::ostream& out);

// Writes primitive-typed variable t to stream out
template<class T>
void read_member(T& t, std::istream& in)
{
    in.read((char*)&t, sizeof(t));
}

// Specialization for std::string
template<>
void read_member<std::string>(std::string& t, std::istream& in);

//! Get the process id of the current process
uint64_t get_pid();

class _id_helper
{
    private:
        static uint64_t id;
    public:
        static uint64_t getId() {
            return id++;
        }
};


//! Get a unique id inside the process
inline uint64_t get_id()
{
    return _id_helper::getId();
}

//! Convert type to string
template<typename T>
std::string to_string(const T& t);

template<typename T>
std::string to_latex_string(const T& t);

std::string to_latex_string(unsigned char c);

//! Delete all files in the file_map in the file system
void delete_all_files(tMSS& file_map);


// thanks to Stefan Arnold for the assign functions
//! Assigns the value x of type T to the value of y of type U.
/*!
 *  \param x	The assigned variable.
 *	\param y	The variable which provides the value that is assigned to x.
 */
template<class T, class U>
void assign(T& x, const U& y)
{
    x = T(y);
}

//! Swaps variables x and y.
/*!
 * \param x Reference to the first variable.
 * \param y Reference to the second variable.
 */
template<class T>
void assign(T& x, T& y)
{
    x.swap(y);
}

//! clear the space used by x
/*!
 * \param x Reference to the data structure.
 */
template<class T>
void clear(T& x)
{
    T y;
    x.swap(y);
}

//! Swap support data structure and assign to new vector
/*! \param s1 First support structure.
 *  \param s2 Second support structure.
 *  \param p1 First supported structure.
 *  \param p2 Second supported structure.
 *  s1 is swapped with s2 and after the execution s1 supports p1 and s2 supports
 *  p2. I.e. if p1 and p2 are members of a complex data structure, we have to
 *  swap p1 and p2 before we use this method.
 */
template<class S, class P>
void swap_support(S& s1, S& s2, const P* p1, const P* p2)
{
    s1.swap(s2);
    s1.set_vector(p1);
    s2.set_vector(p2);
}

//! Initialise support data structure with
/*! \param s Support structure which should be initialized
 *  \param x Pointer to the data structure which should be supported.
 */
template<class S, class X>
void init_support(S& s, const X* x)
{
    S temp(x);			// generate a temporary support object
    s.swap(temp);		// swap its content with the target object
    s.set_vector(x);    // set the support object's  pointer to x
}



}

//==================== Template functions ====================


template<class T>
typename T::size_type util::get_size_in_bytes(const T& t)
{
    if ((&t) == NULL)
        return 0;
    util::nullstream ns;
    return t.serialize(ns);
}

template<class T>
double util::get_size_in_mega_bytes(const T& t)
{
    return get_size_in_bytes(t)/(1024.0*1024.0);
}

template<class T>
bool util::store_to_file(const T& t, const char* file_name)
{
    std::ofstream out;
    out.open(file_name, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out)
        return false;
    t.serialize(out);
    out.close();
    return true;
}

inline bool util::store_to_file(const char* v, const char* file_name)
{
    std::ofstream out;
    out.open(file_name, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out)
        return false;
    uint64_t n = strlen((const char*)v);
    out.write(v, n);
    out.close();
    return true;
}

template<uint8_t fixed_int_width, class size_type_class>
bool util::store_to_file(const int_vector<fixed_int_width, size_type_class>& v, const char* file_name, bool write_fixed_as_variable)
{
    std::ofstream out;
    out.open(file_name, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out)
        return false;
    v.serialize(out, write_fixed_as_variable);
    out.close();
    return true;
}

template<class T>
bool util::load_from_file(T& v, const char* file_name)
{
    std::ifstream in;
    in.open(file_name, std::ios::binary | std::ios::in);
    if (!in)
        return false;
    v.load(in);
    in.close();
    return true;
}


template<class size_type_class>
bool util::load_from_int_vector_buffer(unsigned char*& text, int_vector_file_buffer<8, size_type_class>& text_buf)
{
    text_buf.reset();
    size_type_class n = text_buf.int_vector_size;
    if (text != NULL) {
        delete [] text;
        text = NULL;
    }
    text = new unsigned char[n];
    for (size_type_class i=0, r_sum=0, r=text_buf.load_next_block(); r_sum < n;) {
        for (; i < r_sum+r; ++i) {
            text[i] = text_buf[i-r_sum];
        }
        r_sum += r; r = text_buf.load_next_block();
    }
    return true;
}




template<class int_vector_type>
void util::set_random_bits(int_vector_type& v, int seed)
{
    if (0 == seed) {
        srand48((int)time(NULL));
    } else
        srand48(seed);

    uint64_t* data = v.m_data;
    if (v.empty())
        return;
    *data = (((uint64_t)lrand48()&0xFFFFULL)<<48)
            |(((uint64_t)lrand48()&0xFFFFULL)<<32)
            |(((uint64_t)lrand48()&0xFFFFULL)<<16)
            |((uint64_t)lrand48()&0xFFFFULL);
    for (typename int_vector_type::size_type i=1; i < (v.capacity()>>6); ++i) {
        *(++data) = (((uint64_t)lrand48()&0xFFFFULL)<<48)
                    |(((uint64_t)lrand48()&0xFFFFULL)<<32)
                    |(((uint64_t)lrand48()&0xFFFFULL)<<16)
                    |((uint64_t)lrand48()&0xFFFFULL);
    }
}

// all elements of vector v modulo m
template<class int_vector_type, class size_type_class>
void util::all_elements_mod(int_vector_type& v, size_type_class m)
{
    for (typename int_vector_type::size_type i=0; i < v.size(); ++i) {
        v[i] = v[i] % m;
    }
}

template<class int_vector_type>
void util::set_zero_bits(int_vector_type& v)
{
    uint64_t* data = v.m_data;
    if (v.empty())
        return;
    // TODO: replace by memset() but take care of size_t in the argument!
    *data = 0ULL;
    for (typename int_vector_type::size_type i=1; i < (v.capacity()>>6); ++i) {
        *(++data) = 0ULL;
    }
}

template<class int_vector_type>
void util::set_one_bits(int_vector_type& v)
{
    uint64_t* data = v.m_data;
    if (v.empty())
        return;
    *data = 0xFFFFFFFFFFFFFFFFULL;
    for (typename int_vector_type::size_type i=1; i < (v.capacity()>>6); ++i) {
        *(++data) = 0xFFFFFFFFFFFFFFFFULL;
    }
}

template<class int_vector_type>
void util::bit_compress(int_vector_type& v)
{
    typename int_vector_type::value_type max=0;
    for (typename int_vector_type::size_type i=0; i < v.size(); ++i) {
        if (v[i] > max) {
            max = v[i];
        }
    }
    uint8_t min_width = bit_magic::l1BP(max)+1;
    uint8_t old_width = v.get_int_width();
    if (old_width > min_width) {
        const uint64_t* read_data = v.m_data;
        uint64_t* write_data = v.m_data;
        uint8_t read_offset = 0;
        uint8_t write_offset = 0;
        for (typename int_vector_type::size_type i=0; i < v.size(); ++i) {
            uint64_t x = bit_magic::read_int_and_move(read_data, read_offset, old_width);
            bit_magic::write_int_and_move(write_data,  x, write_offset, min_width);
        }
        v.bit_resize(v.size()*min_width);
        v.set_int_width(min_width);
    }
}


template<class int_vector_type>
void util::set_all_values_to_k(int_vector_type& v, uint64_t k)
{
    uint64_t* data = v.m_data;
    if (v.empty())
        return;
    uint8_t int_width = v.m_int_width;
    if (int_width == 0) {
        throw std::logic_error("util::set_all_values_to_k can not be performed with int_width=0!");
    }
    k = k & (0xFFFFFFFFFFFFFFFFULL >> (64-int_width));
    uint64_t vec[67] = {0}; // allocate memory for the mask and initialize with zeros
    vec[0] = 0;
    uint8_t offset = 0;
    uint64_t n=0, vals=0;
    do { // loop terminates after at most 64 iterations
        vec[n] = vec[n] | (k << offset);
        offset += int_width;
        vals++;
        if (offset >= 64) {
            vec[n+1] = 0;
            vec[++n] = k >> (int_width-(offset-64));
            offset -= 64;
        }
    } while (offset != 0);

    typename int_vector_type::size_type n64 = v.capacity()/64;
    for (typename int_vector_type::size_type i=0; i < n64;) {
        for (uint64_t ii=0; ii < n and i < n64; ++ii,++i) {
            *(data++) = vec[ii];
        }
    }
}


template<class int_vector_type>
void util::set_to_id(int_vector_type& v)
{
    for (typename int_vector_type::size_type i=0; i < v.size(); ++i) {
        v[i] = i;
    }
}

template<class int_vector_type>
typename int_vector_type::size_type util::get_one_bits(const int_vector_type& v)
{
    const uint64_t* data = v.data();
    if (v.empty())
        return 0;
    typename int_vector_type::size_type result = bit_magic::b1Cnt(*data);
    for (typename int_vector_type::size_type i=1; i < (v.capacity()>>6); ++i) {
        result += bit_magic::b1Cnt(*(++data));
    }
    if (v.bit_size()&0x3F) {
        result -= bit_magic::b1Cnt((*data) & (~bit_magic::Li1Mask[v.bit_size()&0x3F]));
    }
    return result;
}


template<class int_vector_type>
typename int_vector_type::size_type util::get_onezero_bits(const int_vector_type& v)
{
    const uint64_t* data = v.data();
    if (v.empty())
        return 0;
    uint64_t carry = 0, oldcarry=0;
    typename int_vector_type::size_type result = bit_magic::b10Cnt(*data, carry);
    for (typename int_vector_type::size_type i=1; i < (v.capacity()>>6); ++i) {
        oldcarry = carry;
        result += bit_magic::b10Cnt(*(++data), carry);
    }
    if (v.bit_size()&0x3F) {// if bit_size is not a multiple of 64, substract the counts of the additional bits
        result -= bit_magic::b1Cnt(bit_magic::b10Map(*data, oldcarry) & bit_magic::Li0Mask[v.bit_size()&0x3F]);
    }
    return result;
}

template<class int_vector_type>
typename int_vector_type::size_type util::get_zeroone_bits(const int_vector_type& v)
{
    const uint64_t* data = v.data();
    if (v.empty())
        return 0;
    uint64_t carry = 1, oldcarry = 1;
    typename int_vector_type::size_type result = bit_magic::b01Cnt(*data, carry);
    for (typename int_vector_type::size_type i=1; i < (v.capacity()>>6); ++i) {
        oldcarry = carry;
        result += bit_magic::b01Cnt(*(++data), carry);
    }
    if (v.bit_size()&0x3F) {// if bit_size is not a multiple of 64, substract the counts of the additional bits
        result -= bit_magic::b1Cnt(bit_magic::b01Map(*data, oldcarry) & bit_magic::Li0Mask[v.bit_size()&0x3F]);
    }
    return result;
}

template<typename T>
std::string util::to_string(const T& t)
{
    std::stringstream ss;
    ss<<t;
    return ss.str();
}

template<typename T>
std::string util::to_latex_string(const T& t)
{
    return to_string(t);
}

}// end namespace sdsl

#endif // end file 
