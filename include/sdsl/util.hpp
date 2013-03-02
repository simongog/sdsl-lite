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
#include "structure_tree.hpp"
#include "config.hpp"  // for constants 
#include <iosfwd>      // forward declaration of ostream
#include <stdint.h>    // for uint64_t uint32_t declaration
#include <cassert>
#include <fstream>     // file stream for storeToFile and loadFromFile
#include <ctime>       // for rand initialization
#include <string>
#include <locale>       // for class_to_hash
#include <string.h>    // for strlen and strdup
#include <libgen.h>    // for basename
#include <cstdlib>
#include <unistd.h>    // for getpid 
#include <sstream>     // for to_string method
#include <stdexcept>   // for std::logic_error
#include <typeinfo>    // for typeid


// macros to transform a defined name to a string
#define SDSL_STR(x) #x
#define SDSL_XSTR(s) SDSL_STR(s)

//! Namespace for the succinct data structure library.
namespace sdsl
{

template<uint8_t>
class int_vector;	 // forward declaration

//! A namespace for helper functions
namespace util
{

static bool verbose = false;

void set_verbose();

//! Get the size of a file in bytes
off_t get_file_size(const char* file_name);

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

//! Expands the integer width to new_width >= v.get_int_width()
template<class int_vector_type>
void expand_width(int_vector_type&v, uint8_t new_width);

//! All elements of v modulo m
template<class int_vector_type>
void all_elements_mod(int_vector_type& v, typename int_vector_type::size_type m);


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

//! Get the smallest position \f$i\geq idx\f$ where a bit is set
/*! \param v The int_vector in which the bit is searched
 *	\param idx The start position for the search \f$ 0\leq idx < v.bit_size()\f$
 *	\return The smallest position greater or equal to idx, where corresponding bit is 1 or v.bit_size() if no such position exists
 *	\par Time complexity
 *		\f$ \Order{n} \f$
 */
template <class int_vector_type>
typename int_vector_type::size_type next_bit(const int_vector_type& v, uint64_t idx);

//! Get the greatest position \f$i\leq idx\f$ where a bit is set
/*! \param v The int_vector in which the bit is searched
 *	\param idx The start position for the search \f$ 0\leq idx < v.bit_size()\f$
 *	\return The greatest position smaller or equal to idx, where corresponding bit is 1 or v.bit_size() if no such position exists
 *	\par Time complexity
 *		\f$ \Order{n} \f$
*/
template <class int_vector_type>
typename int_vector_type::size_type prev_bit(const int_vector_type& v, uint64_t idx);

//! Load a data structure from a file.
/*! The data structure has to provide a load function.
 * \param v Data structure to load.
   \param file_name Name of the serialized file.
 */
template<class T>
bool load_from_file(T& v, const char* file_name);

template<>
bool load_from_file(void*&, const char* file_name);

//! Load an int_vector from a plain array of `num_bytes`-byte integers with X in \{0, 1,2,4,8\} from disk.
// TODO: Remove ENDIAN dependency: currently in BIG_ENDIAN format
template<class int_vector_type>
bool load_vector_from_file(int_vector_type &v, const char* file_name, uint8_t num_bytes=1, uint8_t max_int_width=64){
	if ( (uint8_t)0 == num_bytes ){ // if byte size is variable read int_vector<0> from file
		return load_from_file(v, file_name);
	}else {
		off_t file_size = get_file_size( file_name );
		if ( file_size == 0 ){
			v.resize(0);
			return true;
		}
		if ( file_size % num_bytes != 0 ){
			throw std::logic_error("file size "+to_string(file_size)+" of \""+to_string(file_name)
											   +"\" is not a multiple of "+to_string(num_bytes));
			return false;
		}
		std::ifstream in(file_name);
		if ( in ){
			v.set_int_width( std::min( (int)8*num_bytes, (int)max_int_width ) );
			v.resize( file_size / num_bytes );
			if ( 8 == int_vector_type::fixed_int_width and 1 == num_bytes ){ // if int_vector<8> is created from byte alphabet file 
				in.read((char*)v.m_data, file_size);
			}else{
				size_t idx=0;
				const size_t block_size = constants::SDSL_BLOCK_SIZE*num_bytes;
				uint8_t * buf = new uint8_t[block_size];
				// TODO: check for larger alphabets with num_bytes*8 = v::fixed_int_width

				uint64_t x = 0; // value
				uint8_t  cur_byte = 0;
				do{
					in.read((char*)buf, block_size);
					size_t read = in.gcount();
					uint8_t* begin = buf;
					uint8_t* end   = begin+read;
					while ( begin < end ){
						x |= (*begin) << (cur_byte*8);
						++cur_byte;
						if ( cur_byte == num_bytes ){
							v[idx++] = x;
							cur_byte = 0;
							x = 0ULL;
						}
						++begin;
					}
				} while( idx < v.size() );
				delete [] buf;
				in.close();
			}
			return true;
		}else{
			return false;
		}
	}
}

//! Store an int_vector as plain int_type array to disk 
template<class int_type, class int_vector_type>
bool store_to_plain_array(int_vector_type &v, const char* file_name){
	std::ofstream out(file_name);
	if ( out ){
		for (typename int_vector_type::size_type i=0; i<v.size(); ++i){
			int_type x = v[i];
			out.write((char*)&x, sizeof(int_type));
		}
		return true;
	}else{
		return false;
	}
}

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
template<uint8_t fixed_int_width>
bool store_to_file(const int_vector<fixed_int_width>& v, const char* file_name, bool write_fixed_as_variable=false);

//! Demangle the class name of typeid(...).name()
/*!
 *	\param name A pointer to the the result of typeid(...).name()
 */
std::string demangle(const char* name);

//! Demangle the class name of typeid(...).name() and remove the "sdsl::"-prefix, "unsigned int",...
std::string demangle2(const char* name);

//! Transforms the demangled class name of an object to a hash value.
template<class T>
std::string class_to_hash(const T& t){
	std::locale loc;
	const std::collate<char>& coll = std::use_facet<std::collate<char> >(loc);
	std::string name = sdsl::util::demangle2(typeid(T).name());
	uint64_t my_hash = coll.hash(name.data(),name.data()+name.length());
	return to_string(my_hash);
}

template<class T>
std::string class_name(const T& t)
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

//! Get the size of a data structure in mega bytes (MiB).
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
    nullstream(): std::ios(&m_sbuf), std::ostream(&m_sbuf), m_sbuf() {}
};

// Writes primitive-typed variable t to stream out
template<class T>
size_t write_member(const T& t, std::ostream& out, sdsl::structure_tree_node* v=NULL, std::string name="")
{
    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, util::class_name(t));
    out.write((char*)&t, sizeof(t));
    size_t written_bytes = sizeof(t);
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

// Specialization for std::string
template<>
size_t write_member<std::string>(const std::string& t, std::ostream& out, sdsl::structure_tree_node* v, std::string name);


// Writes primitive-typed variable t to stream out
template<class T>
void read_member(T& t, std::istream& in)
{
    in.read((char*)&t, sizeof(t));
}

// Specialization for std::string
template<>
void read_member<std::string>(std::string& t, std::istream& in);

//! Serialize each element of an std::vector 
/*!
 * \param vec	The vector which should be serialized.
 * \param out	Output stream to which should be written.
 * \param v		Structure tree node. Note: If all elements have the same
 *              structure, then it is tried to combine all elements (i.e.
 *              make one node w with size set to the cumulative sum of all
 *              sizes of the children) 
 */
template<class T>
size_t serialize_vector(const std::vector<T> &vec, std::ostream& out, sdsl::structure_tree_node* v=NULL, std::string name=""){
	if ( vec.size() > 0 ){
		sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, "std::vector<"+util::class_name(vec[0])+">" );
		size_t written_bytes = 0;
		for (typename std::vector<T>::size_type i = 0; i < vec.size(); ++i ){
			written_bytes += vec[i].serialize(out, child, "[]");
		}
		structure_tree::add_size(child, written_bytes);
		sdsl::structure_tree::merge_children(child);
		return written_bytes;
	}else{
		return 0;
	}
}

//! Load all elements of a vector from a input stream
/*! \param vec	Vector whose elements should be loaded.
 *  \param in   Input stream.
 *  \par Note
 *   The vector has to be resized prior the loading 
 *   of its elements.
 */
template<class T>
void load_vector(std::vector<T> &vec, std::istream& in){
	for (typename std::vector<T>::size_type i = 0; i < vec.size(); ++i){
		vec[i].load(in);
	}
}

//! Get the process id of the current process
uint64_t get_pid();

class _id_helper {
    private:
        static uint64_t id;
    public:
        static uint64_t getId() {
            return id++;
        }
};


//! Get a unique id inside the process
inline uint64_t get_id() {
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

template<format_type F, class X>
void write_structure(const X& x, std::ostream& out)
{
    structure_tree_node* v = new structure_tree_node();
	nullstream ns;
    x.serialize(ns, v, "");
    if (v->children.size() > 0) {
        sdsl::write_structure_tree<F>(v->children[0], out);
    }
    delete v;
}

std::string cache_file_name(const char* key, const cache_config &config);

void register_cache_file(const char* key, cache_config &config);

bool cache_file_exists(const char* key, const cache_config &config);

template<class T>
bool load_from_cache(T&v, const char* key, const cache_config &config){
	std::string file_name = cache_file_name(key, config);
	if( load_from_file(v, file_name.c_str()) ){
		return true;
	}
	return false;
}


template<class T>
bool store_to_cache(const T& v, const char* key, cache_config &config){
	std::string file_name = cache_file_name(key, config);
	if ( store_to_file(v, file_name.c_str()) ){
		config.file_map[std::string(key)] = file_name;
		return true;	
	}else{
		return false;
	}
}

} // end namespace util

//==================== Template functions ====================


template<class T>
typename T::size_type util::get_size_in_bytes(const T& t) {
    if ((&t) == NULL)
        return 0;
    util::nullstream ns;
    return t.serialize(ns);
}

template<class T>
double util::get_size_in_mega_bytes(const T& t) {
    return get_size_in_bytes(t)/(1024.0*1024.0);
}

template<class T>
bool util::store_to_file(const T& t, const char* file_name) {
    std::ofstream out;
    out.open(file_name, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out)
        return false;
    t.serialize(out);
    out.close();
    return true;
}

inline bool util::store_to_file(const char* v, const char* file_name) {
    std::ofstream out;
    out.open(file_name, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out)
        return false;
    uint64_t n = strlen((const char*)v);
    out.write(v, n);
    out.close();
    return true;
}

template<uint8_t fixed_int_width>
bool util::store_to_file(const int_vector<fixed_int_width>& v, const char* file_name, bool write_fixed_as_variable) {
    std::ofstream out;
    out.open(file_name, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out)
        return false;
    v.serialize(out, NULL, "", write_fixed_as_variable);
    out.close();
    return true;
}

template<class T>
bool util::load_from_file(T& v, const char* file_name) {
    std::ifstream in;
    in.open(file_name, std::ios::binary | std::ios::in);
    if (!in)
        return false;
    v.load(in);
    in.close();
    return true;
}


template<class int_vector_type>
void util::set_random_bits(int_vector_type& v, int seed) {
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
template<class int_vector_type>
void util::all_elements_mod(int_vector_type& v, typename int_vector_type::size_type m) {
    for (typename int_vector_type::size_type i=0; i < v.size(); ++i) {
        v[i] = v[i] % m;
    }
}

template<class int_vector_type>
void util::set_zero_bits(int_vector_type& v) {
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
void util::set_one_bits(int_vector_type& v) {
    uint64_t* data = v.m_data;
    if (v.empty())
        return;
    *data = 0xFFFFFFFFFFFFFFFFULL;
    for (typename int_vector_type::size_type i=1; i < (v.capacity()>>6); ++i) {
        *(++data) = 0xFFFFFFFFFFFFFFFFULL;
    }
}

template<class int_vector_type>
void util::bit_compress(int_vector_type& v) {
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
void util::expand_width(int_vector_type&v, uint8_t new_width){
	uint8_t old_width = v.get_int_width();
	typename int_vector_type::size_type n = v.size();
	if ( new_width > old_width and n > 0 ){
		typename int_vector_type::size_type i, old_pos, new_pos;
		new_pos = (n-1)*new_width;
		old_pos = (n-1)*old_width;
		v.bit_resize(v.size()*new_width);
    	for (i=0; i < n; ++i, new_pos-=new_width, old_pos-=old_width) {
			v.set_int(new_pos, v.get_int(old_pos, old_width), new_width);
		}
		v.set_int_width(new_width);
	}
}

template<class int_vector_type>
void util::set_all_values_to_k(int_vector_type& v, uint64_t k) {
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

//! Set v[i] = i for i=[0..v.size()-1]
template<class int_vector_type>
void util::set_to_id(int_vector_type& v) {
    for (typename int_vector_type::size_type i=0; i < v.size(); ++i) {
        v[i] = i;
    }
}

template<class int_vector_type>
typename int_vector_type::size_type util::get_one_bits(const int_vector_type& v) {
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
typename int_vector_type::size_type util::get_onezero_bits(const int_vector_type& v) {
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
typename int_vector_type::size_type util::get_zeroone_bits(const int_vector_type& v) {
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

template <class int_vector_type>
typename int_vector_type::size_type util::next_bit(const int_vector_type& v, uint64_t idx)
{
    uint64_t pos = idx>>6;
    uint64_t node = v.data()[pos];
    node >>= (idx&0x3F);
    if(node) {
        return idx+bit_magic::r1BP(node);
    } else {
        ++pos;
        while((pos<<6) < v.bit_size()) {
            if(v.data()[pos]) {
                return (pos<<6)|bit_magic::r1BP(v.data()[pos]);
            }
            ++pos;
        }
        return v.bit_size();
    }
}

template <class int_vector_type>
typename int_vector_type::size_type util::prev_bit(const int_vector_type& v, uint64_t idx)
{
    uint64_t pos = idx>>6;
    uint64_t node = v.data()[pos];
    node <<= 63-(idx&0x3F);
    if(node) {
        return bit_magic::l1BP(node)+(pos<<6)-(63-(idx&0x3F));
    } else {
        --pos;
        while((pos<<6) < v.bit_size() ) {
            if(v.data()[pos]) {
                return (pos<<6)|bit_magic::l1BP(v.data()[pos]);
            }
            --pos;
        }
        return v.bit_size();
    }
}

template<typename T>
std::string util::to_string(const T& t) {
    std::stringstream ss;
    ss<<t;
    return ss.str();
}

template<typename T>
std::string util::to_latex_string(const T& t) {
    return to_string(t);
}

}// end namespace sdsl

#endif // end file 
