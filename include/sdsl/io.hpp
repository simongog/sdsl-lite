/* sdsl - succinct data structures library
    Copyright (C) 2013 Simon Gog

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
/*! \file io.hpp
    \brief io.hpp contains some methods for reading/writing sdsl structures.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_IO
#define INCLUDED_SDSL_IO

#include "util.hpp"
#include "sdsl_concepts.hpp"
#include "structure_tree.hpp"
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include <cctype>

namespace sdsl
{

//! Load sdsl-object v from a file.
/*!
 * \param v sdsl-
 * \param file Name of the serialized file.
 */
template<class T>
bool load_from_file(T& v, const std::string& file);

template<>
bool load_from_file(void*&, const std::string& file);

//! Specialization of load_from_file for a char array
/*  \pre v=nullptr
 */
bool load_from_file(char*& v, const std::string& file);

//! Load an int_vector from a plain array of `num_bytes`-byte integers with X in \{0, 1,2,4,8\} from disk.
// TODO: Remove ENDIAN dependency.
template<class t_int_vec>
bool load_vector_from_file(t_int_vec& v, const std::string& file, uint8_t num_bytes=1, uint8_t max_int_width=64)
{
    if ((uint8_t)0 == num_bytes) {  // if byte size is variable read int_vector<0> from file
        return load_from_file(v, file);
    } else if (num_bytes == 'd') {
        uint64_t x = 0, max_x = 0;
        isfstream in(file);
        if (!in) {
            return false;
        } else {
            std::vector<uint64_t> tmp;
            while (in >> x) {
                tmp.push_back(x);
                max_x = std::max(x, max_x);
            }
            v.width(bits::hi(max_x)+1); v.resize(tmp.size());
            for (size_t i=0; i < tmp.size(); ++i) {
                v[i] = tmp[i];
            }
            return true;
        }
    } else {
        off_t file_size = util::file_size(file);
        if (file_size == 0) {
            v.resize(0);
            return true;
        }
        if (file_size % num_bytes != 0) {
            throw std::logic_error("file size "+util::to_string(file_size)+" of \""+ file
                                   +"\" is not a multiple of "+util::to_string(num_bytes));
            return false;
        }
        isfstream in(file);
        if (in) {
            v.width(std::min((int)8*num_bytes, (int)max_int_width));
            v.resize(file_size / num_bytes);
            if (8 == t_int_vec::fixed_int_width and 1 == num_bytes) {  // if int_vector<8> is created from byte alphabet file
                in.read((char*)v.m_data, file_size);
            } else {
                size_t idx=0;
                const size_t block_size = constants::SDSL_BLOCK_SIZE*num_bytes;
                uint8_t* buf = new uint8_t[block_size];
                // TODO: check for larger alphabets with num_bytes*8 = v::fixed_int_width

                uint64_t x = 0; // value
                uint8_t  cur_byte = 0;
                do {
                    in.read((char*)buf, block_size);
                    size_t read = in.gcount();
                    uint8_t* begin = buf;
                    uint8_t* end   = begin+read;
                    while (begin < end) {
                        x |= ((uint64_t)(*begin)) << (cur_byte*8);
                        ++cur_byte;
                        if (cur_byte == num_bytes) {
                            v[idx++] = x;
                            cur_byte = 0;
                            x = 0ULL;
                        }
                        ++begin;
                    }
                } while (idx < v.size());
                delete [] buf;
                in.close();
            }
            return true;
        } else {
            return false;
        }
    }
}

//! Store a data structure to a file.
/*! The data structure has to provide a serialize function.
 *  \param v Data structure to store.
 *  \param file Name of the file where to store the data structure.
 *  \param Return if the data structure was stored successfully
 */
template<class T>
bool store_to_file(const T& v, const std::string& file);

//! Specialization of store_to_file for a char array
bool store_to_file(const char* v, const std::string& file);

//! Specialization of store_to_file for int_vector
template<uint8_t t_width>
bool store_to_file(const int_vector<t_width>& v, const std::string& file, bool write_fixed_as_variable=false);


//! Store an int_vector as plain int_type array to disk
template<class int_type, class t_int_vec>
bool store_to_plain_array(t_int_vec& v, const std::string& file)
{
    osfstream out(file);
    if (out) {
        for (typename t_int_vec::size_type i=0; i<v.size(); ++i) {
            int_type x = v[i];
            out.write((char*)&x, sizeof(int_type));
        }
        return true;
    } else {
        return false;
    }
}

template<class T>
size_t serialize_empty_object(std::ostream&, structure_tree_node* v=nullptr, std::string name="", const T* t=nullptr)
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*t));
    size_t written_bytes = 0;
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}



//! Get the size of a data structure in bytes.
/*!
 *  \param v A reference to the data structure for which the size in bytes should be calculated.
 */
template<class T>
typename T::size_type size_in_bytes(const T& t);

//! Get the size of a data structure in mega bytes (MiB).
/*!
 *  \param t A reference to the data structure for which the size in bytes should be calculated.
 */
template<class T>
double size_in_mega_bytes(const T& t);

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
size_t write_member(const T& t, std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")
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
 * \param vec The vector which should be serialized.
 * \param out Output stream to which should be written.
 * \param v   Structure tree node. Note: If all elements have the same
 *            structure, then it is tried to combine all elements (i.e.
 *            make one node w with size set to the cumulative sum of all
 *           sizes of the children)
 */
template<class T>
size_t serialize_vector(const std::vector<T>& vec, std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")
{
    if (vec.size() > 0) {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, "std::vector<"+util::class_name(vec[0])+">");
        size_t written_bytes = 0;
        for (typename std::vector<T>::size_type i = 0; i < vec.size(); ++i) {
            written_bytes += vec[i].serialize(out, child, "[]");
        }
        structure_tree::add_size(child, written_bytes);
        sdsl::structure_tree::merge_children(child);
        return written_bytes;
    } else {
        return 0;
    }
}

//! Load all elements of a vector from a input stream
/*! \param vec  Vector whose elements should be loaded.
 *  \param in   Input stream.
 *  \par Note
 *   The vector has to be resized prior the loading
 *   of its elements.
 */
template<class T>
void load_vector(std::vector<T>& vec, std::istream& in)
{
    for (typename std::vector<T>::size_type i = 0; i < vec.size(); ++i) {
        vec[i].load(in);
    }
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

//! Internal function used by csXprintf
uint64_t _parse_number(std::string::const_iterator& c, const std::string::const_iterator& end);

//! Internal function used by csXprintf
template<class t_csa>
const t_csa& _idx_csa(const t_csa& t, csa_tag)
{
    return t;
};

//! Internal function used by csXprintf
template<class t_cst>
const typename t_cst::csa_type& _idx_csa(const t_cst& t, cst_tag)
{
    return t.csa;
};

//! Internal function used by csXprintf
template<class t_csa>
std::string _idx_lcp_val(const t_csa&, uint64_t, uint64_t, csa_tag)
{
    return "";
}

//! Internal function used by csXprintf
template<class t_cst>
std::string _idx_lcp_val(const t_cst& t, uint64_t i, uint64_t w, cst_tag)
{
    return util::to_string(t.lcp[i], w);
}

//! Prints members of CSAs and CSTs
/*! This is a printf like method to write members of CSAs and CSTs into an outstream.
 * \tparam t_idx   Type of the index. Class should be of concept csa_tag or cst_tag.
 * \param out      Output stream.
 * \param format   Format string. See explanation below.
 * \param idx      CSA or CST object.
 * \param sentinel Character which should replace the \0-symbol in BWT/ TEXT.
 *
 * \par Format string
 *   Each line of the output will be formatted according to the format string.
 *   All content, except tokens which start with `%` will be copied. Tokens
 *   which start with `%` will be replaced as follows (let w be a positive
 *    number. setw(w) is used to format single numbers):
 *
 *      Token      |  Replacement | Comment
 *      -----------------------------------------------------------------------
 *       %[w]I     | Row index i.                           |
 *       %[w]S     | SA[i]                                  |
 *       %[w]s     | ISA[i]                                 |
 *       %[w]P     | PSI[i]                                 |
 *       %[w]p     | LF[i]                                  |
 *       %[w]L     | LCP[i]                                 | only for CSTs
 *       %[w]B     | BWT[i]                                 |
 *       %[w[:W]]T | Print min(idx.size(),w) chars of each  |
 *                 | suffix, each char formatted by setw(W).|
 *       %%        | %                                      |
 */
template<class t_idx>
void csXprintf(std::ostream& out, const std::string& format, const t_idx& idx, char sentinel='$')
{
    typename t_idx::index_category cat;
    const typename t_idx::csa_type& csa = _idx_csa(idx, cat);
    vector<std::string> res(csa.size());
    for (std::string::const_iterator c = format.begin(), s=c; c != format.end(); s=c) {
        while (c != format.end() and* c != '%') ++c;   // string before the next `%`
        if (c > s) {  // copy format string part
            vector<std::string> to_copy(csa.size(), std::string(s, c));
            transform(res.begin(), res.end(), to_copy.begin(), res.begin(), std::plus<std::string>());
        }
        if (c == format.end()) break;
        ++c; // skip `%`
        uint64_t w = _parse_number(c, format.end());  // element width
        if (c == format.end()) break;
        uint64_t W = 0; // character width
        if (':' == *c) {
            ++c;
            W = _parse_number(c, format.end());
        }
        if (c == format.end()) break;
        for (uint64_t i=0; i<csa.size(); ++i) {
            switch (*c) {
                case 'I': res[i] += util::to_string(i,w); break;
                case 'S': res[i] += util::to_string(csa[i],w); break;
                case 's': res[i] += util::to_string(csa(i),w); break;
                case 'P': res[i] += util::to_string(csa.psi[i],w); break;
                case 'p': res[i] += util::to_string(csa.psi(i),w); break;
                case 'L': res[i] += _idx_lcp_val(idx,i,w, cat); break;
                case 'B': if (0 == csa.bwt[i]) {
                        res[i] += util::to_string(sentinel,w);
                    } else {
                        res[i] += util::to_string(csa.bwt[i],w);
                    }
                    break;
                case 'T': for (uint64_t k=0; (w>0 and k < w) or(0==w and k < csa.size()); ++k) {
                        if (0 == csa.text[(csa[i]+k)%csa.size()]) {
                            res[i] += util::to_string(sentinel, W);
                        } else {
                            res[i] += util::to_string(csa.text[(csa[i]+k)%csa.size()], W);
                        }
                    }
                    break;
                case '%': res[i] += "%"; break;
            }
        }
        ++c;
    }
    for (size_t i=0; i<res.size(); ++i) out << res[i] << std::endl;
}



//! Returns the file name of the resource.
/*!
 * \param  key        Resource key.
 * \param  config    Cache configuration.
 * \return The file name of the resource.
 */
std::string cache_file_name(const std::string& key, const cache_config& config);

//! Register the existing resource specified by the key to the cache
/*!
 *  \param key        Resource key.
 *  \param config    Cache configuration.
 *
 *  Note: If the resource does not exist under the given key,
 *  it will be not added to the cache configuration.
 */
void register_cache_file(const std::string& key, cache_config& config);

//! Checks if the resource specified by the key exists in the cache.
/*!
  \param key    Resource key.
  \param config Cache configuration.
  \return True, if the file exists, false otherwise.
*/
bool cache_file_exists(const std::string& key, const cache_config& config);

//! Returns a name for a temporary file. I.e. the name was not used before.
std::string tmp_file(const cache_config& config, std::string name_part="");

template<class T>
bool load_from_cache(T& v, const std::string& key, const cache_config& config)
{
    std::string file = cache_file_name(key, config);
    if (load_from_file(v, file)) {
        if (util::verbose) {
            std::cerr << "Load `" << file << std::endl;
        }
        return true;
    } else {
        std::cerr << "WARNING: Could not load file '";
        std::cerr << file << "'" << std::endl;
        return false;
    }
}

//! Stores the object v as a resource in the cache.
/*!
 *  \param
 */
template<class T>
bool store_to_cache(const T& v, const std::string& key, cache_config& config)
{
    std::string file = cache_file_name(key, config);
    if (store_to_file(v, file)) {
        config.file_map[std::string(key)] = file;
        return true;
    } else {
        std::cerr<<"WARNING: store_to_cache: could not store file `"<< file <<"`" << std::endl;
        return false;
    }
}

//==================== Template functions ====================

template<class T>
typename T::size_type size_in_bytes(const T& t)
{
    if ((&t) == nullptr)
        return 0;
    nullstream ns;
    return t.serialize(ns);
}

template<class T>
double size_in_mega_bytes(const T& t)
{
    return size_in_bytes(t)/(1024.0*1024.0);
}

template<class T>
bool store_to_file(const T& t, const std::string& file)
{
    osfstream out(file, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out) {
        if (util::verbose) {
            std::cerr<<"ERROR: store_to_file not successful for: `"<<file<<"`"<<std::endl;
        }
        return false;
    }
    t.serialize(out);
    out.close();
    if (util::verbose) {
        std::cerr<<"INFO: store_to_file: `"<<file<<"`"<<std::endl;
    }
    return true;
}

bool store_to_file(const char* v, const std::string& file);

template<uint8_t t_width>
bool store_to_file(const int_vector<t_width>& v, const std::string& file, bool write_fixed_as_variable)
{
    osfstream out(file, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out) {
        std::cerr<<"ERROR: util::store_to_file:: Could not open file `"<<file<<"`"<<std::endl;
        return false;
    } else {
        if (util::verbose) {
            std::cerr<<"INFO: store_to_file: `"<<file<<"`"<<std::endl;
        }
    }
    v.serialize(out, nullptr, "", write_fixed_as_variable);
    out.close();
    return true;
}

template<class T>
bool load_from_file(T& v, const std::string& file)
{
    isfstream in(file, std::ios::binary | std::ios::in);
    if (!in) {
        if (util::verbose) {
            std::cerr << "Could not load file `" << file << "`" << std::endl;
        }
        return false;
    }
    v.load(in);
    in.close();
    if (util::verbose) {
        std::cerr << "Load file `" << file << "`" << std::endl;
    }
    return true;
}

}
#endif
