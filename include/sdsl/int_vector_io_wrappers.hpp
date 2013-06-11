/*! \author Simon Gog (simon.gog@unimelb.edu.au)
 *  \file int_vector_io_wrappers.hpp
 *  \brief This file contains classes which could be used to encode and decode integer vectors when they are written
 *         to disk.
 *
 *  So far we two coders:
 *    int_vector_(load/serialize)_wrapper is just a dummy class which does not change the writing and reading process
 *    int_vector_(load/serialize)_vbyte_wrapper encodes and decodes each element using vbyte (aka escaping) coding.
 *    int_vector_(load/serialize)_vlen_wrapper encodes and decodes each element using a variable length code (e.g. Delta/Gamma/Fibonacci code)
 *
 *	  TODO: more sophisticated coders which have in contrast to the implemented version an internal state
 *          after the initialization. Like a Huffman coder which calculates the Huffman tree at initialization.
 */
#ifndef INCLUDE_SDSL_INT_VECTOR_IO_WRAPPERS
#define INCLUDE_SDSL_INT_VECTOR_IO_WRAPPERS

#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <sdsl/coder.hpp>

#include <iostream>

namespace sdsl
{

template<uint8_t fixedIntWidth=0>
class int_vector_serialize_vbyte_wrapper
{
    public:
        typedef int_vector<fixedIntWidth> int_vector_type;
        typedef typename int_vector_type::size_type size_type;
        typedef typename int_vector_type::value_type value_type;

    private:
        const int_vector_type& m_vec;

    public:
        int_vector_serialize_vbyte_wrapper(const int_vector_type& vec):m_vec(vec) {}

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            // (1) write size and int_width
            written_bytes += _sdsl_serialize_size_and_int_width(out, fixedIntWidth, m_vec.width(), m_vec.bit_size());
            // (2) write entries in vbyte coding
            for (size_type i=0; i < m_vec.size(); ++i) {
                value_type ww = m_vec[i];
                uint8_t w = ww & 0x7F;
                ww >>= 7;
                while (ww > 0) {
                    w |= 0x80; // mark overflow bit
                    out.write((const char*)&w, sizeof(uint8_t));  // write byte
                    w = ww & 0x7F;
                    ww >>= 7;
                    ++written_bytes;
                }
                out.write((const char*)&w, sizeof(uint8_t));  // write without overflow bit
                ++written_bytes;
            }
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }
};

template<uint8_t fixedIntWidth=0>
class int_vector_load_vbyte_wrapper
{
    public:
        typedef int_vector<fixedIntWidth> int_vector_type;
        typedef typename int_vector_type::size_type size_type;
        typedef typename int_vector_type::value_type value_type;

    private:
        int_vector_type& m_vec;

    public:
        int_vector_load_vbyte_wrapper(int_vector_type& vec):m_vec(vec) {}

        void load(std::istream& in) {
            size_type size;
            typename int_vector_type::int_width_type int_width;
            // (1) read size and int_width
            int_vector_trait<fixedIntWidth>::read_header(size, int_width, in);
            // (2) resize the vector
            m_vec.width(int_width);
            m_vec.bit_resize(size);
            // (3) read vbyte encoded entries an put them into the vector
            size_type i = 0;
            while (i < m_vec.size()) {
                value_type ww=0;
                uint8_t w=0;
                value_type shift=0;
                do {
                    in.read((char*)&w, sizeof(uint8_t));
                    ww |= (((value_type)(w&0x7F))<<shift);
                    shift += 7;
                } while ((w&0x80) > 0);
                m_vec[ i++ ] = ww;
            }
        }
};

template<class coder_type=coder::elias_delta>
class int_vector_serialize_vlen_wrapper
{
    public:
        typedef int_vector<> int_vector_type;
        typedef typename int_vector_type::size_type size_type;
        typedef typename int_vector_type::value_type value_type;

    private:
        const int_vector_type& m_vec;

    public:
        int_vector_serialize_vlen_wrapper(const int_vector_type& vec):m_vec(vec) {}

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            int_vector_type enc_vec;
            coder_type::encode(m_vec, enc_vec);
            written_bytes += enc_vec.serialize(out, child, "enc_vector");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }
};

template<class coder_type=coder::elias_delta>
class int_vector_load_vlen_wrapper
{
    public:
        typedef int_vector<> int_vector_type;
        typedef typename int_vector_type::size_type size_type;
        typedef typename int_vector_type::value_type value_type;

    private:
        int_vector_type& m_vec;

    public:
        int_vector_load_vlen_wrapper(int_vector_type& vec):m_vec(vec) {}

        void load(std::istream& in) {
            int_vector_type enc_vec;
            enc_vec.load(in);
            coder_type::decode(enc_vec, m_vec);
        }
};


template<class int_vector_type=int_vector<> >
class int_vector_serialize_wrapper
{
    public:
        typedef typename int_vector_type::size_type size_type;

    private:
        const int_vector_type& m_vec;

    public:
        int_vector_serialize_wrapper(const int_vector_type& vec):m_vec(vec) {}

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            return m_vec.serialize(out, v, name);
        }
};

template<class int_vector_type=int_vector<> >
class int_vector_load_wrapper
{
    public:
        typedef typename int_vector_type::size_type size_type;
    private:
        int_vector_type& m_vec;

    public:
        int_vector_load_wrapper(int_vector_type& vec):m_vec(vec) {}
        void load(std::istream& in) {
            m_vec.load(in);
        }
};

template<class int_vector_serialize_wrapper_type=int_vector_serialize_wrapper<> >
class int_vector_serialize_min_overhead
{

};

} // end namespace

#endif // end include guard
