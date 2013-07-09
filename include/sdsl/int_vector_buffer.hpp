#ifndef INCLUDED_INT_VECTOR_BUFFER
#define INCLUDED_INT_VECTOR_BUFFER

#include <assert.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sdsl/int_vector.hpp>

using namespace std;
using namespace sdsl;


namespace sdsl
{

template<uint8_t t_width=0>
class int_vector_buffer
{
// TODO: fstream -> ifstream und ofstream ersetzten???
// TODO: iterator (Simon)
// TODO: discuss mit Simon ++ on int_vector_buffer<1>
    private:
        static_assert(t_width <= 64 , "int_vector_buffer: width must be at most 64bits.");
        sdsl::isfstream              m_ifile;
        sdsl::osfstream              m_ofile;
        string                     m_filename;
        int_vector<t_width>        m_buffer;
        bool                       m_need_to_write = false;
        // length of int_vector header in bytes
        static constexpr uint64_t  m_offset        = t_width ? 8 : 9;
        uint64_t                   m_buffersize    = 0;  // in elements! m_buffersize*width() must be a multiple of 8!
        uint64_t                   m_max_elements  = 0;
        uint64_t                   m_begin         = 0; // number in elements
        bool                       m_persistent    = true;
        bool                       m_closed        = true;

        //! Read block containing element v[idx]
        void read_block(const uint64_t idx) {
            m_begin = (idx/m_buffersize)*m_buffersize;
            if (m_begin >= m_max_elements) {
                util::set_to_value(m_buffer, 0);
            } else {
                m_ifile.seekg(m_offset+(m_begin*width())/8);
                assert(m_ifile.good());
                m_ifile.read((char*) m_buffer.data(), (m_buffersize*width())/8);
                if ((uint64_t)m_ifile.gcount() < (m_buffersize*width())/8) {
                    m_ifile.clear();
                }
                assert(m_ifile.good());
                for (uint64_t i=m_max_elements-m_begin; i<m_buffersize; ++i) {
                    m_buffer[i] = 0;
                }
            }
        }

        //! Write current block to file
        void write_block() {
            if (m_need_to_write) {
                m_ofile.seekp(m_offset+(m_begin*width())/8);
                assert(m_ofile.good());
                if (m_begin+m_buffersize >= m_max_elements) {
                    //last block in file
                    uint64_t wb = ((m_max_elements-m_begin)*width()+7)/8;
                    m_ofile.write((char*) m_buffer.data(), wb);
                } else {
                    m_ofile.write((char*) m_buffer.data(), (m_buffersize*width())/8);
                }
                m_ofile.flush();
                assert(m_ofile.good());
                m_need_to_write = false;
            }
        }

    public:
        int_vector_buffer() {
            m_buffer = int_vector<t_width>();
        }

        int_vector_buffer(const std::string _filename, const bool _open_existing_file, const uint64_t _buffersize=1024*1024, const uint8_t _width=t_width, const bool _persistent=true) {
            m_persistent = _persistent;
            m_filename = _filename;
//             if (0 == t_width) {
            m_buffer.width(_width);
//             }
            if (!_open_existing_file) {
                // Create file, if it already exist it will be cleared
                m_ofile.open(_filename.c_str(), std::ios::out|std::ios::binary);
                assert(m_ofile.good());
                m_ofile.close();
                assert(m_ofile.good());
            }
            // Open file for IO
            m_ofile.open(_filename.c_str(), std::ios::in|std::ios::out|std::ios::binary);
            assert(m_ofile.good());
            m_ifile.open(_filename.c_str(), std::ios::in|std::ios::binary);
            assert(m_ifile.good());
            if (_open_existing_file) {
                uint64_t size  = 0;
                uint8_t  width = 0;
                int_vector_trait<t_width>::read_header(size, width, m_ifile);
                assert(m_ifile.good());
//                 if (0 == t_width) {
                m_buffer.width(width);
//                 }
                m_max_elements = size/this->width();
            }
            if (0==(_buffersize*8)%width()) {
                m_buffersize = _buffersize*8/width(); //m_buffersize might not be multiple of 8, but m_buffersize*width is.
            } else {
                uint64_t element_buffersize = (_buffersize*8)/width()+1; //one more element than fits into given buffersize in byte
                m_buffersize = element_buffersize+7 - (element_buffersize+7)%8; // take next multiple of 8
            }
            m_buffer = int_vector<t_width>(m_buffersize, 0, width());
            if (0!=m_buffersize) read_block(0);
            m_need_to_write = false;
            m_closed = false;
        }

        //! Move constructor.
        int_vector_buffer(int_vector_buffer&& ivb) :
            m_filename((string&&)ivb.m_filename),
            m_buffer((int_vector<t_width>&&)ivb.m_buffer),
            m_need_to_write(ivb.m_need_to_write),
            m_buffersize(ivb.m_buffersize),
            m_max_elements(ivb.m_max_elements),
            m_begin(ivb.m_begin),
            m_persistent(ivb.m_persistent),
            m_closed(ivb.m_closed) {
            ivb.m_ifile.close();
            ivb.m_ofile.close();
            m_ifile.open(m_filename.c_str(), std::ios::in|std::ios::binary);
            m_ofile.open(m_filename.c_str(), std::ios::in|std::ios::out|std::ios::binary);
            assert(m_ifile.good());
            assert(m_ofile.good());
            ivb.m_persistent = true;
            ivb.m_closed = true;
        }

        ~int_vector_buffer() {
            if (!m_closed) {
                close();
            }
            if (!m_persistent) {
                sdsl::remove(m_filename);
            }
        }

        //! Move assignment operator.
        int_vector_buffer<t_width>& operator=(int_vector_buffer&& ivb) {
            close();
            ivb.m_ifile.close();
            ivb.m_ofile.close();
            std::swap(m_filename, ivb.m_filename);
            m_ifile.open(m_filename.c_str(), std::ios::in|std::ios::binary);
            m_ofile.open(m_filename.c_str(), std::ios::in|std::ios::out|std::ios::binary);
            assert(m_ifile.good());
            assert(m_ofile.good());
            std::swap(m_closed, ivb.m_closed);
            std::swap(m_buffer, ivb.m_buffer);
            std::swap(m_need_to_write, ivb.m_need_to_write);
            std::swap(m_buffersize, ivb.m_buffersize);
            std::swap(m_max_elements, ivb.m_max_elements);
            std::swap(m_begin, ivb.m_begin);
            std::swap(m_persistent, ivb.m_persistent);
            return *this;
        }

        uint8_t width() const {
            return m_buffer.width();
        }

        uint64_t size() const {
            return m_max_elements;
        }

        bool persistence() const {
            return m_persistent;
        }

        void persistence(bool persistent) {
            m_persistent = persistent;
        }

        uint64_t buffersize() const { //returns buffersize in byte
            assert(m_buffersize*width()%8==0);
            return (m_buffersize*width())/8;
        }

        void buffersize(uint64_t buffersize) { //set buffersize in byte  //remove this?
            write_block();
            if (0==(buffersize*8)%width()) {
                m_buffersize = buffersize*8/width(); //m_buffersize might not be multiple of 8, but m_buffersize*width() is.
            } else {
                uint64_t element_buffersize = (buffersize*8)/width()+1; //one more element than fits into given buffersize in byte
                m_buffersize = element_buffersize+7 - (element_buffersize+7)%8; // take next multiple of 8
            }
            m_buffer.resize(m_buffersize);
            if (0!=m_buffersize) read_block(0);
        }

        void reset() {
            //delete all content
            assert(!m_persistent);
            // reset file
            assert(m_ifile.good());
            assert(m_ofile.good());
            m_ifile.close();
            m_ofile.close();
            m_ofile.open(m_filename.c_str(), std::ios::out|std::ios::binary);
            assert(m_ofile.good());
            m_ofile.close();
            m_ofile.open(m_filename.c_str(), std::ios::in|std::ios::out|std::ios::binary);
            m_ifile.open(m_filename.c_str(), std::ios::in|std::ios::binary);
            assert(m_ifile.good());
            assert(m_ofile.good());
            // reset member variables
            m_need_to_write = false;
            m_max_elements = 0;
            // reset buffer
            read_block(0);
        }

        uint64_t read(const uint64_t idx) {
            assert(!m_closed);
            assert(idx < m_max_elements);
            if (idx < m_begin or m_begin+m_buffersize <= idx) {
                write_block();
                read_block(idx);
            }
            return m_buffer[idx-m_begin];
        }

        void write(const uint64_t idx, const uint64_t value) {
            assert(!m_closed);
            // If idx is not in current block, write current block and load needed block
            if (idx < m_begin or m_begin+m_buffersize <= idx) {
                write_block();
                read_block(idx);
            }
            if (m_max_elements <= idx) {
                m_max_elements = idx+1;
            }
            m_need_to_write = true;
            m_buffer[idx-m_begin] = value;
        }

        class int_vector_buffer_reference;

        int_vector_buffer_reference operator[](uint64_t idx) {
            return int_vector_buffer_reference(this, idx);
        }

        void push_back(const uint64_t value) {
            write(m_max_elements, value);
        }

        void close() {
            if (!m_closed) {
                m_ifile.close();
                assert(m_ifile.good());
                write_block();
                uint64_t size = m_max_elements*width();
                m_ofile.seekp(0, ios::beg);
//                int_vector_trait<t_width>::write_header(size, width, m_file);
                assert(m_ofile.good());
                m_ofile.write((char*) &size, sizeof(size));
                assert(m_ofile.good());
                if (m_offset == 9) {
                    uint8_t width = this->width();
                    m_ofile.write((char*) &width, sizeof(width));
                    assert(m_ofile.good());
                }
                uint64_t wb = (size+7)/8;
                if (wb%8) {
                    m_ofile.seekp(m_offset+wb);
                    assert(m_ofile.good());
                    m_ofile.write("\0\0\0\0\0\0\0\0", 8-wb%8);
                    assert(m_ofile.good());
                }
                m_ofile.close();
                assert(m_ofile.good());
                m_closed = true;
            }
        }

        //old swap
        void swap(int_vector_buffer<t_width>& ivb) {
            if (this != &ivb) {
                m_ifile.close();
                ivb.m_ifile.close();
                m_ofile.close();
                ivb.m_ofile.close();
                std::swap(m_filename, ivb.m_filename);
                m_ifile.open(m_filename.c_str(), std::ios::in|std::ios::binary);
                assert(m_ifile.good());
                m_ofile.open(m_filename.c_str(), std::ios::in|std::ios::out|std::ios::binary);
                assert(m_ofile.good());
                ivb.m_ifile.open(ivb.m_filename.c_str(), std::ios::in|std::ios::binary);
                assert(ivb.m_ifile.good());
                ivb.m_ofile.open(ivb.m_filename.c_str(), std::ios::in|std::ios::out|std::ios::binary);
                assert(ivb.m_ofile.good());
                std::swap(m_buffer, ivb.m_buffer);
                std::swap(m_need_to_write, ivb.m_need_to_write);
                std::swap(m_buffersize, ivb.m_buffersize);
                std::swap(m_max_elements, ivb.m_max_elements);
                std::swap(m_begin, ivb.m_begin);
                std::swap(m_persistent, ivb.m_persistent);
                std::swap(m_closed, ivb.m_closed);
            }
        }

        //new swap ??????
// 		void swap(int_vector_buffer<t_width>& ivb)
// 		{
// 			if (this != &ivb) { // if ivb and _this_ are not the same object
//                 m_file.close();
//                 ivb.m_file.close();
// 				m_filename = ivb.m_filename;
//                 m_file.open(m_filename.c_str(), std::ios::in|std::ios::out|std::ios::binary);
//                 assert(m_file.good());
//                 ivb.m_file.open(ivb.m_filename.c_str(), std::ios::in|std::ios::out|std::ios::binary);
//                 assert(ivb.m_file.good());
//                 m_buffer = ivb.m_buffer;
//                 m_need_to_write = ivb.m_need_to_write;
//                 m_buffersize = ivb.m_buffersize;
//                 m_max_elements = ivb.m_max_elements;
//                 m_begin = ivb.m_begin;
//                 m_persistent = ivb.m_persistent;
//                 m_closed = ivb.m_closed;
// 			}
// 		}


        class int_vector_buffer_reference
        {
                friend class int_vector_buffer<t_width>;
            private:
                int_vector_buffer<t_width>* const m_int_vector_buffer;
                uint64_t m_idx;

                int_vector_buffer_reference() :
                    m_int_vector_buffer(NULL), m_idx(0) {}

                int_vector_buffer_reference(int_vector_buffer<t_width>* _int_vector_buffer, uint64_t _idx) :
                    m_int_vector_buffer(_int_vector_buffer), m_idx(_idx) {}

            public:

                // conversion to int for read operations
                operator uint64_t () {
                    return m_int_vector_buffer->read(m_idx);
                }

                // assignment operator for write operations
                int_vector_buffer_reference& operator = (const uint64_t& val)     {
                    m_int_vector_buffer->write(m_idx, val);
                    return *this;
                }

                //! Prefix increment of the proxy object
                int_vector_buffer_reference& operator++() {
                    uint64_t x = m_int_vector_buffer->read(m_idx);
                    m_int_vector_buffer->write(m_idx, x+1);
                    return *this;
                }

                //! Postfix increment of the proxy object
                uint64_t operator++(int) {
                    uint64_t val = (uint64_t)*this;
                    ++(*this);
                    return val;
                }

                //! Prefix decrement of the proxy object
                int_vector_buffer_reference& operator--() {
                    uint64_t x = m_int_vector_buffer->read(m_idx);
                    m_int_vector_buffer->write(m_idx, x-1);
                    return *this;
                }

                //! Postfix decrement of the proxy object
                uint64_t operator--(int) {
                    uint64_t val = (uint64_t)*this;
                    --(*this);
                    return val;
                }

                //! Add assign from the proxy object
                int_vector_buffer_reference& operator+=(const uint64_t x) {
                    uint64_t w = m_int_vector_buffer->read(m_idx);
                    m_int_vector_buffer->write(m_idx, w+x);
                    return *this;
                }

                //! Subtract assign from the proxy object
                int_vector_buffer_reference& operator-=(const uint64_t x) {
                    uint64_t w = m_int_vector_buffer->read(m_idx);
                    m_int_vector_buffer->write(m_idx, w-x);
                    return *this;
                }

                bool operator==(const int_vector_buffer_reference& x)const {
                    return uint64_t(*this) == uint64_t(x);
                }

                bool operator<(const int_vector_buffer_reference& x)const {
                    return uint64_t(*this) < uint64_t(x);
                }
        };
};

// maybe "delete" for the following:
template<>
int_vector_buffer<1>::int_vector_buffer_reference& int_vector_buffer<1>::int_vector_buffer_reference::operator++()
{
    uint64_t x = m_int_vector_buffer->read(m_idx);
    m_int_vector_buffer->write(m_idx, (x+1)&1ULL);
    return *this;
}

template<>
uint64_t int_vector_buffer<1>::int_vector_buffer_reference::operator++(int)
{
    uint64_t val = (uint64_t)*this;
    ++(*this);
    return val;
}

template<>
int_vector_buffer<1>::int_vector_buffer_reference& int_vector_buffer<1>::int_vector_buffer_reference::operator-=(const uint64_t x)
{
    uint64_t w = m_int_vector_buffer->read(m_idx);
    m_int_vector_buffer->write(m_idx, (w-x)&1ULL);
    return *this;
}

template<>
int_vector_buffer<1>::int_vector_buffer_reference& int_vector_buffer<1>::int_vector_buffer_reference::operator+=(const uint64_t x)
{
    uint64_t w = m_int_vector_buffer->read(m_idx);
    m_int_vector_buffer->write(m_idx, (w+x)&1ULL);
    return *this;
}

} // end of namespace

#endif // include guard
