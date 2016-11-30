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
/*! \file int_vector_buffer.hpp
    \brief int_vector_buffer.hpp contains the sdsl::int_vector_buffer class.
    \author Maike Zwerger, Timo Beller and Simon Gog
*/
#ifndef INCLUDED_INT_VECTOR_BUFFER
#define INCLUDED_INT_VECTOR_BUFFER

#include "int_vector.hpp"
#include "iterators.hpp"
#include <cassert>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>

namespace sdsl
{

template<uint8_t t_width=0>
class int_vector_buffer
{
    public:
        class iterator;
        typedef typename int_vector<t_width>::difference_type difference_type;
        typedef typename int_vector<t_width>::value_type      value_type;

    private:
        static_assert(t_width <= 64 , "int_vector_buffer: width must be at most 64 bits.");
        sdsl::isfstream     m_ifile;
        sdsl::osfstream     m_ofile;
        std::string         m_filename;
        int_vector<t_width> m_buffer;
        bool                m_need_to_write = false;
        // length of int_vector header in bytes: 0 for plain, 8 for int_vector<t_width> (0 < t_width), 9 for int_vector<0>
        uint64_t            m_offset     = 0;
        uint64_t            m_buffersize = 8;    // in elements! m_buffersize*width() must be a multiple of 8!
        uint64_t            m_size       = 0;    // size of int_vector_buffer
        uint64_t            m_begin      = 0;    // number in elements

        //! Read block containing element at index idx.
        void read_block(const uint64_t idx)
        {
            m_begin = (idx/m_buffersize)*m_buffersize;
            if (m_begin >= m_size) {
                util::set_to_value(m_buffer, 0);
            } else {
                m_ifile.seekg(m_offset+(m_begin*width())/8);
                assert(m_ifile.good());
                m_ifile.read((char*) m_buffer.data(), (m_buffersize*width())/8);
                if ((uint64_t)m_ifile.gcount() < (m_buffersize*width())/8) {
                    m_ifile.clear();
                }
                assert(m_ifile.good());
                for (uint64_t i=m_size-m_begin; i<m_buffersize; ++i) {
                    m_buffer[i] = 0;
                }
            }
        }

        //! Write current block to file.
        void write_block()
        {
            if (m_need_to_write) {
                m_ofile.seekp(m_offset+(m_begin*width())/8);
                assert(m_ofile.good());
                if (m_begin+m_buffersize >= m_size) {
                    //last block in file
                    uint64_t wb = ((m_size-m_begin)*width()+7)/8;
                    m_ofile.write((char*) m_buffer.data(), wb);
                } else {
                    m_ofile.write((char*) m_buffer.data(), (m_buffersize*width())/8);
                }
                m_ofile.flush();
                assert(m_ofile.good());
                m_need_to_write = false;
            }
        }

        //! Read value from idx.
        uint64_t read(const uint64_t idx)
        {
            assert(is_open());
            assert(idx < m_size);
            if (idx < m_begin or m_begin+m_buffersize <= idx) {
                write_block();
                read_block(idx);
            }
            return m_buffer[idx-m_begin];
        }

        //! Write value to idx.
        void write(const uint64_t idx, const uint64_t value)
        {
            assert(is_open());
            // If idx is not in current block, write current block and load needed block
            if (idx < m_begin or m_begin+m_buffersize <= idx) {
                write_block();
                read_block(idx);
            }
            if (m_size <= idx) {
                m_size = idx+1;
            }
            m_need_to_write = true;
            m_buffer[idx-m_begin] = value;
        }

    public:

        //! Constructor.
        int_vector_buffer()
        {
            m_buffer = int_vector<t_width>();
        }

        //! Constructor for int_vector_buffer.
        /*! \param filename   File that contains the data read from / written to.
         *  \param mode       Openmode:
         *                    std::ios::in opens an existing file (that must exist already),
         *                    std::ios::out creates a new file (that may exist already).
         *  \param buffersize Buffersize in bytes. This has to be a multiple of 8, if not the next multiple of 8 will be taken
         *  \param int_width  The width of each integer.
         *  \param is_plain   If false (default) the file will be interpreted as int_vector.
         *                    If true the file will be interpreted as plain array with t_width bits per integer.
         *                    In second case (is_plain==true), t_width must be 8, 16, 32 or 64.
         */
        int_vector_buffer(const std::string filename, std::ios::openmode mode=std::ios::in, const uint64_t buffer_size=1024*1024, const uint8_t int_width=t_width, const bool is_plain=false)
        {
            m_filename = filename;
            assert(!(mode&std::ios::app));
            mode &= ~std::ios::app;
            m_buffer.width(int_width);
            if (is_plain) {
                // is_plain is only allowed with width() in {8, 16, 32, 64}
                assert(8==width() or 16==width() or 32==width() or 64==width());
            } else {
                m_offset = t_width ? 8 : 9;
            }

            // Open file for IO
            m_ofile.open(m_filename, mode|std::ios::out|std::ios::binary);
            assert(m_ofile.good());
            m_ifile.open(m_filename, std::ios::in|std::ios::binary);
            assert(m_ifile.good());
            if (mode & std::ios::in) {
                uint64_t size  = 0;
                if (is_plain) {
                    m_ifile.seekg(0, std::ios_base::end);
                    size = m_ifile.tellg()*8;
                } else {
                    uint8_t width = 0;
                    int_vector<t_width>::read_header(size, width, m_ifile);
                    m_buffer.width(width);
                }
                assert(m_ifile.good());
                m_size = size/width();
            }
            buffersize(buffer_size);
        }

        //! Move constructor.
        int_vector_buffer(int_vector_buffer&& ivb) :
            m_filename(std::move(ivb.m_filename)),
            m_buffer(std::move(ivb.m_buffer)),
            m_need_to_write(ivb.m_need_to_write),
            m_offset(ivb.m_offset),
            m_buffersize(ivb.m_buffersize),
            m_size(ivb.m_size),
            m_begin(ivb.m_begin)
        {
            ivb.m_ifile.close();
            ivb.m_ofile.close();
            m_ifile.open(m_filename, std::ios::in|std::ios::binary);
            m_ofile.open(m_filename, std::ios::in|std::ios::out|std::ios::binary);
            assert(m_ifile.good());
            assert(m_ofile.good());
            // set ivb to default-constructor state
            ivb.m_filename = "";
            ivb.m_buffer = int_vector<t_width>();
            ivb.m_need_to_write = false;
            ivb.m_offset = 0;
            ivb.m_buffersize = 8;
            ivb.m_size = 0;
            ivb.m_begin = 0;
        }

        //! Destructor.
        ~int_vector_buffer()
        {
            close();
        }

        //! Move assignment operator.
        int_vector_buffer<t_width>& operator=(int_vector_buffer&& ivb)
        {
            close();
            ivb.m_ifile.close();
            ivb.m_ofile.close();
            m_filename = ivb.m_filename;
            m_ifile.open(m_filename, std::ios::in|std::ios::binary);
            m_ofile.open(m_filename, std::ios::in|std::ios::out|std::ios::binary);
            assert(m_ifile.good());
            assert(m_ofile.good());
            // assign the values of ivb to this
            m_buffer = (int_vector<t_width>&&)ivb.m_buffer;
            m_need_to_write = ivb.m_need_to_write;
            m_offset = ivb.m_offset;
            m_buffersize = ivb.m_buffersize;
            m_size = ivb.m_size;
            m_begin = ivb.m_begin;
            // set ivb to default-constructor state
            ivb.m_filename = "";
            ivb.m_buffer = int_vector<t_width>();
            ivb.m_need_to_write = false;
            ivb.m_offset = 0;
            ivb.m_buffersize = 8;
            ivb.m_size = 0;
            ivb.m_begin = 0;
            return *this;
        }

        //! Returns the width of the integers which are accessed via the [] operator.
        uint8_t width() const
        {
            return m_buffer.width();
        }

        //! Returns the number of elements currently stored.
        uint64_t size() const
        {
            return m_size;
        }

        //! Returns the filename.
        std::string filename() const
        {
            return m_filename;
        }

        //! Returns the buffersize in bytes
        uint64_t buffersize() const
        {
            assert(m_buffersize*width()%8==0);
            return (m_buffersize*width())/8;
        }

        //! Set the buffersize in bytes
        void buffersize(uint64_t buffersize)
        {
            if (0ULL == buffersize)
                buffersize = 8;
            write_block();
            if (0==(buffersize*8)%width()) {
                m_buffersize = buffersize*8/width(); // m_buffersize might not be multiple of 8, but m_buffersize*width() is.
            } else {
                uint64_t element_buffersize = (buffersize*8)/width()+1; // one more element than fits into given buffersize in byte
                m_buffersize = element_buffersize+7 - (element_buffersize+7)%8; // take next multiple of 8
            }
            m_buffer = int_vector<t_width>(m_buffersize, 0, width());
            if (0!=m_buffersize) read_block(0);
        }

        //! Returns whether state of underlying streams are good
        bool good()
        {
            return m_ifile.good() and m_ofile.good();
        }

        //! Returns whether underlying streams are currently associated to a file
        bool is_open()
        {
            return m_ifile.is_open() and m_ofile.is_open();;
        }

        //! Delete all content and set size to 0
        void reset()
        {
            // reset file
            assert(m_ifile.good());
            assert(m_ofile.good());
            m_ifile.close();
            m_ofile.close();
            m_ofile.open(m_filename, std::ios::out|std::ios::binary);
            assert(m_ofile.good());
            m_ifile.open(m_filename, std::ios::in|std::ios::binary);
            assert(m_ifile.good());
            assert(m_ofile.good());
            // reset member variables
            m_need_to_write = false;
            m_size = 0;
            // reset buffer
            read_block(0);
        }

        // Forward declaration
        class reference;

        //! [] operator
        /*! \param i Index the i-th integer of length width().
         *  \return A reference to the i-th integer of length width().
         */
        reference operator[](uint64_t idx)
        {
            return reference(this, idx);
        }

        //! Appends the given element value to the end of the int_vector_buffer
        void push_back(const uint64_t value)
        {
            write(m_size, value);
        }

        //! Close the int_vector_buffer.
        /*! It is not possible to read from / write into the int_vector_buffer after calling this method
         *  \param remove_file If true, the underlying file will be removed on closing.
         */
        void close(bool remove_file=false)
        {
            if (is_open()) {
                if (!remove_file) {
                    write_block();
                    if (0 < m_offset) { // in case of int_vector, write header and trailing zeros
                        uint64_t size = m_size*width();
                        m_ofile.seekp(0, std::ios::beg);
                        int_vector<t_width>::write_header(size, width(), m_ofile);
                        assert(m_ofile.good());
                        uint64_t wb = (size+7)/8;
                        if (wb%8) {
                            m_ofile.seekp(m_offset+wb);
                            assert(m_ofile.good());
                            m_ofile.write("\0\0\0\0\0\0\0\0", 8-wb%8);
                            assert(m_ofile.good());
                        }
                    }
                }
                m_ifile.close();
                assert(m_ifile.good());
                m_ofile.close();
                assert(m_ofile.good());
                if (remove_file) {
                    sdsl::remove(m_filename);
                }
            }
        }

        iterator begin()
        {
            return iterator(*this, 0);
        }

        iterator end()
        {
            return iterator(*this, size());
        }

        //! Swap method for int_vector_buffer.
        void swap(int_vector_buffer<t_width>& ivb)
        {
            if (this != &ivb) {
                m_ifile.close();
                ivb.m_ifile.close();
                m_ofile.close();
                ivb.m_ofile.close();
                std::swap(m_filename, ivb.m_filename);
                m_ifile.open(m_filename, std::ios::in|std::ios::binary);
                assert(m_ifile.good());
                m_ofile.open(m_filename, std::ios::in|std::ios::out|std::ios::binary);
                assert(m_ofile.good());
                ivb.m_ifile.open(ivb.m_filename, std::ios::in|std::ios::binary);
                assert(ivb.m_ifile.good());
                ivb.m_ofile.open(ivb.m_filename, std::ios::in|std::ios::out|std::ios::binary);
                assert(ivb.m_ofile.good());
                std::swap(m_buffer, ivb.m_buffer);
                std::swap(m_need_to_write, ivb.m_need_to_write);
                std::swap(m_offset, ivb.m_offset);
                std::swap(m_buffersize, ivb.m_buffersize);
                std::swap(m_size, ivb.m_size);
                std::swap(m_begin, ivb.m_begin);
            }
        }

        class reference
        {
                friend class int_vector_buffer<t_width>;
            private:
                int_vector_buffer<t_width>* const m_int_vector_buffer = nullptr;
                uint64_t m_idx = 0;

                reference() {}

                reference(int_vector_buffer<t_width>* _int_vector_buffer, uint64_t _idx) :
                    m_int_vector_buffer(_int_vector_buffer), m_idx(_idx) {}

            public:

                //! Conversion to int for read operations
                operator uint64_t ()const
                {
                    return m_int_vector_buffer->read(m_idx);
                }

                //! Assignment operator for write operations
                reference& operator=(const uint64_t& val)
                {
                    m_int_vector_buffer->write(m_idx, val);
                    return *this;
                }

                //! Assignment operator
                reference& operator=(reference& x)
                {
                    return *this = (uint64_t)(x);
                };

                //! Prefix increment of the proxy object
                reference& operator++()
                {
                    uint64_t x = m_int_vector_buffer->read(m_idx);
                    m_int_vector_buffer->write(m_idx, x+1);
                    return *this;
                }

                //! Postfix increment of the proxy object
                uint64_t operator++(int)
                {
                    uint64_t val = (uint64_t)*this;
                    ++(*this);
                    return val;
                }

                //! Prefix decrement of the proxy object
                reference& operator--()
                {
                    uint64_t x = m_int_vector_buffer->read(m_idx);
                    m_int_vector_buffer->write(m_idx, x-1);
                    return *this;
                }

                //! Postfix decrement of the proxy object
                uint64_t operator--(int)
                {
                    uint64_t val = (uint64_t)*this;
                    --(*this);
                    return val;
                }

                //! Add assign from the proxy object
                reference& operator+=(const uint64_t x)
                {
                    uint64_t w = m_int_vector_buffer->read(m_idx);
                    m_int_vector_buffer->write(m_idx, w+x);
                    return *this;
                }

                //! Subtract assign from the proxy object
                reference& operator-=(const uint64_t x)
                {
                    uint64_t w = m_int_vector_buffer->read(m_idx);
                    m_int_vector_buffer->write(m_idx, w-x);
                    return *this;
                }

                bool operator==(const reference& x)const
                {
                    return (uint64_t)*this == (uint64_t)x;
                }

                bool operator<(const reference& x)const
                {
                    return (uint64_t)*this < (uint64_t)x;
                }
        };

        class iterator: public std::iterator<std::random_access_iterator_tag, value_type, difference_type, value_type*, reference>
        {
            private:
                int_vector_buffer<t_width>& m_ivb;
                uint64_t m_idx = 0;
            public:

                iterator() = delete;
                iterator(int_vector_buffer<t_width>& ivb, uint64_t idx=0) : m_ivb(ivb), m_idx(idx) {}

                iterator& operator++()
                {
                    ++m_idx;
                    return *this;
                }

                iterator operator++(int)
                {
                    iterator it = *this;
                    ++(*this);
                    return it;
                }

                iterator& operator--()
                {
                    --m_idx;
                    return *this;
                }

                iterator operator--(int)
                {
                    iterator it = *this;
                    --(*this);
                    return it;
                }

                reference operator*()const
                {
                    return m_ivb[m_idx];
                }

                iterator& operator+=(difference_type i)
                {
                    if (i<0)
                        return *this -= (-i);
                    m_idx += i;
                    return *this;
                }

                iterator& operator-=(difference_type i)
                {
                    if (i<0)
                        return *this += (-i);
                    m_idx -= i;
                    return *this;
                }

                iterator operator+(difference_type i) const
                {
                    iterator it = *this;
                    return it += i;
                }

                iterator& operator-(difference_type i) const
                {
                    iterator it = *this;
                    return it -= i;
                }

                bool operator==(const iterator& it) const
                {
                    return &m_ivb == &(it.m_ivb) and m_idx == it.m_idx;
                }

                bool operator!=(const iterator& it) const
                {
                    return !(*this == it);
                }
                inline difference_type operator-(const iterator& it)
                {
                    return (m_idx - it.m_idx);
                }
        };
};

} // end of namespace

#endif // include guard
