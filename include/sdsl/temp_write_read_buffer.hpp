/*! \file temp_write_read_buffer.hpp contains a helper class
 * for external algorithm implementations
 */

#ifndef INCLUDED_SDSL_TEMP_WRITE_READ_BUFFER
#define INCLUDED_SDSL_TEMP_WRITE_READ_BUFFER
#include "int_vector.hpp"
#include <string>
#include <fstream>

namespace sdsl
{

template<uint8_t int_width=0>
class temp_write_read_buffer
{
    public:
        typedef int_vector<int_width> buffer_type;
        typedef typename buffer_type::size_type size_type;
        typedef typename buffer_type::value_type value_type;

    private:

        buffer_type 	m_buf; // buffer for the data
        size_type   	m_buf_size; // size of the buffer
        std::string 	m_file_name; // filename for the file holding
        // the data which does not fit in the buffer
        size_type 		m_in_buf_idx;// current index in the buffer
        size_type		m_buf_idx;// index of the current buffer
        size_type		m_buf_cnt;// number of buffers written to disk
        std::ofstream	m_out;		 // file out stream
        std::ifstream	m_in;		 // file in stream
        bool			m_output_exists; // if there exists output to the file
        size_type		m_r;// remaining entries in the buffer
        size_type		m_last_block_size; // size of the last block written to disk

    public:

        //! Constructior
        /*! \param buf_size The size of the buffer.
         *	\param width	The width of the integers if template paramter int_width=0.
         *	\param dir		Directory in which the temporary file is stored.
         * */
        temp_write_read_buffer(size_type buf_size, uint8_t width, std::string dir="./") {
            m_buf_size = buf_size;
            m_buf = buffer_type(buf_size, 0, width);    // initialize buffer
            m_in_buf_idx = 0;
            m_buf_cnt = 0;
            m_file_name =  dir + "temp_write_read_buffer_" + util::to_string(util::get_pid())+"_"+util::to_string(util::get_id()).c_str();
            m_output_exists = false;
        }

        // Destructor
        ~temp_write_read_buffer() {
            if (m_out.is_open()) { // if out buffer is still open
                m_out.close();     // close it
            }
            if (m_in.is_open()) { // if in buffer is still open
                m_in.close();	  // close it
            }
            if (m_output_exists)  // if we have written output to a file
                std::remove(m_file_name.c_str());   // delete it
        }

        value_type operator<<(value_type x) {
            if (m_in_buf_idx == m_buf_size) {
                m_in_buf_idx = 0;
                ++m_buf_cnt; ++m_buf_idx;
                if (m_buf_cnt == 1) {
                    m_out.open(m_file_name.c_str(), std::ios::trunc | std::ios::out | std::ios::binary);   // open file buffer
                }
                m_buf.serialize(m_out);   // write buffer to disk
                m_output_exists = true;
            }
            m_buf[m_in_buf_idx++] = x;
            return x;
        }

        void write_close() {
            if (m_buf_cnt > 0) {
                m_buf.serialize(m_out);    // write last buffer to disk
                m_last_block_size = m_in_buf_idx;
                ++m_buf_cnt;
                m_out.close();
                m_r = 0;
                m_in_buf_idx = 0;
                m_buf_idx = 0;
                m_in.open(m_file_name.c_str(), std::ios::in | std::ios::binary);
            } else {
                m_r = m_in_buf_idx;
                m_in_buf_idx = 0;
                m_buf_idx = 0;
                m_last_block_size = 0;
            }
        }

        void reset() {
            m_in_buf_idx = 0;
            m_buf_cnt = 0;
            if (m_out.is_open()) {
                m_out.close();
            }
            if (m_in.is_open()) {
                m_in.close();
            }
        }

        bool operator>>(value_type& x) {
            if (m_in_buf_idx >= m_r) {  // load next buffer
//			std::cout<< "m_in_buf_idx = " << m_in_buf_idx << " m_r "<< m_r <<" m_buf_idx = "<< m_buf_idx<<" m_buf_cnt = "<< m_buf_cnt << std::endl;
                if (m_buf_idx < m_buf_cnt) {
                    ++m_buf_idx; // increase buffer index
                    m_in_buf_idx = 0; // reset in buffer index
                    m_buf.load(m_in);   // load next block
                    if (m_buf_idx == m_buf_cnt)
                        m_r = m_buf_idx;
                    else
                        m_r = m_buf_size;
                } else {
                    x = 0;
                    return false;
                }
            }
            x = m_buf[ m_in_buf_idx++ ];
            return true;
        }
};

} // end namespace sdsl

#endif


