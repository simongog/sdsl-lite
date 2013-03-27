#include "sdsl/construct_lcp_helper.hpp"
#include "sdsl/int_vector.hpp"
#include <algorithm>

namespace sdsl
{

//! Merges a partial LCP array into the LCP array on disk.
/*!
 * \param partial_lcp		Vector containing LCP values for all indexes \f$i\f$ with
 *                      	index_done[i] == 0. Let x=partail_lcp[rank(index_done, i, 0)];
 *                      	LCP[i]=x if x!=0 and index_done[i] == 0
 * \param lcp_file			Path to the LCP array on disk.
 * \param index_done		Entry index_done[i] indicates if LCP[i] is already calculated.
 * \param max_lcp_value 	Maximum known LCP value
 * \param lcp_value_offset	Largest LCP value in lcp_file
 */
void insert_lcp_values(int_vector<>& partial_lcp, bit_vector& index_done, std::string lcp_file, uint64_t max_lcp_value, uint64_t lcp_value_offset)
{
    std::string tmp_lcp_file  = lcp_file+"_TMP";
    const uint64_t buffer_size = 1000000; // has to be a multiple of 64
    typedef int_vector<>::size_type size_type;
    int_vector_file_buffer<> lcp_buffer(lcp_file, buffer_size); // open lcp_file
    uint64_t n = lcp_buffer.int_vector_size;

    // open tmp_lcp_file
    uint8_t int_width = bits::hi(max_lcp_value-1)+1;
    uint64_t bit_size = n*int_width;								// Size of output file
    size_type wb = 0;												// Number of bits that were already written
    int_vector<> out_buf(buffer_size, 0, int_width); 				// Output buffer
    osfstream tmp_lcp_out_buf(tmp_lcp_file, std::ios::binary | std::ios::trunc | std::ios::out);
    tmp_lcp_out_buf.write((char*) &(bit_size), sizeof(bit_size));		// Write length of vector
    tmp_lcp_out_buf.write((char*) &(int_width), sizeof(int_width));	// Write int-width of vector

    // Write values into buffer
    for (size_type i=0, r_sum=0, calc_idx=0, r=0; r_sum < n;) {
        // Copy next r values into buffer
        util::set_to_value(out_buf, 0); // initialize buffer with zeros
        for (; i < r_sum+r; ++i) {
            // If values was already calculated
            if (index_done[i]) {
                out_buf[i-r_sum] = lcp_buffer[i-r_sum]; // Copy value
            } else {
                if (partial_lcp[calc_idx]) {   // If values was now calculated
                    // Insert value
                    out_buf[i-r_sum] = partial_lcp[calc_idx]+lcp_value_offset;
                    index_done[i] = true;
                }
                ++calc_idx;
            }
        }
        // Write next r values from buffer to harddisk
        if (r>0) {
            size_type cur_wb = (r*out_buf.width()+7)/8;
            tmp_lcp_out_buf.write((const char*)out_buf.data(), cur_wb);
            wb += cur_wb;
        }
        // Count how many values were written and how many values will be written next
        r_sum += r;
        r = lcp_buffer.load_next_block();
    }
    // Close file and replace old file with new one
    if (wb%8) {
        tmp_lcp_out_buf.write("\0\0\0\0\0\0\0\0", 8-wb%8);
    }
    tmp_lcp_out_buf.close();
    sdsl::rename(tmp_lcp_file, lcp_file);
}

buffered_char_queue::buffered_char_queue():m_widx(0), m_ridx(0), m_sync(true), m_disk_buffered_blocks(0), m_c('?'),m_rb(0), m_wb(0) {};

void buffered_char_queue::init(const std::string& dir, char c)
{
    m_c = c;
    m_file_name = dir+"buffered_char_queue_"+util::to_string(util::pid());
//		m_stream.rdbuf()->pubsetbuf(0, 0);
}

buffered_char_queue::~buffered_char_queue()
{
    m_stream.close();
    std::remove(m_file_name.c_str());
}

void buffered_char_queue::push_back(uint8_t x)
{
    m_write_buf[m_widx] = x;
    if (m_sync) {
        m_read_buf[m_widx] = x;
    }
    ++m_widx;
    if (m_widx == m_buffer_size) {
        if (!m_sync) { // if not sync, write block to disk
            if (!m_stream.is_open()) {
                m_stream.open(m_file_name.c_str(), std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);
            }
            m_stream.seekp(m_buffer_size * (m_wb++), std::ios::beg);
            m_stream.write((char*) m_write_buf, m_buffer_size);
            ++m_disk_buffered_blocks;
        }
        m_sync = 0;
        m_widx = 0;
    }
}

uint8_t buffered_char_queue::pop_front()
{
    uint8_t x = m_read_buf[m_ridx];
    ++m_ridx;
    if (m_ridx ==  m_buffer_size) {
        if (m_disk_buffered_blocks > 0) {
            m_stream.seekg(m_buffer_size * (m_rb++), std::ios::beg);
            m_stream.read((char*) m_read_buf, m_buffer_size);
            --m_disk_buffered_blocks;
        } else { // m_disk_buffered_blocks == 0
            m_sync = 1;
            memcpy(m_read_buf, m_write_buf, m_widx+1);
        }
        m_ridx = 0;
    }
    return x;
}

void lcp_info(cache_config& config)
{
    typedef int_vector<>::size_type size_type;
    int_vector_file_buffer<> lcp_buf(util::cache_file_name(constants::KEY_LCP, config));
    size_type n = lcp_buf.int_vector_size;

    size_type max_lcp = 0;
    size_type sum_lcp = 0;
    for (size_type i=0, r_sum=0, r=0; i < n;) {
        for (; i < r_sum+r; ++i) {
            if (lcp_buf[i-r_sum] > max_lcp)
                max_lcp = lcp_buf[i-r_sum];
            sum_lcp += lcp_buf[i-r_sum];
        }
        r_sum += r; r = lcp_buf.load_next_block();
    }
    std::cout<<"# max lcp = " << max_lcp << std::endl;
    std::cout<<"# sum lcp = " << sum_lcp << std::endl;
    std::cout<<"# avg lcp = " << sum_lcp/(double)n << std::endl;
}

} // end namespace sdsl
