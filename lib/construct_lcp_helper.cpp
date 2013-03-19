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
    uint8_t int_width = bit_magic::l1BP(max_lcp_value-1)+1;
    uint64_t bit_size = n*int_width;								// Size of output file
    size_type wb = 0;												// Number of bits that were already written
    int_vector<> out_buf(buffer_size, 0, int_width); 				// Output buffer
    osfstream tmp_lcp_out_buf(tmp_lcp_file, std::ios::binary | std::ios::trunc | std::ios::out);
    tmp_lcp_out_buf.write((char*) &(bit_size), sizeof(bit_size));		// Write length of vector
    tmp_lcp_out_buf.write((char*) &(int_width), sizeof(int_width));	// Write int-width of vector

    // Write values into buffer
    for (size_type i=0, r_sum=0, calc_idx=0, r=0; r_sum < n;) {
        // Copy next r values into buffer
        util::set_zero_bits(out_buf); // initialize buffer with zeros
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
    std::rename(tmp_lcp_file.c_str(), lcp_file.c_str());
}

} // end namespace sdsl
