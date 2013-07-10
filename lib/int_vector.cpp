#include "sdsl/int_vector.hpp"

namespace sdsl
{

int_vector_size_type char_array_serialize_wrapper::serialize(std::ostream& out) const
{
    size_type size = m_n*8; // number of bits
    size_type written_bytes = int_vector<8>::write_header(8, size, out);
    const char* p = (const char*)m_cp;
    size_type idx = 0;
    while (idx+constants::SDSL_BLOCK_SIZE < m_n) {
        out.write(p, constants::SDSL_BLOCK_SIZE);
        written_bytes += constants::SDSL_BLOCK_SIZE;
        p += constants::SDSL_BLOCK_SIZE;
        idx += constants::SDSL_BLOCK_SIZE;
    }
    // now: m_n-idx <= SDSL_BLOCK_SIZE
    out.write(p , m_n-idx);
    written_bytes += m_n-idx;
    uint32_t r= (sizeof(uint64_t)-(m_n % sizeof(uint64_t))) % sizeof(uint64_t);
    out.write("\0\0\0\0\0\0\0\0", r);
    written_bytes += r;
    return written_bytes;
}


template<>
void int_vector<0>::width(const uint8_t new_int_width)
{
    if (0 < new_int_width and new_int_width <= 64)
        m_width = new_int_width;
    else
        m_width = 64;
}

template<>
void bit_vector::flip()
{
    if (!empty()) {
        for (uint64_t i=0; i<(capacity()>>6); ++i) {
            m_data[i] = ~m_data[i];
        }
    }
}

template<>
void int_vector_file_buffer<0>::set_width(uint8_t int_width)
{
    if (0 < int_width and int_width <= 64)
        m_width = int_width;
    else
        m_width = 64;
}

}
