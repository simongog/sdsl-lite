#include "sdsl/int_vector.hpp"

namespace sdsl
{

int_vector_size_type char_array_serialize_wrapper::serialize(std::ostream& out) const
{
    size_type size = m_n*8; // number of bits
    size_type written_bytes = 0;
    written_bytes += _sdsl_serialize_size_and_int_width(out, 8, 8, size);
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
void bit_vector::flip()
{
    if (!empty()) {
        for (uint64_t i=0; i<(capacity()>>6); ++i) {
            m_data[i] = ~m_data[i];
        }
    }
}

}
