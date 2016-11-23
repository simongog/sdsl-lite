#include <memory>

#include "sdsl/dac_vector.hpp"
#include "sdsl/rrr_vector.hpp"

namespace sdsl {

uint64_t dac_vector_level::get(size_t idx) const {
    assert(!m_overflow_plain.empty() || !m_overflow_rrr.empty());
    if (!m_overflow_plain.empty())
        return get_impl(m_overflow_plain, m_overflow_rank_plain, idx);
    else
        return get_impl(m_overflow_rrr, m_overflow_rank_rrr, idx);
}

}  // namespace sdsl
