#ifndef DOC_RANK_INDEX
#define DOC_RANK_INDEX

#include <sdsl/construct.hpp>
#include <string>

struct doc_list_tag {};

const std::string KEY_DOCWEIGHT = "docweights";
const std::string KEY_DARRAY = "darray";

template<class t_index>
void
construct(t_index& idx, const std::string& file, sdsl::cache_config& config, uint8_t num_bytes, doc_list_tag)
{
    t_index tmp_idx(file, config, num_bytes);
    idx.swap(tmp_idx);
}

using list_type = std::vector<std::pair<uint64_t,double>>;

class result : public list_type
{
private:
    uint64_t m_range_size;
public:
    uint64_t total_range_sizes() {
        return m_range_size;
    }
    // Constructors for an empty result and for a result in the interval [sp, ep]:
    result(uint64_t rs,list_type&& l) : list_type(l) , m_range_size(rs) {}
    result() : m_range_size(0) {}
    result(uint64_t rs) : m_range_size(rs) {}
    result& operator=(const result& res) {
        if (this != &res) {
            list_type::operator=(res);
            m_range_size = res.m_range_size;
        }
        return *this;
    }
};

#include "doc_rank_greedy.hpp"

#endif
