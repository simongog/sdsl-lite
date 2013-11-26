#ifndef DOC_RANK_INDEX
#define DOC_RANK_INDEX

#include <sdsl/construct.hpp>
#include <string>

struct doc_list_tag {};

const std::string KEY_DOCWEIGHT = "docweights";
const std::string KEY_DARRAY = "darray";
const std::string KEY_DOCPERM = "docperm";
const std::string KEY_SADADF = "sadadf";
const std::string KEY_SADADFSEL = "sadadfsel";
const std::string KEY_WTD = "wtd";
const std::string KEY_C = "C";
const std::string KEY_WTC = "wtc";
const std::string KEY_TMPCST = "tempcst";

template <class t_index>
void construct(t_index& idx, const std::string& file,
               sdsl::cache_config& config, uint8_t num_bytes, doc_list_tag)
{
    t_index tmp_idx(file, config, num_bytes);
    idx.swap(tmp_idx);
}

using list_type = std::vector<std::pair<uint64_t, double>>;

class result : public list_type
{
    public:
        uint64_t range_size = 0;
    public:
        uint64_t total_range_sizes() { return range_size; }
        // Constructors for an empty result and for a result in the interval [sp,
        // ep]:
        result(list_type&& l) : list_type(l), range_size(0) {}
        result() : range_size(0) {}
        result(uint64_t rs) : range_size(rs) {}
        result& operator=(const result& res) {
            if (this != &res) {
                list_type::operator=(res);
                range_size = res.range_size;
            }
            return *this;
        }
};

#include "doc_rank_greedy.hpp"

#endif
