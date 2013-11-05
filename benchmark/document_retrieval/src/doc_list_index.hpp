#ifndef DOC_LIST_INDEX
#define DOC_LIST_INDEX

#include <sdsl/construct.hpp>
#include <string>

struct doc_list_tag {};

template<class t_index>
void
construct(t_index& idx, const std::string& file, sdsl::cache_config& config, uint8_t num_bytes, doc_list_tag)
{
    t_index tmp_idx(file, config, num_bytes);
    idx.swap(tmp_idx);
}

#include "doc_list_index_sada.hpp"
#include "doc_list_index_greedy.hpp"
#include "doc_list_index_qprobing.hpp"
#include "doc_list_index_sort.hpp"

#endif
