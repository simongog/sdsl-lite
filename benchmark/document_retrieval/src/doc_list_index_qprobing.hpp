/*!
 * This file contains a document listing class, which implements
 * strategy QUANTILE in the article:
 * J. S. Culpepperi, G. Navarro, S. J. Puglisi and A. Turpin:
 * ,,Top-k Ranked Document Search in General Text Databases''
 * Proceedings Part II of the 18th Annual European Symposium on
 * Algorithms (ESA 2010)
 */
#ifndef DOCUMENT_LISING_QPROBING_INCLUDED
#define DOCUMENT_LISING_QPROBING_INCLUDED

#include "doc_list_index_greedy.hpp"
#include <string>
#include <vector>
#include <fstream>
#include <list>
#include <utility>

using std::vector;

namespace sdsl
{

template<
class t_csa          = csa_wt<wt_huff<rrr_vector<63>>, 1000000, 1000000>,
      class t_wtd    = wt_int<bit_vector,rank_support_v5<1>,select_support_scan<1>,select_support_scan<0>>,
      typename t_csa::char_type t_doc_delim = 1
      >
class doc_list_index_qprobing : public doc_list_index_greedy<t_csa, t_wtd, t_doc_delim>
{
    private:
        using base_type = doc_list_index_greedy<t_csa, t_wtd, t_doc_delim>;
        using base_type::m_csa;
        using base_type::m_wtd;

    public:
        using size_type = typename base_type::size_type;
        using result = typename base_type::result;

        doc_list_index_qprobing() : base_type() {}
        doc_list_index_qprobing(std::string file_name, sdsl::cache_config& cconfig, uint8_t num_bytes) : base_type(file_name, cconfig, num_bytes) {}

        //! Search for the k documents which contains the search term most frequent
        size_type search(std::string::iterator begin, std::string::iterator end, result& res, size_t k) const {
            size_type sp=1, ep=0;
            if (0 == backward_search(m_csa, 0, m_csa.size()-1, begin, end, sp, ep)) {
                res = result();
                return 0;
            } else {
                auto tmp_res = m_wtd.topk_qprobing(sp,ep,k);
                res = result(sp, ep, std::move(tmp_res));
                return ep-sp+1;
            }
        }

};

} // end namespace

#endif
