/*!
 * This file contains a document listing class, which implements
 * strategy QUANTILE in the article:
 * J. S. Culpepper, G. Navarro, S. J. Puglisi and A. Turpin:
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
          using base_type::m_csa_full;
          using base_type::m_wtd;

          public:
          using size_type = typename base_type::size_type;
          using value_type = typename t_wtd::value_type;
          using result = typename base_type::result;

doc_list_index_qprobing() : base_type() {}
doc_list_index_qprobing(std::string file_name, sdsl::cache_config& cconfig, uint8_t num_bytes) : base_type(file_name, cconfig, num_bytes) {}

//! Search for the k documents which contains the search term most frequent
template<class t_pat_iter>
size_type search(t_pat_iter begin, t_pat_iter end, result& res, size_t k) const {
    size_type sp=1, ep=0;
    if (0 == backward_search(m_csa_full, 0, m_csa_full.size()-1, begin, end, sp, ep)) {
        res = result();
        return 0;
    } else {
        auto tmp_res = topk_qprobing(sp, ep ,k);
        res = result(sp, ep, std::move(tmp_res));
        return ep-sp+1;
    }
}

//! Returns the top-k most frequent documents in m_wtd[lb..rb]
/*!
 *  \param lb left array bound in T
 *  \param rb right array bound in T
 *  \param k the number of documents to return
 *  \returns the top-k items in ascending order.
 */
std::vector< std::pair<value_type,size_type> >
topk_qprobing(size_type lb, size_type rb,size_type k) const {
    using p_t = std::pair<value_type,size_type>;
    std::vector<p_t> results;
    auto comp = [](p_t& a,p_t& b) { return a.second > b.second; };
    std::priority_queue<p_t,std::vector<p_t>,decltype(comp)> heap(comp);
    bit_vector seen(1ULL << m_wtd.max_level); // TODO: better idea?

    /* we start probing using the largest power smaller than len */
    size_type len = rb-lb+1;
    size_type power2greaterlen = 1 << (bits::hi(len)+1);
    size_type probe_interval = power2greaterlen >> 1;

    /* we probe the smallest elem (pos 0 in sorted array) only once */
    auto qf = quantile_freq(m_wtd,lb,rb,0);
    heap.push(qf);
    seen[qf.first] = 1;

    qf = quantile_freq(m_wtd,lb,rb,probe_interval);
    if (!seen[qf.first]) heap.push(qf);
    seen[qf.first] = 1;

    while (probe_interval > 1) {
        size_type probe_pos = probe_interval >> 1;
        while (probe_pos < len) {
            qf = quantile_freq(m_wtd,lb,rb,probe_pos);
            if (!seen[qf.first]) { /* not in heap */
                if (heap.size()<k) {
                    heap.push(qf);
                    seen[qf.first] = 1;
                } else {
                    /* throw out the smallest and add the new one */
                    if (heap.top().second < qf.second) {
                        heap.pop();
                        heap.push(qf);
                        seen[qf.first] = 1;
                    }
                }
            }
            probe_pos += probe_interval;
        }
        probe_interval >>= 1;
        /* we have enough or can't find anything better */
        if (heap.size() == k && probe_interval-1 <= heap.top().second) break;
    }
    /* populate results */
    while (!heap.empty())  {
        results.emplace(results.begin() , heap.top());
        heap.pop();
    }
    return results;
};



      };

} // end namespace

#endif
