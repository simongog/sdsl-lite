/*!
 * This file contains a document listing class, which implements
 * strategy GREEDY in the article:
 * J. S. Culpepper, G. Navarro, S. J. Puglisi and A. Turpin:
 * ,,Top-k Ranked Document Search in General Text Databases''
 * Proceedings Part II of the 18th Annual European Symposium on
 * Algorithms (ESA 2010)
 */
#ifndef DOCUMENT_LISING_GREEDY_INCLUDED
#define DOCUMENT_LISING_GREEDY_INCLUDED

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/rmq_support.hpp>
#include <string>
#include <vector>
#include <fstream>
#include <list>
#include <utility>
#include "doc_list_index.hpp"

using std::vector;

namespace sdsl
{

template<
class t_csa       = csa_wt<wt_huff<rrr_vector<63>>, 1000000, 1000000>,
      class t_wtd = wt_int<bit_vector,rank_support_v5<1>,select_support_scan<1>,select_support_scan<0>>,
      typename t_csa::char_type t_doc_delim = 1
      >
      class doc_list_index_greedy
      {
          public:
          using size_type = typename t_wtd::size_type;
          using value_type = typename t_wtd::value_type;

          typedef t_csa                                       csa_type;
          typedef t_wtd                                       wtd_type;
          typedef std::vector<std::pair<size_type,size_type>> list_type;
          typedef doc_list_tag                                index_category;

          enum { WIDTH = t_csa::alphabet_category::WIDTH };

          class result : public list_type
          {
              private:
              size_type m_sp, m_ep;
              public:
              // Number of occurrences
size_type count() {
    return m_ep-m_sp+1;
}

// Constructors for an empty result and for a result in the interval [sp, ep]:
result(size_type sp, size_type ep,list_type&& l) : list_type(l), m_sp(sp), m_ep(ep) {}
result() : m_sp(1), m_ep(0) {}
result(size_type sp, size_type ep) : m_sp(sp), m_ep(ep) {}
result& operator=(const result& res) {
    if (this != &res) {
        list_type::operator=(res);
        m_sp = res.m_sp;
        m_ep = res.m_ep;
    }
    return *this;
}

          };


struct wt_range_t {
    using node_type = typename wtd_type::node_type;

    node_type v;
    range_type r;

size_t size() const {
    return r.second - r.first + 1;
}

bool operator<(const wt_range_t& x) const {
    if (x.size() != size())
        return size() < x.size();
    return v.sym > x.v.sym;
}

wt_range_t() {}
wt_range_t(const node_type& _v, const range_type& _r):
    v(_v), r(_r) {}
};


protected:
size_type m_doc_cnt; // number of documents in the collection
csa_type  m_csa_full;     // CSA built from the collection text
wtd_type  m_wtd;     // wtd build from the collection text
public:

//! Default constructor
doc_list_index_greedy() { }

doc_list_index_greedy(std::string file_name, sdsl::cache_config& cconfig, uint8_t num_bytes) {
    construct(m_csa_full, file_name, cconfig, num_bytes);

    const char* KEY_TEXT = key_text_trait<WIDTH>::KEY_TEXT;
    std::string text_file = cache_file_name(KEY_TEXT, cconfig);

    bit_vector doc_border;
    construct_doc_border(text_file,doc_border);
    bit_vector::rank_1_type doc_border_rank(&doc_border);
    m_doc_cnt = doc_border_rank(doc_border.size());

    int_vector_buffer<0> sa_buf(cache_file_name(conf::KEY_SA, cconfig));
    {
        int_vector<> D;
        construct_D_array(sa_buf, doc_border_rank, m_doc_cnt, D);
        std::string d_file = cache_file_name("DARRAY", cconfig);
        store_to_file(D, d_file);
        util::clear(D);
        construct(m_wtd, d_file);
        sdsl::remove(d_file);
    }
}

size_type doc_cnt()const {
    return m_wtd.sigma-1; // subtract one, since zero does not count
}

size_type word_cnt()const {
    return m_wtd.size()-doc_cnt();
}

size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += write_member(m_doc_cnt, out, child, "doc_cnt");
    written_bytes += m_csa_full.serialize(out, child, "csa_full");
    written_bytes += m_wtd.serialize(out, child, "wtd");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void load(std::istream& in) {
    read_member(m_doc_cnt, in);
    m_csa_full.load(in);
    m_wtd.load(in);
}

void swap(doc_list_index_greedy& dr) {
    if (this != &dr) {
        std::swap(m_doc_cnt, dr.m_doc_cnt);
        m_csa_full.swap(dr.m_csa_full);
        m_wtd.swap(dr.m_wtd);
    }
}

//! Search for the k documents which contain the search term most frequent
template<class t_pat_iter>
size_type search(t_pat_iter begin, t_pat_iter end, result& res, size_t k) const {
    size_type sp=1, ep=0;
    if (0 == backward_search(m_csa_full, 0, m_csa_full.size()-1, begin, end, sp, ep)) {
        res = result();
        return 0;
    } else {
        auto tmp_res = topk_greedy(sp, ep, k);
        res = result(sp, ep, std::move(tmp_res));
        return ep-sp+1;
    }
}

private:
//! Construct the doc_border bitvector by streaming the text file
void
construct_doc_border(const std::string& text_file, bit_vector& doc_border) {
    int_vector_buffer<WIDTH> text_buf(text_file);
    doc_border = bit_vector(text_buf.size(), 0);
    for (size_type i = 0; i < text_buf.size(); ++i) {
        if (t_doc_delim == text_buf[i]) {
            doc_border[i] = 1;
        }
    }
}

void
construct_D_array(int_vector_buffer<0>& sa_buf,
                  bit_vector::rank_1_type& doc_border_rank,
                  const size_type doc_cnt,
int_vector<>& D) {
    D = int_vector<>(sa_buf.size(), 0, bits::hi(doc_cnt+1)+1);
    for (size_type i = 0; i < sa_buf.size(); ++i) {
        uint64_t d = doc_border_rank(sa_buf[i]+1);
        D[i] = d;
    }
}

//! Returns the top k most frequent documents in D[lb..rb]
/*!
 *  \param lb  Left array border in D.
 *  \param rb  Right array border in D.
 *  \param k   The number of documents to return.
 *  \returns the top-k items in ascending order.
 */
std::vector< std::pair<value_type,size_type> >
topk_greedy(size_type lb, size_type rb, size_type k) const {
    std::vector< std::pair<value_type,size_type> > results;
    std::priority_queue<wt_range_t> heap;

    heap.emplace(wt_range_t(m_wtd.root(), {lb, rb}));

    while (! heap.empty()) {
        wt_range_t e = heap.top(); heap.pop();
        if (m_wtd.is_leaf(e.v)) {
            results.emplace_back(e.v.sym, e.size());
            if (results.size()==k) {
                break;
            }
            continue;
        }

        auto child = m_wtd.expand(e.v);
        auto child_ranges = m_wtd.expand(e.v, e.r);
        auto left_range = std::get<0>(child_ranges);
        auto right_range = std::get<1>(child_ranges);

        if (!empty(left_range)) {
            heap.emplace(wt_range_t(std::get<0>(child), left_range));
        }
        if (!empty(right_range)) {
            heap.emplace(wt_range_t(std::get<1>(child), right_range));
        }
    }
    return results;
};

      };

} // end namespace

#endif
