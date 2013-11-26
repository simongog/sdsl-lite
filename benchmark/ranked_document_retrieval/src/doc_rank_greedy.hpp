/*!
 * This file contains a document listing class, which implements
 * strategy GREEDY in the article:
 * J. S. Culpepper, G. Navarro, S. J. Puglisi and A. Turpin:
 * ,,Top-k Ranked Document Search in General Text Databases''
 * Proceedings Part II of the 18th Annual European Symposium on
 * Algorithms (ESA 2010)
 */
#ifndef DOC_RANK_GREEDY_HPP
#define DOC_RANK_GREEDY_HPP

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>
#include <string>
#include <vector>
#include <fstream>
#include <list>
#include <utility>
#include <stack>
#include "doc_rank_index.hpp"
#include <chrono>

using std::vector;
using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

namespace sdsl
{

template
<
class t_csa = csa_wt<wt_huff<rrr_vector<63>>, 1000000, 1000000>,
      class t_wtd = wt_int<bit_vector_il<1024>>,
      class t_sadadfbv = rrr_vector<63>,
      class t_sadadf_select = typename t_sadadfbv::select_1_type,
      typename t_csa::char_type t_doc_delim = 1
      >
      class doc_rank_greedy
      {
          public:
          typedef t_csa                                                   csa_type;
          typedef t_wtd                                                   wtd_type;
          typedef t_sadadfbv                                              dfbv_type;
          typedef t_sadadf_select                                         dfbv_select_type;
          typedef int_vector<>::size_type                                 size_type;
          typedef doc_list_tag                                            index_category;

          enum { WIDTH = t_csa::alphabet_category::WIDTH };

          protected:
          size_type        m_doc_cnt;            // number of documents in the collection
          csa_type         m_csa_full;           // CSA built from the collection text
          wtd_type         m_wtd;                // wtd build from the collection text
          int_vector<>     m_doc_perm;           // perumtation of the doc_ids
          dfbv_type        m_sada_df;
          dfbv_select_type m_sada_df_select;
          public:

          //! Default constructor
doc_rank_greedy() { }

doc_rank_greedy(const std::string& input_file, sdsl::cache_config& cconfig,uint8_t num_bytes) {
    // (1) check, if the compressed suffix array exists
    //     ->  also creates the bwt and the suffix array
    if (!cache_file_exists(conf::KEY_CSA+std::string("_")+util::class_to_hash(m_csa_full), cconfig)) {
        construct(m_csa_full,input_file, cconfig, num_bytes);
        store_to_cache(m_csa_full,conf::KEY_CSA+std::string("_")+util::class_to_hash(m_csa_full), cconfig);
        util::clear(m_csa_full);
    } else {
        register_cache_file(conf::KEY_CSA+std::string("_")+util::class_to_hash(m_csa_full), cconfig);
        register_cache_file(conf::KEY_SA, cconfig);
        const char* KEY_TEXT = key_text_trait<WIDTH>::KEY_TEXT;
        register_cache_file(KEY_TEXT, cconfig);
    }

    construct_doc_permutation(cconfig);
    m_doc_cnt = m_doc_perm.size();

    if (!cache_file_exists(KEY_DARRAY, cconfig)) {
        const char* KEY_TEXT = key_text_trait<WIDTH>::KEY_TEXT;
        std::string text_file = cache_file_name(KEY_TEXT, cconfig);
        bit_vector doc_border;
        construct_doc_border(text_file,doc_border);
        bit_vector::rank_1_type doc_border_rank(&doc_border);

        construct_D_array(cconfig, doc_border_rank, m_doc_cnt);
        store_to_file(m_doc_perm, cache_file_name(KEY_DOCPERM, cconfig));
        util::clear(m_doc_perm);
    } else {
        register_cache_file(KEY_DARRAY, cconfig);
        register_cache_file(KEY_DOCPERM, cconfig);
        util::clear(m_doc_perm);
    }

    if (!cache_file_exists(KEY_SADADF+std::string("_")+util::class_to_hash(m_sada_df), cconfig)) {
        construct_sada_df(cconfig);
        store_to_file(m_sada_df,cache_file_name(KEY_SADADF+std::string("_")+util::class_to_hash(m_sada_df),cconfig));
        store_to_file(m_sada_df_select,cache_file_name(KEY_SADADFSEL+std::string("_")+util::class_to_hash(m_sada_df_select),cconfig));
    } else {
        load_from_cache(m_sada_df,KEY_SADADF+std::string("_")+util::class_to_hash(m_sada_df), cconfig);
        load_from_cache(m_sada_df_select,KEY_SADADFSEL+std::string("_")+util::class_to_hash(m_sada_df_select),cconfig);
        m_sada_df_select.set_vector(&m_sada_df);
    }

    if (!cache_file_exists(KEY_WTD+std::string("_")+util::class_to_hash(m_wtd), cconfig)) {
        construct(m_wtd, cache_file_name(KEY_DARRAY, cconfig));
        store_to_file(m_wtd,cache_file_name(KEY_WTD+std::string("_")+util::class_to_hash(m_wtd),cconfig));
    } else {
        load_from_cache(m_wtd,KEY_WTD+std::string("_")+util::class_to_hash(m_wtd), cconfig);
    }

    // load the remaining pieces from disk
    load_from_cache(m_csa_full,conf::KEY_CSA+std::string("_")+util::class_to_hash(m_csa_full), cconfig);
    load_from_cache(m_doc_perm,KEY_DOCPERM, cconfig);
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
    written_bytes += m_doc_perm.serialize(out, child, "doc_perm");
    written_bytes += m_sada_df.serialize(out, child, "sada_df");
    written_bytes += m_sada_df_select.serialize(out, child, "sada_df_select");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void load(std::istream& in) {
    read_member(m_doc_cnt, in);
    m_csa_full.load(in);
    m_wtd.load(in);
    m_doc_perm.load(in);
    m_sada_df.load(in);
    m_sada_df_select.load(in, &m_sada_df);
}

void swap(doc_rank_greedy& dr) {
    if (this != &dr) {
        std::swap(m_doc_cnt, dr.m_doc_cnt);
        m_csa_full.swap(dr.m_csa_full);
        m_wtd.swap(dr.m_wtd);
        m_doc_perm.swap(dr.m_doc_perm);
        m_sada_df.swap(dr.m_sada_df);
        util::swap_support(m_sada_df_select, dr.m_sada_df_select, &m_sada_df, &(dr.m_sada_df));
    }
}

//! Search for the k documents which contain the search term most frequent
template<class t_pat>
result search(const t_pat& query, size_t k) const {
    result res;
    std::vector<typename t_wtd::greedy_tfidf_range_t> ranges;
for (const auto& p : query) {
        size_type sp=1, ep=0;
        if (0 != backward_search(m_csa_full, 0, m_csa_full.size()-1, p.begin(), p.end(), sp, ep)) {
            size_t f_t = calculate_df(sp,ep);

        }
    }
    if (ranges.size()) {

    }
    return res;
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
construct_D_array(sdsl::cache_config& cconfig, bit_vector::rank_1_type& doc_border_rank, const size_type doc_cnt) {
    int_vector_buffer<0> sa_buf(cache_file_name(conf::KEY_SA, cconfig));
    std::string d_file = cache_file_name(KEY_DARRAY, cconfig);
    int_vector_buffer<> D(d_file,std::ios::out,1024*1024,bits::hi(doc_cnt)+1,false);
    for (size_type i = 0; i < sa_buf.size(); ++i) {
        uint64_t d = doc_border_rank(sa_buf[i]+1);
        D[i] = m_doc_perm[d];
    }
}

void
construct_doc_permutation(sdsl::cache_config& cconfig) {
    std::string docweight_file = cconfig.file_map[KEY_DOCWEIGHT];
    std::ifstream dw_in(docweight_file);
    if (! dw_in.is_open()) {
        std::cerr << "cannot open docweight file";
        exit(EXIT_FAILURE);
    }
    size_t cur_id = 0;
    using u64p_t = std::pair<uint64_t,double>;
    std::vector< u64p_t > weights;
    char buf[512];
    for (; dw_in.getline(buf, 512);) {
        weights.emplace_back(cur_id,std::strtod(buf,NULL));
        cur_id++;
    }
    auto wsorter = [](const u64p_t& a,const u64p_t& b)
    { return a.second < b.second; };
    std::sort(weights.begin(),weights.end(),wsorter);

    m_doc_perm.width(bits::hi(weights.size())+1);
    m_doc_perm.resize(weights.size());
    for (size_t i=0; i<weights.size(); i++) {
        m_doc_perm[weights[i].first] = i;
    }
}

void
construct_sada_df(sdsl::cache_config& cconfig) {
    // construct cst
    if (!cache_file_exists(conf::KEY_LCP, cconfig)) {
        if (cst_sct3<csa_type>::alphabet_category::WIDTH==8) {
            construct_lcp_semi_extern_PHI(cconfig);
        } else {
            construct_lcp_PHI<cst_sct3<csa_type>::alphabet_category::WIDTH>(cconfig);
        }
    }
    register_cache_file(conf::KEY_LCP, cconfig);

    cst_sct3<csa_type> temp_cst;
    if (!cache_file_exists("tempcst", cconfig)) {
        temp_cst = cst_sct3<csa_type>(cconfig,true);
        store_to_file(temp_cst,cache_file_name("tempcst",cconfig));
    } else {
        load_from_file(temp_cst,cache_file_name("tempcst",cconfig));
    }

    wt_int<> wtc;
    std::string d_file = cache_file_name(KEY_DARRAY, cconfig);
    int_vector_buffer<> D(d_file);
    if (!cache_file_exists("wtc", cconfig)) {
        {
            int_vector<> C(D.size(),0,bits::hi(D.size())+1);
            int_vector<> last_occ(m_doc_cnt,D.size(),bits::hi(D.size())+1);
            for (size_t i=0; i < D.size(); ++i) {
                uint64_t d = D[i];
                C[i] = last_occ[d];
                last_occ[d] = i;
            }
            util::bit_compress(C);
            store_to_file(C, cache_file_name(KEY_C, cconfig));
        }
        construct(wtc, cache_file_name(KEY_C, cconfig),cconfig,0);
        sdsl::remove(cache_file_name(KEY_C, cconfig));
        store_to_file(wtc,cache_file_name("wtc",cconfig));
    } else {
        load_from_file(wtc,cache_file_name("wtc",cconfig));
    }

    // construct the bv
    bit_vector h(2*D.size(), 0);
    size_t h_idx = 0;
    using n_type = std::tuple<cst_sct3<>::node_type, size_t, size_t, bool>;
    std::stack<n_type> s;
    s.emplace(temp_cst.root(), 1, temp_cst.degree(temp_cst.root()), true);
    // invariant: node has two children
    while (!s.empty()) {
        n_type node = s.top(); s.pop();
        auto v = std::get<0>(node);
        auto l_child = std::get<1>(node); auto r_child = std::get<2>(node);
        auto first   = std::get<3>(node);
        if (first) {   // first half
            // recurse down
            std::get<3>(node) = false;
            s.push(node);
            if (r_child == l_child+1) {
                auto w = temp_cst.select_child(v, l_child);
                if (!temp_cst.is_leaf(w))
                    s.emplace(w, 1, temp_cst.degree(w), true);
            } else {
                auto mid = l_child + (r_child-l_child)/2;
                s.emplace(v, l_child, mid, true);
            }
        } else { // second half
            auto lb  = temp_cst.lb(temp_cst.select_child(v, l_child));
            auto rb  = temp_cst.rb(temp_cst.select_child(v, r_child));
            auto mid  = l_child + (r_child-l_child)/2;
            size_t dup_elements = 0;
            if (lb+1 == rb) {
                dup_elements = (wtc[rb] == lb);
            } else {
                auto mid_rb = temp_cst.rb(temp_cst.select_child(v, mid));
                auto mid_lb = mid_rb+1;
                dup_elements = std::get<0>(wtc.range_search_2d(mid_lb, rb, lb, mid_rb, false));
            }
            h_idx+=dup_elements;
            h[h_idx++] = 1;
            if (mid+1 == r_child) {
                auto w = temp_cst.select_child(v, r_child);
                if (!temp_cst.is_leaf(w))
                    s.emplace(w, 1, temp_cst.degree(w), true);
            } else {
                s.emplace(v, mid+1, r_child, true);
            }
        }
    }
    h.resize(h_idx);
    util::clear(temp_cst);
    // convert to proper bv type
    m_sada_df = dfbv_type(h);
    m_sada_df_select = dfbv_select_type(&m_sada_df);
}

size_type
calculate_df(size_t sp,size_t ep) const {
    size_t dup = 0;
    size_t y = m_sada_df_select(ep);
    if (0 == sp) {
        dup = (y+1)-ep;// (#all elements)-#1=#0
    } else {
        size_t x = m_sada_df_select(sp);
        dup = (y+1)-ep - ((x+1)-sp);
    }
    return ep-sp+1-dup;
}
      };

} // end namespace

#endif
