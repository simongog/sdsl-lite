/*! How to code a parametrizable document listing data structure
 *
 * This file contains a document listing class implemented as
 * suggested in Kunihiko Sadakane's article:
 * ,,Succinct Data Structures for Flexible Text Retrieval Systems''
 * Journal of Discrete Algorithms, 2007.
 *
 */
#ifndef DOCUMENT_LISING_SADA_INCLUDED
#define DOCUMENT_LISING_SADA_INCLUDED

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

template<uint8_t t_width>
struct sa_trait {
    typedef uint64_t value_type;
    typedef std::vector<value_type> vec_type;
    enum { num_bytes = 0 };
    template <class t_sa>
    static void calc_sa(t_sa& sa, vec_type& text) {
        qsufsort::construct_sa(sa, text);
    }
};

template<>
struct sa_trait<8> {
    typedef uint8_t value_type;
    typedef std::vector<value_type> vec_type;
    enum { num_bytes = 1 };
    template <class t_sa>
    static void calc_sa(t_sa& sa, vec_type& text) {
        algorithm::calculate_sa(text.data(), text.size(), sa);
    }
};


template<
class t_csa_full                = csa_wt<wt_huff<rrr_vector<63>>, 30, 1000000, text_order_sa_sampling<> >,
      class t_range_min         = rmq_succinct_sct<true>,
      class t_range_max         = rmq_succinct_sct<false>,
      class t_doc_border        = sd_vector<>,
      class t_doc_border_rank   = typename t_doc_border::rank_1_type,
      class t_doc_border_select = typename t_doc_border::select_1_type,
      typename t_csa_full::char_type t_doc_delim = 1
      >
      class doc_list_index_sada
      {
          public:
          typedef t_csa_full                                  csa_full_type;
          typedef t_range_min                                 range_min_type;
          typedef t_range_max                                 range_max_type;
          typedef t_doc_border                                doc_border_type;
          typedef t_doc_border_rank                           doc_border_rank_type;
          typedef t_doc_border_select                         doc_border_select_type;
          typedef int_vector<>::size_type                     size_type;
          typedef std::vector<std::pair<size_type,size_type>> list_type;
          typedef doc_list_tag                                index_category;

          enum { WIDTH = t_csa_full::alphabet_category::WIDTH };

          typedef sa_trait<WIDTH>                             sa_tt;

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


private:
size_type              m_doc_cnt;           // number of documents in the collection
csa_full_type          m_csa_full;          // CSA build from the collection text
vector<int_vector<>>   m_doc_isa;           // array of inverse SAs. m_doc_isa[i] contains the ISA of document i
range_min_type         m_rminq;             // range minimum data structure build over an array Cprev
range_max_type         m_rmaxq;             // range maximum data structure build over an array Cnext
doc_border_type        m_doc_border;        // bitvector indicating the positions of the separators in the collection text
doc_border_rank_type   m_doc_border_rank;   // rank data structure on m_doc_border
doc_border_select_type m_doc_border_select; // select data structure on m_doc_border
size_type              m_doc_max_len;       // maximal length of a document in the collection
mutable bit_vector     m_doc_rmin_marked;   // helper bitvector for search process
mutable bit_vector     m_doc_rmax_marked;   // helper bitvector for search process

public:

//! Default constructor
doc_list_index_sada() { }

doc_list_index_sada(std::string file_name, sdsl::cache_config& cconfig, uint8_t num_bytes) {
    construct(m_csa_full, file_name, cconfig, num_bytes);

    const char* KEY_TEXT = key_text_trait<WIDTH>::KEY_TEXT;
    std::string text_file = cache_file_name(KEY_TEXT, cconfig);

    construct_doc_border(text_file, m_doc_border, m_doc_max_len);
    m_doc_border_rank   = doc_border_rank_type(&m_doc_border);
    m_doc_border_select = doc_border_select_type(&m_doc_border);
    m_doc_cnt = m_doc_border_rank(m_doc_border.size());

    construct_doc_isa(text_file, m_doc_cnt, m_doc_max_len, m_doc_isa);

    int_vector_buffer<0> sa_buf(cache_file_name(conf::KEY_SA, cconfig));
    {
        int_vector<> D;
        construct_D_array(sa_buf, m_doc_border_rank, m_doc_cnt, D);
        {
            int_vector<> Cprev;
            construct_Cprev_array(D, m_doc_cnt, Cprev);
            range_min_type rminq(&Cprev);
            m_rminq = rminq;
        }
        {
            int_vector<> Cnext;
            construct_Cnext_array(D, m_doc_cnt, Cnext);
            range_max_type rmaxq(&Cnext);
            m_rmaxq = rmaxq;
        }
    }
    m_doc_rmin_marked = bit_vector(m_doc_cnt, 0);
    m_doc_rmax_marked = bit_vector(m_doc_cnt, 0);
}

size_type doc_cnt()const {
    return m_doc_cnt;
}

size_type word_cnt()const {
    return m_csa_full.size()-doc_cnt();
}

size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += write_member(m_doc_cnt, out, child, "doc_cnt");
    written_bytes += m_csa_full.serialize(out, child, "csa_full");
    written_bytes += serialize_vector(m_doc_isa, out, child, "doc_isa");
    written_bytes += m_rminq.serialize(out, child, "rminq");
    written_bytes += m_rmaxq.serialize(out, child, "rmaxq");
    written_bytes += m_doc_border.serialize(out, child, "doc_border");
    written_bytes += m_doc_border_rank.serialize(out, child, "doc_border_rank");
    written_bytes += m_doc_border_select.serialize(out, child, "doc_border_select");
    written_bytes += write_member(m_doc_max_len, out, child, "doc_max_len");
    // helper bitvector m_doc_rmin_marked and m_doc_rmax_marked are not serialize
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void load(std::istream& in) {
    read_member(m_doc_cnt, in);
    m_csa_full.load(in);
    m_doc_isa.resize(m_doc_cnt);
    load_vector(m_doc_isa, in);
    m_rminq.load(in);
    m_rmaxq.load(in);
    m_doc_border.load(in);
    m_doc_border_rank.load(in);
    m_doc_border_rank.set_vector(&m_doc_border);
    m_doc_border_select.load(in);
    m_doc_border_select.set_vector(&m_doc_border);
    read_member(m_doc_max_len, in);
    // also initialize the helper bitvectors
    m_doc_rmin_marked = bit_vector(m_doc_cnt);
    m_doc_rmax_marked = bit_vector(m_doc_cnt);
}

void swap(doc_list_index_sada& dr) {
    if (this != &dr) {
        std::swap(m_doc_cnt, dr.m_doc_cnt);
        m_csa_full.swap(dr.m_csa_full);
        m_doc_isa.swap(dr.m_doc_isa);
        m_rminq.swap(dr.m_rminq);
        m_rmaxq.swap(dr.m_rmaxq);
        m_doc_border.swap(dr.m_doc_border);
        util::swap_support(m_doc_border_rank, dr.m_doc_border_rank,
                           &m_doc_border, &(dr.m_doc_border));
        util::swap_support(m_doc_border_select, dr.m_doc_border_select,
                           &m_doc_border, &(dr.m_doc_border));
        std::swap(m_doc_max_len, dr.m_doc_max_len);
        m_doc_rmin_marked.swap(dr.m_doc_rmin_marked);
        m_doc_rmax_marked.swap(dr.m_doc_rmax_marked);
    }
}

//! Search for the k documents which contains the search term most frequent
template<class t_pat_iter>
size_t
search(t_pat_iter begin,
       t_pat_iter end,
       result& res,
size_t k) const {
    size_type sp=1, ep=0;
    if (0 == backward_search(m_csa_full, 0, m_csa_full.size()-1, begin, end, sp, ep)) {
        res = result();
        return 0;
    } else {
        res = result(sp, ep);
        compute_tf_idf(sp, ep, res);
        size_t kprime = std::min(res.size(), k);
        auto comp = [](std::pair<size_type,size_type>& a,std::pair<size_type,size_type>& b) {
            return (a.second != b.second) ? a.second > b.second  : a.first < b.first;
        };
        partial_sort(res.begin(),res.begin()+kprime, res.end(), comp);
        res.resize(kprime);
        return ep-sp+1;
    }
}

private:
void compute_tf_idf(const size_type& sp, const size_type& ep, result& res)const {
    vector<size_type> suffixes;
    get_lex_smallest_suffixes(sp, ep, suffixes);
    get_lex_largest_suffixes(sp, ep, suffixes);
    sort(suffixes.begin(), suffixes.end());

    for (size_type i=0; i < suffixes.size(); i+=2) {
        size_type suffix_1 = suffixes[i];
        size_type suffix_2 = suffixes[i+1];
        size_type doc                 = m_doc_border_rank(suffix_1+1);
        m_doc_rmin_marked[doc]        = 0;  // reset marking, which was set in get_lex_smallest_suffixes
        m_doc_rmax_marked[doc]        = 0;  //                                 get_lex_largest_suffixes

        if (suffix_1 == suffix_2) {  // if pattern occurs exactly once
            res.push_back( {doc,1}); // add the #occurrence
        } else {
            size_type doc_begin = doc ? m_doc_border_select(doc) + 1 : 0;
            size_type doc_sp    = m_doc_isa[doc][ suffix_1 - doc_begin ];
            size_type doc_ep    = m_doc_isa[doc][ suffix_2 - doc_begin ];
            if (doc_sp > doc_ep) {
                std::swap(doc_sp, doc_ep);
            }
            res.push_back( {doc, doc_ep - doc_sp + 1});
        }
    }
}

void get_lex_smallest_suffixes(size_type sp, size_type ep, vector<size_type>& suffixes) const {
    using lex_range_t = std::pair<size_type,size_type>;
    std::stack<lex_range_t> stack;
    stack.emplace(sp,ep);
    while (!stack.empty()) {
        auto range = stack.top();
        stack.pop();
        size_type rsp = std::get<0>(range);
        size_type rep = std::get<1>(range);
        if (rsp <= rep) {
            size_type min_idx = m_rminq(rsp,rep);
            size_type suffix  = m_csa_full[min_idx];
            size_type doc     = m_doc_border_rank(suffix+1);

            if (!m_doc_rmin_marked[doc]) {
                suffixes.push_back(suffix);
                m_doc_rmin_marked[doc] = 1;
                stack.emplace(min_idx+1,rep);
                stack.emplace(rsp,min_idx-1); // min_idx != 0, since `\0` is appended to string
            }
        }
    }
}

void get_lex_largest_suffixes(size_type sp, size_type ep, vector<size_type>& suffixes) const {
    using lex_range_t = std::pair<size_type,size_type>;
    std::stack<lex_range_t> stack;
    stack.emplace(sp,ep);
    while (!stack.empty()) {
        auto range = stack.top();
        stack.pop();
        size_type rsp = std::get<0>(range);
        size_type rep = std::get<1>(range);
        if (rsp <= rep) {
            size_type max_idx = m_rmaxq(rsp,rep);
            size_type suffix  = m_csa_full[max_idx];
            size_type doc     = m_doc_border_rank(suffix+1);

            if (!m_doc_rmax_marked[doc]) {
                suffixes.push_back(suffix);
                m_doc_rmax_marked[doc] = 1;
                stack.emplace(rsp,max_idx - 1); // max_idx != 0, since `\0` is appended to string
                stack.emplace(max_idx+1,rep);
            }
        }
    }
}

//! Construct the doc_border bitvector by streaming the text file
void
construct_doc_border(const std::string& text_file,
                     doc_border_type& doc_border,
size_type& doc_max_len) {
    int_vector_buffer<WIDTH> text_buf(text_file);
    bit_vector tmp_doc_border(text_buf.size(), 0); // create temporary uncompressed vector
    doc_max_len = 0;
    size_type len = 0;
    for (size_type i = 0; i < text_buf.size(); ++i) {
        if (t_doc_delim == text_buf[i]) {
            tmp_doc_border[i] = 1;
            doc_max_len = std::max(doc_max_len, len);
            len = 0;
        } else {
            ++len;
        }
    }
    doc_border = doc_border_type(tmp_doc_border);
}

void
construct_doc_isa(const std::string& text_file,
                  const size_type doc_cnt,
                  SDSL_UNUSED const size_type doc_max_len,
vector<int_vector<> >& doc_isa) {
    doc_isa.resize(doc_cnt);
    typename sa_tt::vec_type doc_buffer;
    int_vector_buffer<WIDTH> text_buf(text_file);
    size_type doc_id = 0;
    for (size_type i = 0; i < text_buf.size(); ++i) {
        if (t_doc_delim == text_buf[i]) {
            if (doc_buffer.size() > 0) {
                doc_buffer.push_back(0);
                construct_doc_isa(doc_buffer, doc_isa[doc_id]);
                ++doc_id;
            }
            doc_buffer.clear();
        } else {
            doc_buffer.push_back(text_buf[i]);
        }
    }
}

void
construct_doc_isa(typename sa_tt::vec_type& doc_buffer,
int_vector<>& doc_isa) {
    int_vector<> sa(doc_buffer.size(), 0, bits::hi(doc_buffer.size())+1);
    sa_tt::calc_sa(sa, doc_buffer);
    util::bit_compress(sa);
    doc_isa = sa;
    for (size_type i = 0; i < doc_buffer.size(); ++i) {
        doc_isa[sa[i]] = i;
    }
}

void
construct_D_array(int_vector_buffer<0>& sa_buf,
                  const doc_border_rank_type& doc_border_rank,
                  const size_type doc_cnt,
int_vector<>& D) {
    D = int_vector<>(sa_buf.size(), 0, bits::hi(doc_cnt+1)+1);
    for (size_type i = 0; i < sa_buf.size(); ++i) {
        D[i] = doc_border_rank(sa_buf[i]+1);
    }
}


void
construct_Cprev_array(const int_vector<>& D,
                      size_type doc_cnt,
int_vector<>& Cprev) {
    Cprev = int_vector<>(D.size(), 0, bits::hi(D.size())+1);
    int_vector<> last_occ(doc_cnt+1, 0, bits::hi(D.size())+1);
    for (size_type i = 0; i < D.size(); ++i) {
        size_type doc = D[i];
        Cprev[i]      = last_occ[doc];
        last_occ[doc] = i;
    }
}

void
construct_Cnext_array(const int_vector<>& D,
                      size_type doc_cnt,
int_vector<>& Cnext) {
    Cnext = int_vector<>(D.size(), 0, bits::hi(D.size())+1);
    int_vector<> last_occ(doc_cnt+1, D.size(), bits::hi(D.size())+1);
    for (size_type i = 0, j = D.size()-1; i < D.size(); ++i, --j) {
        size_type doc = D[j];
        Cnext[j]      = last_occ[doc];
        last_occ[doc] = j;
    }
}
      };

} // end namespace

#endif
