/*!
 * this file contains a simple SORT baseline
 */
#ifndef DOCUMENT_LISING_SORT
#define DOCUMENT_LISING_SORT

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
      typename t_csa::char_type t_doc_delim = 1
      >
class doc_list_index_sort
{
    public:
        typedef t_csa                                       csa_type;
        typedef int_vector<>                                d_type;
        typedef int_vector<>::size_type                     size_type;
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

                result(size_type sp, size_type ep,list_type&& l) : list_type(l), m_sp(1), m_ep(0) {}
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

    protected:
        size_type m_doc_cnt; // number of documents in the collection
        csa_type  m_csa_full;     // CSA built from the collection text
        d_type    m_d;     // wtd build from the collection text
    public:

        //! Default constructor
        doc_list_index_sort() { }

        doc_list_index_sort(std::string file_name, sdsl::cache_config& cconfig, uint8_t num_bytes) {
            construct(m_csa_full, file_name, cconfig, num_bytes);

            const char* KEY_TEXT = key_text_trait<WIDTH>::KEY_TEXT;
            std::string text_file = cache_file_name(KEY_TEXT, cconfig);

            bit_vector doc_border;
            construct_doc_border(text_file,doc_border);
            bit_vector::rank_1_type doc_border_rank(&doc_border);
            m_doc_cnt = doc_border_rank(doc_border.size());

            int_vector_buffer<0> sa_buf(cache_file_name(conf::KEY_SA, cconfig));
            construct_D_array(sa_buf, doc_border_rank, m_doc_cnt, m_d);
        }

        size_type doc_cnt()const {
            return m_doc_cnt; // subtract one, since zero does not count
        }

        size_type word_cnt()const {
            return m_d.size()-doc_cnt();
        }

        size_type sigma()const {
            return m_csa_full.sigma;
        }


        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_doc_cnt, out, child, "doc_cnt");
            written_bytes += m_csa_full.serialize(out, child, "csa_full");
            written_bytes += m_d.serialize(out, child, "D");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            read_member(m_doc_cnt, in);
            m_csa_full.load(in);
            m_d.load(in);
        }

        void swap(doc_list_index_sort& dr) {
            if (this != &dr) {
                std::swap(m_doc_cnt, dr.m_doc_cnt);
                m_csa_full.swap(dr.m_csa_full);
                m_d.swap(dr.m_d);
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
                res = result(sp, ep);
                size_t n = ep-sp+1;
                std::vector<uint64_t> tmp(n);
                std::copy(m_d.begin()+sp,m_d.begin()+ep+1,tmp.begin());
                std::sort(tmp.begin(),tmp.end());
                size_t last = tmp[0];
                size_t f_dt = 1;
                for (size_t i=1; i<n; i++) {
                    if (tmp[i] != last) {
                        res.emplace_back(last,f_dt);
                        last = tmp[i];
                        f_dt = 1;
                    } else {
                        f_dt++;
                    }
                }
                res.emplace_back(last,f_dt);
                if (res.size() < k) k = res.size();
                static auto freq_cmp = [](const std::pair<size_type,size_type>& a,
                const std::pair<size_type,size_type>& b) {
                    return a.second > b.second;
                };
                std::partial_sort(res.begin(),res.begin()+k,res.end(),freq_cmp);
                res.resize(k);

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
};

} // end namespace

#endif
