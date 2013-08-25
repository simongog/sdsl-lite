/*!
 * This file contains a document listing class, which implements
 * strategy GREEDY in the article:
 * J. S. Culpepperi, G. Navarro, S. J. Puglisi and A. Turpin:
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
class t_csa          = csa_wt<wt_huff<rrr_vector<63>>, 1000000, 1000000>,
      class t_wtd    = wt_int<bit_vector,rank_support_v5<1>,select_support_scan<1>,select_support_scan<0>>,
      typename t_csa::char_type t_doc_delim = 1
      >
class doc_list_index_greedy
{
    public:
        typedef t_csa                                       csa_type;
        typedef t_wtd                                       wtd_type;
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

                // Constructors for an empty result and for a result in the interval [sp, ep]:
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
        size_type                   m_doc_cnt;              // number of documents in the collection
        csa_type                    m_csa;                  // CSA built from the collection text
        wtd_type                    m_wtd;                  // wtd build from the collection text
    public:

        //! Default constructor
        doc_list_index_greedy() { }

        doc_list_index_greedy(std::string file_name, sdsl::cache_config& cconfig, uint8_t num_bytes) {
            construct(m_csa, file_name, cconfig, num_bytes);

            const char* KEY_TEXT = key_text_trait<WIDTH>::KEY_TEXT;
            std::string text_file = cache_file_name(KEY_TEXT, cconfig);

            bit_vector doc_border;
            construct_doc_border(text_file,doc_border);
            bit_vector::rank_1_type doc_border_rank(&doc_border);
            m_doc_cnt = doc_border_rank(doc_border.size());

            int_vector_buffer<0> sa_buf(cache_file_name(constants::KEY_SA, cconfig));
            {
                int_vector<> D;
                construct_D_array(sa_buf, doc_border_rank, m_doc_cnt, D);
                construct_im(m_wtd,D);
            }
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_doc_cnt, out, child, "doc_cnt");
            written_bytes += m_csa.serialize(out, child, "csa");
            written_bytes += m_wtd.serialize(out, child, "wtd");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            read_member(m_doc_cnt, in);
            m_csa.load(in);
            m_wtd.load(in);
        }

        void swap(doc_list_index_greedy& dr) {
            if (this != &dr) {
                std::swap(m_doc_cnt, dr.m_doc_cnt);
                m_csa.swap(dr.m_csa);
                m_wtd.swap(dr.m_wtd);
            }
        }

        //! Search for the k documents which contains the search term most frequent
        size_type search(std::string::iterator begin, std::string::iterator end, result& res,size_t k) const {
            size_type sp=1, ep=0;
            if (0 == backward_search(m_csa, 0, m_csa.size()-1, begin, end, sp, ep)) {
                res = result();
                return 0;
            } else {
                auto tmp_res = m_wtd.topk_greedy(sp,ep,k);
                res = result(sp, ep, std::move(tmp_res));
                return ep-sp+1;
            }
        }

    private:
        //! Construct the doc_border bitvector by streaming the text file
        void
        construct_doc_border(const std::string& text_file, bit_vector& doc_border) {
            int_vector_buffer<8> text_buf(text_file);
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
                D[i] = doc_border_rank(sa_buf[i]+1);
            }
        }
};

} // end namespace

#endif
