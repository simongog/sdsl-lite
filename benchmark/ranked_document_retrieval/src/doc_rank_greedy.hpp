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
#include <sdsl/rmq_support.hpp>
#include <string>
#include <vector>
#include <fstream>
#include <list>
#include <utility>
#include "doc_rank_index.hpp"

using std::vector;

namespace sdsl
{

template
    <
      class t_csa = csa_wt<wt_huff<rrr_vector<63>>, 1000000, 1000000>,
      class t_wtd = wt_int<bit_vector_il<1024>>,
      typename t_csa::char_type t_doc_delim = 1
    >
class doc_rank_greedy
{
    public:
        typedef t_csa                                                   csa_type;
        typedef t_wtd                                                   wtd_type;
        typedef int_vector<>::size_type                                 size_type;
        typedef doc_list_tag                                            index_category;

        enum { WIDTH = t_csa::alphabet_category::WIDTH };

    protected:
        size_type       m_doc_cnt;            // number of documents in the collection
        csa_type        m_csa_full;           // CSA built from the collection text
        wtd_type        m_wtd;                // wtd build from the collection text
        int_vector<>    m_doc_perm;           // perumtation of the doc_ids
    public:

        //! Default constructor
        doc_rank_greedy() { }

        doc_rank_greedy(const std::string& input_file, sdsl::cache_config& cconfig,uint8_t num_bytes) {

            construct(m_csa_full,input_file, cconfig, num_bytes);

            construct_doc_permutation(cconfig);

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
                std::string d_file = cache_file_name(KEY_DARRAY, cconfig);
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
            written_bytes += m_doc_perm.serialize(out, child, "doc_perm");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            read_member(m_doc_cnt, in);
            m_csa_full.load(in);
            m_wtd.load(in);
            m_doc_perm.load(in);
        }

        void swap(doc_rank_greedy& dr) {
            if (this != &dr) {
                std::swap(m_doc_cnt, dr.m_doc_cnt);
                m_csa_full.swap(dr.m_csa_full);
                m_wtd.swap(dr.m_wtd);
                m_doc_perm.swap(dr.m_doc_perm);
            }
        }

        //! Search for the k documents which contain the search term most frequent
        template<class t_pat>
        result search(const t_pat& query, size_t k) const {
            result res;
            /*
            size_type sp=1, ep=0;
            if (0 == backward_search(m_csa_full, 0, m_csa_full.size()-1, begin, end, sp, ep)) {
                res = result();
                return 0;
            } else {
                auto tmp_res = m_wtd.topk_greedy(sp,ep,k);
                res = result(sp, ep, std::move(tmp_res));
                return ep-sp+1;
            }*/
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

        void
        construct_doc_permutation(const sdsl::cache_config& cconfig)
        {
            std::string docweight_file = cache_file_name(KEY_DOCWEIGHT, cconfig);
            std::ifstream dw_in(docweight_file);
            if(! dw_in.is_open() ) {
                std::cerr << "cannot open docweight file";
                exit(EXIT_FAILURE);
            }
            size_t cur_id = 0;
            using u64p_t = std::pair<uint64_t,double>;
            std::vector< u64p_t > weights;
            char buf[512];
            for(;dw_in.getline(buf, 512); ) {
                weights.emplace_back(cur_id,std::strtod(buf,NULL));
            }
            auto wsorter = [](const u64p_t& a,const u64p_t& b)
                { return a.second < b.second; };
            std::sort(weights.begin(),weights.end(),wsorter);

            m_doc_perm.width(bits::hi(weights.size())+1);
            m_doc_perm.resize(weights.size());
            for(size_t i=0;i<weights.size();i++) {
                m_doc_perm[weights[i].first] = i;
            }
        }

};

} // end namespace

#endif
