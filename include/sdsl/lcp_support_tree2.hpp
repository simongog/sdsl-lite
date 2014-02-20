#ifndef INCLUDED_SDSL_SUPPORT_TREE2
#define INCLUDED_SDSL_SUPPORT_TREE2

#include "lcp.hpp"
#include "util.hpp"
#include "rank_support_v.hpp"
#include "wt_huff.hpp"
#include "sorted_multi_stack_support.hpp"
#include <iostream>
#include <string>

namespace sdsl
{

// Forward declaration of helper method
template<uint32_t t_dens, uint8_t t_bwt_width>
void construct_first_child_and_lf_lcp(int_vector_buffer<>&,
                                      int_vector_buffer<t_bwt_width>&,
                                      const std::string&,
                                      const std::string&, int_vector<>&);


/*! An lcp array class for cst_sct3 and cst_sada.
 *    The time of the []-operator depends on:
 *    - The time of the []-operation of the wt_huff
 *    - The time of the LF calculation of the underlying CSA of the CST
 *    - The time of the tlcp_idx function of the CST
 *
 *  \tparam t_dens Sample density in the CST.
 *  \tparam t_cst  Underlying CST.
 */
template<uint32_t t_dens, class t_cst>
class _lcp_support_tree2
{
    public:
        typedef int_vector<>::value_type                         value_type;
        typedef random_access_const_iterator<_lcp_support_tree2> const_iterator;
        typedef const_iterator                                   iterator;
        typedef const value_type                                 const_reference;
        typedef const_reference                                  reference;
        typedef const_reference*                                 pointer;
        typedef const pointer                                    const_pointer;
        typedef int_vector<>::size_type                          size_type;
        typedef int_vector<>::difference_type                    difference_type;
        typedef t_cst                                            cst_type;
        typedef wt_huff<bit_vector, rank_support_v5<>,
                select_support_scan<1>,
                select_support_scan<0> >                         small_lcp_type;

        typedef lcp_tree_and_lf_compressed_tag                   lcp_category;

        enum { fast_access = 0,
               text_order = 0,
               sa_order = 0
             };

        template<class CST>
        struct type {
            typedef _lcp_support_tree2 lcp_type;
        };

    private:
        const cst_type*    m_cst;
        small_lcp_type  m_small_lcp; // vector for lcp values < 254
        int_vector<> m_big_lcp;      // vector for lcp values >= 254

    public:

        //! Default constructor
        _lcp_support_tree2() {}

        //! Copy / Move constructor
        _lcp_support_tree2(const _lcp_support_tree2&)  = default;
        _lcp_support_tree2(_lcp_support_tree2&&)  = default;
        _lcp_support_tree2& operator=(const _lcp_support_tree2&)  = default;
        _lcp_support_tree2& operator=(_lcp_support_tree2&&) = default;


        //! Constructor
        /*! \param config Cache configuration.

         */
        _lcp_support_tree2(cache_config& config, const cst_type* cst = nullptr) {
            m_cst = cst;

            int_vector_buffer<> lcp_buf(cache_file_name(conf::KEY_LCP, config));
            std::string bwt_file = cache_file_name(key_trait<t_cst::csa_type::alphabet_type::int_width>::KEY_BWT, config);
            int_vector_buffer<t_cst::csa_type::alphabet_type::int_width> bwt_buf(bwt_file);

            std::string sml_lcp_file = tmp_file(config, "_fc_lf_lcp_sml");
            std::string big_lcp_file = tmp_file(config, "_fc_lf_lcp_big");

            construct_first_child_and_lf_lcp<t_dens>(lcp_buf, bwt_buf, sml_lcp_file, big_lcp_file, m_big_lcp);
            int_vector_buffer<8> sml_lcp_buf(sml_lcp_file);

            {
                small_lcp_type tmp_small_lcp(sml_lcp_buf, sml_lcp_buf.size());
                m_small_lcp.swap(tmp_small_lcp);
            }
            sml_lcp_buf.close(true);
            sdsl::remove(big_lcp_file);
        }

        void set_cst(const cst_type* cst) {
            m_cst = cst;
        }

        size_type size()const {
            return m_cst->size();
        }

        static size_type max_size() {
            return int_vector<>::max_size();
        }

        size_type empty()const {
            return m_small_lcp.empty();
        }

        void swap(_lcp_support_tree2& lcp_c) {
            m_small_lcp.swap(lcp_c.m_small_lcp);
            m_big_lcp.swap(lcp_c.m_big_lcp);
        }

        //! Returns a const_iterator to the first element.
        const_iterator begin()const {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        const_iterator end()const {
            return const_iterator(this, size());
        }

        //! []-operator
        /*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
         * \par Time complexity
         *     \f$ \Order{t_{find\_close} + t_{rank}} \f$
         */
        inline value_type operator[](size_type i)const {
            size_type idx, offset=0;
            uint8_t val;
start:
            idx = m_cst->tlcp_idx(i);
            val = m_small_lcp[idx];
            if (val < 254) {
                return val;// - offset;
            } else if (val == 254) { // if lcp value is >= 254 and position i is reducible
                i = m_cst->csa.lf[i]; // i = LF[i]    // (*m_psi)(i);
                ++offset; // goto lcp value, which is one bigger
                goto start;
            } else { // if lcp value is >= 254 and (not reducable or sampled)
                return m_big_lcp[m_small_lcp.rank(idx ,255)] - offset;
            }
        }

        //! Serialize to a stream.
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_small_lcp.serialize(out, child, "small_lcp");
            written_bytes += m_big_lcp.serialize(out, child, "large_lcp");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load from a stream.
        void load(std::istream& in, const t_cst* cst=nullptr) {
            m_small_lcp.load(in);
            m_big_lcp.load(in);
            m_cst = cst;
        }
};

//! Helper class which provides _lcp_support_tree2 the context of a CST.
template<uint32_t t_dens=16>
struct lcp_support_tree2 {
    template<class t_cst>
    using type = _lcp_support_tree2<t_dens, t_cst>;
};


/*!
 * \tparam t_dens       Sample an LCP value x if x modulo t_dens == 0
 * \tparam t_bwt_width  Width of the integers of the streamed BWT array.
 * \tparam
 */
template<uint32_t t_dens, uint8_t t_bwt_width>
void construct_first_child_and_lf_lcp(int_vector_buffer<>& lcp_buf,
                                      int_vector_buffer<t_bwt_width>& bwt_buf,
                                      const std::string& small_lcp_file,
                                      const std::string& big_lcp_file,
                                      int_vector<>& big_lcp)
{
    typedef int_vector<>::size_type size_type;
    const size_type M = 255;	// limit for values represented in the small LCP part
    size_type buf_len = 1000000;
    lcp_buf.buffersize(buf_len);
    bwt_buf.buffersize(buf_len);
    size_type n = lcp_buf.size();

    osfstream sml_lcp_out(small_lcp_file, std::ios::out | std::ios::trunc);
    uint64_t bit_size = 8*n;
    sml_lcp_out.write((char*) &bit_size, sizeof(bit_size));

    osfstream big_lcp_out(big_lcp_file, std::ios::out | std::ios::trunc);

    size_type fc_cnt = 0; // number of lcp values at the first child r
    size_type fc_cnt_big = 0; // number of lcp values at the first child which are big and not reducible
    size_type fc_cnt_big2 = 0;
    sorted_multi_stack_support vec_stack(n); // occupies 2n bits
    bit_vector is_big_and_not_reducable(n, 0); // initialized with 0s
    bool is_one_big_and_not_reducable = false; // all positions have to be reducible

    size_type y, max_lcp=0;
    uint64_t last_bwti=0, val;
    for (size_type i=0, x; i < n; ++i) {
        x = lcp_buf[i];
        is_one_big_and_not_reducable = false;

        while (!vec_stack.empty() and x < vec_stack.top()) {
            y = vec_stack.top();
            is_one_big_and_not_reducable |= is_big_and_not_reducable[vec_stack.size()-1];
            if (vec_stack.pop()) { // if y was the last copy of y on the stack
                if (y > M-2) {
                    if (is_one_big_and_not_reducable) {
                        val = M;
                        big_lcp_out.write((char*)&y, sizeof(y));
                        ++fc_cnt_big;
                        if (y > max_lcp) max_lcp = y;
                    } else {
                        val = M-1;
                        ++fc_cnt_big2;
                    }
                } else {
                    val = y;
                }
                sml_lcp_out.write((const char*)&val, 1);
                ++fc_cnt;
                is_one_big_and_not_reducable = false;
            }
        }
        if (x > M-2 and(0 == i or last_bwti != bwt_buf[i] or x % t_dens == 0)) {
            is_big_and_not_reducable[vec_stack.size()] = 1;
        } else {
            is_big_and_not_reducable[vec_stack.size()] = 0;
        }
        vec_stack.push(x);
        last_bwti = bwt_buf[i];
    }

    while (!vec_stack.empty()) {
        y = vec_stack.top();
        if (vec_stack.pop()) {
            if (y > M-2) {
                if (is_big_and_not_reducable[vec_stack.size()]) {
                    val = M;
                    big_lcp_out.write((char*)&y, sizeof(y));
                    ++fc_cnt_big;
                    if (y > max_lcp) max_lcp = y;
                } else {
                    val = M-1;
                    ++fc_cnt_big2;
                }
            } else {
                val = y;
            }
            sml_lcp_out.write((const char*)&val, 1);
            ++fc_cnt;
        }
    }
    // write number of elements of sml_lcp into the out file stream
    sml_lcp_out.seekp(0);
    bit_size = 8*fc_cnt;
    sml_lcp_out.write((char*) &bit_size, sizeof(bit_size));
    sml_lcp_out.close();

    big_lcp_out.close();
    isfstream big_lcp_in(big_lcp_file);
    big_lcp.width(bits::hi(max_lcp)+1);
    big_lcp.resize(fc_cnt_big);

    for (size_type i=0; i<fc_cnt_big; ++i) {
        big_lcp_in.read((char*)&y, sizeof(y));
        big_lcp[i] = y;
    }
}

} // end namespace
#endif
