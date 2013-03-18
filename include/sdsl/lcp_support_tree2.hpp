#ifndef INCLUDED_SDSL_SUPPORT_TREE2
#define INCLUDED_SDSL_SUPPORT_TREE2

#include "lcp.hpp"
#include "util.hpp"
#include "algorithms_for_compressed_suffix_trees.hpp"
#include "rank_support_v.hpp"
#include "wt_huff.hpp"
#include <iostream>
#include <string>

namespace sdsl
{


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
        typedef int_vector<>::value_type                         value_type;      // STL Container requirement
        typedef random_access_const_iterator<_lcp_support_tree2> const_iterator;  // STL Container requirement
        typedef const_iterator                                   iterator;         // STL Container requirement
        typedef const value_type                                 const_reference;
        typedef const_reference                                  reference;
        typedef const_reference*                                 pointer;
        typedef const pointer                                    const_pointer;
        typedef int_vector<>::size_type                          size_type;         // STL Container requirement
        typedef int_vector<>::difference_type                    difference_type; // STL Container requirement
        typedef t_cst                                            cst_type;
        typedef wt_huff<bit_vector, rank_support_v5<>,
                select_support_scan<1>,
                select_support_scan<0> >                 small_lcp_type;

        typedef lcp_tree_and_lf_compressed_tag                   lcp_category;

        enum { fast_access = 0,
               text_order = 0,
               sa_order = 0
             };

        template<class CST>  // template inner class which is used in CSTs to parametrize lcp classes
        class type           // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef _lcp_support_tree2 lcp_type;
        };

    private:
        const cst_type*    m_cst;
        small_lcp_type  m_small_lcp; // vector for lcp values < 254
        int_vector<> m_big_lcp;      // vector for lcp values >= 254

        void copy(const _lcp_support_tree2& lcp_c) {
            m_small_lcp = lcp_c.m_small_lcp;
            m_big_lcp = lcp_c.m_big_lcp;
            m_cst = lcp_c.m_cst;
        }

    public:

        //! Default constructor
        _lcp_support_tree2() {}

        //! Copy constructor
        _lcp_support_tree2(const _lcp_support_tree2& lcp) {
            copy(lcp);
        }

        //! Construct the lcp array from an lcp array
        /*! \param lcp_buf Buffer to the uncompressed lcp array
         *  \param sa_buf
         */
        template<uint8_t int_width>
        _lcp_support_tree2(int_vector_file_buffer<int_width>& lcp_buf,
                           int_vector_file_buffer<t_cst::csa_type::alphabet_type::int_width>& bwt_buf,
                           const cst_type* cst = NULL) {
            m_cst = cst;
            std::string small_lcp_file_name =  util::to_string(util::pid())+"_"+util::to_string(util::id()).c_str() + "_fc_lf_lcp_sml";
            std::string big_lcp_file_name =  util::to_string(util::pid())+"_"+util::to_string(util::id()).c_str() + "_fc_lf_lcp_big";

            algorithm::construct_first_child_and_lf_lcp<t_dens>(lcp_buf, bwt_buf, small_lcp_file_name, big_lcp_file_name, m_big_lcp);
            // construct wavelet tree huffman from file buffer
            int_vector_file_buffer<8> small_lcp_buf(small_lcp_file_name);
            util::assign(m_small_lcp, small_lcp_type(small_lcp_buf, small_lcp_buf.int_vector_size));
            std::remove(small_lcp_file_name.c_str());
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
         * \par Note
         *  Required for the STL Random Access Container Concept.
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
                i = m_cst->csa.psi(i); // i = LF[i]    // (*m_psi)(i);
                ++offset; // goto lcp value, which is one bigger
                goto start;
            } else { // if lcp value is >= 254 and (not reducable or sampled)
                return m_big_lcp[m_small_lcp.rank(idx ,255)] - offset;
            }
        }

        //! Assignment Operator.
        /*!
         *    Required for the Assignable Concept of the STL.
         */
        _lcp_support_tree2& operator=(const _lcp_support_tree2& lcp_c) {
            if (this != &lcp_c) {
                copy(lcp_c);
            }
            return *this;
        }

        //! Serialize to a stream.
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_small_lcp.serialize(out, child, "small_lcp");
            written_bytes += m_big_lcp.serialize(out, child, "large_lcp");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load from a stream.
        void load(std::istream& in, const t_cst* cst=NULL) {
            m_small_lcp.load(in);
            m_big_lcp.load(in);
            m_cst = cst;
        }
};

//! Helper class which provides _lcp_support_tree2 the context of a CST.
template<uint32_t t_dens=16>
class lcp_support_tree2
{
    public:
        template<class t_cst> // template inner class which is used in CSTs to parametrize lcp classes
        class type            // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef _lcp_support_tree2<t_dens, t_cst> lcp_type;
        };
};



} // end namespace

#endif
