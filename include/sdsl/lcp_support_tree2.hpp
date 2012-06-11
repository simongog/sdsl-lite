#ifndef INCLUDED_SDSL_SUPPORT_TREE2
#define INCLUDED_SDSL_SUPPORT_TREE2

#include "lcp.hpp"
#include "util.hpp"
#include "algorithms_for_compressed_suffix_trees.hpp"
#include "rank_support_v.hpp"
#include "select_support_dummy.hpp"
#include "wt_huff.hpp"
#include <iostream>
#include <string>

namespace sdsl
{


/*! An lcp array class for cst_sct3 and cst_sada.
 *	The time consumption of the []-operator depends on:
 *    - The time consumption of the []-operation of the wt_huff
 *    - The time consumption of the LF calculation of the underlying CSA of the CST
 *    - The time consumption of the tlcp_idx function of the CST
 *
 */
template<uint32_t SampleDens, class Cst>
class _lcp_support_tree2
{
    public:
        typedef int_vector<>::value_type				value_type; 	 // STL Container requirement
        typedef random_access_const_iterator<_lcp_support_tree2>	const_iterator;  // STL Container requirement
        typedef const_iterator									iterator;		 // STL Container requirement
        typedef const value_type								const_reference;
        typedef const_reference									reference;
        typedef const_reference*								pointer;
        typedef const pointer									const_pointer;
        typedef int_vector<>::size_type							size_type;		 // STL Container requirement
        typedef int_vector<>::difference_type					difference_type; // STL Container requirement
        typedef Cst												cst_type;

        typedef lcp_tree_and_lf_compressed_tag					lcp_category;

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
        const cst_type*	m_cst;
        typedef select_support_dummy tDummySS;
        wt_huff<bit_vector, rank_support_v5<>, tDummySS, tDummySS > m_small_lcp; // vector for lcp values < 254
        int_vector<> m_big_lcp;  // vector for lcp values >= 254

        void copy(const _lcp_support_tree2& lcp_c) {
            m_small_lcp = lcp_c.m_small_lcp;
            m_big_lcp = lcp_c.m_big_lcp;
            m_cst = lcp_c.m_cst;
        }

    public:

        //! Default constructor
        _lcp_support_tree2() {}

        // Destructor
        ~_lcp_support_tree2() {}

        //! Copy constructor
        _lcp_support_tree2(const _lcp_support_tree2& lcp) {
            copy(lcp);
        }

        /*!
         * \see _lcp_support_tree2
         */
        template<class Text, class Sa>
        void construct(const Text& text, const Sa& sa, const cst_type* cst);


        //! Construct the lcp array from an lcp array
        /*! \param lcp_buf Buffer to the uncompressed lcp array
         *  \param sa_buf
         */
        template<uint8_t int_width, class size_type_class>
        void construct(int_vector_file_buffer<int_width, size_type_class>& lcp_buf,
                       int_vector_file_buffer<8, size_type_class>& bwt_buf,
                       const cst_type* cst = NULL) {
            m_cst = cst;
            std::string small_lcp_file_name =  util::to_string(util::get_pid())+"_"+util::to_string(util::get_id()).c_str() + "_fc_lf_lcp_sml";
            std::string big_lcp_file_name =  util::to_string(util::get_pid())+"_"+util::to_string(util::get_id()).c_str() + "_fc_lf_lcp_big";

            algorithm::construct_first_child_and_lf_lcp<SampleDens>(lcp_buf, bwt_buf,
                    small_lcp_file_name, big_lcp_file_name, m_big_lcp);
            // construct wavelet tree huffman from file buffer
            int_vector_file_buffer<8> small_lcp_buf(small_lcp_file_name.c_str());
            m_small_lcp.construct(small_lcp_buf, small_lcp_buf.int_vector_size);
            std::remove(small_lcp_file_name.c_str());
        }

        void set_cst(const cst_type* cst) {
            m_cst = cst;
        }

        size_type size()const {
            return m_cst->size(); // corresponds to the length of the
            // original lcp array
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
        /*! Required for the STL Container Concept.
         *  \sa end
         */
        const_iterator begin()const {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        /*! Required for the STL Container Concept.
         *  \sa begin.
         */
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
         *	Required for the Assignable Concept of the STL.
         */
        _lcp_support_tree2& operator=(const _lcp_support_tree2& lcp_c) {
            if (this != &lcp_c) {
                copy(lcp_c);
            }
            return *this;
        }

        //! Equality Operator
        /*! Two Instances of _lcp_support_tree2 are equal if
         *  all their members are equal.
         *  \par Required for the Equality Comparable Concept of the STL.
         *  \sa operator!=
         */
        bool operator==(const _lcp_support_tree2& lcp_c)const {
            if (this == &lcp_c)
                return true;
            return m_cst == lcp_c.m_cst and
                   m_small_lcp == lcp_c.m_small_lcp and
                   m_big_lcp == lcp_c.m_big_lcp;
        }

        //! Unequality Operator
        /*! Two Instances of _lcp_support_tree2 are equal if
         *  not all their members are equal.
         *  \par Required for the Equality Comparable Concept of the STL.
         *  \sa operator==
         */
        bool operator!=(const _lcp_support_tree2& lcp_c)const {
            return !(*this == lcp_c);
        }

        //! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_small_lcp.serialize(out, child, "small_lcp");
            written_bytes += m_big_lcp.serialize(out, child, "large_lcp");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load from a stream.
        /*! \param in Inputstream to load the data structure from.
         *	\param fc_bpss A pointer to the balanced parentheses support for the first child bit_vector
         *	\param fcr	A pointer to the rank support for the first child bit_vector
         */
        void load(std::istream& in, const Cst* cst=NULL) {
            m_small_lcp.load(in);
            m_big_lcp.load(in);
            m_cst = cst;
        }
};

//! Helper class which provides _lcp_support_tree2 the context of a CST.
template<uint32_t SampleDens=16>
class lcp_support_tree2
{
    public:
        template<class Cst>  // template inner class which is used in CSTs to parametrize lcp classes
        class type           // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef _lcp_support_tree2<SampleDens, Cst> lcp_type;
        };
};



} // end namespace

#endif
