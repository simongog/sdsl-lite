#ifndef INCLUDED_SDSL_SUPPORT_LCP_TREE
#define INCLUDED_SDSL_SUPPORT_LCP_TREE

#include "lcp.hpp"
#include "util.hpp"
#include "algorithms_for_compressed_suffix_trees.hpp"
#include "lcp_wt.hpp"
#include <iostream>
#include <string>

namespace sdsl
{


/*! This class composes a virtual LCP array from a LCP arrays which is in suffix array order
 * (e.g. lcp_byte or lcp_bitcompressed) and a CST.
 *    The time consumption of the []-operator depends on:
 *    - The time consumption of the tlcp_idx function of the CST
 *    - The access time to the suffix array ordered LCP array
 *
 * \tparam t_lcp Type of the underlying LCP array. Must be an suffix array ordered one.
 * \tparam t_cst Type of the underlying CST.
 */
template<class t_lcp, class t_cst>
class _lcp_support_tree
{
    public:
        typedef typename t_lcp::value_type                      value_type;      // STL Container requirement
        typedef random_access_const_iterator<_lcp_support_tree> const_iterator;  // STL Container requirement
        typedef const_iterator                                  iterator;         // STL Container requirement
        typedef const value_type                                const_reference;
        typedef const_reference                                 reference;
        typedef const_reference*                                pointer;
        typedef const pointer                                   const_pointer;
        typedef typename t_lcp::size_type                       size_type;         // STL Container requirement
        typedef typename t_lcp::difference_type                 difference_type; // STL Container requirement

        typedef lcp_tree_compressed_tag                         lcp_category;

        enum { fast_access = 0,
               text_order  = t_lcp::text_order,
               sa_order    = t_lcp::sa_order
             };

        template<class CST>  // template inner class which is used in CSTs to parametrize lcp classes
        class type           // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef _lcp_support_tree lcp_type;
        };

    private:

        const t_cst* m_cst;
        t_lcp        m_lcp;

        void copy(const _lcp_support_tree& lcp_c) {
            m_cst = lcp_c.m_cst;
            m_lcp = lcp_c.m_lcp; // works for lcp_bitcompressed and lcp_kurtz
        }

    public:

        //! Default constructor
        _lcp_support_tree() {}

        // Destructor
        ~_lcp_support_tree() {}

        //! Copy constructor
        _lcp_support_tree(const _lcp_support_tree& lcp) {
            m_cst = lcp.m_cst;
            m_lcp = lcp.m_lcp;
        }

        //! Construct the lcp array from an int_vector_file_buffer of the lcp array
        /*!
         *  \param lcp_buf An int_vector_file_buf of the lcp array
         *  \param fc_bpss A pointer to the CST
        */
        template<uint8_t int_width>
        _lcp_support_tree(int_vector_file_buffer<int_width>& lcp_buf, const t_cst* cst = NULL) {
            m_cst = cst;
            std::string id =  util::to_string(util::pid())+"_"+util::to_string(util::id()).c_str() + "_fc_lcp";
            {
                int_vector<int_width> temp_lcp;
                algorithm::construct_first_child_lcp(lcp_buf, temp_lcp, (int_vector_size_type) 0);
                // TODO: store lcp values directly to disk
                util::store_to_file(temp_lcp, id);
            }
            {
                int_vector_file_buffer<int_width> temp_lcp_buf(id);
                m_lcp = t_lcp(temp_lcp_buf); // works for lcp_kurtz, lcp_wt and lcp_bitcompressed
            }
            std::remove(id.c_str());
        }

        size_type size()const {
            return m_cst->size();
        }

        void set_cst(const t_cst* cst) {
            m_cst = cst;
        }

        static size_type max_size() {
            return t_lcp::max_size();
        }

        size_type empty()const {
            return m_lcp.empty();
        }

        void swap(_lcp_support_tree& lcp_c) {
            m_lcp.swap(lcp_c.m_lcp);
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
            return m_lcp[ m_cst->tlcp_idx(i) ];
        }

        //! Assignment Operator.
        _lcp_support_tree& operator=(const _lcp_support_tree& lcp_c) {
            if (this != &lcp_c) {
                copy(lcp_c);
            }
            return *this;
        }

        //! Serialize to a stream.
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            size_type written_bytes = 0;
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            written_bytes += m_lcp.serialize(out, child, "lcp");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load from a stream.
        void load(std::istream& in, const t_cst* cst=NULL) {
            m_lcp.load(in); // works for lcp_byte and lcp_bitcompressed
            m_cst = cst;
        }
};

//! Helper class which provides _lcp_support_tree the context of a CST.
template<class t_lcp = lcp_wt<> >
class lcp_support_tree
{
    public:
        template<class t_cst>  // template inner class which is used in CSTs to parametrize lcp classes
        class type           // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef _lcp_support_tree<t_lcp, t_cst> lcp_type;
        };
};

} // end namespace

#endif
