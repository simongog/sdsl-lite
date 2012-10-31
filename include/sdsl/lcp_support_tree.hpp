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


/*! This class composes a virtual lcp array from a lcp arrays which is ordered in suffix array order (e.g. lcp_kurtz or lcp_bitcompressed) and a compressed suffix tree like cst_sct3 or cst_sada.
 *	The time consumption of the []-operator depends on:
 *    - The time consumption of the tlcp_idx function of the CST
 *    - The access time to the suffix array ordered LCP array
 *
 */
template<class Lcp, class Cst>
class _lcp_support_tree
{
    public:
        typedef typename Lcp::value_type						value_type; 	 // STL Container requirement
        typedef random_access_const_iterator<_lcp_support_tree>	const_iterator;  // STL Container requirement
        typedef const_iterator									iterator;		 // STL Container requirement
        typedef const value_type								const_reference;
        typedef const_reference									reference;
        typedef const_reference*								pointer;
        typedef const pointer									const_pointer;
        typedef typename Lcp::size_type							size_type;		 // STL Container requirement
        typedef typename Lcp::difference_type					difference_type; // STL Container requirement

        typedef lcp_tree_compressed_tag					lcp_category;

        enum { fast_access = 0,
               text_order = Lcp::text_order,
               sa_order = Lcp::sa_order
             };

        template<class CST>  // template inner class which is used in CSTs to parametrize lcp classes
        class type           // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef _lcp_support_tree lcp_type;
        };

    private:

        const Cst* m_cst;
        Lcp        m_lcp;

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

        //! Construct the lcp array  from a text and the corresponding suffix array.
        template<class Text, class Sa>
        _lcp_support_tree(const Text& text, const Sa& sa, const Cst* cst);
        // TODO: implement



        //! Construct the lcp array from an int_vector_file_buffer of the lcp array
        /*!
         *	\param lcp_buf 	An int_vector_file_buf of the lcp array
         *	\param fc_bpss 	A pointer to balanced parentheses support of the balanced parentheses representation of the Super-Cartesian Tree of the lcp array
         *  \param fcr		A pointer to the rank support of the first child bit_vector.
        */
        template<uint8_t int_width>
        _lcp_support_tree(int_vector_file_buffer<int_width>& lcp_buf, const Cst* cst = NULL) {
            m_cst = cst;
            std::string id =  util::to_string(util::get_pid())+"_"+util::to_string(util::get_id()).c_str() + "_fc_lcp";
            {
                int_vector<int_width> temp_lcp;
                algorithm::construct_first_child_lcp(lcp_buf, temp_lcp, (int_vector_size_type) 0);
                // TODO: store lcp values directly to disk
                util::store_to_file(temp_lcp, id.c_str());
            }
            {
                int_vector_file_buffer<int_width> temp_lcp_buf(id.c_str());
                m_lcp = Lcp(temp_lcp_buf); // works for lcp_kurtz, lcp_wt and lcp_bitcompressed
            }
            std::remove(id.c_str());
        }

        size_type size()const {
            return m_cst->size(); // corresponds to the length of the
            // original lcp array
        }

        void set_cst(const Cst* cst) {
            m_cst = cst;
        }

        static size_type max_size() {
            return Lcp::max_size();
        }

        size_type empty()const {
            return m_lcp.empty();
        }

        void swap(_lcp_support_tree& lcp_c) {
            m_lcp.swap(lcp_c.m_lcp);
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
            return m_lcp[ m_cst->tlcp_idx(i) ];
        }

        //! Assignment Operator.
        /*!
         *	Required for the Assignable Concept of the STL.
         */
        _lcp_support_tree& operator=(const _lcp_support_tree& lcp_c) {
            if (this != &lcp_c) {
                copy(lcp_c);
            }
            return *this;
        }

        //! Serialize to a stream.
        /*! \param out Outstream to write the data structure.
         *  \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const {
            size_type written_bytes = 0;
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            written_bytes += m_lcp.serialize(out, child, "lcp");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load from a stream.
        /*! \param in Inputstream to load the data structure from.
         *	\param fc_bpss A pointer to the balanced parentheses support for the first child bit_vector
         *	\param fcr	A pointer to the rank support for the first child bit_vector
         */
        void load(std::istream& in, const Cst* cst=NULL) {
            m_lcp.load(in); // works for lcp_kurtz and lcp_bitcompressed
            m_cst = cst;
        }
};

//! Helper class which provides _lcp_support_tree the context of a CST.
template<class Lcp = lcp_wt<> >
class lcp_support_tree
{
    public:
        template<class Cst>  // template inner class which is used in CSTs to parametrize lcp classes
        class type           // with information about the CST. Thanks Stefan Arnold! (2011-03-02)
        {
            public:
                typedef _lcp_support_tree<Lcp, Cst> lcp_type;
        };
};



} // end namespace

#endif
