#ifndef INCLUDED_SDSL_SUFFIX_TREE_HELPER
#define INCLUDED_SDSL_SUFFIX_TREE_HELPER

#include <stdint.h>
#include <cstdlib>
#include <cassert>
#include <stack>
#include "sorted_multi_stack_support.hpp"
#include "sorted_stack_support.hpp"
#include "iterators.hpp"

namespace sdsl
{


template <class t_cst>
class cst_node_child_proxy_iterator : public std::iterator<std::forward_iterator_tag, typename t_cst::node_type>
{
    public:
        using node_type = typename t_cst::node_type;
        using value_type = node_type;
        using const_reference = const node_type;
        using iterator_type =  cst_node_child_proxy_iterator<t_cst>;
    private:
        const t_cst* m_cst;
        node_type m_cur_node;
    public:
        cst_node_child_proxy_iterator() : m_cst(nullptr) {};
        cst_node_child_proxy_iterator(const t_cst* cst,node_type v) : m_cst(cst) , m_cur_node(v) {}
        cst_node_child_proxy_iterator(const iterator_type& it): m_cst(it.m_cst), m_cur_node(it.m_cur_node) {}
    public:
        const_reference operator*() const
        {
            return m_cur_node;
        }
        iterator_type& operator++()
        {
            m_cur_node = m_cst->sibling(m_cur_node);
            return *this;
        }
        iterator_type operator++(int)
        {
            iterator_type it = *this;
            ++(*this);
            return it;
        }
        bool operator==(const iterator_type& it)const
        {
            return it.m_cur_node == m_cur_node;
        }
        bool operator!=(const iterator_type& it)const
        {
            return !(*this==it);
        }
};

template <class t_cst>
class cst_node_child_proxy
{
    public: // types
        using iterator_type = cst_node_child_proxy_iterator<t_cst>;
        using node_type = typename t_cst::node_type;
        using size_type = typename t_cst::size_type;
    private: // data
        node_type m_parent;
        const t_cst* m_cst;
    public: // constructors
        cst_node_child_proxy() = delete;
        explicit cst_node_child_proxy(const t_cst* cst,node_type v) : m_parent(v) , m_cst(cst) {};
        cst_node_child_proxy(const cst_node_child_proxy& p) : m_parent(p.m_parent) , m_cst(p.m_cst) {};
    public: // methods
        node_type operator[](size_type i) const { return m_cst->select_child(m_parent,i+1); } // enumeration starts with 1 not 0
        size_type size() { return m_cst->degree(m_parent); }
        iterator_type begin() const { return iterator_type(m_cst,m_cst->select_child(m_parent,1)); }
        iterator_type end() const { return iterator_type(m_cst,m_cst->root()); }
};

//! Calculate the balanced parentheses of the Super-Cartesian tree, described in Ohlebusch and Gog (SPIRE 2009).
/*! \param vec Random access container for which the Super-Cartesian tree representation should be calculated.
 *             The value_type of vec should be an unsigned integer type.
 *  \param bp Reference to the balanced parentheses sequence which represents the Super-Cartesian tree.
 *  \param minimum Specifies if the higher levels contains minima or maxima. Default is maxima.
 *  \par Time complexity
 *       \f$ \Order{2n} \f$, where \f$ n=\f$vec.size()
 *  \par Space complexity
 *       \f$ \Order{n \cdot \log n } \f$ bits.
 */
template<class t_rac>
void construct_supercartesian_tree_bp(const t_rac& vec, bit_vector& bp, const bool minimum=true)
{
    typedef typename t_rac::size_type size_type;
    bp.resize(2*vec.size());      // resize bit vector for balanaced parantheses to 2 n bits
    util::set_to_value(bp, 0);
    std::stack<typename t_rac::value_type> vec_stack;

    size_type k=0;
    for (size_type i=0; i < vec.size(); ++i) {
        typename t_rac::value_type l = vec[i];
        if (minimum) {
            while (vec_stack.size() > 0 and l < vec_stack.top()) {
                vec_stack.pop(); ++k; /*bp[k++] = 0; bp is already initialized to zero*/ // writing a closing parenthesis
            }

        } else {
            while (vec_stack.size() > 0 and l > vec_stack.top()) {
                vec_stack.pop(); ++k; /*bp[k++] = 0; bp is already initialized to zero*/ // writing a closing parenthesis
            }
        }
        vec_stack.push(l);
        bp[k++] = 1; // writing an opening  parenthesis
    }
    while (vec_stack.size() > 0) {
        vec_stack.pop();
        bp[k++] = 0; // writing a closing parenthesis
    }
    assert(k == 2*vec.size());
}

//! Calculate the balanced parentheses of the Super-Cartesian tree, described in Ohlebusch and Gog (SPIRE 2009).
/*! \param vec Random access container for which the Super-Cartesian tree representation should be calculated.
 *             The value_type of vec should be an unsigned integer type.
 *  \param minimum Specifies if the higher levels contains minima or maxima. Default is maxima.
 *  \return The balanced parentheses sequence representing the Super-Cartesian tree.
 *  \par Time complexity
 *       \f$ \Order{2n} \f$, where \f$ n=\f$vec.size()
 *  \par Space complexity
 *       \f$\Order{n}\f$ bits
 */
template<class t_rac>
bit_vector
construct_supercartesian_tree_bp_succinct(const t_rac& vec, const bool minimum=true)
{
    typedef typename t_rac::size_type size_type;
    bit_vector bp(2*vec.size(), 0); // initialize result
    if (vec.size() > 0) {
        sorted_stack_support vec_stack(vec.size());

        size_type k=0;
        if (minimum) {
            bp[k++] = 1;
            for (size_type i=1; i < vec.size(); ++i) {
                if (vec[i] < vec[i-1]) {
                    ++k;
                    while (vec_stack.size() > 0 and vec[i] < vec[vec_stack.top()]) {
                        vec_stack.pop(); ++k; // writing a closing parenthesis, bp is already initialized to zero
                    }
                } else {
                    vec_stack.push(i-1); // "lazy stack" trick: speed-up approx. 25%
                }
                bp[k++] = 1; // writing an opening  parenthesis
            }
        } else {
            // no "lazy stack" trick used here
            for (size_type i=0; i < vec.size(); ++i) {
                while (vec_stack.size() > 0 and vec[i] > vec[vec_stack.top()]) {
                    vec_stack.pop(); ++k; /*bp[k++] = 0; bp is already initialized to zero*/ // writing a closing parenthesis
                }
                vec_stack.push(i);
                bp[k++] = 1; // writing an opening  parenthesis
            }
        }
    }
    return bp;
}

//! Calculate the balanced parentheses of the Super-Cartesian tree, described in Ohlebusch and Gog (SPIRE 2009).
/*! \param lcp_buf int_vector_buffer of the LCP Array for which the Super-Cartesian tree representation should be calculated.
 *             The value_type of vec should be an unsigned integer type.
 *  \param minimum Specifies if the higher levels contains minima or maxima. Default is maxima.
 *  \return The balanced parentheses sequence representing the Super-Cartesian tree.
 *  \par Time complexity
 *       \f$ \Order{2n} \f$, where \f$ n=\f$vec.size()
 *  \par Space complexity
 *       \f$\Order{2n}\f$ bits, by the multi_stack_support
 *  \pre
 *  The largest value in lcp_buf has to be smaller than lcp_buf.size().
 */
template<uint8_t t_width>
bit_vector
construct_supercartesian_tree_bp_succinct(int_vector_buffer<t_width>& lcp_buf, const bool minimum=true)
{
    typedef bit_vector::size_type size_type;
    bit_vector bp(2*lcp_buf.size(), 0); // initialize result
    if (lcp_buf.size() > 0)	{
        sorted_multi_stack_support vec_stack(lcp_buf.size());

        size_type k=0;
        if (minimum) {
            bp[k++] = 1;
            size_type last = lcp_buf[0];
            for (size_type i=1, x; i < lcp_buf.size(); ++i) {
                x = lcp_buf[i];
                if (x < last) {
                    ++k; // writing a closing parenthesis for last
                    while (!vec_stack.empty() and x < vec_stack.top()) {
                        vec_stack.pop(); ++k; // writing a closing parenthesis, bp is already initialized to zeros
                    }
                } else {
                    vec_stack.push(last); // "lazy stack" trick: speed-up about 25 %
                }
                bp[k++] = 1; // writing an opening parenthesis
                last = x;
            }
        } else {
            // no "lazy stack" trick use here
            for (size_type i=0, x; i < lcp_buf.size(); ++i) {
                x = lcp_buf[i];
                while (!vec_stack.empty() and x > vec_stack.top()) {
                    vec_stack.pop(); ++k; // writing a closing parenthesis, bp is already initialized to zeros
                }
                vec_stack.push(x);
                bp[k++] = 1; // writing an opening parenthesis
            }
        }
    }
    return bp;
}

//! Calculate the balanced parentheses of the Super-Cartesian tree, described in Ohlebusch and Gog (SPIRE 2009) and the first_child bit_vector
/*! \param lcp_buf int_vector_buffer for the lcp array for which the Super-Cartesian tree representation should be calculated.
 *             The value_type of vec should be an unsigned integer type.
 *  \param bp Reference to the balanced parentheses sequence which represents the Super-Cartesian tree.
 *  \param bp_fc Reference to the first child bit_vector of bp.
 *  \param minimum Specifies if the higher levels contains minima or maxima. Default is maxima.
 *  \par Time complexity
 *       \f$ \Order{2n} \f$, where \f$ n=\f$vec.size()
 *  \par Space complexity
 *       \f$\Order{2n}\f$ bits, by the multi_stack_support
 */
template<uint8_t t_width>
bit_vector::size_type
construct_supercartesian_tree_bp_succinct_and_first_child(int_vector_buffer<t_width>& lcp_buf, bit_vector& bp, bit_vector& bp_fc, const bool minimum=true)
{
    typedef bit_vector::size_type size_type;
    size_type n = lcp_buf.size();
    bp.resize(2*n);      // resize bit vector for balanced parentheses to 2 n bits
    bp_fc.resize(n);
    if (n == 0)	// if n == 0 we are done
        return 0;
    size_type fc_cnt=0; // first child counter
    util::set_to_value(bp, 0);
    util::set_to_value(bp_fc, 0);
    sorted_multi_stack_support vec_stack(n);

    size_type k=0;
    size_type k_fc=0; // first child index
    if (minimum) {
        // no "lazy stack" trick used here
        for (size_type i=0, x; i < n; ++i) {
            x = lcp_buf[i];
            while (!vec_stack.empty() and x < vec_stack.top()) {
                if (vec_stack.pop()) {
                    bp_fc[k_fc] = 1;
                    ++fc_cnt;
                }
                ++k; // writing a closing parenthesis, bp is already initialized to zeros
                ++k_fc; // write a bit in first_child
            }
            vec_stack.push(x);
            bp[k++] = 1; // writing an opening parenthesis
        }

    } else {
        // no "lazy stack" trick used here
        for (size_type i=0, x; i < n; ++i) {
            x = lcp_buf[i];
            while (!vec_stack.empty() and x > vec_stack.top()) {
                if (vec_stack.pop()) {
                    bp_fc[k_fc] = 1;
                    ++fc_cnt;
                }
                ++k; // writing a closing parenthesis, bp is already initialized to zeros
                ++k_fc; // write a bit in first_child
            }
            vec_stack.push(x);
            bp[k++] = 1; // writing an opening parenthesis
        }
    }
    while (!vec_stack.empty()) {
        if (vec_stack.pop()) {
            bp_fc[k_fc] = 1;
            ++fc_cnt;
        }
        // writing a closing parenthesis in bp, not necessary as bp is initialized with zeros
        ++k;
        ++k_fc;
    }
    return fc_cnt;
}


// Gets ISA[SA[idx]+d]
// d = depth of the character 0 = first position
template<class t_csa>
typename t_csa::size_type get_char_pos(typename t_csa::size_type idx, typename t_csa::size_type d, const t_csa& csa)
{
    if (d == 0)
        return idx;
    // if we have to apply \f$\LF\f$ or \f$\Phi\f$ more
    // than 2*d times to calc csa(csa[idx]+d), we opt to
    // apply \f$ \Phi \f$ d times
    if (csa.sa_sample_dens + csa.isa_sample_dens > 2*d+2) {
        for (typename t_csa::size_type i=0; i < d; ++i)
            idx = csa.psi[idx];
        return idx;
    }
    return csa.isa[csa[idx] + d];
}


// has_id<X>::value is true if class X has
// implement method id
// Adapted solution from jrok's proposal:
// http://stackoverflow.com/questions/87372/check-if-a-class-has-a-member-function-of-a-given-signature
template<typename t_wt>
struct has_id {
    template<typename T>
    static constexpr auto check(T*)
    -> typename
    std::is_same<
    decltype(std::declval<T>().id(
                 std::declval<typename T::node_type&>()
             )),
             typename T::size_type>::type {return std::true_type();}
             template<typename>
    static constexpr std::false_type check(...) {return std::false_type();}
    typedef decltype(check<t_wt>(nullptr)) type;
    static constexpr bool value = type::value;
};




}

#endif
