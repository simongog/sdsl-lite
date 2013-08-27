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
class cst_node_child_proxy_iterator
{
    public:
        using node_type = typename t_cst::node_type;
        using const_reference = const node_type;
        using iterator_type =  cst_node_child_proxy_iterator<t_cst>;
    private:
        const t_cst& m_cst;
        node_type current_node;
    public:
        cst_node_child_proxy_iterator() = delete;
        cst_node_child_proxy_iterator(const t_cst& cst,const node_type& v) : m_cst(cst) , current_node(v) {}
    public:
        const_reference operator*() const {
            return current_node;
        }
        iterator_type& operator++() {
            current_node = m_cst.sibling(current_node);
            return *this;
        }
        bool operator==(const iterator_type& it)const {
            return it.current_node == current_node;
        }
        bool operator!=(const iterator_type& it)const {
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
        const node_type& m_parent;
        const t_cst& m_cst;
    public: // constructors
        cst_node_child_proxy() = delete;
        explicit cst_node_child_proxy(const t_cst& cst,const node_type& v) : m_parent(v) , m_cst(cst) {};
    public: // methods
        node_type operator[](size_type i) const { return m_cst.select_child(m_parent,i+1); } // enumeration starts with 1 not 0
        size_type size() { return m_cst.degree(m_parent); }
        iterator_type begin() const { return iterator_type(m_cst,m_cst.select_child(m_parent,1)); }
        iterator_type end() const { return iterator_type(m_cst,m_cst.root()); }
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
template<class RandomAccessContainer>
void construct_supercartesian_tree_bp(const RandomAccessContainer& vec, bit_vector& bp, const bool minimum=true)
{
    typedef typename RandomAccessContainer::size_type size_type;
    bp.resize(2*vec.size());      // resize bit vector for balanaced parantheses to 2 n bits
    util::set_to_value(bp, 0);
    std::stack<typename RandomAccessContainer::value_type> vec_stack;

    size_type k=0;
    for (size_type i=0; i < vec.size(); ++i) {
        typename RandomAccessContainer::value_type l = vec[i];
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
 *  \param bp Reference to the balanced parentheses sequence which represents the Super-Cartesian tree.
 *  \param minimum Specifies if the higher levels contains minima or maxima. Default is maxima.
 *  \par Time complexity
 *       \f$ \Order{2n} \f$, where \f$ n=\f$vec.size()
 *  \par Space complexity
 *       \f$\Order{n}\f$ bits, by the stack_support described in the paper "Optimal Succinctness For Range Minimum Queries" of Johannes Fischer.
 */
// TODO: sorted_multi_stack_support einbauen, RandomAccessContainer durch int_vector_buffer ersetzen
template<class RandomAccessContainer>
void construct_supercartesian_tree_bp_succinct(const RandomAccessContainer& vec, bit_vector& bp, const bool minimum=true)
{
    typedef typename RandomAccessContainer::size_type size_type;
    bp.resize(2*vec.size());      // resize bit vector for balanced parentheses to 2 n bits
    if (vec.size() > 0) {
        util::set_to_value(bp, 0);
        sorted_stack_support vec_stack(vec.size()); // <- ist das ein Problem fuer int_vector_buffer

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
                    vec_stack.push(i-1); // "lazy stack" trick: speed-up ca. 25%
                }
                bp[k++] = 1; // writing an opening  parenthesis
            }
            /*
            vec_stack.push(0);
            bp[k++] = 1;
            for(size_type i=1,j, start_run=1; i < vec.size(); ++i){
            	if( vec[i] < vec[i-1] ){
            		j = i;
            		while( --j >= start_run and vec[i] < vec[j]) ++k;
            		while(start_run <= j){	// auf den stack pushen
            			vec_stack.push(start_run++);
            		}
            		while( vec_stack.size() > 0 and vec[i] < vec[vec_stack.top()] ){
            			vec_stack.pop(); ++k;
            		}
            		start_run = i;
            	}
            	bp[k++] = 1;
            }
            */
        } else {
            // hier noch ohne "lazy stack" trick
            for (size_type i=0; i < vec.size(); ++i) {
                while (vec_stack.size() > 0 and vec[i] > vec[vec_stack.top()]) {
                    vec_stack.pop(); ++k; /*bp[k++] = 0; bp is already initialized to zero*/ // writing a closing parenthesis
                }
                vec_stack.push(i);
                bp[k++] = 1; // writing an opening  parenthesis
            }
        }
#ifdef SDSL_DEBUG
        // not necessary as bp is already initialized to zero
        while (!vec_stack.empty()) {
            vec_stack.pop();
            bp[k++] = 0; // writing a closing parenthesis
        }
        assert(k == 2*vec.size());
#endif
    }
}

//! Calculate the balanced parentheses of the Super-Cartesian tree, described in Ohlebusch and Gog (SPIRE 2009).
/*! \param lcp_buf int_vector_buffer of the LCP Array for which the Super-Cartesian tree representation should be calculated.
 *             The value_type of vec should be an unsigned integer type.
 *  \param bp Reference to the balanced parentheses sequence which represents the Super-Cartesian tree.
 *  \param minimum Specifies if the higher levels contains minima or maxima. Default is maxima.
 *  \par Time complexity
 *       \f$ \Order{2n} \f$, where \f$ n=\f$vec.size()
 *  \par Space complexity
 *       \f$\Order{2n}\f$ bits, by the multi_stack_support
 */
template<uint8_t fixedIntWidth>
void construct_supercartesian_tree_bp_succinct(int_vector_buffer<fixedIntWidth>& lcp_buf, bit_vector& bp, const bool minimum=true)
{
    typedef int_vector_size_type size_type;
    size_type n = lcp_buf.size();
    bp.resize(2*n);      // resize bit vector for balanced parentheses to 2 n bits
    if (n == 0)	// if n == 0 we are done
        return;
    util::set_to_value(bp, 0);
    sorted_multi_stack_support vec_stack(n);

    size_type k=0;
    if (minimum) {
        bp[k++] = 1;
        size_type last = lcp_buf[0];
        for (size_type i=1, x; i < n; ++i) {
            x = lcp_buf[i];
            if (x < last) {
                ++k; // writing a closing parenthesis for last
                while (!vec_stack.empty() and x < vec_stack.top()) {
                    vec_stack.pop(); ++k; // writing a closing parenthesis, bp is already initialized to zeros
                }
            } else {
                vec_stack.push(last); // "lazy stack" trick: Beschleunigung: ca 25 %
            }
            bp[k++] = 1; // writing an opening parenthesis
            last = x;
        }
    } else {
        // hier noch ohne "lazy stack" trick
        for (size_type i=0, x; i < n; ++i) {
            x = lcp_buf[i];
            while (!vec_stack.empty() and x > vec_stack.top()) {
                vec_stack.pop(); ++k; // writing a closing parenthesis, bp is already initialized to zeros
            }
            vec_stack.push(x);
            bp[k++] = 1; // writing an opening parenthesis
        }
    }
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
template<uint8_t fixedIntWidth>
int_vector_size_type construct_supercartesian_tree_bp_succinct_and_first_child(int_vector_buffer<fixedIntWidth>& lcp_buf, bit_vector& bp, bit_vector& bp_fc, const bool minimum=true)
{
    typedef int_vector_size_type size_type;
    size_type n = lcp_buf.size();
    bp.resize(2*n);      // resize bit vector for balanaced parantheses to 2 n bits
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
        // hier noch ohne "lazy stack" trick
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
        // hier noch ohne "lazy stack" trick
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
        // writing a closing parenthesis in bp, not necessary as bp is initalized with zeros
        ++k;
        ++k_fc;
    }
//	assert( k == 2*vec.size() );
    return fc_cnt;
}


template<class RandomAccessContainer>
void construct_supercartesian_tree_bp_succinct2(const RandomAccessContainer& vec, bit_vector& bp,
        SDSL_UNUSED const bool minimum=true)
{
    typedef typename RandomAccessContainer::size_type size_type;
    bp.resize(2*vec.size());      // resize bit vector for balanced parentheses to 2 n bits
    util::set_to_value(bp, 0);
    sorted_stack_support vec_stack(vec.size()); // <- ist das ein Problem fuer int_vector_buffer

    size_type k=0;
//	uint64_t wbuf=0;
    for (size_type i=0/*, cnt64=0*/; i < vec.size(); ++i) {
        while (vec_stack.size() > 0 and vec[i] < vec[vec_stack.top()]) {
            vec_stack.pop(); ++k; /*bp[k++] = 0; bp is already initialized to zero*/ // writing a closing parenthesis
        }
        vec_stack.push(i);
        bp[k++] = 1; // writing an opening  parenthesis
        while (i+1 < vec.size() and vec[i+1] >= vec[i]) {
            vec_stack.push(++i);
            bp[k++];
        }
    }
#ifdef SDSL_DEBUG
// not neccessary as bp is already initialized to zero
    while (vec_stack.size() > 0) {
        vec_stack.pop();
        bp[k++] = 0; // writing a closing parenthesis
    }
    assert(k == 2*vec.size());
#endif
}

}

#endif
