#ifndef INCLUDED_SDSL_SUFFIX_TREE_HELPER
#define INCLUDED_SDSL_SUFFIX_TREE_HELPER

#include <stdint.h>
#include <cstdlib>
#include <cassert>
#include "iterators.hpp"

namespace sdsl
{


template <class t_cst>
class cst_node_child_proxy
{
    public: // types
        using iterator_type = random_access_const_iterator<cst_node_child_proxy>;
        using node_type = typename t_cst::node_type;
        using size_type = typename t_cst::size_type;
        using difference_type = size_type;
        using value_type = node_type;
    private: // data
        const node_type& m_parent;
        const t_cst& m_cst;
    public: // constructors
        cst_node_child_proxy() = delete;
        explicit cst_node_child_proxy(const t_cst& cst,const node_type& v) : m_parent(v) , m_cst(cst) {};
    public: // methods
        value_type operator[](size_type i) const { return m_cst.select_child(m_parent,i+1); } // enumeration starts with 1 not 0
        iterator_type begin() const { return iterator_type(this); }
        iterator_type end() const { return iterator_type(this,m_cst.degree(m_parent)); }
};


}

#endif
