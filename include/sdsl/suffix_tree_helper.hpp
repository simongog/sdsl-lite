#ifndef INCLUDED_SDSL_SUFFIX_TREE_HELPER
#define INCLUDED_SDSL_SUFFIX_TREE_HELPER

#include <stdint.h>
#include <cstdlib>
#include <cassert>
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

}

#endif
