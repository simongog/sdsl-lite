#ifndef SDSL_TEST_CST_HELPER
#define SDSL_TEST_CST_HELPER

#include <sstream>
#include <iostream>
#include "gtest/gtest.h"
#ifdef WIN32
#include "iso646.h"
#endif

template<class Cst>
std::string format_node(const Cst& cst, const typename Cst::node_type& v)
{
    std::stringstream ss;
    ss << cst.depth(v) << "-["<<cst.lb(v)<<","<<cst.rb(v)<<"]";
    return ss.str();
}

template<class tCst>
void check_node_method(const tCst& cst)
{
    typedef typename tCst::const_iterator const_iterator;
    typedef typename tCst::node_type node_type;
    typedef typename tCst::size_type size_type;
    for (const_iterator it = cst.begin(), end = cst.end(); it != end; ++it) {
        if (it.visit() == 1) {
            node_type v = *it;
            size_type lb = cst.lb(v), rb = cst.rb(v);
            ASSERT_EQ(v, cst.node(lb, rb));
        }
    }
}

template<class Cst>
typename Cst::node_type naive_lca(const Cst& cst, typename Cst::node_type v, typename Cst::node_type w, bool output=false)
{
    typedef typename Cst::size_type size_type;
    size_type steps = 0;
    while (v != w  and steps < cst.csa.size()) {
        if (cst.depth(v) > cst.depth(w)) {
            v = cst.parent(v);
            if (output) {
                std::cout << "v="<<format_node(cst, v) << std::endl;
            }
        } else {
            w = cst.parent(w);
            if (output) {
                std::cout << "w="<<format_node(cst, v) << std::endl;
            }
        }
        steps++;
    }
    return v;
}



#endif
