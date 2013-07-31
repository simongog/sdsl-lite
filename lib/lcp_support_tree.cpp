#include "sdsl/lcp_support_tree.hpp"

namespace sdsl
{

void construct_first_child_lcp(int_vector_buffer<>& lcp_buf, int_vector<>& fc_lcp)
{
    typedef int_vector_size_type size_type;
    size_type n = lcp_buf.size();
    if (n == 0) {	// if n == 0 we are done
        fc_lcp = int_vector<>(0);
    }
    {
        int_vector<> tmp(n, 0, bits::hi(n)+1);
        fc_lcp.swap(tmp);
    }

    size_type fc_cnt=0; // first child counter
    sorted_multi_stack_support vec_stack(n);
    size_type y;
    for (size_type i=0, x; i < n; ++i) {
        x = lcp_buf[i];
        while (!vec_stack.empty() and x < vec_stack.top()) {
            y = vec_stack.top();
            if (vec_stack.pop()) {
                fc_lcp[fc_cnt++] = y;
            }
        }
        vec_stack.push(x);
    }

    while (!vec_stack.empty()) {
        y = vec_stack.top();
        if (vec_stack.pop()) {
            fc_lcp[fc_cnt++] = y;
        }
    }
    if (fc_cnt < fc_lcp.size()) {
        fc_lcp.resize(fc_cnt);
    }
}

}
