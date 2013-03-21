#include "sdsl/lcp_support_tree.hpp"

namespace sdsl
{

void construct_first_child_lcp(int_vector_file_buffer<>& lcp_buf, int_vector<>& fc_lcp)
{
    typedef int_vector_size_type size_type;
    lcp_buf.reset();
    size_type n = lcp_buf.int_vector_size;
    if (n == 0) {	// if n == 0 we are done
        fc_lcp = int_vector<>(0);
    }
    {
        int_vector<> tmp(n, 0, bits::l1BP(n)+1);
        fc_lcp.swap(tmp);
    }

    size_type fc_cnt=0; // first child counter
    sorted_multi_stack_support vec_stack(n);
    size_type y;
    for (size_type i=0, r_sum = 0, r = lcp_buf.load_next_block(), x; r_sum < n;) {
        for (; i < r_sum +r; ++i) {
            x = lcp_buf[i-r_sum];
            while (!vec_stack.empty() and x < vec_stack.top()) {
                y = vec_stack.top();
                if (vec_stack.pop()) {
                    fc_lcp[fc_cnt++] = y;
                }
            }
            vec_stack.push(x);
        }
        r_sum += r; r = lcp_buf.load_next_block();
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
