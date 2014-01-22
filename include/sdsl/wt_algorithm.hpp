#ifndef INCLUDED_SDSL_WT_ALGORITHM
#define INCLUDED_SDSL_WT_ALGORITHM

namespace sdsl
{

//! Intersection of elements in WT[s_0,e_0], WT[s_1,e_1],...,WT[s_k,e_k]
/*! \param wt     The wavelet tree object.
 *  \param ranges The ranges.
 *  \param t      Threshold in how many distinct ranges the value has to be
 *                present. Default: t=ranges.size()
 */
template<class t_wt>
std::vector< std::pair<typename t_wt::value_type, typename t_wt::size_type> >
intersect(const t_wt& wt, const std::vector<typename t_wt::range_type>& ranges, typename t_wt::size_type t=0)
{
    using size_type      = typename t_wt::size_type;
    using value_type     = typename t_wt::value_type;
    using node_type      = typename t_wt::node_type;
    using range_type     = typename t_wt::range_type;
    using range_vec_type = typename t_wt::range_vec_type;
    using pnvr_type      = std::pair<node_type, range_vec_type>;

    using p_t = std::pair<value_type,size_type>;
    std::vector<p_t> res;

    if (ranges.empty())
        return res;

    t = (t==0) ? ranges.size() : t;

    std::stack<pnvr_type> stack;
    size_type n = wt.size();
    stack.emplace(pnvr_type(wt.root(), ranges));

    while (!stack.empty()) {
        pnvr_type x = stack.top(); stack.pop();

        if (wt.is_leaf(x.first)) {
            if (t <= x.second.size()) {
                size_type freq = 0;
                for (auto& r : x.second) freq += (r.second - r.first + 1);
                res.emplace_back(x.first.sym,freq);
            }
        } else {
            auto child        = wt.expand(x.first);
            auto child_ranges = wt.expand(x.first, x.second);

            auto push_node = [&](node_type& child, range_vec_type& child_range) {
                auto end = std::remove_if(child_range.begin(), child_range.end(), [&](const range_type& x) { return wt.empty(x);});
                if (end > child_range.begin() + t - 1) {
                    stack.emplace(pnvr_type(child, range_vec_type(child_range.begin(), end)));
                }
            };
            push_node(std::get<0>(child), std::get<0>(child_ranges));
            push_node(std::get<1>(child), std::get<1>(child_ranges));
        }
    }
    return res;
}
}

#endif
