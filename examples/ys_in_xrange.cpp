#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_algorithm.hpp>
#include <iostream>
#include <stack>

using namespace sdsl;
using namespace std;

//! Returns if range1 overlaps with range2
/*!
 * \param range1 First range [r1, l1]
 * \param range2 Second range [r2, l2].
 * \return If the two ranges overlap.
 */
template<typename r_t1, typename r_t2>
bool
overlaps(const r_t1& range1, const r_t2& range2)
{
    return std::get<1>(range1) >= std::get<0>(range2)
           and std::get<0>(range1) <= std::get<1>(range2);
}

//! Returns if range1 is contained in range2
/*!
 * \param range1 First range [r1, l1]
 * \param range2 Second range [r2, l2].
 * \return If range1 is contained in range2.
 */
template<typename r_t1, typename r_t2>
bool
is_contained(const r_t1& range1, const r_t2& range2)
{
    return std::get<0>(range1) >= std::get<0>(range2)
           and std::get<1>(range1) <= std::get<1>(range2);
}



//! Count how many points are in [x_0,x_1] x [y_0, y_1]
/*!
 * \param wt      Wavelet tree representing a sequence.
 * \param x_range x-range [x_0,x_1] in wt.
 * \param y_range y-range [y_0,y_1] in wt.
 *
 * \return
 */
template<class t_wt>
typename t_wt::size_type
count(const t_wt& wt, const range_type& x_range, const range_type& y_range)
{
    typedef typename t_wt::size_type size_type;
    typedef typename t_wt::node_type node_type;

    if (x_range.first >= wt.size())
        return 0;
    size_type res = 0;
    stack<std::pair<node_type,range_type>> s;

    auto push_node = [&wt, &s, &y_range, &res](const node_type& v, const range_type& x_range) {
        if (size(x_range) > 0) {
            auto value_range = wt.value_range(v);
            if (is_contained(value_range, y_range)) {
                res += size(x_range);
            } else if (overlaps(y_range, value_range)) {
                s.emplace(v, x_range);
            }
        }
    };

    push_node(wt.root(), x_range);
    while (!s.empty()) {
        auto v       = std::get<0>(s.top());
        auto x_range = std::get<1>(s.top());
        s.pop();
        if (!wt.is_leaf(v)) {
            auto child_v = wt.expand(v);
            auto child_r = wt.expand(v, x_range);
            for (size_t i=0, j=child_v.size()-1; i < child_v.size(); ++i,--j) {
                push_node(child_v[j], child_r[j]);
            }
        }
    }
    return res;
}


int main()
{
    typedef wm_int<> t_wt;
    t_wt wt;
    construct_im(wt, "9 4 3 2 1 4 6 3 1 4 6 5 3 2 1 3 5 3 2 3 4",'d');
    cout << wt << endl;
    auto y_it = ys_in_x_range(wt, 0, wt.size());
    while (y_it) {
        auto y = *y_it;
        cout << get<0>(y) << " ("<< get<1>(y) << "," << get<2>(y) << ")" << endl;
        ++y_it;
    }

    stack<typename t_wt::node_type> s;
    s.push(wt.root());
    while (!s.empty()) {
        auto v = s.top();
        s.pop();
        auto range = wt.value_range(v);
        cout << range[0] << "," << range[1] << endl;
        if (!wt.is_leaf(v)) {
            auto child_v = wt.expand(v);
            for (auto it = child_v.crbegin(); it != child_v.crend(); ++it) {
                s.push(*it);
            }
        }
    }

    cout << "count[0,"<<wt.size()-1<<"][1,16] = " << count(wt, {0, wt.size()-1}, {1,16}) << endl;
    cout << "count[0,"<<wt.size()-1<<"][1,8] = " << count(wt, {0, wt.size()-1}, {1,8}) << endl;
    cout << "count[0,"<<wt.size()-1<<"][2,3] = " << count(wt, {0, wt.size()-1}, {2,3}) << endl;
}

