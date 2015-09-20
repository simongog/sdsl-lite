#ifndef INCLUDED_SDSL_WT_ALGORITHM
#define INCLUDED_SDSL_WT_ALGORITHM

#include <algorithm>
#include <utility>

namespace sdsl
{

template<typename t_wt>
struct has_interval_symbols;

template<typename t_wt, bool t_has_interval_symbols>
struct _interval_symbols_wt;

template<typename, typename T>
struct has_expand;

//! Intersection of elements in WT[s_0,e_0], WT[s_1,e_1],...,WT[s_k,e_k]
/*! \param wt     The wavelet tree object.
 *  \param ranges The ranges.
 *  \param t      Threshold in how many distinct ranges the value has to be
 *                present. Default: t=ranges.size()
 *  \return       A vector containing (value, frequency) - of value which are
 *                contained in t different ranges. Frequency = accumulated
 *                frequencies in all ranges. The tuples are ordered according
 *                to value, if t_wt::lex_ordered=1.
 */
template<class t_wt>
std::vector< std::pair<typename t_wt::value_type, typename t_wt::size_type> >
intersect(const t_wt& wt, const std::vector<range_type>& ranges, typename t_wt::size_type t=0)
{
    using std::get;
    using size_type      = typename t_wt::size_type;
    using value_type     = typename t_wt::value_type;
    using node_type      = typename t_wt::node_type;
    using pnvr_type      = std::pair<node_type, range_vec_type>;
    typedef std::stack<pnvr_type> stack_type;

    static_assert(has_expand<t_wt, std::pair<node_type,node_type>(const node_type&)>::value,
                  "intersect requires t_wt to have expand(const node_type&)");

    using p_t = std::pair<value_type,size_type>;
    std::vector<p_t> res;

    auto push_node = [&wt,&t](stack_type& s, node_type& child,
    range_vec_type& child_range) {
        auto end = std::remove_if(child_range.begin(), child_range.end(),
        [&](const range_type& x) { return empty(x);});
        if (end > child_range.begin() + t - 1) {
            s.emplace(pnvr_type(child, range_vec_type(child_range.begin(),
                                end)));
        }
    };

    if (ranges.empty())
        return res;

    t = (t==0) ? ranges.size() : t;

    std::stack<pnvr_type> stack;
    stack.emplace(pnvr_type(wt.root(), ranges));

    while (!stack.empty()) {
        pnvr_type x = stack.top(); stack.pop();

        if (wt.is_leaf(x.first)) {
            const auto& iv = x.second;
            if (t <= iv.size()) {
                auto freq = std::accumulate(iv.begin(), iv.end(), 0ULL,
                [](size_type acc, const range_type& r) {
                    return acc+(r.second-r.first+1);
                });
                res.emplace_back(wt.sym(x.first),freq);
            }
        } else {
            auto child        = wt.expand(x.first);
            auto child_ranges = wt.expand(x.first, x.second);

            push_node(stack, get<1>(child), get<1>(child_ranges));
            push_node(stack, get<0>(child), get<0>(child_ranges));
        }
    }
    return res;
}


//! Returns the q-th smallest element and its frequency in wt[lb..rb].
/*! \param wt The wavelet tree.
 *  \param lb Left array bound in T
 *  \param rb Right array bound in T
 *  \param q q-th largest element ('quantile'), 0-based indexed.
 */
template<class t_wt>
std::pair<typename t_wt::value_type, typename t_wt::size_type>
quantile_freq(const t_wt& wt, typename t_wt::size_type lb,
              typename t_wt::size_type rb, typename t_wt::size_type q)
{
    static_assert(t_wt::lex_ordered,
                  "quantile_freq requires a lex_ordered WT");
    using std::get;
    using node_type      = typename t_wt::node_type;
    static_assert(has_expand<t_wt, std::pair<node_type,node_type>(const node_type&)>::value,
                  "quantile_freq requires t_wt to have expand(const node_type&)");

    node_type v = wt.root();
    range_type r(lb,rb);

    while (!wt.is_leaf(v)) {
        auto child        = wt.expand(v);
        auto child_ranges = wt.expand(v, r);
        auto num_zeros    = size(get<0>(child_ranges));

        if (q >= num_zeros) {
            q -= num_zeros;
            v = get<1>(child);
            r = get<1>(child_ranges);
        } else {
            v = get<0>(child);
            r = get<0>(child_ranges);
        }
    }
    return {wt.sym(v), size(r)};
};


template<class t_wt>
void
_interval_symbols_rec(const t_wt& wt, range_type r,
                      typename t_wt::size_type& k,
                      std::vector<typename t_wt::value_type>& cs,
                      std::vector<typename t_wt::size_type>& rank_c_i,
                      std::vector<typename t_wt::size_type>& rank_c_j,
                      const typename t_wt::node_type& v)
{
    using std::get;
    if (wt.is_leaf(v)) {
        rank_c_i[k] = r.first;
        rank_c_j[k] = r.second+1;
        cs[k++] = wt.sym(v);
    } else {
        auto child        = wt.expand(v);
        auto child_ranges = wt.expand(v, r);
        if (!empty(get<0>(child_ranges))) {
            _interval_symbols_rec(wt, get<0>(child_ranges), k, cs, rank_c_i,
                                  rank_c_j, get<0>(child));
        }
        if (!empty(get<1>(child_ranges))) {
            _interval_symbols_rec(wt, get<1>(child_ranges), k, cs, rank_c_i,
                                  rank_c_j, get<1>(child));
        }
    }
}

template<class t_wt>
void
_interval_symbols(const t_wt& wt, typename t_wt::size_type i,
                  typename t_wt::size_type j,
                  typename t_wt::size_type& k,
                  std::vector<typename t_wt::value_type>& cs,
                  std::vector<typename t_wt::size_type>& rank_c_i,
                  std::vector<typename t_wt::size_type>& rank_c_j)
{

    assert(i <= j and j <= wt.size());
    k=0;
    if ((i+1)==j) {
        auto res = wt.inverse_select(i);
        cs[0]=res.second;
        rank_c_i[0]=res.first;
        rank_c_j[0]=res.first+1;
        k=1;
        return;
    } else if (j>i) {
        _interval_symbols_rec(wt, range_type(i,j-1), k, cs,
                              rank_c_i, rank_c_j, wt.root());
    }
}

//! For each symbol c in wt[i..j-1] get rank(i,c) and rank(j,c).
/*!
 * \param i        The start index (inclusive) of the interval.
 * \param j        The end index (exclusive) of the interval.
 * \param k        Reference for number of different symbols in [i..j-1].
 * \param cs       Reference to a vector that will contain in
 *                 cs[0..k-1] all symbols that occur in [i..j-1] in
 *                 ascending order.
 * \param rank_c_i Reference to a vector which equals
 *                 rank_c_i[p] = rank(i,cs[p]), for \f$ 0 \leq p < k \f$.
 * \param rank_c_j Reference to a vector which equals
 *                 rank_c_j[p] = rank(j,cs[p]), for \f$ 0 \leq p < k \f$.
 * \par Time complexity
 *      \f$ \Order{\min{\sigma, k \log \sigma}} \f$
 *
 * \par Precondition
 *      \f$ i \leq j \leq size() \f$
 *      \f$ cs.size() \geq \sigma \f$
 *      \f$ rank_{c_i}.size() \geq \sigma \f$
 *      \f$ rank_{c_j}.size() \geq \sigma \f$
 */
template<class t_wt>
void
interval_symbols(const t_wt& wt, typename t_wt::size_type i,
                 typename t_wt::size_type j,
                 typename t_wt::size_type& k,
                 std::vector<typename t_wt::value_type>& cs,
                 std::vector<typename t_wt::size_type>& rank_c_i,
                 std::vector<typename t_wt::size_type>& rank_c_j)
{
    // check if wt has a built-in interval_symbols method
    constexpr bool has_own = has_interval_symbols<t_wt>::value;
    if (has_own) {  // if yes, call it
        _interval_symbols_wt<t_wt, has_own>::call(wt, i, j, k,
                cs, rank_c_i, rank_c_j);
    } else { // otherwise use generic implementation based on expand
        _interval_symbols(wt, i,j, k, cs, rank_c_i, rank_c_j);
    }
}



// has_interval_symbols<X>::value is true if class X has
// implement method interval_symbols
// Adapted solution from jrok's proposal:
// http://stackoverflow.com/questions/87372/check-if-a-class-has-a-member-function-of-a-given-signature
template<typename t_wt>
struct has_interval_symbols {
    template<typename T>
    static constexpr auto check(T*)
    -> typename
    std::is_same<
    decltype(std::declval<T>().interval_symbols(
                 std::declval<typename T::size_type>(),
                 std::declval<typename T::size_type>(),
                 std::declval<typename T::size_type&>(),
                 std::declval<std::vector<typename T::value_type>&>(),
                 std::declval<std::vector<typename T::size_type>&>(),
                 std::declval<std::vector<typename T::size_type>&>()
             )),
             void>::type {return std::true_type();}
             template<typename>
    static constexpr std::false_type check(...) {return std::false_type();}
    typedef decltype(check<t_wt>(nullptr)) type;
    static constexpr bool value = type::value;
};

template<typename t_wt, bool t_has_interval_symbols>
struct _interval_symbols_wt {
    typedef typename t_wt::size_type  size_type;
    typedef typename t_wt::value_type value_type;

    static void call(const t_wt& wt, size_type i, size_type j, size_type& k,
                     std::vector<value_type>& cs, std::vector<size_type>& rank_c_i,
                     std::vector<size_type>& rank_c_j)
    {
        wt.interval_symbols(i,j,k,cs,rank_c_i,rank_c_j);
    }
};


template<typename t_wt>
struct _interval_symbols_wt<t_wt, false> {
    typedef typename t_wt::size_type  size_type;
    typedef typename t_wt::value_type value_type;

    static void call(const t_wt&, size_type, size_type, size_type&,
                     std::vector<value_type>&, std::vector<size_type>&,
                     std::vector<size_type>&)
    {
    }
};

template<typename, typename T>
struct has_expand {
    static_assert(std::integral_constant<T, false>::value,
                  "Second template parameter needs to be of function type.");
};

template<typename t_wt, typename t_ret, typename... t_args>
struct has_expand<t_wt, t_ret(t_args...)> {
    template<typename T>
    static constexpr auto check(T*)
    -> typename
    std::is_same<
    decltype(std::declval<T>().expand(std::declval<t_args>()...)),
    t_ret>::type { return std::true_type();}
    template<typename>
static constexpr std::false_type check(...) { return std::false_type();}
typedef decltype(check<t_wt>(nullptr)) type;
static constexpr bool value = type::value;
            };

template<typename t_wt>
struct has_range_search_2d {
    template<typename T>
    static constexpr auto check(T*)
    -> typename
    std::is_same<
    decltype(std::declval<T>().range_search_2d(//
                 std::declval<typename T::size_type>(),
                 std::declval<typename T::size_type>(),
                 std::declval<typename T::value_type>(),
                 std::declval<typename T::value_type>(),
                 false
             )),
             std::pair<typename T::size_type,
             std::vector<std::pair<typename T::value_type,
             typename T::size_type>>>>::type {return std::true_type();}

             template<typename>
    static constexpr std::false_type check(...) {return std::false_type();}
    typedef decltype(check<t_wt>(nullptr)) type;
    static constexpr bool value = type::value;
};


//! Returns for a symbol c the previous smaller or equal symbol in the WT.
/*! \param c the symbol
 *  \return A pair. The first element of the pair consititues if
 *          a valid answer was found (true) or no valid answer (false)
 *          could be found. The second element contains the found symbol.
 */
template<class t_wt>
std::pair<bool,typename t_wt::value_type>
_symbol_lte(const t_wt& wt,typename t_wt::value_type c)
{
    if (((1ULL) << (wt.max_level)) <= c) {
        // c is greater than any symbol in wt. return the largest symbol!
        c = sdsl::bits::lo_set[wt.max_level];
    }
    auto node = wt.root();
    auto predecessor_subtree = node;
    uint64_t mask = (1ULL) << (wt.max_level - 1);
    while (!wt.is_leaf(node)) {
        auto children = wt.expand(node);
        auto left_child = std::get<0>(children);
        auto right_child = std::get<1>(children);
        if (c & (mask >> node.level)) { // go right
            if (right_child.size) {
                node = right_child;
                if (left_child.size) { // potential predecessor subtree?
                    predecessor_subtree = left_child;
                }
            } else { // dead end
                // left child can't be empty if left child is
                node = left_child;
                c = sdsl::bits::all_set;
            }
        } else { // go left
            if (left_child.size) {
                node = left_child;
            } else { // dead end
                if (predecessor_subtree == wt.root()) {
                    // there is no valid predecessor. symbol must be
                    // smaller than the smallest symbol in the wt.
                    return {false, 0};
                }
                node = predecessor_subtree;
                c = sdsl::bits::all_set;
            }
        }
    }
    return {true, node.sym};
}


//! Returns for a symbol c the next larger or equal symbol in the WT.
/*! \param c the symbol
 *  \return A pair. The first element of the pair consititues if
 *          a valid answer was found (true) or no valid answer (false)
 *          could be found. The second element contains the found symbol.
 */
template<class t_wt>
std::pair<bool,typename t_wt::value_type>
_symbol_gte(const t_wt& wt,typename t_wt::value_type c)
{
    if (((1ULL) << (wt.max_level)) <= c) {
        // c is greater than any symbol in wt
        return {false, 0};
    }
    auto node = wt.root();
    auto successor_subtree = node;
    uint64_t mask = (1ULL) << (wt.max_level - 1);
    while (!wt.is_leaf(node)) {
        auto children = wt.expand(node);
        auto left_child = std::get<0>(children);
        auto right_child = std::get<1>(children);
        if (c & (mask >> node.level)) { // go right
            if (right_child.size) {
                node = right_child;
            } else { // dead end
                if (successor_subtree == wt.root()) {
                    // there is no valid successor. symbol must be
                    // bigger than the largest symbol in the wt.
                    return {false, 0};
                }
                node = successor_subtree;
                c = 0;
            }
        } else { // go left
            if (left_child.size) {
                node = left_child;
                if (right_child.size) { // potential successor subtree?
                    successor_subtree = right_child;
                }
            } else { // dead end
                // right child can't be empty if left child is
                node = right_child;
                c = 0;
            }
        }
    }
    return {true, node.sym};
}

template<class t_wt, bool t_has_interval_symbols>
struct _symbols_calls_wt {
    typedef typename t_wt::value_type value_type;

    static std::pair<bool, value_type>
    call_symbol_gte(const t_wt& wt,value_type c)
    {
        return wt.symbol_gte(c);
    }

    static std::pair<bool,value_type>
    call_symbol_lte(const t_wt& wt,value_type c)
    {
        return wt.symbol_lte(c);
    }
};


template<class t_wt>
struct _symbols_calls_wt<t_wt, false> {
    typedef typename t_wt::value_type value_type;

    static std::pair<bool,value_type>
    call_symbol_gte(const t_wt& wt,value_type c)
    {
        return _symbol_gte(wt,c);
    }

    static std::pair<bool,value_type>
    call_symbol_lte(const t_wt& wt,value_type c)
    {
        return _symbol_lte(wt,c);
    }
};

template<typename t_wt>
struct has_symbols_wt {
    template<typename T>
    static constexpr auto check(T*)
    -> typename
    std::is_same<
    decltype(std::declval<T>().symbol_gte(std::declval<typename T::value_type>())),
             std::pair<bool,typename T::value_type>
             >::type {return std::true_type();}

             template<typename>
    static constexpr std::false_type check(...) {return std::false_type();}
    typedef decltype(check<t_wt>(nullptr)) type;
    static constexpr bool value = type::value;
};

//! Returns for a symbol c the previous smaller or equal symbol in the WT.
/*! \param c the symbol
 *  \return A pair. The first element of the pair consititues if
 *          a valid answer was found (true) or no valid answer (false)
 *          could be found. The second element contains the found symbol.
 */
template<class t_wt>
std::pair<bool,typename t_wt::value_type>
symbol_lte(const t_wt& wt, typename t_wt::value_type c)
{
    static_assert(t_wt::lex_ordered, "symbols_lte requires a lex_ordered WT");
    // check if wt has a built-in interval_symbols method
    constexpr bool has_own = has_symbols_wt<t_wt>::value;
    return _symbols_calls_wt<t_wt, has_own>::call_symbol_lte(wt,c);
}

//! Returns for a symbol c the next larger or equal symbol in the WT.
/*! \param c the symbol
 *  \return A pair. The first element of the pair consititues if
 *          a valid answer was found (true) or no valid answer (false)
 *          could be found. The second element contains the found symbol.
 */
template<class t_wt>
std::pair<bool,typename t_wt::value_type>
symbol_gte(const t_wt& wt, typename t_wt::value_type c)
{
    static_assert(t_wt::lex_ordered, "symbols_gte requires a lex_ordered WT");
    // check if wt has a built-in interval_symbols method
    constexpr bool has_own = has_symbols_wt<t_wt>::value;
    return _symbols_calls_wt<t_wt, has_own>::call_symbol_gte(wt,c);
}

//! Returns for a x range [x_i,x_j] and a value range [y_i,y_j] all unique y
//! values occuring in [x_i,x_j] in ascending order.
/*! \param x_i lower bound of the x range
 *  \param x_j upper bound of the x range
 *  \param y_i lower bound of the y range
 *  \param y_j upper bound of the y range
 *  \return a vector of increasing y values occuring in the range [x_i,x_j]
 */
template <class t_wt>
std::vector<typename t_wt::value_type>
restricted_unique_range_values(const t_wt& wt,
                               typename t_wt::size_type x_i,
                               typename t_wt::size_type x_j,
                               typename t_wt::value_type y_i,
                               typename t_wt::value_type y_j)
{
    static_assert(t_wt::lex_ordered, "restricted_unique_range_values requires a lex_ordered WT");

    std::vector<typename t_wt::value_type> unique_values;

    // make sure things are within bounds
    if (x_j > wt.size()-1) x_j = wt.size()-1;
    if ((x_i > x_j) || (y_i > y_j)) {
        return unique_values;
    }
    auto lower_y_bound = symbol_gte(wt,y_i);
    auto upper_y_bound = symbol_lte(wt,y_j);
    // is the y range valid?
    if (!lower_y_bound.first || !upper_y_bound.first
        || (lower_y_bound.second > upper_y_bound.second))  {
        return unique_values;
    }

    auto lower_y_bound_path = wt.path(lower_y_bound.second);
    auto upper_y_bound_path = wt.path(upper_y_bound.second);

    auto compare_path = [](uint64_t node_path,uint64_t node_path_len,
    std::pair<uint64_t,uint64_t> bound_path) -> int {
        auto bound_path_len = bound_path.first;
        auto bound_path_val = bound_path.second;
        /* align to same length */
        if (bound_path_len > node_path_len)
            bound_path_val = bound_path_val >> (bound_path_len-node_path_len);
        if (bound_path_len < node_path_len)
            bound_path_val = bound_path_val << (node_path_len-bound_path_len);
        /* cmp */
        if (node_path < bound_path_val) return -1;
        if (node_path > bound_path_val) return 1;
        return 0;
    };

    std::stack<std::tuple<typename t_wt::node_type,sdsl::range_type,uint64_t,uint64_t>> stack;
    sdsl::range_type initial_range = {x_i,x_j};
    stack.emplace(wt.root(),initial_range,0,0);
    while (!stack.empty()) {
        auto node_data = stack.top(); stack.pop();
        auto node = std::get<0>(node_data);
        auto range = std::get<1>(node_data);
        auto node_path = std::get<2>(node_data);
        auto node_level = std::get<3>(node_data);
        if (wt.is_leaf(node)) {
            unique_values.emplace_back(wt.sym(node));
        } else {
            auto children = wt.expand(node);
            auto left_path = node_path<<1ULL;
            auto right_path = (node_path<<1ULL)|1ULL;
            auto child_ranges = wt.expand(node,range);
            if (compare_path(right_path,node_level+1,upper_y_bound_path) < 1) {
                auto right_child = std::get<1>(children);
                auto right_range = std::get<1>(child_ranges);
                if (!sdsl::empty(right_range))
                    stack.emplace(right_child,right_range,right_path,node_level+1);
            }
            if (compare_path(left_path,node_level+1,lower_y_bound_path) > -1) {
                auto left_child = std::get<0>(children);
                auto left_range = std::get<0>(child_ranges);
                if (!sdsl::empty(left_range))
                    stack.emplace(left_child,left_range,left_path,node_level+1);
            }
        }
    }

    return unique_values;
}



// Check for node_type of wavelet_tree
// http://stackoverflow.com/questions/7834226/detecting-typedef-at-compile-time-template-metaprogramming
// + trick to make it work for VC++
// https://connect.microsoft.com/VisualStudio/feedback/details/790269/compile-error-with-explicitly-specified-template-arguments
template<typename T>
struct void_ { typedef void type; };

template<typename t_wt, typename T = void>
struct has_node_type {
    typedef std::false_type t_expr;
    enum { value = t_expr::value };
};

template<typename t_wt>
struct has_node_type<t_wt, typename void_<typename t_wt::node_type>::type> {
    typedef std::true_type t_expr;
    enum { value = t_expr::value };
};

} // end namespace

#endif
