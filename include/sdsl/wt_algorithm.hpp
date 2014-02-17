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
    decltype(std::declval<T>().interval_symbols(//
                 std::declval<typename T::size_type>(),
                 std::declval<typename T::size_type>(),
                 std::declval<typename T::size_type&>(),
                 std::declval<std::vector<typename T::value_type>&>(),
                 std::declval<std::vector<typename T::size_type>&>(),
                 std::declval<std::vector<typename T::size_type>&>()
             )),
             void>::type;
    template<typename>
    static constexpr std::false_type check(...);
    typedef decltype(check<t_wt>(nullptr)) type;
    static constexpr bool value = type::value;
};

template<typename t_wt, bool t_has_interval_symbols>
struct _interval_symbols_wt {
    typedef typename t_wt::size_type  size_type;
    typedef typename t_wt::value_type value_type;

    static void call(const t_wt& wt, size_type i, size_type j, size_type& k,
                     std::vector<value_type>& cs, std::vector<size_type>& rank_c_i,
                     std::vector<size_type>& rank_c_j) {
        wt.interval_symbols(i,j,k,cs,rank_c_i,rank_c_j);
    }
};


template<typename t_wt>
struct _interval_symbols_wt<t_wt, false> {
    typedef typename t_wt::size_type  size_type;
    typedef typename t_wt::value_type value_type;

    static void call(const t_wt&, size_type, size_type, size_type&,
                     std::vector<value_type>&, std::vector<size_type>&,
                     std::vector<size_type>&) {
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
    t_ret>::type;
    template<typename>
    static constexpr std::false_type check(...);
    typedef decltype(check<t_wt>(nullptr)) type;
    static constexpr bool value = type::value;
};

} // end namespace

#endif
