#ifndef INCLUDED_SDSL_WT_ALGORITHM
#define INCLUDED_SDSL_WT_ALGORITHM

#include <algorithm>
#include <stack>
#include <utility>

namespace sdsl
{

template<typename, typename T>
struct has_expand;

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
    std::stack<std::pair<node_type,range_type>> s;

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


template<typename t_wt>
class map_to_sorted_iterator
{
        static_assert(t_wt::traversable, "map_to_sorted_sequence requires t_wt to be traversable.");
        static_assert(t_wt::lex_ordered, "map_to_sorted_sequence requires t_wt to be lex_ordered.");
    public:
        typedef void(*t_mfptr)();
        typedef typename t_wt::value_type value_type;
        typedef typename t_wt::size_type  size_type;
        typedef std::tuple<value_type, range_type> ret_type;
    private:
        typedef typename t_wt::node_type node_type;
        typedef std::tuple<node_type,range_type,size_type> state_type;

        const t_wt*             m_wt = nullptr;
        range_type              m_y_range;
        std::stack<state_type>  m_stack;
        ret_type                m_ret;
        bool                    m_valid = false;

        void cond_push(const node_type& v, const range_type& x_range, size_type lex_sml)
        {
            if (size(x_range) > 0) {
                auto value_range = m_wt->value_range(v);
                if (overlaps(m_y_range, value_range)) {
                    m_stack.emplace(v, x_range, lex_sml);
                }
            }
        }

    public:
        map_to_sorted_iterator() = default;
        map_to_sorted_iterator(const map_to_sorted_iterator&) = default;
        map_to_sorted_iterator(map_to_sorted_iterator&&) = default;
        map_to_sorted_iterator& operator=(const map_to_sorted_iterator&) = default;
        map_to_sorted_iterator& operator=(map_to_sorted_iterator&&) = default;

        map_to_sorted_iterator(const t_wt* wt, const range_type& x_range,
                               const range_type& y_range) : m_wt(wt),
            m_y_range(y_range)
        {
            if (wt!=nullptr and size(x_range) > 0 and std::get<0>(x_range) < wt->size()) {
                if (overlaps(y_range, m_wt->value_range(wt->root()))) {
                    m_stack.emplace(m_wt->root(), x_range, 0);
                    ++(*this);
                }
            }
        }

        map_to_sorted_iterator& operator++()
        {
            m_valid = false;
            while (!m_stack.empty()) {
                auto v       = std::get<0>(m_stack.top());
                auto x_range = std::get<1>(m_stack.top());
                auto lex_sml = std::get<2>(m_stack.top());
                m_stack.pop();
                if (!m_wt->is_leaf(v)) {
                    auto child_v = m_wt->expand(v);
                    auto child_r = m_wt->expand(v, x_range);
                    for (int i=1; i >= 0; --i) {
                        if (size(child_r[i]) > 0) {
                            if (i==1)
                                cond_push(child_v[i], child_r[i], lex_sml+m_wt->size(child_v[0]));
                            else
                                cond_push(child_v[i], child_r[i], lex_sml);
                        }
                    }
                } else {
                    m_ret = ret_type(m_wt->sym(v), range_type(std::get<0>(x_range)+lex_sml,std::get<1>(x_range)+lex_sml));
                    m_valid = true;
                    break;
                }
            }
            return *this;
        }

        //! Postfix increment of the iterator
        map_to_sorted_iterator operator++(int)
        {
            map_to_sorted_iterator it = *this;
            ++(*this);
            return it;
        }

        ret_type operator*() const
        {
            return m_ret;
        }

        operator t_mfptr() const
        {
            return (t_mfptr)(m_valid);
        }
};

template<typename t_wt>
map_to_sorted_iterator<t_wt>
map_to_sorted_sequence(const t_wt& wt, const range_type& x_range, const range_type& y_range)
{
    static_assert(t_wt::traversable, "map_to_sorted_sequence requires t_wt to be traversable.");
    static_assert(t_wt::lex_ordered, "map_to_sorted_sequence requires t_wt to be lex_ordered.");
    return map_to_sorted_iterator<t_wt>(&wt, x_range, y_range);
}



template<typename t_wt>
class y_iterator
{
        static_assert(t_wt::traversable, "y_iterator requires t_wt to be traversable.");
    public:
        typedef void(*t_mfptr)();
        typedef typename t_wt::value_type value_type;
        typedef typename t_wt::size_type  size_type;
        typedef std::tuple<value_type, size_type, size_type> t_ret;

    private:
        typedef typename t_wt::node_type node_type;
        typedef std::pair<node_type, range_type> t_state;

        const t_wt*         m_wt = nullptr;
        std::stack<t_state> m_stack;
        t_ret               m_ret;
        bool                m_valid = false;

    public:
        y_iterator() = default;
        y_iterator(const y_iterator&) = default;
        y_iterator(y_iterator&&) = default;
        y_iterator& operator=(const y_iterator&) = default;
        y_iterator& operator=(y_iterator&&) = default;
        // wt wavelet tree
        // lb inclusive
        // rb exclusive
        y_iterator(const t_wt& wt, size_type lb, size_type rb) :
            m_wt(&wt), m_valid(wt.size()>0)
        {
            if (m_wt->size() > 0) {
                if ((lb+1) == rb) {
                    auto res = m_wt->inverse_select(lb);
                    m_ret = t_ret(res.second, res.first, res.first+1);
                } else if (rb > lb) {
                    m_stack.emplace(wt.root(), range_type(lb, rb-1));
                    ++(*this);
                }
            }
        }

        //! Prefix increment of the iterator
        y_iterator& operator++()
        {
            m_valid = false;
            while (!m_stack.empty()) {
                auto v = std::get<0>(m_stack.top());
                auto r = std::get<1>(m_stack.top());
                m_stack.pop();
                if (m_wt->is_leaf(v)) {
                    m_ret = t_ret(m_wt->sym(v), r.first, r.second+1);
                    m_valid = true;
                    break;
                } else {
                    auto child_v = m_wt->expand(v);
                    auto child_r = m_wt->expand(v, r);
                    if (!sdsl::empty(std::get<1>(child_r))) {
                        m_stack.emplace(std::get<1>(child_v), std::get<1>(child_r));
                    }
                    if (!sdsl::empty(std::get<0>(child_r))) {
                        m_stack.emplace(std::get<0>(child_v), std::get<0>(child_r));
                    }
                }
            }
            return *this;
        }

        //! Postfix increment of the iterator
        y_iterator operator++(int)
        {
            y_iterator it = *this;
            ++(*this);
            return it;
        }

        t_ret operator*() const
        {
            return m_ret;
        }

        operator t_mfptr() const
        {
            return (t_mfptr)(m_valid);
        }
};

//! Returns an iterator over all values in wt[i..j-1]
/*!
 * \param wt The wavelet tree.
 * \param i  The start index (inclusive) of the interval.
 * \param j  The end index (exclusive) of the interval.
 * \return   Iterator to the result. The iterator points to
 *           triples (y-value, sp, ep), where [sp,ep) is the
 *           range in the leaf. I.e. ep-sp is the number of
 *           ys in the x-range [i..j).
 *
 * \par Time complexity
 *      Iterating over all k values in wt[i..j] takes \f$\Order{k\log\sigma} time.
 *
 * \par Precondition
 *      \f$ i \leq j \leq size() \f$
 */
template<typename t_wt>
y_iterator<t_wt>
ys_in_x_range(const t_wt& wt, typename t_wt::size_type i,
              typename t_wt::size_type j)
{
    static_assert(t_wt::traversable, "ys_in_x_range requires t_wt to be traversable.");
    return y_iterator<t_wt>(wt, i, j);
}


template<typename t_wt>
class yoff_iterator
{
        static_assert(t_wt::traversable, "yoff_iterator requires t_wt to be traversable.");
        static_assert(t_wt::lex_ordered, "yoff_iterator requires t_wt to be lex_ordered.");
    public:
        typedef void(*t_mfptr)();
        typedef typename t_wt::value_type value_type;
        typedef typename t_wt::size_type  size_type;
        typedef std::tuple<value_type, size_type, size_type, size_type> t_ret;

    private:
        typedef typename t_wt::node_type node_type;
        typedef std::tuple<node_type, range_type, size_type> t_state;

        const t_wt*         m_wt = nullptr;
        std::stack<t_state> m_stack;
        t_ret               m_ret;
        bool                m_valid = false;

    public:
        yoff_iterator() = default;
        yoff_iterator(const yoff_iterator&) = default;
        yoff_iterator(yoff_iterator&&) = default;
        yoff_iterator& operator=(const yoff_iterator&) = default;
        yoff_iterator& operator=(yoff_iterator&&) = default;
        // wt wavelet tree
        // lb inclusive
        // rb exclusive
        yoff_iterator(const t_wt& wt, size_type lb, size_type rb) :
            m_wt(&wt), m_valid(wt.size()>0)
        {
            if (m_wt->size() > 0) {
                if (rb > lb) {
                    m_stack.emplace(wt.root(), range_type(lb, rb-1), 0);
                    ++(*this);
                }
            }
        }

        //! Prefix increment of the iterator
        yoff_iterator& operator++()
        {
            m_valid = false;
            while (!m_stack.empty()) {
                auto v = std::get<0>(m_stack.top());
                auto r = std::get<1>(m_stack.top());
                auto lex_smaller = std::get<2>(m_stack.top());
                m_stack.pop();
                if (m_wt->is_leaf(v)) {
                    m_ret = t_ret(m_wt->sym(v), r.first, r.second+1, lex_smaller);
                    m_valid = true;
                    break;
                } else {
                    auto child_v = m_wt->expand(v);
                    auto child_r = m_wt->expand(v, r);
                    if (!sdsl::empty(std::get<1>(child_r))) {
                        m_stack.emplace(std::get<1>(child_v), std::get<1>(child_r),
                                        lex_smaller + m_wt->size(std::get<0>(child_v)));
                    }
                    if (!sdsl::empty(std::get<0>(child_r))) {
                        m_stack.emplace(std::get<0>(child_v), std::get<0>(child_r), lex_smaller);
                    }
                }
            }
            return *this;
        }

        //! Postfix increment of the iterator
        yoff_iterator operator++(int)
        {
            yoff_iterator it = *this;
            ++(*this);
            return it;
        }

        t_ret operator*() const
        {
            return m_ret;
        }

        operator t_mfptr() const
        {
            return (t_mfptr)(m_valid);
        }
};

//! Returns an iterator over all values in wt[i..j-1]
/*!
 * \param wt The wavelet tree.
 * \param i  The start index (inclusive) of the interval.
 * \param j  The end index (exclusive) of the interval.
 * \return   Iterator to the result. The iterator points to
 *           triples (y-value, sp, ep), where [sp,ep) is the
 *           range in the leaf. I.e. ep-sp is the number of
 *           ys in the x-range [i..j).
 *
 * \par Time complexity
 *      Iterating over all k values in wt[i..j] takes \f$\Order{k\log\sigma} time.
 *
 * \par Precondition
 *      \f$ i \leq j \leq size() \f$
 */
template<typename t_wt>
yoff_iterator<t_wt>
ys_and_off_in_x_range(const t_wt& wt, typename t_wt::size_type i,
                      typename t_wt::size_type j)
{
    static_assert(t_wt::traversable, "ys_in_x_range requires t_wt to be traversable.");
    return yoff_iterator<t_wt>(wt, i, j);
}

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
 *  \return A pair. The first element of the pair indicates if
 *          a valid answer was found (true) or no valid answer (false)
 *          could be found. The second element contains the found symbol.
 */
template<class t_wt>
std::pair<bool,typename t_wt::value_type>
symbol_lte(const t_wt& wt, typename t_wt::value_type c)
{
    static_assert(t_wt::lex_ordered, "symbols_lte requires a lex_ordered WT");
    // check if wt has a built-in symbols_gte method
    constexpr bool has_own = has_symbols_wt<t_wt>::value;
    return _symbols_calls_wt<t_wt, has_own>::call_symbol_lte(wt,c);
}

//! Returns for a symbol c the next larger or equal symbol in the WT.
/*! \param c the symbol
 *  \return A pair. The first element of the pair indicates if
 *          a valid answer was found (true) or no valid answer (false)
 *          could be found. The second element contains the found symbol.
 */
template<class t_wt>
std::pair<bool,typename t_wt::value_type>
symbol_gte(const t_wt& wt, typename t_wt::value_type c)
{
    static_assert(t_wt::lex_ordered, "symbols_gte requires a lex_ordered WT");
    // check if wt has a built-in symbol_gte_method
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

template<typename T>
struct void_ { typedef void type; };

template<typename t_wt, typename T = void>
struct has_node_type {
    static constexpr std::false_type value = std::false_type();
};

template<typename t_wt>
struct has_node_type<t_wt, typename void_<typename t_wt::node_type>::type> {
    static constexpr std::true_type value = std::true_type();
};



} // end namespace

#endif
