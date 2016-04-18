#ifndef INCLUDED_SDSL_VLG_INDEX
#define INCLUDED_SDSL_VLG_INDEX

#include "suffix_arrays.hpp"
#include <vector>

namespace sdsl
{

template<typename alphabet_tag>
struct gapped_pattern {
    typedef int_vector<alphabet_tag::WIDTH> string_type;
    
    std::vector<string_type> subpatterns;
    std::vector<std::pair<uint64_t,uint64_t>> gaps;
    gapped_pattern() { }
    gapped_pattern(const std::string& raw_regexp)
    {
        bool string_patterns = alphabet_tag::WIDTH == 8;
        auto parse_subpattern = [&](std::string s) {
            string_type subpattern;
            if (string_patterns) {
                subpattern.resize(s.size());
                std::copy(s.begin(), s.end(), subpattern.begin());
            } else {
                std::istringstream symbol_stream(s);
                uint64_t symbol;
                while (symbol_stream >> symbol) {
                    auto index = subpattern.size();
                    subpattern.resize(index + 1);
                    subpattern[index] = symbol;
                }
            }
            return subpattern;
        };

        int64_t last_gap_end = -1;
        size_t gap_pos;
        while ((gap_pos = raw_regexp.find(".{", last_gap_end + 1)) != std::string::npos) {
            // extract subpattern
            auto subptrlen = gap_pos - (last_gap_end+1);
            auto subpattern = raw_regexp.substr(last_gap_end+1,subptrlen);
            auto subptr = parse_subpattern(subpattern);
            subpatterns.push_back(subptr);
            
            // extract gap
            auto gap_end = raw_regexp.find("}", gap_pos) + 1;
            std::string gap_str = raw_regexp.substr(gap_pos, gap_end - gap_pos);
            /* try parsing the gap  description */
            auto first_num_sep = gap_str.find(",");
            auto first_num_str = gap_str.substr(2,first_num_sep-2);
            auto second_num_str = gap_str.substr(first_num_sep+1);
            second_num_str.pop_back();
            uint64_t first_gap_num = std::stoull(first_num_str);
            uint64_t second_gap_num = std::stoull(second_num_str);
            if (first_gap_num > second_gap_num) {
                throw std::runtime_error("invalid gap description: min-gap > max-gap");
            }
            gaps.emplace_back(first_gap_num + subptr.size(),second_gap_num + subptr.size());
            // validate lazy semantics
            if (gap_end == raw_regexp.size() || raw_regexp[gap_end] != '?') {
                throw std::runtime_error("invalid gap description: expected '?' (lazy semantics)");
            }
            last_gap_end = gap_end;
        }
        // extract remaining subpattern
        auto last_subpattern = raw_regexp.substr(last_gap_end+1);
        subpatterns.push_back(parse_subpattern(last_subpattern));
    };
};

template<typename type_index>
struct node_cache;

template<typename alphabet_tag=byte_alphabet_tag,
         typename t_wt=wt_int<
             bit_vector_il<>,
             rank_support_il<>>>
class vlg_index
{
        static_assert(std::is_same<typename index_tag<t_wt>::type, wt_tag>::value,
                      "Second template argument has to be a wavelet tree.");
        
    private:
        typedef vlg_index<alphabet_tag, t_wt>   index_type;
    public:
        typedef alphabet_tag                    alphabet_category;
        typedef t_wt                          wt_type;
        typedef typename wt_type::node_type   node_type;
        typedef typename wt_type::size_type   size_type;

        typedef int_vector<alphabet_tag::WIDTH>  text_type;


    private:
        text_type m_text;
        wt_type   m_wt;
        
    public:
        const text_type& text = m_text;
        const wt_type&   wt  = m_wt;

        //! Default constructor
        vlg_index() = default;

        //! Copy constructor
        vlg_index(const vlg_index& idx)
            : m_text(idx.m_text), m_wt(idx.m_wt)
        {
        }

        //! Copy constructor
        vlg_index(vlg_index&& idx)
        {
            *this = std::move(idx);
        }

        vlg_index(text_type text, wt_type wt)
            : m_text(text), m_wt(wt)
        { }

        //! Assignment move operator
        vlg_index& operator=(vlg_index&& idx)
        {
            if (this != &idx) {
                m_text     = std::move(idx.m_text);
                m_wt       = std::move(idx.m_wt);
            }
            return *this;
        }

        //! Swap operation
        void swap(vlg_index& idx)
        {
            if (this != &idx) {
                m_text.swap(idx.m_text);
                m_wt.swap(idx.m_wt);
            }
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_text.serialize(out, child, "text");
            written_bytes += m_wt.serialize(out, child, "wt");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in)
        {
            m_text.load(in);
            m_wt.load(in);
        }
};

template<typename type_index>
struct node_cache {
    typedef typename type_index::node_type node_type;
    typedef typename type_index::size_type size_type;

    const type_index& index;
    node_type node;
    size_type range_begin;
    size_type range_end;
    bool is_leaf;

    size_type range_size()
    {
        return range_end - range_begin + 1;
    }

    node_cache(
        node_type node,
        const type_index& index)
        : index(index)
    {
        this->node = node;
        auto range = index.wt.value_range(node);
        this->range_begin = std::get<0>(range);
        this->range_end = std::get<1>(range);
        this->is_leaf = index.wt.is_leaf(node);
    }
};

template<typename type_index>
class wavelet_tree_range_walker
{
    private:
        typedef node_cache<type_index> node_type;
        const type_index& index;
        std::vector<std::pair<range_type,node_type>> dfs_stack;

    public:
        typedef decltype(dfs_stack) state_type;
        wavelet_tree_range_walker(const type_index& index, range_type initial_range, node_type root_node)
            : index(index)
        {
            dfs_stack.reserve(64);
            dfs_stack.emplace_back(initial_range, root_node);
        }

        inline bool has_more() const
        {
            return !dfs_stack.empty();
        }

        inline node_type current_node() const
        {
            return dfs_stack.back().second;
        }

        inline void skip_subtree()
        {
            dfs_stack.pop_back();
        }

        inline void expand()
        {
            auto top = dfs_stack.back(); dfs_stack.pop_back();
            auto& node = top.second;
            auto children = index.wt.expand(node.node);
            auto exp_range = index.wt.expand(node.node, top.first);
            if (!empty(exp_range[1]))
                dfs_stack.emplace_back(exp_range[1], node_type(children[1], index));
            if (!empty(exp_range[0]))
                dfs_stack.emplace_back(exp_range[0], node_type(children[0], index));
        }

        inline bool move_next_leaf()
        {
            if (has_more() && current_node().is_leaf)
                skip_subtree();
            while (has_more() && !current_node().is_leaf)
                expand();
            return has_more();
        }
};

template<typename type_index>
class vlg_iterator : public std::iterator<std::forward_iterator_tag, typename type_index::size_type>
{
    private:
        typedef typename type_index::node_type node_type;
        typedef typename type_index::size_type size_type;
        typedef typename type_index::wt_type   wt_type;
        typedef typename type_index::size_type result_type;

        // (lex_range, node)
        std::vector<wavelet_tree_range_walker<type_index>> lex_ranges;

        std::vector<std::pair<uint64_t,uint64_t>> gaps;
        size_type last_subpattern_size;

        bool finished;

        bool relax() {
            bool redo = true;
            while (redo)
            {
                redo = false;
                for (size_t i = 1; i < size(); ++i) {
                    if (lex_ranges[i - 1].current_node().range_end + gaps[i - 1].second < lex_ranges[i].current_node().range_begin) {
                        lex_ranges[i - 1].skip_subtree();
                        redo = true;
                        if (!lex_ranges[i - 1].has_more())
                            return false;
                    }
                    if (lex_ranges[i - 1].current_node().range_begin + gaps[i - 1].first > lex_ranges[i].current_node().range_end) {
                        lex_ranges[i].skip_subtree();
                        redo = true;
                        if (!lex_ranges[i].has_more())
                            return false;
                    }
                }
            }
            return true;
        }

        bool pull_forward()
        {
            // pull a forward
            auto x = lex_ranges[lex_ranges.size() - 1].current_node().range_begin;

            while (lex_ranges[0].has_more() && lex_ranges[0].current_node().range_end <= x) lex_ranges[0].skip_subtree();
            while (lex_ranges[0].move_next_leaf() && lex_ranges[0].current_node().range_begin < x + last_subpattern_size) ;
            
            return lex_ranges[0].has_more();
        }
        void next()
        {
            while (relax()) {
                size_t r = 1;
                size_t j = 0;
                bool skip = false;
                for (size_t i = 0; i < size(); ++i) {
                    auto lr = lex_ranges[i].current_node().range_size();
                    if (lr > r) {
                        r = lr;
                        j = i;
                        skip = true;
                    }
                }

                if (skip)
                    lex_ranges[j].expand();
                else
                    return;
            }
            
            finished = true;
        }

    public:
        vlg_iterator()
        {
            finished = true;
        }
        vlg_iterator(const type_index& index,
                     const gapped_pattern<typename type_index::alphabet_category>& pattern)
            : gaps(pattern.gaps), last_subpattern_size(pattern.subpatterns[pattern.subpatterns.size() - 1].size())
        {
            finished = false;
            
            auto root_node = node_cache<type_index>(index.wt.root(), index);
            size_type sp = 1, ep = 0;

            for (auto sx : pattern.subpatterns) {
                forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, sx.begin(), sx.end(), sp, ep);
                lex_ranges.emplace_back(index, range_type(sp, ep), root_node);
            }
            if (valid())
                next();
        }

        bool valid() const
        {
            for (auto lr : lex_ranges)
                if (!lr.has_more())
                    return false;
            return true;
        }

        size_t size() const
        {
            return lex_ranges.size();
        }
        result_type operator[](int subpattern_index) const
        {
            return lex_ranges[subpattern_index].current_node().range_begin;
        }
        result_type const operator*() const
        {
            return this->operator[](0);
        }
        result_type* operator->()
        {
            return &this->operator[](0);
        }

        vlg_iterator& operator++()
        {
            if (pull_forward())
                next();
            else
                finished = true;
            return *this;
        }

        friend bool operator==(
            const vlg_iterator& a,
            const vlg_iterator& b)
        {
            return a.finished && b.finished;
        }

        friend bool operator!=(
            const vlg_iterator& a,
            const vlg_iterator& b)
        {
            return !(a == b);
        }
};

template<typename alphabet_tag, typename t_wt>
void construct(vlg_index<alphabet_tag, t_wt>& idx, const std::string& file, cache_config& config, uint8_t num_bytes)
{
    config = sdsl::cache_config(false,".","WCSEARCH_TMP"); // TODO: make better!!! keep in mind: wts needs SA!
                    
    int_vector<alphabet_tag::WIDTH> text;
    load_vector_from_file(text, file, num_bytes);

    csa_wt<wt_int<>> csa;
    construct(csa, file, config, num_bytes);

    t_wt wts;
    construct(wts, cache_file_name(conf::KEY_SA, config));

    util::delete_all_files(config.file_map);

    idx = std::move(vlg_index<alphabet_tag, t_wt>(text, wts));
}

template<typename iter>
class container
{
    private:
        const iter m_it_begin;
        const iter m_it_end;

    public:
        container(const iter it_begin, const iter it_end)
            : m_it_begin(it_begin), m_it_end(it_end) { }

        iter begin() const {
            return m_it_begin;
        }
        iter end() const {
            return m_it_end;
        }
};

template<typename type_index>
container<vlg_iterator<type_index>> locate(type_index idx, std::string pattern) {
    return container<vlg_iterator<type_index>>(
        vlg_iterator<type_index>(idx, pattern),
        vlg_iterator<type_index>()
    );
}

template<typename type_index>
typename type_index::size_type count(type_index idx, std::string pattern) {
    typename type_index::size_type result = 0;
    auto cont = locate(idx, pattern);
    for (auto it = cont.begin(); it != cont.end(); ++it)
        ++result;
    return result;
}

} // end namespace sdsl
#endif
