#pragma once

#include <vector>
#include "collection.hpp"
#include "utils.hpp"
#include "sdsl/suffix_arrays.hpp"

template<class type_index>
struct node_cache;

template<class t_wt=sdsl::wt_int<
             sdsl::bit_vector_il<>,
             sdsl::rank_support_il<>>,
         class t_bv=sdsl::rrr_vector<>>
class matching_index
{
    private:
        typedef matching_index<t_wt, t_bv> index_type;
    public:
        typedef t_wt                          wt_type;
        typedef t_bv                          bv_type;
        typedef typename bv_type::rank_1_type rank_type;
        typedef typename wt_type::node_type   node_type;
        typedef typename wt_type::size_type   size_type;

        typedef sdsl::int_vector<0>           text_type;


    private:
        text_type m_text;
        wt_type   m_wt;
        
    public:
        const text_type& text = m_text;
        const wt_type&   wt  = m_wt;

        //! Default constructor
        matching_index() = default;

        //! Copy constructor
        matching_index(const matching_index& idx)
            : m_text(idx.m_text), m_wt(idx.m_wt)
        {
        }

        //! Copy constructor
        matching_index(matching_index&& idx)
        {
            *this = std::move(idx);
        }

        matching_index(text_type text, wt_type wt)
            : m_text(text), m_wt(wt)
        { }

        //! Assignment move operator
        matching_index& operator=(matching_index&& idx)
        {
            if (this != &idx) {
                m_text     = std::move(idx.m_text);
                m_wt       = std::move(idx.m_wt);
            }
            return *this;
        }

        //! Swap operation
        void swap(matching_index& idx)
        {
            if (this != &idx) {
                m_text.swap(idx.m_text);
                m_wt.swap(idx.m_wt);
            }
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_text.serialize(out, child, "text");
            written_bytes += m_wt.serialize(out, child, "wt");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in)
        {
            m_text.load(in);
            m_wt.load(in);
        }
};

template<class t_wt, class t_bv>
void construct(matching_index<t_wt, t_bv>& idx, const std::string& file, sdsl::cache_config& config, uint8_t num_bytes)
{
    sdsl::int_vector<0> text;
    {
        //auto event = memory_monitor::event("text");
        load_vector_from_file(text, file, num_bytes);
    }
    sdsl::csa_wt<sdsl::wt_int<>> csa;
    {
        //auto event = memory_monitor::event("csa");
        construct(csa, file, config, num_bytes);
    }
    t_wt wts;
    {
        //auto event = memory_monitor::event("wt");
        construct(wts, cache_file_name(sdsl::conf::KEY_SA, config));
    }

    sdsl::util::delete_all_files(config.file_map);

    {
        //auto event = memory_monitor::event("compose"); // contains rank support initialization
        idx = std::move(matching_index<t_wt, t_bv>(text, wts));
    }
}

template<class type_index>
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

template<class type_index>
class wavelet_tree_range_walker
{
    private:
        typedef node_cache<type_index> node_type;
        const type_index& index;
        std::vector<std::pair<sdsl::range_type,node_type>> dfs_stack;

    public:
        typedef decltype(dfs_stack) state_type;
        wavelet_tree_range_walker(const type_index& index, sdsl::range_type initial_range, node_type root_node)
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
            if (!sdsl::empty(exp_range[1]))
                dfs_stack.emplace_back(exp_range[1], node_type(children[1], index));
            if (!sdsl::empty(exp_range[0]))
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

template<class type_index>
class wild_card_match_iterator3 : public std::iterator<std::forward_iterator_tag, typename type_index::size_type>
{
    private:
        typedef typename type_index::node_type node_type;
        typedef typename type_index::size_type size_type;
        typedef typename type_index::wt_type   wt_type;
        typedef typename type_index::size_type result_type;

        // (lex_range, node)
        std::vector<wavelet_tree_range_walker<type_index>> lex_ranges;

        size_t min_gap;
        size_t max_gap;

        size_t size3;

        result_type current;

        bool relax() {
            bool redo = true;
            while (redo)
            {
                redo = false;
                for (size_t i = 1; i < lex_ranges.size(); ++i) {
                    if (lex_ranges[i - 1].current_node().range_end + max_gap < lex_ranges[i].current_node().range_begin) {
                        lex_ranges[i - 1].skip_subtree();
                        redo = true;
                        if (!lex_ranges[i - 1].has_more())
                            return false;
                    }
                    if (lex_ranges[i - 1].current_node().range_begin + min_gap > lex_ranges[i].current_node().range_end) {
                        lex_ranges[i].skip_subtree();
                        redo = true;
                        if (!lex_ranges[i].has_more())
                            return false;
                    }
                }
            }
            return true;
        }

        bool next()
        {
            //if (valid())
            // find next independent match
            while (relax()) {
                size_t r = 1;
                size_t j;
                bool skip = false;
                for (size_t i = 0; i < lex_ranges.size(); ++i) {
                    auto lr = lex_ranges[i].current_node().range_size();
                    if (lr > r) {
                        r = lr;
                        j = i;
                        skip = true;
                    }
                }

                if (skip)
                    lex_ranges[j].expand();
                else {
                    current = lex_ranges[0].current_node().range_begin;
                    //std::cerr << "HIT" << std::endl;
                    //for (size_t i = 0; i < lex_ranges.size(); ++i)
                    //    std::cerr << "\t" << lex_ranges[i].current_node()->range_begin << std::endl;

                    auto x = lex_ranges[lex_ranges.size() - 1].current_node().range_begin;

                    // pull a forward
                    while (lex_ranges[0].has_more() && lex_ranges[0].current_node().range_end <= x) lex_ranges[0].skip_subtree();
                    while (lex_ranges[0].move_next_leaf() && lex_ranges[0].current_node().range_begin < x + size3) ;

                    return true;
                }
            }

            current = (result_type)-1;
            return false;
        }

    public:
        wild_card_match_iterator3()
        {
            current = (result_type)-1;
        }
        wild_card_match_iterator3(const type_index& index,
                                  const std::vector<string_type>& s,
                                  size_t min_gap,
                                  size_t max_gap)
            : min_gap(min_gap), max_gap(max_gap), size3(s[s.size() - 1].size())
        {
            
            auto root_node = node_cache<type_index>(index.wt.root(), index);
            size_type sp = 1, ep = 0;

            for (auto sx : s) {
                forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, sx.begin(), sx.end(), sp, ep);
                lex_ranges.emplace_back(index, sdsl::range_type(sp, ep), root_node);
                //std::cerr << std::string(sx.begin(), sx.end()) << ": " << sp << " " << ep << std::endl;
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

        result_type operator*() const
        {
            return current;
        }
        result_type* operator->()
        {
            return &current;
        }

        wild_card_match_iterator3& operator++()
        {
            next();
            return *this;
        }

        friend bool operator==(
            const wild_card_match_iterator3& a,
            const wild_card_match_iterator3& b)
        {
            return a.current == b.current;
        }

        friend bool operator!=(
            const wild_card_match_iterator3& a,
            const wild_card_match_iterator3& b)
        {
            return !(a == b);
        }
};



class index_wcsearch3
{
    private:
        typedef matching_index<> index_type;
        index_type index;

    public:
        typedef sdsl::int_vector<0>::size_type size_type;
        std::string name() const
        {
            std::string index_name = IDXNAME;
            return "WCSEARCH-"+index_name;
        }

    public:
        index_wcsearch3() { }
        index_wcsearch3(collection& col)
        {
            sdsl::cache_config cc(false,".","WCSEARCH_TMP");
            construct(index, col.file_map[consts::KEY_TEXT], cc, 0);
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=NULL, std::string name="")const
        {
            return index.serialize(out, v, name);
        }

        void load(std::istream& in)
        {
            index.load(in);
        }

        void swap(index_wcsearch3& ir)
        {
            if (this != &ir) {
                index.swap(ir.index);
            }
        }

        std::string info(const gapped_pattern& pat) const
        {
            // output SA-ranges (gives a good estimation about potential matches)
            index_type::size_type total_range = 0, sp = 0, ep = 0;

            for (size_t i = 0; i < pat.subpatterns.size(); ++i)
                total_range += forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, pat.subpatterns[i].begin(), pat.subpatterns[i].end(), sp, ep);

            return std::to_string(total_range);
        }
        void prepare(const gapped_pattern& pat)
        {
            (void)pat;
        }

        gapped_search_result
        search(const gapped_pattern& pat) const
        {
            gapped_search_result res;
            size_type min_gap;
            size_type max_gap;

            std::cerr << "REGEX ::: " << pat.raw_regexp << std::endl;

            min_gap = pat.subpatterns[0].size() + pat.gaps[0].first;
            max_gap = pat.subpatterns[0].size() + pat.gaps[0].second;

            // smart scan
            auto it = wild_card_match_iterator3<index_type>(index, pat.subpatterns, min_gap, max_gap);
            auto ite = wild_card_match_iterator3<index_type>();
            for (; it != ite; ++it) {
                res.positions.push_back(*it);
            }
            return res;
        }
};
