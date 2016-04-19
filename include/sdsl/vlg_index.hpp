/* sdsl - succinct data structures library
    Copyright (C) 2016 Simon Gog, Johannes Bader

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file vlg_index.hpp
    \brief vlg_index.hpp contains an index supporting count and locate queries
            conaining gaps of variable length. This is an implementation of the
            algorithm described in the SEA paper 
            ,,Practical Variable Length Gap Pattern Matching using Wavelet Trees''
            of Johannes Bader, Simon Gog and Matthias Petri.
    \author Johannes Bader, Matthias Petri, Simon Gog
*/
#ifndef INCLUDED_SDSL_VLG_INDEX
#define INCLUDED_SDSL_VLG_INDEX

#include "suffix_arrays.hpp"
#include <vector>

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! Struct representing a query with variable length gaps.
/*!
 * \tparam alphabet_tag   Type of alphabet used by the query (and text to be searched within).
 *
 * This struct provides access to the individual subpatterns and gaps between them.
 * Each gap is represented by its minimum and maximum length.
 * The struct also provides functionality for parsing regex-style query strings.
 */
template<typename alphabet_tag>
struct gapped_pattern_query {
    typedef int_vector<alphabet_tag::WIDTH> string_type;

    std::vector<string_type> subpatterns;
    std::vector<std::pair<uint64_t,uint64_t>> gaps;

    //! Default constructor
    gapped_pattern_query() { }

    //! Constructor parsing the provided regex query
    gapped_pattern_query(const std::string& raw_regexp)
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
            if (gap_end == raw_regexp.size() or raw_regexp[gap_end] != '?') {
                throw std::runtime_error("invalid gap description: expected '?' (lazy semantics)");
            }
            last_gap_end = gap_end;
        }
        // extract remaining subpattern
        auto last_subpattern = raw_regexp.substr(last_gap_end+1);
        subpatterns.push_back(parse_subpattern(last_subpattern));
    };
};

//! An index supporting variable length gap pattern matching.
/*!
 * \tparam alphabet_tag   Type of alphabet used by the indexed text and thus also the index.
 * \tparam t_wt           Wavelet tree used for storing the suffix array.
 *
 * This class provides the datastructures required for variable length gap pattern matching.
 * Specifically, it contains the original text as well as a wavelet tree over the suffix array.
 */
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
        typedef alphabet_tag                  alphabet_category;
        typedef t_wt                          wt_type;
        typedef typename wt_type::node_type   node_type;
        typedef typename wt_type::size_type   size_type;

        typedef int_vector<alphabet_tag::WIDTH>    text_type;
        typedef gapped_pattern_query<alphabet_tag> query_type;


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
        { }

        //! Copy constructor
        vlg_index(vlg_index&& idx)
        {
            *this = std::move(idx);
        }

        //! Constructor
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

//! Stores a wavelet tree node together with some precomputed properties.
/*!
 * \tparam wt_type   Type of wavelet tree this node is from.
 */
template<typename wt_type>
struct wt_node_cache {
    typedef typename wt_type::node_type node_type;
    typedef typename wt_type::size_type size_type;

    const node_type node;
    size_type range_begin;
    size_type range_end;
    bool is_leaf;

    //! Constructor, precalculates frequently used values.
    wt_node_cache(
        const node_type& node,
        const wt_type& wt)
            : node(node)
    {
        auto range = wt.value_range(node);
        this->range_begin = std::get<0>(range);
        this->range_end = std::get<1>(range);
        this->is_leaf = wt.is_leaf(node);
    }
};

//! Provides a convenient way of traversing a wavelet tree from left to right.
/*!
 * \tparam wt_type   Type of wavelet tree to traverse.
 *
 * This class keeps track of the state required to perform a depth-first traversal 
 * of a wavelet tree. Furthermore it provides methods to traverse the tree in 
 * various ways.
 */
template<typename wt_type>
class wt_range_walker
{
    private:
        typedef wt_node_cache<wt_type> node_type;
        const wt_type& wt;
        std::vector<std::pair<range_type,node_type>> dfs_stack;

    public:
        //! Constructor
        wt_range_walker(const wt_type& wt, range_type initial_range, node_type root_node)
            : wt(wt)
        {
            dfs_stack.reserve(64); // TODO: magic number? or rather something with wt_type::size_type? or is there something like max_depth(wt)?
            dfs_stack.emplace_back(initial_range, root_node);
        }

        //! Returns whether the traversal has not yet reached the end of the wavelet tree.
        inline bool has_more() const
        {
            return !dfs_stack.empty();
        }

        //! Returns the wavelet tree node currently pointed at by the walker.
        inline node_type current_node() const
        {
            return dfs_stack.back().second;
        }

        //! Traverse to the next node, discarding any child nodes of the current node.
        inline void next_right()
        {
            dfs_stack.pop_back();
        }

        //! Traverse to the first non-empty child node of the current node.
        inline void next_down()
        {
            auto top = dfs_stack.back(); dfs_stack.pop_back();
            auto& node = top.second;
            auto children = wt.expand(node.node);
            auto exp_range = wt.expand(node.node, top.first);
            if (!empty(exp_range[1]))
                dfs_stack.emplace_back(exp_range[1], node_type(children[1], wt));
            if (!empty(exp_range[0]))
                dfs_stack.emplace_back(exp_range[0], node_type(children[0], wt));
        }

        //! Traverse to the next leaf. Returns false if there is no more, i.e. the traversal has finished.
        inline bool next_leaf()
        {
            if (has_more() and current_node().is_leaf)
                next_right();
            while (has_more() and !current_node().is_leaf)
                next_down();
            return has_more();
        }
};

//! An iterator implementing the variable length gap pattern search as described in the paper.
/*!
 * \tparam type_index   Type of index to use for the search.
 *
 * This class provides the variable length gap pattern search functionality in an on-demand fashion.
 * As long as the iterator is valid (i.e. !is_end()), 
 * it provides vector-like access to the subpattern positions of the current match.
 */
template<typename type_index>
class vlg_iterator : public std::iterator<std::forward_iterator_tag, typename type_index::size_type>
{
    private:
        typedef typename type_index::node_type node_type;
        typedef typename type_index::size_type size_type;
        typedef typename type_index::wt_type   wt_type;

        // current state of iteration
        std::vector<wt_range_walker<wt_type>> lex_ranges;
        bool finished = true;

        // required query information
        const std::vector<std::pair<uint64_t,uint64_t>> gaps;
        const size_type last_subpattern_size;

        // Conservatively enforces gap constraints.
        // Returns false if the iteration has finished due to this operation.
        bool relax()
        {
            bool redo = true;
            while (redo) {
                redo = false;
                for (size_t i = 1; i < size(); ++i) {
                    if (lex_ranges[i - 1].current_node().range_end + gaps[i - 1].second < lex_ranges[i].current_node().range_begin) {
                        lex_ranges[i - 1].next_right();
                        redo = true;
                        if (!lex_ranges[i - 1].has_more())
                            return false;
                    }
                    if (lex_ranges[i - 1].current_node().range_begin + gaps[i - 1].first > lex_ranges[i].current_node().range_end) {
                        lex_ranges[i].next_right();
                        redo = true;
                        if (!lex_ranges[i].has_more())
                            return false;
                    }
                }
            }
            return true;
        }

        // Pulls first subpattern position behind last subpattern position.
        // => enforces non-overlapping match semantics
        // Returns false if the iteration has finished due to this operation.
        bool pull_forward()
        {
            auto last_pos = lex_ranges[lex_ranges.size() - 1].current_node().range_begin;

            // skip entire subtrees as long as save
            while (lex_ranges[0].has_more() and 
                lex_ranges[0].current_node().range_end <= last_pos) lex_ranges[0].next_right();
            // finds first leaf with required position
            while (lex_ranges[0].next_leaf() and 
                lex_ranges[0].current_node().range_begin < last_pos + last_subpattern_size) ;
            return lex_ranges[0].has_more();
        }

        // Finds the next match of the query.
        void next()
        {
            // While relaxation has not reached the end of the wavelet tree...
            while (relax()) {
                // ...determine the largest wavelet tree node
                size_t r = 1;
                size_t j = 0;
                bool found = false;
                for (size_t i = 0; i < size(); ++i) {
                    auto lr = lex_ranges[i].current_node().node.size;
                    if (lr > r) {
                        r = lr;
                        j = i;
                        found = true;
                    }
                }

                if (found) // if there is one, expand it
                    lex_ranges[j].next_down();
                else       // otherwise we found a match!
                    return;
            }

            finished = true;
        }

    public:
        // Type of subpattern positions pointed to by this itercr
        
        //! Default constructor.
        vlg_iterator()
            : gaps()
            , last_subpattern_size(0) { }

        //! Constructor.
        vlg_iterator(const type_index& index,
                     const typename type_index::query_type& query)
            : gaps(query.gaps)
            , last_subpattern_size(query.subpatterns[query.subpatterns.size() - 1].size())
        {
            // initialize wavelet tree iterators using the SA range of each subpattern
            auto root_node = wt_node_cache<wt_type>(index.wt.root(), index.wt);
            for (auto sx : query.subpatterns) {
                size_type sp, ep;
                forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, sx.begin(), sx.end(), sp, ep);
                lex_ranges.emplace_back(index.wt, range_type(sp, ep), root_node);

                // shortcut on empty range
                if (sp > ep) return;
            }

            // find first match
            finished = false;
            next();
        }

        //! Returns whether this iterator has ended, i.e. does not point to a match anymore.
        bool is_end() const
        {
            return finished;
        }

        //! Returns the number of subpattern positions this iterator points to.
        size_t size() const
        {
            return lex_ranges.size();
        }

        //! Returns a subpatterns text position of the current match.
        position_type operator[](int subpattern_index) const
        {
            return lex_ranges[subpattern_index].current_node().range_begin;
        }

        //! Returns the starting position of the current match (which is the first subpattern's position).
        position_type const operator*() const
        {
            return this->operator[](0);
        }

        //! Advances the iterator.
        vlg_iterator& operator++()
        {
            if (pull_forward())
                next();
            else
                finished = true;
            return *this;
        }

        //! Checks to iterators for equality (which only holds for two end-iterators).
        friend bool operator==(
            const vlg_iterator& a,
            const vlg_iterator& b)
        {
            return a.is_end() and b.is_end();
        }

        //! Compares to iterators for inequality .
        friend bool operator!=(
            const vlg_iterator& a,
            const vlg_iterator& b)
        {
            return !(a == b);
        }
};

// Pseudo-container encapsulating begin and end iterator.
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

// Retrieves a container representing all occurrences of the provided pattern.
template<typename type_index>
container<vlg_iterator<type_index>> locate(const type_index& idx, const typename type_index::query_type& pattern) {
    return container<vlg_iterator<type_index>>(
        vlg_iterator<type_index>(idx, pattern),
        vlg_iterator<type_index>()
    );
}

// Retrieves the number of occurrences of the provided pattern.
template<typename type_index>
typename type_index::size_type count(const type_index& idx, const typename type_index::query_type& pattern) {
    typename type_index::size_type result = 0;
    auto cont = locate(idx, pattern);
    for (auto it = cont.begin(); it != cont.end(); ++it)
        ++result;
    return result;
}

} // end namespace sdsl
#endif
