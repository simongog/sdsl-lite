#pragma once

#include "utils.hpp"
#include "collection.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/matching.hpp"

// NO DOCUMENT BOUNDATY SUPPORT!
class index_wcsearch_full
{
    private:
        typedef sdsl::matching_index<sdsl::wt_int<>, sdsl::rrr_vector<>> index_type;
        index_type index;

    public:
        typedef typename index_type::node_type node_type;
        typedef typename index_type::size_type size_type;
        typedef typename index_type::wt_type   wt_type;
        typedef pair<size_type, size_type>     result_type;
        std::string name() const
        {
            std::string index_name = IDXNAME;
            return "WCSEARCH-FULL-"+index_name;
        }

    public:
        index_wcsearch_full() { }
        index_wcsearch_full(collection& col)
        {
            sdsl::cache_config cc(false,".","WCSEARCH_FULL_TMP");
            sdsl::construct(index, col.file_map[consts::KEY_TEXT], cc, 0);
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=NULL, std::string name="")const
        {
            return index.serialize(out, v, name);
        }

        void load(std::istream& in)
        {
            index.load(in);
        }

        void swap(index_wcsearch_full& ir)
        {
            if (this != &ir) {
                index.swap(ir.index);
            }
        }

        std::string info(const gapped_pattern& pat) const { (void)pat; return ""; }
        void prepare(const gapped_pattern& pat) { (void)pat; }

        gapped_search_result
        search(const gapped_pattern& pat) const
        {
            sdsl::PERFCTR_NUM_PROCESSED_WT_NODES = 1;
            
            gapped_search_result res;
            string_type s1;
            string_type s2;
            size_type min_gap;
            size_type max_gap;

            if (pat.subpatterns.size() == 1) {
                s1 = pat.subpatterns[0];
                s2 = pat.subpatterns[0];
                min_gap = 0;
                max_gap = 0;
            } else {
                s1 = pat.subpatterns[0];
                s2 = pat.subpatterns[1];
                min_gap = s1.size() + pat.gaps[0].first;
                max_gap = s1.size() + pat.gaps[0].second;
            }

            // SA ranges
            auto root_node = make_shared<sdsl::node_cache<index_type>>(index.wt.root(), index, nullptr, nullptr);
            size_type sp, ep;

            forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, s1.begin(), s1.end(), sp, ep);
            sdsl::wavelet_tree_range_walker<index_type> lex_range0(index, sdsl::range_type(sp, ep),root_node);
            
            forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, s2.begin(), s2.end(), sp, ep);
            sdsl::wavelet_tree_range_walker<index_type> lex_range1(index, sdsl::range_type(sp, ep),root_node);
            
            // iterate
            size_type a;
            deque<size_type> bs;
            while (lex_range0.has_more()) {
                const auto& top0 = lex_range0.current_node();
                if (lex_range1.has_more() && top0->range_begin + min_gap > lex_range1.current_node()->range_end)
                    lex_range1.skip_subtree();
                else if (!bs.empty() && top0->range_end + max_gap < bs.front())
                    lex_range0.skip_subtree();
                else if (top0->is_leaf) {                    
                    a = top0->range_begin;
                    lex_range0.skip_subtree();

                    // push b forward
                    while (lex_range1.has_more()
                           && a + max_gap >= lex_range1.current_node()->range_begin)
                    {
                        auto leaf = lex_range1.retrieve_leaf_and_traverse();
                        if (leaf != nullptr)
                            bs.push_back(leaf->range_begin);
                    }

                    // shrink bs
                    while (!bs.empty() && a + min_gap > bs.front())
                        bs.pop_front();
                           
                    // output
                    for (size_t i = 0; i < bs.size(); ++i)
                        res.positions.push_back(a);
                    
                } else
                    lex_range0.expand();
            }
            
            cout << "Processed nodes: " << sdsl::PERFCTR_NUM_PROCESSED_WT_NODES << endl;
            return res;
        }
};
