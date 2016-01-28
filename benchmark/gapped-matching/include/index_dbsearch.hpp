#pragma once

#include "utils.hpp"
#include "collection.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/matching.hpp"

// implements double binary search
class index_dbsearch
{
    private:
        typedef sdsl::matching_index<sdsl::wt_int<>, sdsl::rrr_vector<>> index_type;
        index_type index;
        typedef std::pair<index_type::size_type, index_type::size_type> range_type;
        typedef std::pair<range_type, range_type> double_range_type;

    public:
        typedef sdsl::int_vector<0>::size_type size_type;
        std::string name() const
        {
            std::string index_name = IDXNAME;
            return "DBSEARCH-"+index_name;
        }

    public:
        index_dbsearch() { }
        index_dbsearch(collection& col)
        {
            sdsl::cache_config cc(false,".","DBSEARCH_TMP");
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

        void swap(index_dbsearch& ir)
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
            gapped_search_result res;
            string_type s1;
            string_type s2;
            size_type min_gap = 0;
            size_type max_gap = 0;

            if (pat.subpatterns.size() == 1) {
                s1 = pat.subpatterns[0];
                s2 = pat.subpatterns[0];
            } else {
                s1 = pat.subpatterns[0];
                s2 = pat.subpatterns[1];
                min_gap = pat.gaps[0].first + s1.size();
                max_gap = pat.gaps[0].second + s1.size();
            }

            // get ranges
            size_type sp1, ep1;
            size_type sp2, ep2;
            vector<size_type> range_a(forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, s1.begin(), s1.end(), sp1, ep1));
            vector<size_type> range_b(forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, s2.begin(), s2.end(), sp2, ep2));
            std::copy(index.wt.begin() + sp1, index.wt.begin() + ep1 + 1, range_a.begin());
            std::copy(index.wt.begin() + sp2, index.wt.begin() + ep2 + 1, range_b.begin());
            std::sort(range_a.begin(), range_a.end());
            std::sort(range_b.begin(), range_b.end());
            
            // double binary search
            std::stack<double_range_type> todo;
            todo.emplace(
                std::make_pair(0, range_a.size()), 
                std::make_pair(0, range_b.size()));
            while (!todo.empty())
            {                
                double_range_type top = todo.top(); todo.pop();
                range_type& subrange_a = top.first;
                range_type& subrange_b = top.second;
                auto size_a = subrange_a.second - subrange_a.first;
                auto size_b = subrange_b.second - subrange_b.first;
                if (size_a > 0 && size_b > 0)
                {
                    // always perform binary search on larger interval
                    if (size_a > size_b)
                    {
                        size_type median_index = (subrange_a.first + subrange_a.second) / 2;
                        size_type median_element = range_a[median_index];
                        auto lb = std::lower_bound(
                            range_b.begin() + subrange_b.first, 
                            range_b.begin() + subrange_b.second,
                            median_element + min_gap);
                        auto rb = std::upper_bound(
                            range_b.begin() + subrange_b.first, 
                            range_b.begin() + subrange_b.second,
                            median_element + max_gap);
                        size_type lb_idx = lb - range_b.begin();
                        size_type rb_idx = rb - range_b.begin();
                        
                        for (auto i = lb; i < rb; ++i)
                            res.positions.push_back(median_element);
                        
                        todo.emplace(
                            std::make_pair(subrange_a.first, median_index),
                            std::make_pair(subrange_b.first, rb_idx));
                        todo.emplace(
                            std::make_pair(median_index + 1, subrange_a.second),
                            std::make_pair(lb_idx, subrange_b.second));
                    }
                    else
                    {
                        size_type median_index = (subrange_b.first + subrange_b.second) / 2;
                        size_type median_element = range_b[median_index];
                        auto lb = std::lower_bound(
                            range_a.begin() + subrange_a.first, 
                            range_a.begin() + subrange_a.second,
                            median_element < max_gap ? 0 : (median_element - max_gap));
                        auto rb = std::upper_bound(
                            range_a.begin() + subrange_a.first, 
                            range_a.begin() + subrange_a.second,
                            median_element < min_gap ? 0 : (median_element - min_gap));
                        size_type lb_idx = lb - range_a.begin();
                        size_type rb_idx = rb - range_a.begin();
                        
                        for (auto i = lb; i < rb; ++i)
                            res.positions.push_back(*i);
                        
                        todo.emplace(
                            std::make_pair(subrange_a.first, rb_idx),
                            std::make_pair(subrange_b.first, median_index));
                        todo.emplace(
                            std::make_pair(lb_idx, subrange_a.second),
                            std::make_pair(median_index + 1, subrange_b.second));
                    }
                }
            }
            
            return res;
        }
};
