#pragma once

#include "collection.hpp"
#include "utils.hpp"
#include "sdsl/matching.hpp"
#include "sdsl/suffix_arrays.hpp"

#include <boost/regex.hpp>

#define CONST_LINEAR_THRESH 35

class index_wcsearch
{
    private:
        typedef sdsl::matching_index<sdsl::wt_int<>, sdsl::rrr_vector<>> index_type;
        index_type index;
        string text = "";

    public:
        typedef sdsl::int_vector<0>::size_type size_type;
        std::string name() const
        {
            std::string index_name = IDXNAME;
            return "WCSEARCH-"+index_name;
        }

    public:
        index_wcsearch() { }
        index_wcsearch(collection& col)
        {
            sdsl::cache_config cc(false,".","WCSEARCH_TMP");
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

        void swap(index_wcsearch& ir)
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
            if (index.text.width() <= 8)
                text = string(index.text.begin(), index.text.end());
        }

        gapped_search_result
        search(const gapped_pattern& pat) const
        {
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
                // add "s1.size()" because "match2" currently requires word-beginning-relative gaps
                // (this is an important concept, as it allows single-term matching by setting min/max_gap=0)
            }

            // linear scan?
            if (text.size() > 0)
            {
                index_type::size_type total_range = 0, sp = 0, ep = 0;
                for (size_t i = 0; i < pat.subpatterns.size(); ++i)
                    total_range += forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, pat.subpatterns[i].begin(), pat.subpatterns[i].end(), sp, ep);
                
                // check for: |range| * log n > n
                LOG(INFO) << "range = " << total_range << "; |text| = " << text.size();
                total_range *= sdsl::bits::hi(text.size());
                total_range *= CONST_LINEAR_THRESH;
                LOG(INFO) << total_range << " > " << text.size() << " ==> " << (total_range > text.size());

                if (total_range > text.size())
                {
                    // linear scan
                    auto rx = boost::regex(pat.raw_regexp.begin(),pat.raw_regexp.end(),std::regex::ECMAScript);
                    auto matches_begin = boost::sregex_iterator(
                        text.begin(),
                        text.end(),
                        rx,
                        boost::regex_constants::match_flag_type::match_not_dot_newline);
                    auto matches_end = boost::sregex_iterator();

                    for (boost::sregex_iterator it = matches_begin; it != matches_end; ++it) {
                        res.positions.push_back(it->position());
                    }
                    return res;
                }
            }

            // smart scan
            sdsl::PERFCTR_NUM_PROCESSED_WT_NODES = 1;
            for (auto hit : index.match(s1, s2, min_gap, max_gap)) {
                res.positions.push_back(hit.first);
            }
            cout << "Processed nodes: " << sdsl::PERFCTR_NUM_PROCESSED_WT_NODES << endl;
            return res;
        }
};
