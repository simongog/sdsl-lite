#pragma once

#include "utils.hpp"
#include "collection.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/matching.hpp"

class index_sasearch
{
    private:
        sdsl::int_vector<> m_text;
        sdsl::int_vector<> m_sa;

    public:
        typedef sdsl::int_vector<0>::size_type size_type;
        std::string name() const
        {
            std::string index_name = IDXNAME;
            return "SASEARCH-"+index_name;
        }

    public:
        index_sasearch() { }
        index_sasearch(collection& col)
        {
            sdsl::load_from_file(m_text, col.file_map[consts::KEY_TEXT]);

            sdsl::csa_wt<sdsl::wt_huff<>, 1> csa;
            sdsl::construct(csa, col.file_map[consts::KEY_TEXT], 0);
            m_sa = sdsl::int_vector<>(csa.size(), 0, sdsl::bits::hi(csa.size() - 1) + 1);
            std::copy(csa.begin(), csa.end(), m_sa.begin());
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=NULL, std::string name="")const
        {
            (void)name;
            auto size = m_sa.serialize(out, v, "sa");
            size += m_text.serialize(out, v, "text");
            return size;
        }

        void load(std::istream& in)
        {
            m_sa.load(in);
            m_text.load(in);
        }

        void swap(index_sasearch& ir)
        {
            if (this != &ir) {
                m_sa.swap(ir.m_sa);
                m_text.swap(ir.m_text);
            }
        }

        std::string info(const gapped_pattern& pat) const { (void)pat; return ""; }
        void prepare(const gapped_pattern& pat) { (void)pat; }

        gapped_search_result
        search(const gapped_pattern& pat) const
        {
            gapped_search_result res;
            vector<string_type> s;
            size_type min_gap;
            size_type max_gap;

            std::cerr << "REGEX ::: " << pat.raw_regexp << std::endl;

            s.push_back(pat.subpatterns[0]);
            s.push_back(pat.subpatterns[1]);
            for (size_t i = 2; i < NUM_PATTERNS; ++i)
                s.push_back(pat.subpatterns[1]);

            min_gap = s[0].size() + pat.gaps[0].first;
            max_gap = s[0].size() + pat.gaps[0].second;

            // get ranges
            vector<vector<size_type>> ranges;
            vector<vector<size_type>::iterator> its;

            for (auto sx : s) {
                size_type sp, ep;
                ranges.emplace_back(forward_search(m_text.begin(), m_text.end(), m_sa, 0, m_sa.size()-1, sx.begin(), sx.end(), sp, ep));
                vector<size_type>& range = ranges[ranges.size() - 1];
                std::copy(m_sa.begin() + sp, m_sa.begin() + ep + 1, range.begin());
                std::sort(range.begin(), range.end());
                its.push_back(range.begin());
            }


            // linear search

            while (its[0] != ranges[0].end()) {
                auto prev_pos = *its[0];
                bool b_break = false;
                bool b_continue = false;
                for (size_t i = 1; i < ranges.size(); ++i) {
                    // enforcing min_gap ab
                    bool valid;
                    while ((valid = (its[i] != ranges[i].end())) && prev_pos + min_gap > *its[i]) ++its[i];
                    if (!valid) { b_break = true; break; }

                    // check whether within max_gap ab
                    auto new_pos = *its[i];
                    if (prev_pos + max_gap < new_pos) { ++its[i-1]; b_continue = true; break; }
                    prev_pos = new_pos;
                }
                if (b_break) break;
                if (b_continue) continue;

                // situation: VALID match, but LAZY

                res.positions.push_back(*its[0]);
                //std::cerr << "HIT" << std::endl;
                //for (size_t i = 0; i < its.size(); ++i)
                //    std::cerr << "\t" << *its[i] << std::endl;

                // pull a beyond previous c (non-overlapping)
                auto posx = *its[its.size() - 1] + s[0].size();
                while (its[0] != ranges[0].end() && *(its[0]) < posx)
                    ++its[0];
            }
            return res;
        }
};
