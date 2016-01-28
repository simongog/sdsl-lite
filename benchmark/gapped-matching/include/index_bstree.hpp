#pragma once

#include "utils.hpp"
#include "collection.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/matching.hpp"

class index_bstree
{
    private:
        typedef std::pair<size_t, size_t> pii;
        sdsl::int_vector<> m_text;
        sdsl::int_vector<> m_sa;

    public:
        typedef sdsl::int_vector<0>::size_type size_type;
        std::string name() const
        {
            std::string index_name = IDXNAME;
            return "BSTREE-"+index_name;
        }

        index_bstree() { }
        index_bstree(collection& col)
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

        void swap(index_bstree& ir)
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

            auto last_size = s[s.size() - 1].size();

            // get ranges
            vector<sdsl::int_vector_const_iterator<sdsl::int_vector<0>>::const_iterator> its;
            vector<sdsl::int_vector_const_iterator<sdsl::int_vector<0>>::const_iterator> its2;

            for (auto sx : s) {
                size_type sp, ep;
                forward_search(m_text.begin(), m_text.end(), m_sa, 0, m_sa.size()-1, sx.begin(), sx.end(), sp, ep);
                its.push_back(m_sa.begin() + sp);
                its2.push_back(m_sa.begin() + ep + 1);
            }

            vector<pii> spans1, spans2;
            for (auto it = its[0]; it != its2[0]; ++it)
		spans1.emplace_back(*it, *it);

            auto pspan1 = &spans1;
            auto pspan2 = &spans2;

            // incremental search
            for (size_t i = 1; i < its.size(); ++i) {
                std::sort(pspan1->begin(), pspan1->end());
                pspan2->clear();
                for (auto it = its[i]; it != its2[i]; ++it) {
                    auto pos = *it;
                    for (auto pot_match = lower_bound(pspan1->begin(), pspan1->end(), make_pair(pos - max_gap, pos));
                         pot_match != pspan1->end() && pot_match->first <= pos - min_gap;
                         ++pot_match)
                        pspan2->emplace_back(pos, pot_match->second);
                }
                std::swap(pspan1, pspan2);
            }
            std::sort(pspan1->begin(), pspan1->end(), [](const pii &left, const pii &right) { return left.second < right.second; });
            
            for (auto it = pspan1->begin(); it != pspan1->end(); ++it) {
                res.positions.push_back(it->second);
                auto end = it->first + last_size;
                while (it != pspan1->end() && it->second < end)
                    ++it;
            }

            return res;
        }
};
