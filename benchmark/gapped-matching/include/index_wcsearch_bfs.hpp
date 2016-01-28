#pragma once

#include "utils.hpp"
#include "collection.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/matching.hpp"

class index_wcsearch_bfs
{
    private:
        typedef sdsl::matching_index<sdsl::wt_int<>, sdsl::rrr_vector<>> index_type;
        index_type index;

    public:
        typedef sdsl::int_vector<0>::size_type size_type;
        std::string name() const
        {
            std::string index_name = IDXNAME;
            return "WCSEARCH-BFS-"+index_name;
        }

    public:
        index_wcsearch_bfs() { }
        index_wcsearch_bfs(collection& col)
        {
            sdsl::cache_config cc(false,".","WCSEARCH_BFS_TMP");
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

        void swap(index_wcsearch_bfs& ir)
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
            const auto& wts = index.wt;
            using node_type = remove_reference<decltype(wts)>::type::node_type;
            using size_type = remove_reference<decltype(wts)>::type::size_type;
            using range_type = sdsl::range_type;

            gapped_search_result res;
            string_type s1;
            string_type s2;
            size_t min_gap = 0;
            size_t max_gap = 0;

            if (pat.subpatterns.size() == 1) {
                s1 = pat.subpatterns[0];
                s2 = pat.subpatterns[0];
            } else {
                s1 = pat.subpatterns[0];
                s2 = pat.subpatterns[1];
                min_gap = pat.gaps[0].first;
                max_gap = pat.gaps[0].second;
            }

            size_t cnt_nodes = 0;
            array<vector<pair<range_type,size_t>>, 2> lex_ranges;
            vector<node_type> nodes;
            nodes.push_back(wts.root());
            size_type sp = 1, ep = 0;
            if (forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, s1.begin(), s1.end(), sp, ep) > 0)
                lex_ranges[0].emplace_back(range_type(sp, ep),0);
            else
                nodes.clear();
            if (forward_search(index.text.begin(), index.text.end(), index.wt, 0, index.wt.size()-1, s2.begin(), s2.end(), sp, ep) > 0)
                lex_ranges[1].emplace_back(range_type(sp,ep),0);
            else
                nodes.clear();

            cnt_nodes += nodes.size();

            auto _lb = [&lex_ranges, &nodes, &wts](size_type k, size_type i) {
                return get<0>(wts.value_range(nodes[lex_ranges[k][i].second]));
            };

            auto _rb = [&lex_ranges, &nodes, &wts](size_type k, size_type i) {
                return get<1>(wts.value_range(nodes[lex_ranges[k][i].second]));
            };

            while (!nodes.empty()) {
                if (wts.is_leaf(nodes[0]))
                    break;
                vector<node_type> _nodes;
                array<vector<pair<range_type,size_t>>, 2> _lex_ranges;
                std::array<decltype(lex_ranges[0].begin()),2> lex_range_it = {lex_ranges[0].begin(), lex_ranges[1].begin()};
                for (size_t i=0; i<nodes.size(); ++i) {
                    auto exp_node = wts.expand(nodes[i]);
                    _nodes.push_back(exp_node[0]);
                    _nodes.push_back(exp_node[1]);
                    for (size_t j=0; j<2; ++j) {
                        if (lex_range_it[j] != lex_ranges[j].end() and lex_range_it[j]->second == i) {
                            auto exp_range = wts.expand(nodes[i], lex_range_it[j]->first);
                            ++lex_range_it[j];
                            for (size_t k=0; k<2; ++k) {
                                if (!sdsl::empty(exp_range[k])) {
                                    _lex_ranges[j].emplace_back(exp_range[k], _nodes.size()-2+k);
                                }
                            }
                        }
                    }
                }
                cnt_nodes += _nodes.size();
                nodes = std::move(_nodes);
                lex_ranges = std::move(_lex_ranges);

                // filtering
                size_t n0 = 0, n1 = 0;
                for (size_t i=0,j=0; i<lex_ranges[0].size(); ++i) {
                    while (j < lex_ranges[1].size() and _rb(1,j) < _lb(0,i)+s1.size()+min_gap) {
                        ++j;
                    }
                    if (j < lex_ranges[1].size() and _lb(1,j) <= _rb(0,i)+s1.size()+max_gap) {
                        // lex_ranges[0][i] overlaps lex_ranges[1][j]
                        lex_ranges[0][n0++] = lex_ranges[0][i];
                    }
                }
                lex_ranges[0].resize(n0);
                for (size_t i=0,j=0; j<lex_ranges[1].size(); ++j) {
                    while (i < lex_ranges[0].size() and _lb(0,n0-i-1)+s1.size()+min_gap > _rb(1,lex_ranges[1].size()-j-1)) {
                        ++i;
                    }
                    if (i < lex_ranges[0].size() and _lb(1,lex_ranges[1].size()-j-1) <= _rb(0,n0-i-1)+s1.size()+max_gap) {
                        lex_ranges[1][lex_ranges[1].size()-n1-1] = lex_ranges[1][lex_ranges[1].size()-j-1];
                        ++n1;
                    }
                }
                for (size_t j=0; j<n1; ++j) {
                    lex_ranges[1][j] = lex_ranges[1][lex_ranges[1].size()-n1+j];
                }
                lex_ranges[1].resize(n1);
                // remove unused nodes
                n0 = 0;
                for (size_t k=0, i=0, j=0; k<nodes.size(); ++k) {
                    bool occurs = false;
                    if (i < lex_ranges[0].size() and lex_ranges[0][i].second == k) {
                        occurs = true;
                        lex_ranges[0][i++].second = n0;
                    }
                    if (j < lex_ranges[1].size() and lex_ranges[1][j].second == k) {
                        occurs = true;
                        lex_ranges[1][j++].second = n0;
                    }
                    if (occurs) {
                        nodes[n0++] = nodes[k];
                    }
                }
                nodes.resize(n0);
            }

            // GREEDY
            /*            min_gap += s1.size();
                        max_gap += s1.size();

                        vector<size_type> range_a;
                        vector<size_type> range_b;
                        for (size_t i = 0; i < lex_ranges[0].size(); ++i)
                            range_a.push_back(get<0>(wts.value_range(nodes[lex_ranges[0][i].second])));
                        for (size_t i = 0; i < lex_ranges[1].size(); ++i)
                            range_b.push_back(get<0>(wts.value_range(nodes[lex_ranges[1][i].second])));

                        std::sort(range_a.begin(), range_a.end());
                        std::sort(range_b.begin(), range_b.end());

                        // linear search
                        auto a_it = range_a.begin();
                        auto b_it = range_b.begin();

                        while (a_it != range_a.end()) {
                            auto a_pos = *a_it;

                            // enforcing min_gap
                            bool b_valid;
                            while ((b_valid = (b_it != range_b.end())) && a_pos + min_gap > *b_it)
                                ++b_it;
                            if (!b_valid)
                                break;

                            // check whether within max_gap
                            auto b_pos = *b_it;
                            if (a_pos + max_gap < b_pos) {
                                ++a_it;
                                continue;
                            }

                            // push greedy beyond max_gap
                            ++b_it;
                            while (b_it != range_b.end()) {
                                auto b_pos2 = *b_it;
                                if (a_pos + max_gap >= b_pos2)
                                    b_pos = b_pos2;
                                else
                                    break;
                                ++b_it;
                            }

                            res.positions.push_back(a_pos);

                            // pull a beyond previous b (non-overlapping)
                            b_pos += s2.size();
                            while (a_it != range_a.end() && *a_it < b_pos)
                                ++a_it;
                        }
            */
            // ALL
            for (size_t i = 0; i < lex_ranges[0].size(); ++i) {
                auto pos0 = _lb(0,i);
                for (size_t j = 0; j < lex_ranges[1].size(); ++j) {
                    auto pos1 = _lb(1,j);
                    if (pos0 + s1.size() + min_gap <= pos1
                        && pos0 + s1.size() + max_gap >= pos1)
                        res.positions.push_back(pos0);
                }
            }
            cout << "Processed nodes: " << cnt_nodes << " (" <<wts.size() << ")" << endl;

            return res;
        }
};
