/* sdsl - succinct data structures library
    Copyright (C) 2014 Simon Gog

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
/*! \file hybrid_k2_tree.hpp
    \brief hybrid_k2_tree.hpp contains a compact hybrid k^2-tree.
    \author Jan Broß, based on the k2 treap code of Simon Gog, leaf compression is based on the libk2tree implemenetation which uses DACs impelemented by Brisaboa, Ladra et.al.
*/
#ifndef INCLUDED_SDSL_HYBRID_K2_TREE
#define INCLUDED_SDSL_HYBRID_K2_TREE

#include "k2_tree_base.hpp"
#include <gtest/gtest_prod.h>

//! Namespace for the succinct data structure library.
namespace sdsl {

/*! A hybrid k2 tree implementation
 *  \par References
 *       [1] Nieves R. Brisaboa, Susana Ladra, and Gonzalo Navarro:
 *           Compact representation of Web graphs with extended functionality.
 *           Inf. Syst. 39 (January 2014), 152-174. DOI=http://dx.doi.org/10.1016/j.is.2013.08.003
 */
    template<uint8_t t_k_l_1,
            uint8_t t_k_l_1_size,
            uint8_t t_k_l_2,
            uint8_t t_k_leaves,
            typename t_lev=bit_vector,
            typename t_leaf=bit_vector,
            typename t_rank=typename t_lev::rank_1_type>

    class k2_tree_hybrid : public k2_tree_base<t_k_l_1, t_lev, t_leaf, t_rank> {
        static_assert(t_k_l_1 == 2    || t_k_l_1 == 4    || t_k_l_1 == 8    || t_k_l_1 == 16,    "t_k_l_1 has to one of 2,4,8,16");
        static_assert(t_k_l_2 == 2    || t_k_l_2 == 4    || t_k_l_2 == 8    || t_k_l_2 == 16,    "t_k_l_2 has to one of 2,4,8,16");
        static_assert(t_k_leaves == 2 || t_k_leaves == 4 || t_k_leaves == 8 || t_k_leaves == 16, "t_k_leaves has to one of 2,4,8,16");
        static_assert(t_k_leaves >= t_k_l_1,
                      "t_k_leaves has to be larger than t_k_l_1,  otherwise this could lead to different word sizes and thus to a problem for the k2part approach"); //if smaller than t_k_l_1 it could be that t_k_leaves is not used
        static_assert(t_k_leaves <= 8, "t_k can at most be 8 because of the current dac compression implementation.");

    private:

        std::vector<uint8_t> m_k_for_level;
        std::vector<uint8_t> m_shift_table; //for fast div and mod operations by submatrix size

    public:
        typedef stxxl::VECTOR_GENERATOR<std::pair<uint32_t, uint32_t>>::result stxxl_32bit_pair_vector;
        typedef stxxl::VECTOR_GENERATOR<std::pair<uint64_t, uint64_t>>::result stxxl_64bit_pair_vector;
        typedef int_vector<>::size_type size_type;

        k2_tree_hybrid() = default;

        k2_tree_hybrid(const k2_tree_hybrid &tr)
                : k2_tree_base<t_k_l_1, t_lev, t_leaf, t_rank>(tr) {
            *this = tr;
        }

        k2_tree_hybrid(k2_tree_hybrid &&tr)
                : k2_tree_base<t_k_l_1, t_lev, t_leaf, t_rank>(tr) {
            *this = std::move(tr);
        }

        template<typename t_vector>
        k2_tree_hybrid(std::string temp_file_prefix, bool use_counting_sort, t_vector &v, uint64_t max_hint = 0) {

            using namespace k2_treap_ns;
            if (v.size() > 0) {
                if (max_hint == 0) {
                    max_hint = get_maximum(v);
                }
                this->m_tree_height = get_tree_height(max_hint);
                initialize_shift_table();

                if (use_counting_sort) {
                    construct_counting_sort(v, temp_file_prefix);
                    //construct_bottom_up(v, temp_file_prefix);
                } else {
                    this->construct(v, temp_file_prefix);
                }

                this->post_init();
            }
        }

        k2_tree_hybrid(int_vector_buffer<> &buf_x,
                       int_vector_buffer<> &buf_y, bool use_counting_sort = false, uint64_t max_hint = 0)  {
            using namespace k2_treap_ns;
            typedef int_vector_buffer<> *t_buf_p;

            if (buf_x.size() == 0){
                return;
            }

            std::vector<t_buf_p> bufs = {&buf_x, &buf_y};

            uint64_t max;
            if (max_hint != 0) {
                max = max_hint;//temporarily set
            } else {
                auto max_element = [](int_vector_buffer<> &buf) {
                    uint64_t max_val = 0;
                    for (auto val : buf) {
                        max_val = std::max((uint64_t) val, max_val);
                    }
                    return max_val;
                };

                auto max_buf_element = [&]() {
                    uint64_t max_v = 0;
                    for (auto buf : bufs) {
                        uint64_t _max_v = max_element(*buf);
                        max_v = std::max(max_v, _max_v);
                    }
                    return max_v;
                };

                max = max_buf_element();
            };

            this->m_tree_height = get_tree_height(max);
            initialize_shift_table();

            if (this->m_max_element <= std::numeric_limits<uint32_t>::max()) {
                auto v = read<uint32_t, uint32_t>(
                        bufs);
                if (use_counting_sort) {
                    construct_counting_sort(v, buf_x.filename());
                } else {
                    this->construct(v, buf_x.filename());
                }

            } else {
                auto v = read<uint64_t, uint64_t>(
                        bufs);
                if (use_counting_sort) {
                    construct_counting_sort(v, buf_x.filename());
                } else {
                    this->construct(v, buf_x.filename());
                }
            }

            this->post_init();
        }

        size_type serialize(std::ostream &out, structure_tree_node *v, std::string name) const override {
            return k2_tree_base<t_k_l_1, t_lev, t_leaf, t_rank>::serialize(out, v, name);
        }

        uint word_size() const override {
            return div_ceil((uint) t_k_leaves * t_k_leaves, kUcharBits);
        }


        size_t words_count() const override {
            return this->m_leaves.size() / t_k_leaves / t_k_leaves;
        }

        void load(std::istream &in) override {
            k2_tree_base<t_k_l_1, t_lev, t_leaf, t_rank>::load(in);
            if (this->m_tree_height > 0) {
                for (int i = 1; i <= std::min(t_k_l_1_size, (uint8_t) (this->m_tree_height - 1)); ++i) {
                    m_k_for_level.push_back(t_k_l_1);
                }

                for (int j = t_k_l_1_size + 1; j <= (this->m_tree_height - 1); ++j) {
                    m_k_for_level.push_back(t_k_l_2);
                }

                m_k_for_level.push_back(t_k_leaves);

                initialize_shift_table();

                if (this->m_access_shortcut_size > 0) {
                    this->perform_access_shortcut_precomputations();
                }
            }
        }

        //hack a the moment, because construct cannot be virtual
        void load_from_ladrabin(std::string fileName, bool use_counting_sort = false, uint8_t access_shortcut_size = 0, std::string temp_file_prefix = "") {
            using namespace k2_treap_ns;
            if (!has_ending(fileName, ".ladrabin")) {
                fileName.append(".ladrabin");
                std::cout << "Appending .graph-txt to filename as file has to be in .ladrabin format" << std::endl;
            }

            std::fstream fileStream(fileName, std::ios_base::in);

            if (fileStream.is_open()) {
                uint number_of_nodes;
                ulong number_of_edges;

                read_member(number_of_nodes, fileStream);
                read_member(number_of_edges, fileStream);

                this->m_max_element = number_of_nodes - 1;
                this->m_size = number_of_edges;
                this->m_tree_height = get_tree_height(this->m_max_element);
                initialize_shift_table();

                uint nodes_read = 0;
                uint source_id;
                int target_id;

                std::vector<std::pair<uint, uint>> coords;
                coords.reserve(number_of_edges);
                /*typedef stxxl::VECTOR_GENERATOR<pair<uint32_t, uint32_t>>::result stxxl_pair_vector;
                stxxl_pair_vector coords(number_of_nodes);*/
                for (uint64_t i = 0; i < number_of_nodes + number_of_edges; i++) {
                    read_member(target_id, fileStream);
                    if (target_id < 0) {
                        nodes_read++;
                    } else {
                        source_id = nodes_read - 1;
                        coords.push_back(std::make_pair(source_id, target_id));
                    }
                }
                fileStream.close();

                if (coords.size() > 0) {
                    if (use_counting_sort) {
                        construct_counting_sort(coords, temp_file_prefix);
                        //construct_bottom_up(v, temp_file_prefix);
                    } else {
                        this->construct(coords, temp_file_prefix);
                    }

                    std::cout << "Finished Construction" << std::endl;

                    coords.clear();

                    this->m_access_shortcut_size = access_shortcut_size;
                    this->post_init();
                }
            } else {
                throw std::runtime_error("Could not load ladrabin file");
            }
        }

        k2_tree_hybrid &operator=(k2_tree_hybrid &&tr) {
            k2_tree_base<t_k_l_1, t_lev, t_leaf, t_rank>::operator=(tr);
            if (this != &tr) {
                m_k_for_level = tr.m_k_for_level;
                m_shift_table = tr.m_shift_table;
            }
            return *this;
        }

        k2_tree_hybrid &operator=(const k2_tree_hybrid &tr) {
            k2_tree_base<t_k_l_1, t_lev, t_leaf, t_rank>::operator=(tr);
            if (this != &tr) {
                m_k_for_level = tr.m_k_for_level;
                m_shift_table = tr.m_shift_table;
            }
            return *this;
        }

        bool operator==(const k2_tree_hybrid &tr) const {
            if (!k2_tree_base<t_k_l_1, t_lev, t_leaf, t_rank>::operator==(tr)) {
                return false;
            }

            if (m_k_for_level.size() !=
                tr.m_k_for_level.size()) {//must be the same for the same template parameters and data
                return false;
            }

            if (m_shift_table.size() !=
                tr.m_shift_table.size()) {//must be the same for the same template parameters and data
                return false;
            }

            return true;
        }

        void swap(k2_tree_hybrid &tr) {
            k2_tree_base<t_k_l_1, t_lev, t_leaf, t_rank>::swap(tr);
            std::swap(m_k_for_level, tr.m_k_for_level);
            std::swap(m_shift_table, tr.m_shift_table);
        }

        virtual void construct_access_shortcut(uint8_t access_shortcut_size) override {
            if (this->m_access_shortcut_size <= t_k_l_1_size){
                k2_tree_base<t_k_l_1, t_lev, t_leaf, t_rank>::construct_access_shortcut(access_shortcut_size);
            } else {
                throw std::runtime_error("access shortcut is only supported over levels with the same k values, therefore t_k_l_1_size must be >= access_shortcut_size");
            }
        }

        std::string get_type_string_without_compression() const {
            return "k2_tree_hybrid<"+std::to_string(t_k_l_1)+","+std::to_string(t_k_l_1_size)+","+std::to_string(t_k_l_2)+","+std::to_string(t_k_leaves)+">";
        }

        std::string get_type_string() const {
            return "k2_tree_hybrid<"+std::to_string(t_k_l_1)+","+std::to_string(t_k_l_1_size)+","+std::to_string(t_k_l_2)+","+std::to_string(t_k_leaves)+","+get_compression_name(this->m_used_compression)+">";
        }

    private:

        /**
        * Constructs the tree corresponding to the points in the links vector by using counting
        * sort with k² buckets and rearranging the input on every level of the tree
        * @param links
        * @param temp_file_prefix
        */
        template<typename t_vector>
        void construct_counting_sort(t_vector &links, std::string temp_file_prefix = "") {
            using namespace k2_treap_ns;
            typedef decltype(links[0].first) t_x;
            typedef decltype(links[0].second) t_y;
            using t_e = std::pair<t_x, t_y>;

            this->m_size = links.size();

            if (this->m_size == 0) {
                return;
            }

            //std::cout << "Using counting sort" << std::endl;

            std::string id_part = util::to_string(util::pid())
                                  + "_" + util::to_string(util::id());

            {
                //                  upper left          lower right                 interval in links   level
                typedef std::tuple<std::pair<t_x, t_y>, std::pair<uint64_t, uint64_t>, t_e, uint8_t> t_queue;
                std::queue<t_queue> queue;

                //partition recursively until reaching the leaves
                queue.push(t_queue(std::make_pair<t_x, t_y>(0, 0),
                                   std::make_pair(this->m_max_element - 1, this->m_max_element - 1),
                                   t_e(0, links.size()), 0));

                uint64_t number_of_bits = 0; //for speed comparison purposes of different k
                uint8_t previous_level = 0;
                uint64_t ctr = 0;
                uint64_t set_bits_in_level = 0;
                bit_vector buffer(get_k(0) * get_k(0));
                this->m_levels.resize(this->m_tree_height - 1);
                while (!queue.empty()) {
                    auto upper_left = std::get<0>(queue.front());
                    auto lower_right = std::get<1>(queue.front());
                    auto links_interval = std::get<2>(queue.front());
                    auto current_level = std::get<3>(queue.front());
                    auto inv_level = this->m_tree_height - current_level - 1; //inverse level for accessing precomps

                    const uint8_t k = get_k(current_level);
                    if (current_level > previous_level) {
                        previous_level = current_level;
                        this->m_levels[current_level - 1] = buffer;
                        bit_vector tmp(set_bits_in_level * k * k);
                        tmp.swap(buffer);
                        ctr = 0;
                        set_bits_in_level = 0;
                    }

                    std::vector<t_x> intervals(k * k + 1);

                    queue.pop();

                    //do counting sort
                    auto x1 = upper_left.first;
                    auto y1 = upper_left.second;
                    for (uint64_t j = links_interval.first; j < links_interval.second; ++j) {
                        auto x = links[j].first;
                        auto y = links[j].second;
                        auto p1 = divexp(x - x1, inv_level);
                        auto p2 = divexp(y - y1, inv_level);

                        //FIXME: replace by bitshift for k
                        auto corresponding_matrix = p1 * k + p2;
                        intervals[corresponding_matrix +
                                  1]++;//offset corresponding matrix by one to allow for more efficient in interval comparision
                    }

                    intervals[0] = 0;

                    //append bits to level_vectors[level] based on result
                    for (uint i = 1; i < intervals.size(); ++i) {
//                        item_count[current_level].push_back(intervals[i])
//                        utilization[current_level].push_back(
//                                ((double) intervals[i]) / (submatrix_size * submatrix_size) * 100);

                        if (intervals[i] > 0) {
                            buffer[ctr] = 1;
                            //std::cout << "1";
                            set_bits_in_level++;
                        } /*else {
                            buffer[ctr] = 0;
                        }*/
                        number_of_bits++;
                        ctr++;
                        intervals[i] += intervals[i - 1]; //build prefix sum
                    }

                    //leaves not reached yet --> append to level_vector & reorder
                    if (inv_level > 0) {
                        std::vector<t_x> offset(k * k);
                        offset[0] = 0;
                        for (size_t l = 1; l < offset.size(); ++l) {
                            offset[l] = intervals[l] + 1;
                        }

                        auto begin = links.begin() + links_interval.first;
                        auto it = begin;

                        //reorder links based on counting sort offsets
                        uint64_t index = 0;
                        while (it != links.begin() + links_interval.second) {
                            uint corresponding_matrix =
                                    divexp(((*it).first - x1), inv_level) * k +
                                    divexp((*it).second - y1, inv_level);

                            if (index >= intervals[corresponding_matrix] &&
                                index < intervals[corresponding_matrix + 1]) {
                                //element is at correct position
                                offset[corresponding_matrix]++;

                                it++;
                                index++;
                            } else {
                                //swap to correct position either swapping a match or a non match back
                                //there are at most m swaps if all m elements are at the wrong position
                                //every value is at most read twice
                                std::iter_swap(it, begin + offset[corresponding_matrix] - 1);
                                offset[corresponding_matrix]++;
                            }
                        }

                        //enqueue submatrixes
                        auto new_submatrix_size = exp(inv_level);
                        for (uint x = 0; x < k; ++x) {
                            for (uint y = 0; y < k; ++y) {
                                auto new_interval = std::make_pair(intervals[x * k + y] + links_interval.first,
                                                                   intervals[x * k + y + 1] + links_interval.first);
                                if (new_interval.first != new_interval.second) {
                                    auto new_upper_left = std::make_pair<t_x, t_y>(
                                            x * new_submatrix_size + upper_left.first,
                                            y * new_submatrix_size + upper_left.second);
                                    auto new_lower_right = std::make_pair<t_x, t_y>(
                                            (x + 1) * new_submatrix_size - 1 + upper_left.first,
                                            (y + 1) * new_submatrix_size - 1 + upper_left.second);
                                    queue.push(std::make_tuple(new_upper_left, new_lower_right, new_interval,
                                                               current_level + 1));
                                }
                            }
                        }
                    }
                }

                this->m_leaves = t_leaf(buffer);

                this->m_levels_rank.resize(this->m_levels.size());
                for (uint64_t i = 0; i < this->m_levels.size(); ++i) {
                    util::init_support(this->m_levels_rank[i], &this->m_levels[i]);
                }
            }
        }

        /**
        * Constructs the tree corresponding to the points in the links vector by partitioning the input multiple times
        * @param links
        * @param temp_file_prefix
        */
        template<typename t_vector>
        void construct(t_vector &links, std::string temp_file_prefix = "") {
            using namespace k2_treap_ns;
            typedef decltype(links[0].first) t_x;
            typedef decltype(links[0].second) t_y;
            using t_e = std::pair<t_x, t_y>;

            this->m_size = links.size();

            if (this->m_tree_height == 0) {//might occur in k2part
                return;
            }

            std::string id_part = util::to_string(util::pid())
                                  + "_" + util::to_string(util::id());

            {
                std::vector<int_vector_buffer<1>> level_buffers = this->create_level_buffers(temp_file_prefix, id_part);

                auto end = std::end(links);
                //uint64_t last_level_bits = 0;
                //uint64_t level_bits = 0;

                //recursively partition that stuff
                uint64_t submatrix_size = this->m_max_element;
                for (int l = this->m_tree_height; l + 1 > 0; --l) {

                    //std::cout << "Processing Level " << l << std::endl;
                    //level_bits = 0;

                    uint8_t k = 0;
                    uint64_t current_submatrix_size = 0;
                    if (l > 0) {
                        k = get_k(this->m_tree_height - l);
                            current_submatrix_size = submatrix_size / k;
                    }

                    auto sp = std::begin(links);
                    for (auto ep = sp; ep != end;) {

                        //Iterator which only returns the nodes within a certain subtree
                        ep = std::find_if(sp, end, [=, &sp, &l](const t_e &e) {
                            auto x1 = std::get<0>(*sp);
                            auto y1 = std::get<1>(*sp);
                            auto x2 = std::get<0>(e);
                            auto y2 = std::get<1>(e);
                            bool in_sub_tree = divexp(x1, l) != divexp(x2, l)
                                               or divexp(y1, l) != divexp(y2, l);

                            return in_sub_tree;
                        });


                        if (l > 0) {
                            auto _sp = sp;

                            for (uint8_t i = 0; i < k; ++i) {
                                auto _ep = ep;
                                if (i + 1 < k) {  //partition t_k -1 times vertically (1 in the case of k=2)
                                    _ep = std::partition(_sp, _ep, [=, &i, &l](const t_e &e) {
                                        return (divexp(std::get<0>(e), l - 1) & (k-1))<= i;
                                    });
                                }
                                auto __sp = _sp;

                                for (uint8_t j = 0;
                                     j < k; ++j) { //partition the t_k vertical partitions t_k -1 times horizontally
                                    auto __ep = _ep;
                                    if (j + 1 < k) {
                                        __ep = std::partition(__sp, _ep, [=, &j, &l](const t_e &e) {
                                            return (divexp(std::get<1>(e), l - 1) & (k-1)) <= j;
                                        });
                                    }
                                    bool not_empty = __ep > __sp;
                                    level_buffers[this->m_tree_height - l].push_back(not_empty);
                                    //level_bits++;
                                    __sp = __ep;
                                }
                                _sp = _ep;
                            }
                        }
                        if (ep != end) {
                            //++ep;
                            sp = ep;
                        }
                    }
                    //last_level_bits = level_bits;
                    submatrix_size = current_submatrix_size;
                }
            }
            this->load_vectors_from_file(temp_file_prefix, id_part);

            /*
            std::cout << "Fallen Levels: " << std::endl;
            for (uint i = 0; i < this->m_levels.size(); i++){
                std::cout << "Level " << i << "\n";
                for (uint j = 0; j < this->m_levels[i].size(); j++){
                    std::cout << this->m_levels[i][j];
                }
                std::cout << "\n";
            }
            std::cout << std::endl;


            std::cout << "Fallen Leaves: " << std::endl;
            for (uint i = 0; i < this->m_leaves.size(); i++){
                std::cout << this->m_leaves[i];
            }
            std::cout << std::endl;
            */
            //std::cout << "Leaves size" << this->m_leaves.size() << std::endl;
        }

        inline uint64_t exp(uint8_t l) {
            assert(l >= 0 && l <= m_tree_height);
            return 1Ull << m_shift_table[l];
        }

        inline uint64_t divexp(uint64_t x, uint8_t l) {
            assert(l >= 0 && l <= m_tree_height);
            return x >> m_shift_table[l];
        }

        inline uint64_t modexp(uint64_t x, uint8_t l) {
            assert(l >= 0 && l <= m_tree_height);
            return x & bits::lo_set[m_shift_table[l]];
        }

        void initialize_shift_table() {
            m_shift_table.resize(this->m_tree_height + 1);

            m_shift_table[0] = 0;
            for (uint i = 0; i < m_k_for_level.size(); ++i) {
                m_shift_table[i + 1] = m_shift_table[i] + bits::hi(m_k_for_level[this->m_tree_height-i-1]);
            }
        }

        inline uint8_t get_k(uint8_t level) const {
            return m_k_for_level[level];
        }

        uint8_t get_tree_height(const uint64_t max) {
            uint8_t res = 1;
            if (t_k_l_1 <= max) {
                this->m_max_element = t_k_l_1;
                m_k_for_level.push_back(t_k_l_1);
            } else {
                //in case only one level use k_leaves (has to be higher than t_k_l_1 (compile time-checked)
                this->m_max_element = t_k_leaves;
                m_k_for_level.push_back(t_k_leaves);
            }

            while (this->m_max_element <= max) {
                if ((this->m_max_element * t_k_leaves) > max) {
                    this->m_max_element = this->m_max_element * t_k_leaves;
                    m_k_for_level.push_back(t_k_leaves);
                } else if (res < t_k_l_1_size) {
                    this->m_max_element = this->m_max_element * t_k_l_1;
                    m_k_for_level.push_back(t_k_l_1);
                } else {
                    this->m_max_element = this->m_max_element * t_k_l_2;
                    m_k_for_level.push_back(t_k_l_2);
                }
                ++res;
            }

            m_k_for_level.resize(res);
            if (res <= t_k_l_1_size) {

                std::cerr
                        << "The tree height is smaller than t_k_l_1_size using t_k_l_1 up to m_tree_height-1 and then t_k_leaves"
                        << std::endl;
            } else if (res == t_k_l_1_size + 1) {
                std::cerr << "The tree equals t_k_l_1_size+1, only using t_k_l_1 and t_k_leaves" << std::endl;
            }

            return res;
        }
    };
}
#endif

