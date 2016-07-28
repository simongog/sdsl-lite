#ifndef INCLUDED_SDSL_K2_TREE_PARTITONED
#define INCLUDED_SDSL_K2_TREE_PARTITONED

#include "vectors.hpp"
#include "bits.hpp"
#include <tuple>
#include <algorithm>
#include <climits>
#include <vector>
#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/rank_support_v.hpp>
#include <gtest/gtest_prod.h>
//#include <stxxl/vector>

namespace sdsl {

    template<uint8_t t_k0,
            typename subk2_tree>
    class k2_tree_partitioned {
        static_assert(t_k0 > 1, "t_k has to be larger than 1.");
        static_assert(t_k0 <= 16, "t_k has to be smaller than 17.");

    public:

        typedef int_vector<>::size_type size_type;

        enum {
            k0 = t_k0
        };

        size_type m_size = 0;
        uint64_t m_matrix_size;

    private:
        //typedef k2_tree<t_k,t_bv,t_rank> k2;
        std::vector<subk2_tree> m_k2trees;

    public:
        k2_tree_partitioned() = default;

        k2_tree_partitioned(const k2_tree_partitioned &tr) {
            *this = tr;
        }

        k2_tree_partitioned(k2_tree_partitioned &&tr) {
            *this = std::move(tr);
        }

        template<typename t_vector>
        k2_tree_partitioned(std::string temp_file_prefix, bool bottom_up,
                            uint access_shortcut_size, t_vector &v) {
            //FIXME: build multiple k trees
            //partition into k02 parts
            build_k2_trees(v, temp_file_prefix, bottom_up, access_shortcut_size);
            //build tree
        }

        k2_tree_partitioned(int_vector_buffer<> &buf_x,
                            int_vector_buffer<> &buf_y, bool bottom_up=false, uint access_shortcut_size = 0) {
            using namespace k2_treap_ns;
            typedef int_vector_buffer<> *t_buf_p;
            std::vector<t_buf_p> bufs = {&buf_x, &buf_y};

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

            uint64_t x = max_buf_element();
            uint8_t res = 0;
            while (res <= 64 and precomp<t_k0>::exp(res) <= x) { ++res; }
            if (res == 65) {
                throw std::logic_error("Maximal element of input is too big.");
            }

            if (precomp<t_k0>::exp(res) <= std::numeric_limits<uint32_t>::max()) {
                auto v = read < uint32_t, uint32_t>(bufs);
                build_k2_trees(v, buf_x.filename(), bottom_up, access_shortcut_size);
            } else {
                auto v = read < uint64_t, uint64_t>(bufs);
                build_k2_trees(v, buf_x.filename(), bottom_up, access_shortcut_size);
            }
        }

        template<typename t_x=uint64_t, typename t_y=uint64_t>
        std::vector<std::pair<t_x, t_y>> read(std::vector<int_vector_buffer<> * > &bufs) {
            typedef std::vector<std::pair<t_x, t_y>> t_tuple_vec;
            t_tuple_vec v = t_tuple_vec(bufs[0]->size());
            for (uint64_t j = 0; j < v.size(); ++j) {
                std::get<0>(v[j]) = (*(bufs[0]))[j];
            }
            for (uint64_t j = 0; j < v.size(); ++j) {
                std::get<1>(v[j]) = (*(bufs[1]))[j];
            }

            return v;
        }

        template<typename t_x>
        void direct_links(t_x source_id, std::vector<t_x> &result) const {
            result.clear();
            uint y = source_id/(m_matrix_size/t_k0);
            //uint submatrix_size = m_matrix_size/k0;

            //TODO: might be slow to use extra vector as reference, maybe it's better to remove clear() in k2treap and add directly to endresult
            std::vector<t_x> tmp_result;
            for (int j = 0; j < t_k0; ++j) {
                    m_k2trees[y*t_k0+j].direct_links(source_id, tmp_result);
                    result.insert(result.end(), tmp_result.begin(), tmp_result.end());
            }
        }

        template<typename t_x>
        void direct_links2(t_x source_id, std::vector<t_x> &result) const {
            result.clear();
            uint y = source_id/(m_matrix_size/t_k0);
            //uint submatrix_size = m_matrix_size/k0;

            //TODO: might be slow to use extra vector as reference, maybe it's better to remove clear() in k2treap and add directly to endresult
            std::vector<t_x> tmp_result;
            for (int j = 0; j < t_k0; ++j) {
                m_k2trees[y*t_k0+j].direct_links2(source_id, tmp_result);
                result.insert(result.end(), tmp_result.begin(), tmp_result.end());
            }
        }

        template<typename t_x>
        void inverse_links(t_x source_id, std::vector<t_x> &result) const {
            result.clear();
            uint x = source_id/(m_matrix_size/t_k0);
            //uint submatrix_size = m_matrix_size/k0;

            //TODO: might be slow to use extra vector as reference, maybe it's better to remove clear() in k2treap and add directly to endresult
            std::vector<t_x> tmp_result;
            for (int j = 0; j < t_k0; ++j) {
                m_k2trees[j*t_k0+x].inverse_links(source_id, tmp_result);
                result.insert(result.end(), tmp_result.begin(), tmp_result.end());
            }
        }

        template<typename t_x>
        void inverse_links2(t_x source_id, std::vector<t_x> &result) const {
            result.clear();
            uint x = source_id/(m_matrix_size/t_k0);
            //uint submatrix_size = m_matrix_size/k0;

            //TODO: might be slow to use extra vector as reference, maybe it's better to remove clear() in k2treap and add directly to endresult
            std::vector<t_x> tmp_result;
            for (int j = 0; j < t_k0; ++j) {
                m_k2trees[j*t_k0+x].inverse_links(source_id, tmp_result);
                result.insert(result.end(), tmp_result.begin(), tmp_result.end());
            }
        }

        /**
         * Checks whether link from p = link.first to q = link.second is present i.e. matrix entry a_pq = 1
         */
        template<typename t_x, typename t_y>
        bool check_link(std::pair<t_x, t_y> link) const {
            uint x = link.first/(m_matrix_size/t_k0);
            uint y = link.second/(m_matrix_size/t_k0);

            uint corresponding_tree = x*t_k0+y;
            return m_k2trees[corresponding_tree].check_link(link);
        }

        //! Move assignment operator
        k2_tree_partitioned &operator=(k2_tree_partitioned &&tr) {
            if (this != &tr) {
                m_size = tr.m_size;
                m_k2trees = std::move(tr.m_k2trees);
                m_matrix_size = std::move(tr.m_matrix_size);
            }
            return *this;
        }

        //! Assignment operator
        k2_tree_partitioned &operator=(k2_tree_partitioned &tr) {
            if (this != &tr) {
                m_size = tr.m_size;
                m_k2trees = tr.m_k2trees;
                m_matrix_size = tr.m_matrix_size;
            }
            return *this;
        }

        //! Equals operator
        bool operator==(const k2_tree_partitioned &tr) const {
            if (m_size != tr.m_size)
                return false;

            if (m_k2trees != tr.m_k2trees){
                return false;
            }

            if (m_matrix_size != tr.m_matrix_size){
                return false;
            }

            return true;
        }

        //! Number of points in the 2^k treap
        size_type
        size() const {
            return m_size;
        }

        //! Swap operator
        void swap(k2_tree_partitioned &tr) {
            if (this != &tr) {
                std::swap(m_size, tr.m_size);
                std::swap(m_k2trees, tr.m_k2trees);
                std::swap(m_matrix_size, tr.m_matrix_size);
            }
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream &out, structure_tree_node *v = nullptr,
                            std::string name = "") const {

            structure_tree_node *child = structure_tree::add_child(
                    v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "s");
            written_bytes += write_member(m_matrix_size, out, child, "matrix_size");

            for (uint j = 0; j < m_k2trees.size(); ++j) {
                written_bytes += m_k2trees[j].serialize(out, child, "k2_tree_"+j);
            }

            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream &in) {
            read_member(m_size, in);
            read_member(m_matrix_size, in);
            m_k2trees.resize(t_k0*t_k0);
            for (uint j = 0; j < t_k0*t_k0; ++j) {
                m_k2trees[j].load(in);
            }

        }

    private:
        template<typename t_tv>
        uint8_t get_tree_height(const t_tv &v) {
            using namespace k2_treap_ns;
            if (v.size() == 0) {
                return 0;
            }
            using t_e = typename t_tv::value_type;
            auto tupmax = [](t_e a) {
                return std::max(a.first, a.second);
            };
            auto max_it = std::max_element(std::begin(v), std::end(v), [&](t_e a, t_e b) {
                return tupmax(a) < tupmax(b);
            });
            uint64_t x = tupmax(*max_it);

            uint8_t res = 0;
            while (precomp<t_k0>::exp(res) <= x) { ++res; }
            return res;
        }

        template<typename t_vector>
        void build_k2_trees(t_vector &links, std::string temp_file_prefix = "", bool bottom_up = false, uint access_shortcut_size = 0) {
            using namespace k2_treap_ns;
            typedef decltype(links[0].first) t_x;
            typedef decltype(links[0].second) t_y;
            using t_e = std::pair<t_x, t_y>;

            m_size = links.size();

            //typedef typename stxxl::VECTOR_GENERATOR<t_e>::result stxxl_pair_vector;

            //FIXME replace with stxxl vector
            std::vector<std::vector<t_e>> buffers;
            buffers.resize(t_k0*t_k0);

            uint tree_height = get_tree_height(links);
            {
                m_matrix_size = precomp<t_k0>::exp(tree_height); //could be bigger than 32 bit although biggest value in links is 32 bit

                auto upper_left = t_e(0, 0);
                auto lower_right = std::make_pair(m_matrix_size - 1, m_matrix_size - 1);
                auto links_interval = std::make_pair<t_x, t_x>(0, links.size());

                auto submatrix_size =
                        lower_right.first - upper_left.first + 1;//precomp<k>::exp(m_tree_height-level);
                std::vector<t_x> intervals(t_k0 * t_k0 + 1);

                //do counting sort
                auto x1 = upper_left.first;
                auto y1 = upper_left.second;
                auto subDivK = (submatrix_size / t_k0);
                for (uint64_t j = links_interval.first; j < links_interval.second; ++j) {
                    auto x = links[j].first;
                    auto y = links[j].second;
                    uint p1 = (x - x1) / subDivK;
                    uint p2 = (y - y1) / subDivK;
                    uint corresponding_matrix = p1 * t_k0 + p2;
                    intervals[corresponding_matrix + 1]++;//offset corresponding matrix by one to allow for more efficient in interval comparision
                    buffers[corresponding_matrix].push_back(links[j]);
                }
            }

            m_k2trees.resize(t_k0*t_k0);
            for (int l = 0; l < t_k0*t_k0; ++l) {
                const subk2_tree k2tree(temp_file_prefix, bottom_up, access_shortcut_size, buffers[l]);
                m_k2trees[l] = k2tree;//buffers[l]);//, temp_file_prefix, bottom_up, access_shortcut_size);
            }
        }
    };
}

#endif