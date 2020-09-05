#ifndef INCLUDED_SDSL_K2_TREE_ITERATOR
#define INCLUDED_SDSL_K2_TREE_ITERATOR

#include <tuple>
#include <iterator>
#include <iostream>
#include <memory>
#include <deque>
#include <vector>

#include <sdsl/k2_tree.hpp>

using namespace std;

namespace sdsl
{
    using idx_type = k2_tree_ns::idx_type;
    using size_type = k2_tree_ns::size_type;
    using edge = std::tuple<idx_type, idx_type>;


    // typedef struct tree_node
    // {
    //     int node;
    //     size_type n;
    //     unsigned row, col;
    //     idx_type level;
    //     unsigned j;
    //     size_type y;
    // } tree_node;

    // template <class k_tree>
    // class edge_iterator
    // {
    // public:
    //     using value_type = edge;
    //     using pointer = shared_ptr<edge>;
    //     using reference = edge &;
    //     using iterator_category = std::forward_iterator_tag;

    //     edge_iterator() {}

    //     edge_iterator(const k_tree *tree)
    //     {
    //         this->tree = tree;
    //         this->k = tree->k_();
    //         _initialize();
    //     }

    //     idx_type x() {
    //         return get<0>(_ptr);
    //     }

    //     idx_type y() {
    //         return get<1>(_ptr);
    //     }

    //     value_type operator*()
    //     {
    //         return _ptr;
    //     }

    //     bool operator==(const edge_iterator<k_tree> &rhs) const
    //     {
    //         return equal_edge(rhs._ptr, _ptr);
    //     }

    //     bool operator!=(const edge_iterator<k_tree> &rhs) const
    //     {
    //         return !(*this == rhs);
    //     }

    //     void print()
    //     {
    //         cout << "_ptr " << std::get<0>(_ptr) << ", " << std::get<1>(_ptr) << endl;
    //         cout << "_node " << _node << endl;
    //         cout << "_level " << _level << endl;
    //         cout << "_row " << _row << endl;
    //         cout << "st " << endl;

    //         for (tree_node st_i : st)
    //         {

    //             cout << "   node: " << st_i.node << endl;
    //             cout << "   row: " << st_i.row << endl;
    //             cout << "   n: " << st_i.n << endl;
    //             cout << "   level: " << st_i.level << endl;
    //             cout << "   j: " << st_i.j << endl;
    //             cout << "   y: " << st_i.y << endl;
    //             cout << endl;
    //         }
    //     }

    //     edge_iterator<k_tree> &operator++(int)
    //     {
    //         operator++();
    //         return *this;
    //     }

    //     edge_iterator<k_tree> &operator++()
    //     {
    //         while (!st.empty()) // did not go until the end of the subtree
    //         {
    //             tree_node &last_found = st.back();
    //             last_found.j++;
    //             idx_type neigh;
    //             if (last_found.n == _n / k && last_found.j < k)
    //             {
    //                 if (_find_next_recursive(last_found.n / k,
    //                                          last_found.row % last_found.n,
    //                                          last_found.col + last_found.n * last_found.j,
    //                                          last_found.y + last_found.j,
    //                                          neigh, 0))
    //                 {
    //                     _ptr = edge(last_found.node, neigh);
    //                     return *this;
    //                 }
    //             }
    //             if (last_found.j < k)
    //             {
    //                 if (_find_next_recursive(last_found.n / k,
    //                                          last_found.row % last_found.n,
    //                                          last_found.col + last_found.n * last_found.j,
    //                                          last_found.y + last_found.j,
    //                                          neigh, last_found.j))
    //                 {
    //                     _ptr = edge(last_found.node, neigh);
    //                     return *this;
    //                 }
    //             }
    //             st.pop_back();

    //         } // move to another subtree
    //         _row++;
    //         if (_row >= k)
    //         {
    //             _row = 0;
    //             _node++;
    //             _level = k * std::floor(_node / static_cast<double>(_n));
    //         }
    //         _ptr = _find_next();
    //         return *this;
    //     }

    //     edge_iterator<k_tree> end()
    //     {
    //         edge_iterator<k_tree> it = *this;
    //         it._ptr = edge(size, size); //end node

    //         it._node = size - 1;
    //         it._level = k * std::floor(it._node / static_cast<double>(it._n));
    //         it._row = k;

    //         return it;
    //     }

    //     friend void swap(edge_iterator<k_tree> &rhs, edge_iterator<k_tree> &lhs)
    //     {
    //         if (&rhs != &lhs)
    //         {
    //             std::swap(rhs.tree, lhs.tree);
    //             std::swap(rhs.k, lhs.k);

    //             std::swap(rhs._ptr, lhs._ptr);
    //             std::swap(rhs.size, lhs.size);
    //             std::swap(rhs._node, lhs._node);
    //             std::swap(rhs._level, lhs._level);
    //             std::swap(rhs._row, lhs._row);
    //         }
    //     }

    // private:
    //     bool equal_edge(const edge &e1, const edge &e2) const
    //     {
    //         idx_type e1_x = std::get<0>(e1);
    //         idx_type e1_y = std::get<1>(e1);

    //         idx_type e2_x = std::get<0>(e2);
    //         idx_type e2_y = std::get<1>(e2);

    //         return e1_x == e2_x && e1_y == e2_y;
    //     }

    //     value_type _ptr;
    //     // container
    //     const k_tree *tree;
    //     uint8_t k = 2;
    //     //
    //     // iterator state //
    //     deque<tree_node> st;
    //     size_type _n;
    //     idx_type _level;
    //     int _row, _node, size;
    //     //

    //     void
    //     _initialize()
    //     {
    //         _n = static_cast<size_type>(std::pow(k, tree->height())) / k;
    //         _node = 0;
    //         _row = 0;
    //         _level = k * std::floor(_node / static_cast<double>(_n));
    //         size = std::pow(k, tree->height());

    //         if (tree->l_size() > 0)
    //         {
    //             _ptr = _find_next();
    //         }
    //         else
    //         {
    //             // if its empty the begin == end
    //             _ptr = edge(size, size); //end node
    //         }
    //     }

    //     edge _find_next()
    //     {
    //         idx_type neigh;
    //         if (_node < size)
    //         {
    //             for (; _row < k; _row++)
    //             {
    //                 neigh = size;
    //                 _find_next_recursive(_n / k, _node % _n, _n * _row, _level + _row, neigh, 0);
    //                 if (neigh < (idx_type)size)
    //                 {
    //                     return edge(_node, neigh);
    //                 }
    //                 assert(st.empty());
    //             }
    //             _row = 0;
    //             _node++;
    //             _level = k * std::floor(_node / static_cast<double>(_n));
    //             return _find_next();
    //         }
    //         return edge(size, size); // end node
    //     }

    //     bool _find_next_recursive(size_type n, unsigned row, unsigned col, size_type level,
    //                               idx_type &neigh, unsigned initial_j)
    //     {
    //         if (level >= tree->t().size())
    //         { // Last level
    //             if (tree->l()[level - tree->t().size()] == 1)
    //             {
    //                 neigh = col;
    //                 return true;
    //             }
    //             return false;
    //         }

    //         if (tree->t()[level] == 1)
    //         {
    //             size_type y = tree->rank_t()(level + 1) * k * k +
    //                           k * std::floor(row / static_cast<double>(n));

    //             for (unsigned j = initial_j; j < k; j++)
    //             {
    //                 tree_node next_node = {_node, n, row, col, level, j, y};
    //                 st.push_back(next_node);
    //                 if (_find_next_recursive(n / k, row % n, col + n * j, y + j, neigh, 0))
    //                     return true;
    //                 st.pop_back();
    //             }
    //         }
    //         return false;
    //     }
    // };
    ///////////////////////////////////////////////////////////////////////////////////////////////
    template <class k_tree>
    class edge_iterator
    {
    public:
        using value_type = edge;
        using pointer = shared_ptr<edge>;
        using reference = edge &;
        using iterator_category = std::forward_iterator_tag;

        edge_iterator() = default;


        edge_iterator(const k_tree *tree)
        {
            this->tree = tree;
            this->k = tree->k_();
            nodes = tree->get_number_nodes();
            _initialize();
        }

        idx_type x() {
            return get<0>(_ptr);
        }

        idx_type y() {
            return get<1>(_ptr);
        }

        value_type operator*()
        {
            return _ptr;
        }

        bool operator==(const edge_iterator<k_tree> &rhs) const
        {
            return equal_edge(rhs._ptr, _ptr);
        }

        bool operator!=(const edge_iterator<k_tree> &rhs) const
        {
            return !(*this == rhs);
        }

        edge_iterator<k_tree> &operator++(int)
        {
            operator++();
            return *this;
        }

        edge_iterator<k_tree> &operator++()
        {
            edge e(nodes, nodes);
            while (!st.empty()) // did not go until the end of the subtree
            {
                struct edge_node last_found = st.back();
                st.pop_back();

                int type = last_found.type;
                uint i = last_found.i;
                uint j = last_found.j;

                uint64_t dp = last_found.dp;
                uint64_t dq = last_found.dq;
                int64_t y = last_found.y;
                int l = last_found.l;
                //

                j++;
                if(j < k) {
                    if(_find_next_rec(dp, dq, y, l, type, i, j, e))
                        break;
                } else {
                    j = 0;
                    i++;
                    if(i < k) {
                        if(_find_next_rec(dp, dq, y, l, type, i, j, e))
                            break;
                    }
                }
            }
            _ptr = e;
            return *this;
        }

        edge_iterator<k_tree> end()
        {
            edge_iterator<k_tree> it;
            it._ptr = edge(nodes, nodes); //end node
            return it;
        }

        friend void swap(edge_iterator<k_tree> &rhs, edge_iterator<k_tree> &lhs)
        {
            if (&rhs != &lhs)
            {
                std::swap(rhs.tree, lhs.tree);
                std::swap(rhs.k, lhs.k);

                std::swap(rhs._ptr, lhs._ptr);
                std::swap(rhs.nodes, lhs.nodes);
                std::swap(rhs.max_level, lhs.max_level);
                std::swap(rhs.div_level_table, lhs.div_level_table);
                std::swap(rhs.st, lhs.st);
            }
        }

    private:
        bool equal_edge(const edge &e1, const edge &e2) const
        {
            idx_type e1_x = std::get<0>(e1);
            idx_type e1_y = std::get<1>(e1);

            idx_type e2_x = std::get<0>(e2);
            idx_type e2_y = std::get<1>(e2);

            return e1_x == e2_x && e1_y == e2_y;
        }

        value_type _ptr;
        // container
        const k_tree *tree;
        uint8_t k = 2;
        //
        // iterator state //
        class edge_node {
        public:
            uint64_t dp, dq, i, j;
            int64_t y;
            int type, l;
            
            edge_node() {}
            edge_node(uint64_t dp, uint64_t dq, int64_t y, uint64_t i, uint64_t j, int type, int l) : dp(dp), dq(dq),
            i(i), j(j), y(y), type(type), l(l) {}
        };
        std::deque<edge_node> st;
        uint32_t max_level;
        uint64_t nodes;
        std::vector<uint64_t> div_level_table;
        std::vector<uint64_t> pointerL;
        //

        uint exp_pow(uint base, uint pow)
        {
            uint i, result = 1;
            for (i = 0; i < pow; i++)
                result *= base;

            return result;
        }

        void
        _initialize()
        {
            if (tree->l_size() > 0)
            {
                max_level = floor(log(nodes)/log(k));
                if(floor(log(nodes)/log(k)) == (log(nodes)/log(k)))
                    max_level = max_level-1;

                div_level_table = std::vector<uint64_t>(max_level+1);
                for(uint i = 0; i <= max_level; i++)
                    div_level_table[i] = exp_pow(k, max_level-i);

                pointerL = std::vector<uint64_t>(max_level+1);

                pointerL[0] = 0;
                pointerL[1] = k*k;

                for(uint64_t i = 2; i < max_level; i++) {
                    pointerL[i] = (tree->rank_t(pointerL[i-1])+1)*k*k;
                }
                pointerL[max_level] = 0;

                edge e = edge(nodes, nodes);
                _find_next(0, 0, -1, -1, e);
                _ptr = e;
            }
            else
            {
                // if its empty the begin == end
                _ptr = edge(nodes, nodes); //end node
            }
        }

        bool _find_next(uint64_t dp, uint64_t dq, int64_t x, int l, edge &e)
        {
            if(l == (int64_t)max_level){
                if(tree->is_bit_set_l(x)) {
                    e = edge(dp, dq);
                    return true;
                }
                return false;
            }

            if((l == (int64_t)max_level-1) && tree->is_bit_set_t(x)) {
                int64_t y = pointerL[l+1];
                pointerL[l+1] += k*k;

                for(uint64_t i = 0; i < k; i++) {
                    for(uint64_t j = 0; j < k; j++) {
                        st.push_back(edge_node(dp, dq, y, i, j, 1, l));
                        if(_find_next(dp+i, dq+j, y+k*i+j, l+1, e))
                            return true;
                        st.pop_back();
                    }
                }
                return false;
            }
            if((x == -1) || ((l < (int64_t)max_level-1) && tree->is_bit_set_t(x))) {
                int64_t y = pointerL[l+1];
                pointerL[l+1] += k*k;
                
                uint64_t div_level = div_level_table[l+1];

                for(uint64_t i = 0; i < k; i++) {
                    for(uint64_t j = 0; j < k; j++) {
                        st.push_back(edge_node(dp, dq, y, i, j, 2, l));
                        if(_find_next(dp+div_level*i, dq+div_level*j, y+k*i+j, l+1, e))
                            return true;
                        st.pop_back();
                    }
                }
                return false;
            }

            return false; 
        }


        bool _find_next_rec(uint64_t dp, uint64_t dq, int64_t y, int l, int type, uint64_t i, uint64_t j, edge &e) {
            if(type == 1) {
                st.push_back(edge_node(dp, dq, y, i, j, 1, l ));
                return _find_next(dp+i, dq +j, y+k*i+j, l+1, e);
            } else if(type == 2){
                uint64_t div_level = div_level_table[l+1];
                st.push_back(edge_node(dp, dq, y, i, j, 2, l));
                return _find_next(dp+div_level*i, dq+div_level*j, y+k*i+j, l+1, e);
            }
            return false;
        }
    };
    //////////////////////////////////////////////////////////////////////////////////////////////

    template <class k2_tree>
    class node_iterator
    {
        using node_type = k2_tree_ns::idx_type;
        using size_type = k2_tree_ns::size_type;

    public:
        using value_type = node_type;
        using pointer = shared_ptr<node_type>;
        using reference = node_type &;
        using iterator_category = std::forward_iterator_tag;

        node_iterator() {}

        node_iterator(const k2_tree *k_tree)
        {
            this->tree = k_tree;

            if (tree->get_number_edges() == 0)
            {
                value_type num_nodes = tree->get_number_nodes();
                _ptr = num_nodes;
                curr_i = num_nodes;
                curr_j = num_nodes;
                return;
            }

            for (value_type i = 0; i < tree->get_number_nodes(); i++)
            {
                for (value_type j = i; j < tree->get_number_nodes(); j++)
                    if (tree->adj(i, j))
                    {
                        _ptr = i;
                        curr_i = i;
                        curr_j = j;
                        return;
                    }
                for (value_type j = i; j < tree->get_number_nodes(); j++)
                    if (tree->adj(j, i))
                    {
                        _ptr = i;
                        curr_i = i + 1;
                        curr_j = i;
                        return;
                    }
            }
        }

        node_iterator(node_iterator<k2_tree> *other)
        {
            this->_ptr = other->_ptr;
            this->tree = other->tree;
            this->curr_i = other->curr_i;
            this->curr_j = other->curr_j;
        }

        ~node_iterator() {}

        value_type operator*() const
        {
            return _ptr;
        }

        node_iterator<k2_tree> &operator=(const node_iterator<k2_tree> &other)
        {
            if (this != &other)
            {
                tree = other.tree;
                _ptr = other._ptr;
                curr_i = other.curr_i;
                curr_j = other.curr_j;
            }
            return *this;
        }

        bool operator==(const node_iterator<k2_tree> &rhs) const
        {
            return rhs._ptr == _ptr;
        }

        bool operator!=(const node_iterator<k2_tree> &rhs) const
        {
            return !(*this == rhs);
        }

        node_iterator<k2_tree> &operator++()
        {
            for (value_type i = curr_i; i < tree->get_number_nodes(); i++)
            {
                for (value_type j = curr_j; j < tree->get_number_nodes(); j++)
                    if (tree->adj(i, j))
                        return update_state(i);
                for (value_type j = 0; j < tree->get_number_nodes(); j++)
                {
                    if (j == i)
                        continue;
                    if (tree->adj(j, i))
                        return update_state(i);
                }
                curr_j = 0;
            }
            return *this;
        }

        node_iterator<k2_tree> &operator++(int)
        {
            if (_ptr == tree->get_number_nodes() - 1)
            {
                *this = end();
                return *this;
            }
            value_type prev = _ptr;
            operator++();
            if (_ptr == prev)
                this->end();
            return *this;
        }

        node_iterator<k2_tree> end()
        {
            if (tree == NULL)
                return *this;
            node_iterator<k2_tree> it = *this;
            value_type num_nodes = it.tree->get_number_nodes();
            it._ptr = num_nodes; //end node
            it.curr_i = num_nodes;
            it.curr_j = num_nodes;
            return it;
        }

    private:
        value_type _ptr;
        const k2_tree *tree = NULL;
        value_type curr_i, curr_j;

        node_iterator<k2_tree> &update_state(node_type i)
        {
            _ptr = i;
            curr_i = i + 1;
            curr_j = 0;

            return *this;
        }
    };

    template <class k_tree>
    class neighbour_iterator
    {
        using idx_type = k2_tree_ns::idx_type;
        using size_type = k2_tree_ns::size_type;

    public:
        using value_type = int;
        using pointer = shared_ptr<int>;
        using reference = int &;
        using iterator_category = std::forward_iterator_tag;

        neighbour_iterator() {}

        neighbour_iterator(const k_tree *tree, value_type node)
        {
            this->tree = tree;
            this->k = tree->k_();
            this->_node = node;
            _initialize();
        }

        neighbour_iterator(const k_tree *tree)
        {
            this->tree = tree;
            this->k = tree->k_();
            this->_node = 0;
            _initialize();
        }

        value_type operator*()
        {
            return _ptr;
        }

        neighbour_iterator<k_tree> &operator++(int)
        {
            if (_ptr != -1)
                operator++();
            return *this;
        }

        neighbour_iterator<k_tree> &operator++()
        {
            while (!st.empty()) // did not go until the end of the subtree
            {
                tree_node2 &last_found = st.back();
                last_found.j++;
                value_type neigh;

                if (last_found.j < k)
                {
                    if (_find_next_recursive(last_found.n / k,
                                             last_found.row % last_found.n,
                                             last_found.col + last_found.n * last_found.j,
                                             last_found.y + last_found.j,
                                             neigh, 0))
                    {
                        _ptr = neigh;
                        return *this;
                    }
                }
                st.pop_back();

            } // move to another subtree
            _row++;
            if (_row >= k)
            {
                _ptr = -1;
                return *this;
            }
            _ptr = _find_next();
            return *this;
        }

        bool operator==(const neighbour_iterator<k_tree> &rhs) const
        {
            return rhs._ptr == this->_ptr;
        }

        bool operator!=(const neighbour_iterator<k_tree> &rhs) const
        {
            return !(*this == rhs);
        }

        neighbour_iterator<k_tree> end()
        {
            _ptr = -1;
            return *this;
        }

        neighbour_iterator<k_tree> &operator=(const neighbour_iterator<k_tree> &other) {
            if (this != &other) {
                _ptr = other._ptr;
                tree = other.tree;
                k = other.k;
                st = other.st;
                size = other.size;
                _n = other._n; 
                _node = other._node;
                _level = other._level;
                _row = other._row;
            }
            return *this;
        }

        friend void swap(neighbour_iterator<k_tree> &rhs, neighbour_iterator<k_tree> &lhs) {
            std::swap(rhs._ptr, lhs._ptr);
            std::swap(rhs.tree, lhs.tree);
            std::swap(rhs.k, lhs.k);
            std::swap(rhs.st, lhs.st);
            std::swap(rhs.size, lhs.size);
            std::swap(rhs._n, lhs._n);
            std::swap(rhs._node, lhs._node);
            std::swap(rhs._level, lhs._level);
            std::swap(rhs._row, lhs._row);
        }

    private:
        void
        _initialize()
        {
            _n = static_cast<size_type>(std::pow(k, tree->height())) / k;
            _row = 0;
            _level = k * std::floor(_node / static_cast<double>(_n));
            size = std::pow(k, tree->height());

            if (tree->l_size() > 0)
            {
                _ptr = _find_next();
            }
            else
            {
                // if its empty the begin == end
                _ptr = -1; //end node
            }
        }

        value_type _find_next()
        {
            value_type neigh;
            if (_node < size)
            {
                for (; _row < k; _row++)
                {
                    neigh = size;
                    _find_next_recursive(_n / k, _node % _n, _n * _row, _level + _row, neigh, 0);
                    if (neigh < size)
                    {
                        return neigh;
                    }
                }
            }
            return -1; // end node
        }

        bool _find_next_recursive(value_type n, value_type row, value_type col, uint level, value_type &neigh, unsigned initial_j)
        {
            if (level >= tree->t_size()) // Last level
            {
                if (tree->is_bit_set_l(level - tree->t_size()) == 1)
                {
                    neigh = col;
                    return true;
                }
                return false;
            }

            if (tree->is_bit_set_t(level) == 1)
            {
                value_type y = tree->rank_t(level + 1) * k * k +
                              k * std::floor(row / static_cast<double>(n));

                for (unsigned j = initial_j; j < k; j++)
                {
                    tree_node2 next_node = {n, row, col, level, j, y};
                    st.push_back(next_node);
                    if (_find_next_recursive(n / k, row % n, col + n * j, y + j, neigh, 0))
                        return true;
                    st.pop_back();
                }
                return false;
            }
            return false;
        }

        typedef struct tree_node2
        {
            value_type n, row, col;
            uint level, j;
            value_type y;
        } tree_node2;

        value_type _ptr;
        const k_tree *tree;
        uint8_t k = 2;

        // iterator state //
        deque<tree_node2> st;
        value_type size, _n, _node, _level, _row;
        //
    };
} // namespace sdsl

#endif