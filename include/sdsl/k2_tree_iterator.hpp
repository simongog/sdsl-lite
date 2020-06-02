#ifndef INCLUDED_SDSL_K2_TREE_ITERATOR
#define INCLUDED_SDSL_K2_TREE_ITERATOR

#include <tuple>
#include <iterator>
#include <iostream>

#include <sdsl/k2_tree.hpp>

using namespace std;

namespace sdsl
{

    template <uint8_t k,
              typename t_bv = bit_vector,
              typename t_rank = typename t_bv::rank_1_type>
    class edge_iterator
    {
        using idx_type = k2_tree_ns::idx_type;
        using size_type = k2_tree_ns::size_type;
        using edge = std::tuple<idx_type, idx_type>;
        // using k2_tree = k2_tree<k, t_bv, t_rank, l_rank>;

    public:
        using value_type = edge;
        using pointer = edge *;
        using reference = edge &;
        using iterator_category = std::forward_iterator_tag;

        edge_iterator() {}
        edge_iterator(const t_bv &k_t, const t_bv &k_l, const t_rank &k_t_rank, const uint16_t &k_height) : k_t(&k_t), k_l(&k_l), k_t_rank(&k_t_rank), k_height(k_height)
        {
            _initialize();
        }

        value_type operator*()
        {
            return *_ptr;
        }

        edge_iterator<k, t_bv, t_rank> &operator=(const edge_iterator<k, t_bv, t_rank> &other)
        {
            if (this != &other)
            {
                this->k_t = other.k_t;
                this->k_l = other.k_l;
                this->k_t_rank = other.k_t_rank;
                this->k_height = other.k_height;

                this->_ptr = other._ptr;
                this->size = other.size;
                this->_n = other._n;
                this->_node = other._node;
                this->_level = other._level;
                this->_row = other._row;
                this->_col = other._col;
            }
            return *this;
        }

        bool operator==(const edge_iterator<k, t_bv, t_rank> &rhs) const
        {
            if (rhs._ptr != NULL && this->_ptr != NULL)
                return equal_edge(*(rhs._ptr), *(this->_ptr));
            else if (rhs._ptr == NULL && this->_ptr == NULL)
                return true;
            return false;
        }
        bool operator!=(const edge_iterator<k, t_bv, t_rank> &rhs) const
        {
            return !(*this == rhs);
        }

        friend void swap(edge_iterator<k, t_bv, t_rank> &rhs, edge_iterator<k, t_bv, t_rank> &lhs)
        {
            if (lhs != rhs)
            {
                std::swap(lhs._ptr, rhs._ptr);
                std::swap(lhs.k_t, rhs.k_t);
                std::swap(lhs.k_l, rhs.k_l);
                std::swap(lhs.k_t_rank, rhs.k_t_rank);
                std::swap(lhs.size, rhs.size);
                std::swap(lhs._node, rhs._node);
                std::swap(lhs._level, rhs._level);
                std::swap(lhs._row, rhs._row);
                std::swap(lhs._col, rhs._col);
            }
        }

        edge_iterator<k, t_bv, t_rank> &operator++()
        {
            if (st.valid)
            {
                st.j++;
                if (st.j < k)
                {
                    edge e;
                    e = _find_next_middle();
                    _ptr = new edge(std::get<0>(e), std::get<1>(e));
                    return *this;
                }
            }
            st.valid = false;
            _col++;
            if (_col >= k)
            {
                _col = 0;
                _row++;
            }
            edge e;
            e = _find_next();
            _ptr = new edge(std::get<0>(e), std::get<1>(e));
            return *this;
        }
        edge_iterator<k, t_bv, t_rank> &operator++(int)
        {
            edge_iterator<k, t_bv, t_rank> *tmp;
            tmp = new edge_iterator<k, t_bv, t_rank>(*(this->k_t), *(this->k_l), *(this->k_t_rank), this->k_height);
            operator++();
            return *tmp;
        }

        // edge_iterator<k, t_bv, t_rank> &operator--()
        // {
        //     if (!equal_edge(*_ptr, *first_edge))
        //     {
        //         _col--;
        //         if (_col >= 0)
        //         {
        //             _col = k - 1;
        //             _row = _row >= k ? 0 : (_row-1);
        //         }
        //         edge e = _find_prev();
        //         _ptr = new edge(std::get<0>(e), std::get<1>(e));
        //     }
        //     return *this;
        // }
        // edge_iterator<k, t_bv, t_rank> &operator--(int)
        // {
        //     edge_iterator<k, t_bv, t_rank> *tmp;
        //     tmp = new edge_iterator<k, t_bv, t_rank>(*(this->k_t), *(this->k_l), *(this->k_t_rank), this->k_height);
        //     operator--();
        //     return *tmp;
        // }

        ~edge_iterator()
        {
            // if(_ptr != NULL) delete _ptr;
        }

        edge_iterator<k, t_bv, t_rank> end()
        {
            edge_iterator<k, t_bv, t_rank> it = *this;
            it._ptr = new edge(size, size); //end node

            it._node = size - 1;
            it._level = k * std::floor(it._node / static_cast<double>(it._n));
            it._row = k;
            it._col = k;

            return it;
        }

        friend ostream &operator<<(ostream &os, const edge_iterator<k, t_bv, t_rank> &edg)
        {
            os << " ==== ktree Edge Iterator ==== " << endl;
            edge e = *edg._ptr;
            if (get<0>(e) == edg.size && get<1>(e) == edg.size)
            {
                os << " END NODE" << endl;
                //     os << " ============================= " << endl;
                //     return os;
            }

            os << " ptr (" << get<0>(e) << ", " << get<1>(e) << ")" << endl;
            os << " CONTAINER " << endl;
            os << "     k_t ";
            if (edg.k_t != NULL)
            {
                for (uint i = 0; i < edg.k_t->size(); i++)
                    os << (*edg.k_t)[i];
                os << endl;
            }
            else
            {
                os << "null" << endl;
            }
            os << "     k_l ";
            if (edg.k_l != NULL)
            {
                for (uint i = 0; i < edg.k_l->size(); i++)
                    os << (*edg.k_l)[i];
                os << endl;
            }
            else
            {
                os << "null" << endl;
            }
            os << "     k_height " << edg.k_height << endl;
            os << " STATE " << endl;
            os << "     size " << edg.size << endl;
            os << "     _n " << edg._n << endl;
            os << "     _node " << edg._node << endl;
            os << "     _level " << edg._level << endl;
            os << "     _row " << edg._row << endl;
            os << "     _col " << edg._col << endl;
            os << " ============================= " << endl;
            return os;
        }

    protected:
        pointer _ptr = NULL;

        // container
        const t_bv *k_t = NULL;
        const t_bv *k_l = NULL;
        const t_rank *k_t_rank = NULL;
        uint16_t k_height;
        //

        typedef struct state
        {
            unsigned node, row, col, n, level, j, y;
            bool valid = false;

            void set(unsigned _node, unsigned _n, unsigned _row, unsigned _col, unsigned _level, unsigned _j, unsigned _y)
            {
                node = _node;
                row = _row;
                col = _col;
                n = _n;
                level = _level;
                j = _j;
                y = _y;
                valid = true;
            }
        } state;

        // iterator state //
        state st;
        size_t size;
        size_type _n;
        idx_type _node, _level;
        unsigned _col;
        int _row;
        // size_t last_neigh; //marks if there are still other children at the leaf to search
        //

        void
        _initialize()
        {
            _n = static_cast<size_type>(std::pow(k, k_height)) / k;
            _node = 0;
            _row = 0;
            _col = 0;
            _level = k * std::floor(_node / static_cast<double>(_n));
            size = std::pow(k, k_height);

            if (k_l->size() > 0)
            {
                edge first = _find_next();
                _ptr = new edge(std::get<0>(first), std::get<1>(first));
                ;
            }
            else
            {
                // if its empty the begin == end
                _ptr = new edge(size, size); //end node
            }
        }

        edge _find_next_middle()
        {
            if (st.valid)
            {
                idx_type neigh = size;
                _find_next_recursive(st.n / k, st.row % st.n, st.col + st.n * st.j, st.y + st.j, neigh, st.j);
                if (neigh != size)
                {
                    return edge(_node, neigh);
                }
            }

            st.valid = false;
            _col++;
            if (_col >= k)
            {
                _col = 0;
                _row++;
            }
            return _find_next();
        }

        edge _find_next()
        {
            idx_type neigh;
            if (_node < size)
            {
                for (; _row < k; _row++)
                {
                    neigh = size;
                    _find_next_recursive(_n / k, _node % _n, _n * _row, _level + _row, neigh, _col);
                    if (neigh < size)
                    {
                        // cout << "x: " << _node << " y: " << neigh << endl;
                        return edge(_node, neigh);
                    }
                }
                _row = 0;
                _col = 0;
                _node++;
                _level = k * std::floor(_node / static_cast<double>(_n));
                st.valid = false;
                return _find_next();
            }
            return edge(size, size); // end node
        }

        bool _find_next_recursive(size_type n, idx_type row, idx_type col, size_type level, idx_type &neigh, unsigned initial_j)
        {
            if (level >= k_t->size()) // Last level
            {
                if ((*k_l)[level - k_t->size()] == 1)
                {
                    neigh = col;
                    return true;
                }
                return false;
            }

            if ((*k_t)[level] == 1)
            {
                size_type y = (*k_t_rank)(level + 1) * k * k +
                              k * std::floor(row / static_cast<double>(n));

                for (unsigned j = initial_j; j < k; j++)
                {
                    st.set(_node, n, row, col, level, j, y);
                    if (_find_next_recursive(n / k, row % n, col + n * j, y + j, neigh, 0))
                    {
                        _col = j;
                        return true;
                    }
                }
            }
            return false;
        }

        // bool _find_prev_recursive(size_type n, idx_type row, idx_type col, size_type level, idx_type &neigh, unsigned &col_state, unsigned initial_j)
        // {
        //     if (level >= k_t->size()) // Last level
        //     {
        //         if ((*k_l)[level - k_t->size()] == 1)
        //         {
        //             neigh = col;
        //             return true;
        //         }
        //         return false;
        //     }

        //     if ((*k_t)[level] == 1 && n > 0)
        //     {
        //         size_type y = (*k_t_rank)(level + 1) * std::pow(k, 2) +
        //                       k * std::floor(row / static_cast<double>(n));
        //         for (int j = initial_j; j >= 0; j--) {
        //             if (_find_prev_recursive(n / k, row % n, col + n * j, y + j, neigh, col_state, k-1))
        //             {
        //                 col_state = j;
        //                 return true;
        //             }
        //         }
        //     }
        //     return false;
        // }
        // value_type _find_prev()
        // {
        //     idx_type neigh;
        //     if (_node >= 0)
        //     {
        //         for (; _row >= 0; _row--)
        //         {
        //             neigh = size;
        //             _find_prev_recursive(_n / k, _node % _n, _n * _row, _level + _row, neigh, _col, _col);
        //             if (neigh < size)
        //             {
        //                 // cout << "x: " << _node << " y: " << neigh << endl;
        //                 return edge(_node, neigh);
        //             }
        //         }
        //         _row = k - 1;
        //         _col = k - 1;
        //         _node--;
        //         _level = k * std::floor(_node / static_cast<double>(_n));

        //         return _find_prev();
        //     } else {
        //         return *first_edge;
        //     }
        // }

    private:
        bool equal_edge(const edge &e1, const edge &e2) const
        {
            idx_type e1_x = std::get<0>(e1);
            idx_type e1_y = std::get<1>(e1);

            idx_type e2_x = std::get<0>(e2);
            idx_type e2_y = std::get<1>(e2);

            return e1_x == e2_x && e1_y == e2_y;
        }
    };

    template <class k2_tree>
    class node_iterator
    {
        using node_type = k2_tree_ns::idx_type;
        using size_type = k2_tree_ns::size_type;

    public:
        using value_type = node_type;
        using pointer = node_type *;
        using reference = node_type &;
        using iterator_category = std::forward_iterator_tag;

        node_iterator() {}
        node_iterator(const k2_tree *k_tree)
        {
            this->tree = k_tree;
            for (value_type i = 0; i < tree->get_number_nodes(); i++)
            {
                for (value_type j = i; j < tree->get_number_nodes(); j++)
                    if (tree->adj(i, j))
                    {
                        _ptr = &i;
                        curr_i = i;
                        curr_j = j;
                        return;
                    }
                for (value_type j = i; j < tree->get_number_nodes(); j++)
                    if (tree->adj(j, i))
                    {
                        _ptr = &i;
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

        ~node_iterator()
        {
        }

        value_type operator*() const
        {
            return *_ptr;
        }

        bool operator==(const node_iterator<k2_tree> &rhs) const
        {
            if (rhs._ptr != NULL && this->_ptr != NULL)
                return rhs._ptr == this->_ptr;
            else if (rhs._ptr == NULL && this->_ptr == NULL)
                return true;
            return false;
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
                    {
                        return update_state(i);
                    }
                for (value_type j = 0; j < tree->get_number_nodes(); j++)
                {
                    if (j == i)
                        continue;
                    if (tree->adj(j, i))
                    {
                        return update_state(i);
                    }
                }
                curr_j = 0;
            }
            return *this;
        }

        node_iterator<k2_tree> &operator++(int)
        {
            node_iterator<k2_tree> *tmp;
            tmp = new node_iterator<k2_tree>(this);
            operator++();
            return *tmp;
        }

        node_iterator<k2_tree> end()
        {
            node_iterator<k2_tree> it = *this;
            value_type num_nodes = it.tree->get_number_nodes();
            *(it._ptr) = it.tree->get_number_nodes(); //end node
            it.curr_i = num_nodes;
            it.curr_j = num_nodes;
            return it;
        }

        friend void swap(node_iterator<k2_tree> &rhs, node_iterator<k2_tree> &lhs)
        {
            std::swap(rhs._ptr, lhs._ptr);
            std::swap(rhs.tree, lhs.tree);
            std::swap(rhs.curr_i, lhs.curr_i);
            std::swap(rhs.curr_j, lhs.curr_j);
        }

    private:
        pointer _ptr = NULL;
        const k2_tree *tree = NULL;
        value_type curr_i, curr_j;

        node_iterator<k2_tree> &update_state(int i)
        {
            *_ptr = i;
            curr_i = i + 1;
            curr_j = 0;

            return *this;
        }
    };

} // namespace sdsl

#endif