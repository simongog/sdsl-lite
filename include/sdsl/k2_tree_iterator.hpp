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

        //TODO: Change edge to value_type
        edge operator*()
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
                this->curr_node = other.curr_node;
                this->curr_neigh = other.curr_neigh;
                this->curr_row = other.curr_row;
                this->curr_col = other.curr_col;
            }
            return *this;
        }

        bool operator==(const edge_iterator<k, t_bv, t_rank> &rhs) const
        {
            return equal_edge(*(rhs._ptr), *(this->_ptr));
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
                std::swap(lhs.curr_node, rhs.curr_node);
                std::swap(lhs.curr_neigh, rhs.curr_neigh);
                std::swap(lhs.curr_row, rhs.curr_row);
                std::swap(lhs.curr_col, rhs.curr_col);
            }
        }

        edge_iterator<k, t_bv, t_rank> &operator++()
        {
            curr_col++;
            if (curr_col >= k)
            {
                curr_col = 0;
                curr_row++;
            }
            edge e = _find_next();
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
        //     // curr_col++;
        //     // if (curr_col >= k)
        //     // {
        //     //     curr_col = 0;
        //     //     curr_row++;
        //     // }
        //     // edge e = _find_next();
        //     // _ptr = new edge(std::get<0>(e), std::get<1>(e));
        //     // return *this;
        // }
        // edge_iterator<k, t_bv, t_rank> &operator--(int)
        // {
        //     edge_iterator<k, t_bv, t_rank> *tmp;
        //     tmp = new edge_iterator<k, t_bv, t_rank>(*(this->k_t), *(this->k_l), *(this->k_t_rank), this->k_height);
        //     operator--();
        //     return *tmp;
        // }

        ~edge_iterator() {}

        edge_iterator<k, t_bv, t_rank> end()
        { //TODO: find out about backtrack the tree

            // if(k_l->size() == 0) {
            //     return *this;
            // }
            // idx_type old_curr_node = curr_node;
            // idx_type old_curr_neigh = curr_neigh;
            // unsigned old_curr_row = curr_row;
            // unsigned old_curr_col = curr_col;
            // pointer old_ptr = _ptr;

            // edge last_edge = edge(size, size);
            // edge_iterator<k, t_bv, t_rank> prev = *this;
            // while (!equal_edge(*_ptr, last_edge))
            // {
            //     prev = *this;
            //     operator++();
            // }

            // curr_node = old_curr_node;
            // curr_neigh = old_curr_neigh;
            // curr_col = old_curr_col;
            // curr_row = old_curr_row;
            // _ptr = old_ptr;
            edge_iterator<k, t_bv, t_rank> it;
            it._ptr = new edge(size, size); //end node
            return it;
        }
        friend ostream &operator<<(ostream &os, const edge_iterator<k, t_bv, t_rank> &edg) {
            os << " ==== ktree Edge Iterator ==== " << endl;
            edge e = *edg._ptr;
            os << " ptr (" << get<0>(e) << ", " << get<1>(e) << ")" << endl;
            os << " CONTAINER " << endl;
            os << "     k_t ";
            for(uint i = 0; i < edg.k_t->size(); i++)
                os << (*edg.k_t)[i];
            os << endl;
            os << "     k_l ";
            for(uint i = 0; i < edg.k_l->size(); i++)
                os << (*edg.k_l)[i];
            os << endl;
            os << "     k_height " << edg.k_height << endl;
            os << " STATE " << endl;
            os << "     size " << edg.size << endl;
            os << "     _n " << edg._n << endl;
            os << "     curr_node " << edg.curr_node << endl;
            os << "     curr_neigh " << edg.curr_neigh << endl;
            os << "     curr_row " << edg.curr_row << endl;
            os << "     curr_col " << edg.curr_col << endl;
            os << " ============================ " << endl;
            return os;
        }

    protected:
        pointer _ptr;

        // container TODO: pass the tree (?)
        const t_bv *k_t;
        const t_bv *k_l;
        const t_rank *k_t_rank;
        uint16_t k_height;
        //

        // iterator state //
        size_t size;
        size_type _n;
        idx_type curr_node, curr_neigh;
        unsigned curr_row, curr_col;
        //

        void
        _initialize()
        {
            _n = static_cast<size_type>(std::pow(k, k_height)) / k;
            curr_node = 0;
            curr_neigh = k * std::floor(curr_node / static_cast<double>(_n));
            curr_row = 0;
            curr_col = 0;
            size = std::pow(k, k_height);
            if (k_l->size() > 0)
            {
                edge first = _find_next();
                first_edge = new edge(std::get<0>(first), std::get<1>(first));
                _ptr = first_edge;
            }
            else
            {
                // if its empty the begin == end
                first_edge = new edge(size, size);
                _ptr = new edge(size, size); //end node
            }
        }

        edge _find_next()
        {
            idx_type neigh;
            if (curr_node < size)
            {
                for (; curr_row < k; curr_row++)
                {
                    neigh = size;
                    _find_next_recursive(_n / k, curr_node % _n, _n * curr_row, curr_neigh + curr_row, neigh, curr_col, curr_col);
                    if (neigh < size)
                    {
                        // cout << "x: " << curr_node << " y: " << neigh << endl;
                        return edge(curr_node, neigh);
                    }
                }
                curr_row = 0;
                curr_col = 0;
                curr_node++;
                curr_neigh = k * std::floor(curr_node / static_cast<double>(_n));
                return _find_next();
            }
            return edge(size, size); // end node
        }

        unsigned _find_next_recursive(size_type n, idx_type row, idx_type col, size_type level, idx_type &neigh, unsigned &col_state, unsigned initial_j)
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

            if ((*k_t)[level] == 1 && n > 0)
            {
                size_type y = (*k_t_rank)(level + 1) * std::pow(k, 2) +
                              k * std::floor(row / static_cast<double>(n));
                for (unsigned j = initial_j; j < k; j++)
                    if (_find_next_recursive(n / k, row % n, col + n * j, y + j, neigh, col_state, 0))
                    {
                        col_state = j;
                        return true;
                    }
            }
            return false;
        }

        edge _find_prev()
        {
            idx_type neigh;
            if (curr_node > 0)
            {
                for (; curr_row > 0; curr_row--)
                {
                    neigh = size;
                    _find_prev_recursive(_n / k, curr_node % _n, _n * curr_row, curr_neigh + curr_row, neigh, curr_col, curr_col);
                    if (neigh < size)
                    {
                        // cout << "x: " << curr_node << " y: " << neigh << endl;
                        return edge(curr_node, neigh);
                    }
                }
                curr_row = 0;
                curr_col = 0;
                curr_node--;
                curr_neigh = k * std::floor(curr_node / static_cast<double>(_n));
                return _find_prev();
            }
            return edge(size, size); // end node
        }

        unsigned _find_prev_recursive(size_type n, idx_type row, idx_type col, size_type level, idx_type &neigh, unsigned &col_state, unsigned initial_j)
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

            if ((*k_t)[level] == 1 && n > 0)
            {
                size_type y = (*k_t_rank)(level + 1) * std::pow(k, 2) +
                              k * std::floor(row / static_cast<double>(n));
                for (unsigned j = initial_j; j > 0; j--)
                    if (_find_next_recursive(n / k, row % n, col + n * j, y + j, neigh, col_state, 0))
                    {
                        col_state = j;
                        return true;
                    }
            }
            return false;
        }

    private:
        edge *first_edge;

        bool equal_edge(const edge &e1, const edge &e2) const
        {
            idx_type e1_x = std::get<0>(e1);
            idx_type e1_y = std::get<1>(e1);

            idx_type e2_x = std::get<0>(e2);
            idx_type e2_y = std::get<1>(e2);

            return e1_x == e2_x && e1_y == e2_y;
        }
    };
} // namespace sdsl

#endif