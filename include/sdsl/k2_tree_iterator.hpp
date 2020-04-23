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
    // using difference_type = std::ptrdiff_t;
    using pointer = edge *;
    using reference = edge &;
    using iterator_category = std::forward_iterator_tag;

    edge_iterator() {}
    edge_iterator(const t_bv &k_t, const t_bv &k_l, const t_rank &k_t_rank, const uint16_t &k_height) : k_t(&k_t), k_l(&k_l), k_t_rank(&k_t_rank), k_height(k_height)
    {
        _initialize();
    }

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
            this->curr_node = other.curr_node;
            this->curr_neigh = other.curr_neigh;
            this->curr_row = other.curr_row;
            this->curr_col = other.curr_col;
            this->_n = other._n;
        }
        return *this;
    }

    ~edge_iterator() {}

protected:
    pointer _ptr;

    // container TODO: pass the tree
    const t_bv *k_t;
    const t_bv *k_l;
    const t_rank *k_t_rank;
    uint16_t k_height;
    //

    // iterator state //
    size_t size;
    idx_type curr_node, curr_neigh;
    unsigned curr_row, curr_col;
    size_type _n;
    //

    void _initialize()
    {
        // if (k_l.size() == 0 && k_t.size() == 0)
        //     return acc;
        //TODO: Take care of this edge case
        _n = static_cast<size_type>(std::pow(k, k_height)) / k;
        curr_node = 0;
        curr_neigh = k * std::floor(curr_node / static_cast<double>(_n));
        curr_row = 0;
        curr_col = 0;

        size = std::pow(k, k_height);
        edge first = _find_next();

        _ptr = new edge(std::get<0>(first), std::get<1>(first));
    }

    edge _find_next()
    {
        idx_type neigh;
        if (curr_node < size)
        {
            for (; curr_row < k; curr_row++)
            {
                neigh = size;
                _find_next_recursive(_n / k, curr_node % _n, _n * curr_row, curr_neigh + curr_row, neigh);
                if (neigh < size)
                {
                    // cout << "x: " << curr_node << " y: " << neigh << endl;
                    return edge(curr_node, neigh);
                }
            }
            if (curr_row >= k)
                curr_row = 0;
            curr_node++;
            curr_neigh = k * std::floor(curr_node / static_cast<double>(_n));
            curr_col = 0;
            return _find_next();
        }
        return edge(size, size);
    }

    bool _find_next_recursive(size_type n, idx_type row, idx_type col, size_type level, idx_type &neigh)
    {
        if (level >= k_t->size())
        { // Last level
            if ((*k_l)[level - k_t->size()] == 1)
            {
                neigh = col;
                return false;
            }
        }

        if ((*k_t)[level] == 1)
        {
            curr_neigh = (*k_t_rank)(level + 1) * std::pow(k, 2) +
                         k * std::floor(row / static_cast<double>(n));
            while (curr_col < k)
            {
                if (_find_next_recursive(n / k, row % n, col + n * curr_col, curr_neigh + curr_col, neigh))
                {
                    curr_col++;
                }
                else
                    break;
            }
            if (curr_col < k)
                curr_col = 0;
        }
        // did not found
        return true;
    }
};
} // namespace sdsl

#endif