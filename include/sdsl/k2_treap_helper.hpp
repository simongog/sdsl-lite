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
/*! \file k2_treap_helper.hpp
    \brief k2_treap_helper.hpp contains helper functions and definitions for a k^2-treap implementation.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_K2_TREAP_HELPER
#define INCLUDED_SDSL_K2_TREAP_HELPER

#include "vectors.hpp"
#include "bits.hpp"
#include <tuple>
#include <algorithm>
#include <iterator>
#include <vector>
#include <complex>
#include <queue>
#include <array>

//! Namespace for the succinct data structure library.
namespace sdsl
{

namespace k2_treap_ns
{

// Precomputed value for fast k^2 treap operations
template<uint8_t t_k>
struct precomp {
    static struct impl {
        uint64_t exp[65];
        impl()
        {
            exp[0] = 1;
            for (uint8_t i=1; i<65; ++i) {
                exp[i] = t_k * exp[i-1];
            }
        }
    } data;

    static uint64_t exp(uint8_t l)
    {
        return data.exp[l];
    }

    static uint64_t divexp(uint64_t x, uint8_t l)
    {
        return x/data.exp[l];
    }

    static uint64_t modexp(uint64_t x, uint8_t l)
    {
        return x%data.exp[l];
    }
};

template<>
struct precomp<2> {
    static uint64_t exp(uint8_t l)
    {
        return 1ULL<<l;
    }

    static uint64_t divexp(uint64_t x, uint8_t l)
    {
        return x>>l;
    }

    static uint64_t modexp(uint64_t x, uint8_t l)
    {
        return x & bits::lo_set[l];
    }
};

template<>
struct precomp<4> {
    static uint64_t exp(uint8_t l)
    {
        return 1ULL<<(2*l);
    }

    static uint64_t divexp(uint64_t x, uint8_t l)
    {
        return x>>(2*l);
    }

    static uint64_t modexp(uint64_t x, uint8_t l)
    {
        return x & bits::lo_set[2*l];
    }
};

template<>
struct precomp<8> {
    static uint64_t exp(uint8_t l)
    {
        return 1ULL<<(3*l);
    }

    static uint64_t divexp(uint64_t x, uint8_t l)
    {
        return x>>(3*l);
    }

    static uint64_t modexp(uint64_t x, uint8_t l)
    {
        return x & bits::lo_set[3*l];
    }
};

template<>
struct precomp<16> {
    static uint64_t exp(uint8_t l)
    {
        return 1ULL<<(4*l);
    }

    static uint64_t divexp(uint64_t x, uint8_t l)
    {
        return x>>(4*l);
    }

    static uint64_t modexp(uint64_t x, uint8_t l)
    {
        return x & bits::lo_set[4*l];
    }
};


template<uint8_t t_k>
typename precomp<t_k>::impl precomp<t_k>::data;

typedef std::complex<uint64_t> t_p;
typedef t_p                    point_type;

struct node_type {
    uint8_t  t;   // level; size of node 1<<t
    t_p      p;   // lower left corner; upper left in matrix!
    uint64_t idx; // index in bp

    node_type() = default;
    node_type(uint8_t _t, t_p _p, uint64_t _idx) : t(_t), p(_p), idx(_idx)
    {}
    node_type(node_type&&) = default;
    node_type(const node_type&) = default;
    node_type& operator=(node_type&&) = default;
    node_type& operator=(const node_type&) = default;
};

} // end namepsace k2_treap_ns

} // end nomespace sdsl
#endif
