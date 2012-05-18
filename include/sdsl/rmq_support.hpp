/* sdsl - succinct data structures library
    Copyright (C) 2009 Simon Gog

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
/*! \file rmq_support.hpp
    \brief rmq_support.hpp contains different range minimum support data structures.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_RMQ_SUPPORT
#define INCLUDED_SDSL_RMQ_SUPPORT

/** \defgroup rmq_group Range Minimum/Maximum Support (RMS) */

template<class RandomAccessContainer, bool Minimum>	 // for range minimum queries
struct min_max_trait {
    static inline bool strict_compare(const typename RandomAccessContainer::value_type v1, const typename RandomAccessContainer::value_type v2) {
        return v1 < v2;
    }
    static inline bool compare(const typename RandomAccessContainer::value_type v1, const typename RandomAccessContainer::value_type v2) {
        return v1 <= v2;
    }
};

template<class RandomAccessContainer> // for range maximum queries
struct min_max_trait<RandomAccessContainer, false> {
    static inline bool strict_compare(const typename RandomAccessContainer::value_type v1, const typename RandomAccessContainer::value_type v2) {
        return v1 > v2;
    }
    static inline bool compare(const typename RandomAccessContainer::value_type v1, const typename RandomAccessContainer::value_type v2) {
        return v1 >= v2;
    }
};

#include "rmq_support_sparse_table.hpp"
#include "rmq_succinct_sct.hpp"
#include "rmq_succinct_sada.hpp"

#endif
