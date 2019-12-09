/* sdsl - succinct data structures library
    Copyright (C) 2012 Simon Gog, Matthias Petri

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
/*! \file uint128_t.hpp
   \brief uint128_t.hpp contains contains the definition of a 128-bit unsigned integer type.
   \author Simon Gog, Matthias Petri
*/
#ifndef INCLUDED_SDSL_UINT128
#define INCLUDED_SDSL_UINT128

#include <iostream>
#include "bits.hpp"

namespace sdsl
{

#ifdef MODE_TI
typedef unsigned int uint128_t __attribute__((mode(TI)));
#else

class uint128_t
{
    public:
        friend std::ostream& operator << (std::ostream&, const uint128_t&);
    private:
        uint64_t m_lo;
        uint64_t m_high;

    public:
        inline uint128_t(uint64_t lo = 0, uint64_t high = 0) : m_lo(lo), m_high(high) {}

        inline uint8_t popcount() const
        {
            return bits::cnt(m_lo) + bits::cnt(m_high);
        }

        inline uint16_t hi() const
        {
            if (m_high) {
                return bits::hi(m_high) + 64;
            } else {
                return bits::hi(m_lo);
            }
        }

        inline uint16_t select(uint32_t i) const
        {
            uint16_t x = 0;
            if ((x = bits::cnt(m_lo)) >= i) {
                return bits::sel(m_lo, i);
            }
            i -= x;
            return bits::sel(m_high, i) + 64;
        }

        inline uint128_t& operator+=(const uint128_t& x)
        {
            m_high += x.m_high + ((m_lo + x.m_lo) < m_lo);
            m_lo += x.m_lo;
            return *this;
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint128_t& operator+=(T x)
        {
            m_high += (m_lo + x) < m_lo;
            m_lo += x;
            return *this;
        }

        inline uint128_t operator+(const uint128_t& x) const
        {
            return { m_lo + x.m_lo, m_high + x.m_high + ((m_lo + x.m_lo) < m_lo) };
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint128_t operator+(T x) const
        {
            return { m_lo + x, m_high + ((m_lo + x) < m_lo) };
        }

        inline uint128_t& operator-=(const uint128_t& x)
        {
            m_high -= x.m_high + (m_lo < x.m_lo);
            m_lo -= x.m_lo;
            return *this;
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint128_t& operator-=(T x)
        {
            m_high -= (m_lo < x);
            m_lo -= x;
            return *this;
        }

        inline uint128_t operator-(const uint128_t& x) const
        {
            return { m_lo - x.m_lo, m_high - x.m_high - (m_lo < x.m_lo) };
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint128_t operator-(T x) const
        {
            return { m_lo - x, m_high - (m_lo < x) };
        }

        inline uint128_t operator~() const
        {
            return { ~m_lo, ~m_high };
        }

        inline uint128_t& operator|=(const uint128_t& x)
        {
            m_lo |= x.m_lo;
            m_high |= x.m_high;
            return *this;
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint128_t& operator|=(T x)
        {
            m_lo |= x;
            return *this;
        }

        inline uint128_t operator|(const uint128_t& x) const
        {
            return { m_lo | x.m_lo, m_high | x.m_high };
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint128_t operator|(T x) const
        {
            return { m_lo | x, m_high };
        }

        inline uint128_t& operator&=(const uint128_t& x)
        {
            m_lo &= x.m_lo;
            m_high &= x.m_high;
            return *this;
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint128_t& operator&=(T x)
        {
            m_lo &= x;
            m_high = 0ULL;
            return *this;
        }

        inline uint128_t operator&(const uint128_t& x) const
        {
            return { m_lo & x.m_lo, m_high & x.m_high };
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint128_t operator&(T x) const
        {
            return { m_lo & x, 0ULL };
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint128_t& operator<<=(T x)
        {
            if (x >= 64) {
                m_high = x >= 128 ? 0ULL : m_lo << (x - 64);
                m_lo = 0;
            } else {
                // FYI: avoids UB (shifting by the word size)
                m_high = (m_high << x) | ((m_lo >> (63 - x)) >> 1);
                m_lo <<= x;
            }
            return *this;
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint128_t operator<<(T x) const
        {
            if (x >= 64) {
                return { 0ULL, x >= 128 ? 0ULL : m_lo << (x - 64) };
            } else {
                // FYI: avoids UB (shifting by the word size)
                return { m_lo << x, (m_high << x) | ((m_lo >> (63 - x)) >> 1) };
            }
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint128_t& operator>>=(T x)
        {
            if (x >= 64) {
                m_lo = x >= 128 ? 0ULL : m_high >> (x - 64);
                m_high = 0;
            } else {
                // FYI: avoids UB (shifting by the word size)
                m_lo = (m_lo >> x) | ((m_high << (63 - x)) << 1);
                m_high >>= x;
            }
            return *this;
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint128_t operator>>(T x) const
        {
            if (x >= 64) {
                return { x >= 128 ? 0ULL : m_high >> (x - 64), 0ULL };
            } else {
                // FYI: avoids UB (shifting by the word size)
                return { (m_lo >> x) | ((m_high << (63 - x)) << 1), m_high >> x };
            }
        }

        inline bool operator==(const uint128_t& x) const
        {
            return (m_high == x.m_high) and (m_lo == x.m_lo);
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline bool operator==(T x) const
        {
            return !m_high and (m_lo == x);
        }

        inline bool operator!=(const uint128_t& x) const
        {
            return (m_high != x.m_high) or (m_lo != x.m_lo);
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline bool operator!=(T x) const
        {
            return m_high or (m_lo != x);
        }

        inline bool operator>=(const uint128_t& x) const
        {
            if (m_high != x.m_high) {
                return m_high > x.m_high;
            } else {
                return m_lo >= x.m_lo;
            }
        }

        inline bool operator<=(const uint128_t& x) const
        {
            if (m_high != x.m_high) {
                return m_high < x.m_high;
            } else {
                return m_lo <= x.m_lo;
            }
        }

        inline bool operator>(const uint128_t& x) const
        {
            if (m_high != x.m_high) {
                return m_high > x.m_high;
            } else {
                return m_lo > x.m_lo;
            }
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline bool operator>(T x) const
        {
            return m_high || m_lo > x;
        }

        inline bool operator<(const uint128_t& x) const
        {
            if (m_high != x.m_high) {
                return m_high < x.m_high;
            } else {
                return m_lo < x.m_lo;
            }
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline bool operator<(T x) const
        {
            return !m_high && m_lo < x;
        }

        inline operator bool() const { return static_cast<bool>(m_lo | m_high); }
        inline operator uint8_t() const { return static_cast<uint8_t>(m_lo); }
        inline operator uint16_t() const { return static_cast<uint8_t>(m_lo); }
        inline operator uint32_t() const { return static_cast<uint8_t>(m_lo); }
        inline operator uint64_t() const { return m_lo; }
};

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator==(T number, const uint128_t& x) { return x.operator==(number); }

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator!=(T number, const uint128_t& x) { return x.operator!=(number); }

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator<(T number, const uint128_t& x) { return x.operator<(number); }

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator<=(T number, const uint128_t& x) { return x.operator<=(number); }

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator>(T number, const uint128_t& x) { return x.operator>(number); }

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator>=(T number, const uint128_t& x) { return x.operator>=(number); }

#endif

std::ostream& operator<<(std::ostream& os, const uint128_t& x);

} // end namespace

#endif
