/* sdsl - succinct data structures library
    Copyright (C) 2012 Simon Gog

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
/*! \file uint256_t.hpp
   \brief uint256_t.hpp contains a class for 256-bit unsigned integers.
   \author Simon Gog
*/
#ifndef INCLUDED_SDSL_UINT256
#define INCLUDED_SDSL_UINT256

#include <iostream>
#include "bits.hpp"
#include "uint128_t.hpp"

namespace sdsl
{

class uint256_t
{
    public:
        friend std::ostream& operator << (std::ostream&, const uint256_t&);
    private:
        uint128_t m_lo;
        uint128_t m_high;

    public:
        uint256_t() : m_lo(0), m_high(0) {}

        inline uint256_t(const uint128_t& lo, const uint128_t& high = 0)
            : m_lo(lo), m_high(high) {}

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t(T lo, T mid = 0, uint128_t high = 0)
#ifndef MODE_TI
            : m_lo(lo, mid), m_high(high) {}
#else
            : m_lo((uint128_t(mid) << 64) | lo), m_high(high) {}
#endif

        inline uint16_t popcount() const
        {
#ifndef MODE_TI
            return m_lo.popcount() + m_high.popcount();
#else
            return bits::cnt(m_lo) + bits::cnt(m_lo >> 64)
                    + bits::cnt(m_high) + bits::cnt(m_high >> 64);
#endif
        }

        inline uint16_t hi() const
        {
#ifndef MODE_TI
            if (m_high) {
                return m_high.hi() + 128;
            } else {
                return m_lo.hi();
            }
#else
            if (m_high) {
                uint64_t hh = (m_high >> 64);
                if (hh) {
                    return bits::hi(hh) + 192;
                } else {
                    return bits::hi(m_high) + 128;
                }
            } else {
                uint64_t hh = (m_lo >> 64);
                if (hh) {
                    return bits::hi(hh) + 64;
                } else {
                    return bits::hi(m_lo);
                }
            }
#endif
        }

        inline uint16_t select(uint32_t i) const
        {
#ifndef MODE_TI
            uint16_t x = 0;
            if ((x = m_lo.popcount()) >= i) {
                return m_lo.select(i);
            }
            i -= x;
            return m_high.select(i) + 128;
#else
            uint16_t x = 0;
            uint64_t v = m_lo;
            if ((x = bits::cnt(v)) >= i)
                return bits::sel(v, i);

            i -= x;
            v = m_lo >> 64;
            if ((x = bits::cnt(v)) >= i) {
                return bits::sel(v, i) + 64;
            }
            i -= x;
            v = m_high;
            if ((x= bits::cnt(v)) >= i) {
                return bits::sel(v, i) + 128;
            }
            i -= x;
            return bits::sel(m_high >> 64, i) + 192;
#endif
        }

        inline uint256_t& operator+=(const uint256_t& x)
        {
            m_high += x.m_high + ((m_lo + x.m_lo) < m_lo);
            m_lo += x.m_lo;
            return *this;
        }

        inline uint256_t& operator+=(const uint128_t& x)
        {
            m_high += ((m_lo + x) < m_lo);
            m_lo += x;
            return *this;
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t& operator+=(T x)
        {
            m_high += ((m_lo + x) < m_lo);
            m_lo += x;
            return *this;
        }

        inline uint256_t operator+(const uint256_t& x) const
        {
            return { m_lo + x.m_lo,
                m_high + x.m_high + ((m_lo + x.m_lo) < m_lo) };
        }

        inline uint256_t operator+(const uint128_t& x) const
        {
            return { m_lo + x, m_high + ((m_lo + x) < m_lo) };
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t operator+(T x) const
        {
            return { m_lo + x, m_high + ((m_lo + x) < m_lo) };
        }

        inline uint256_t& operator-=(const uint256_t& x)
        {
            m_high -= x.m_high + (m_lo < x.m_lo);
            m_lo -= x.m_lo;
            return *this;
        }

        inline uint256_t& operator-=(const uint128_t& x)
        {
            m_high -= (m_lo < x);
            m_lo -= x;
            return *this;
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t& operator-=(T x)
        {
            m_high -= (m_lo < x);
            m_lo -= x;
            return *this;
        }

        inline uint256_t operator-(const uint256_t& x) const
        {
            return { m_lo - x.m_lo,
                m_high - x.m_high - (m_lo < x.m_lo) };
        }

        inline uint256_t operator-(const uint128_t& x) const
        {
            return { m_lo - x, m_high - (m_lo < x) };
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t operator-(T x) const
        {
            return { m_lo - x, m_high - (m_lo < x) };
        }

        inline uint256_t operator~() const
        {
            return { ~m_lo, ~m_high };
        }

        inline uint256_t& operator|=(const uint256_t& x)
        {
            m_lo |= x.m_lo;
            m_high |= x.m_high;
            return *this;
        }

        inline uint256_t& operator|=(const uint128_t& x)
        {
            m_lo |= x;
            return *this;
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t& operator|=(T x)
        {
            m_lo |= x;
            return *this;
        }

        inline uint256_t operator|(const uint256_t& x) const
        {
            return { m_lo | x.m_lo, m_high | x.m_high };
        }

        inline uint256_t operator|(const uint128_t& x) const
        {
            return { m_lo | x, m_high };
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t operator|(T x) const
        {
            return { m_lo | x, m_high };
        }

        inline uint256_t& operator&=(const uint256_t& x)
        {
            m_lo &= x.m_lo;
            m_high &= x.m_high;
            return *this;
        }

        inline uint256_t& operator&=(const uint128_t& x)
        {
            m_lo &= x;
            m_high = 0ULL;
            return *this;
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t& operator&=(T x)
        {
            m_lo &= x;
            m_high = 0ULL;
            return *this;
        }

        inline uint256_t operator&(const uint256_t& x) const
        {
            return { m_lo & x.m_lo, m_high & x.m_high };
        }

        inline uint256_t operator&(const uint128_t& x) const
        {
            return { m_lo & x, uint128_t(0) };
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t operator&(T x) const
        {
            return { m_lo & x, uint128_t(0) };
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t& operator<<=(T x)
        {
            if (x >= 128) {
                m_lo <<= (x - 128);
                m_high = m_lo;
                m_lo = 0;
            } else {
                m_high <<= x;
                m_high |= (m_lo >> (128 - x));
                m_lo <<= x;
            }
            return *this;
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t operator<<(T x) const
        {
            if (x >= 128) {
                return { uint128_t(0), m_lo << (x - 128),  };
            } else {
                return { m_lo << x, (m_high << x) | (m_lo >> (128 - x)) };
            }
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t& operator>>=(T x)
        {
            if (x >= 128) {
                m_high >>= (x - 128);
                m_lo = m_high;
                m_high = 0;
            } else {
                m_lo >>= x;
                m_lo |= (m_high << (128 - x));
                m_high >>= x;
            }
            return *this;
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t operator>>(T x) const
        {
            if (x >= 128) {
                return { m_high >> (x - 128), uint128_t(0) };
            } else {
                return { (m_lo >> x) | (m_high << (128 - x)), m_high >> x };
            }
        }

        inline bool operator==(const uint256_t& x) const
        {
            return (m_high == x.m_high) and (m_lo == x.m_lo);
        }

        inline bool operator==(const uint128_t& x) const
        {
            return !m_high and (m_lo == x);
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline bool operator==(T x) const
        {
            return !m_high and (m_lo == x);
        }

        inline bool operator!=(const uint256_t& x) const
        {
            return (m_high != x.m_high) or (m_lo != x.m_lo);
        }

        inline bool operator!=(const uint128_t& x) const
        {
            return m_high or (m_lo != x);
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline bool operator!=(T x) const
        {
            return m_high or (m_lo != x);
        }

        inline bool operator>=(const uint256_t& x) const
        {
            if (m_high != x.m_high) {
                return m_high > x.m_high;
            } else {
                return m_lo >= x.m_lo;
            }
        }

        inline bool operator>=(const uint128_t& x) const
        {
            return m_high || m_lo >= x;
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline bool operator>=(T x) const
        {
            return m_high || m_lo >= x;
        }

        inline bool operator<=(const uint256_t& x) const
        {
            if (m_high != x.m_high) {
                return m_high < x.m_high;
            } else {
                return m_lo <= x.m_lo;
            }
        }

        inline bool operator<=(const uint128_t& x) const
        {
            return !m_high && m_lo <= x;
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline bool operator<=(T x) const
        {
            return !m_high && m_lo <= x;
        }

        inline bool operator>(const uint256_t& x) const
        {
            if (m_high != x.m_high) {
                return m_high > x.m_high;
            } else {
                return m_lo > x.m_lo;
            }
        }

        inline bool operator>(const uint128_t& x) const
        {
            return m_high || m_lo > x;
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline bool operator>(T x) const
        {
            return m_high || m_lo > x;
        }

        inline bool operator<(const uint256_t& x) const
        {
            if (m_high != x.m_high) {
                return m_high < x.m_high;
            } else {
                return m_lo < x.m_lo;
            }
        }

        inline bool operator<(const uint128_t& x) const
        {
            return !m_high && m_lo < x;
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline bool operator<(T x) const
        {
            return !m_high && m_lo < x;
        }

        inline operator bool() const { return static_cast<bool>(m_lo) || static_cast<bool>(m_high); }
        inline operator uint8_t() const { return static_cast<uint8_t>(m_lo); }
        inline operator uint16_t() const { return static_cast<uint16_t>(m_lo); }
        inline operator uint32_t() const { return static_cast<uint32_t>(m_lo); }
        inline operator uint64_t() const { return static_cast<uint64_t>(m_lo); }
        inline operator uint128_t() const { return m_lo; }
};

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator==(T number, const uint256_t& x) { return x.operator==(number); }
inline bool operator==(const uint128_t& number, const uint256_t& x) { return x.operator==(number); }

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator!=(T number, const uint256_t& x) { return x.operator!=(number); }
inline bool operator!=(const uint128_t& number, const uint256_t& x) { return x.operator!=(number); }

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator<(T number, const uint256_t& x) { return x.operator>(number); }
inline bool operator<(const uint128_t& number, const uint256_t& x) { return x.operator>(number); }

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator<=(T number, const uint256_t& x) { return x.operator>=(number); }
inline bool operator<=(const uint128_t& number, const uint256_t& x) { return x.operator>=(number); }

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator>(T number, const uint256_t& x) { return x.operator<(number); }
inline bool operator>(const uint128_t& number, const uint256_t& x) { return x.operator<(number); }

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator>=(T number, const uint256_t& x) { return x.operator<=(number); }
inline bool operator>=(const uint128_t& number, const uint256_t& x) { return x.operator<=(number); }

std::ostream& operator<<(std::ostream& os, const uint256_t& x);

} // end namespace

#endif
