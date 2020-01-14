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
        friend std::ostream& operator<<(std::ostream&, const uint256_t&);
    private:
        uint128_t m_lo;
        uint128_t m_high;

    public:
        inline uint256_t() : m_lo(0), m_high(0) {}

        inline uint256_t(const uint256_t& x) : m_lo(x.m_lo), m_high(x.m_high) {}
        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t(T lo) : m_lo(lo), m_high(0) {}

        inline uint256_t(uint128_t lo, uint128_t high) : m_lo(lo), m_high(high) {}

        inline uint256_t(uint64_t lo, uint64_t mid, uint128_t high)
#ifndef MODE_TI
            : m_lo(lo, mid), m_high(high) {}
#else
            : m_lo((uint128_t(mid) << 64) | lo), m_high(high) {}
#endif

        inline uint256_t& operator=(const uint256_t& x) { m_lo = x.m_lo; m_high = x.m_high; return *this; }
        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t& operator=(T lo) { m_lo = lo; m_high = 0; return *this; }

        inline uint16_t popcount() const
        {
#ifndef MODE_TI
            return m_lo.popcount() + m_high.popcount();
#else
            return bits::cnt(static_cast<uint64_t>(m_lo))
                    + bits::cnt(static_cast<uint64_t>(m_lo >> 64))
                    + bits::cnt(static_cast<uint64_t>(m_high))
                    + bits::cnt(static_cast<uint64_t>(m_high >> 64));
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

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t operator&(T x) const
        {
            return { m_lo & x, uint128_t(0) };
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t& operator<<=(T x)
        {
            if (x >= 128) {
                m_high = m_lo << (x - 128);
                m_lo = 0;
            } else {
                // FYI: avoids UB (shifting by the word size)
                m_high = (m_high << x) | ((m_lo >> (127 - x)) >> 1);
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
                // FYI: avoids UB (shifting by the word size)
                return { m_lo << x, (m_high << x) | ((m_lo >> (127 - x)) >> 1) };
            }
        }

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline uint256_t& operator>>=(T x)
        {
            if (x >= 128) {
                m_lo = m_high >> (x - 128);
                m_high = 0;
            } else {
                // FYI: avoids UB (shifting by the word size)
                m_lo = (m_lo >> x) | ((m_high << (127 - x)) << 1);
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
                // FYI: avoids UB (shifting by the word size)
                return { (m_lo >> x) | ((m_high << (127 - x)) << 1), m_high >> x };
            }
        }

        inline bool operator==(const uint256_t& x) const
        {
            return (m_high == x.m_high) and (m_lo == x.m_lo);
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

        template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
        inline bool operator<(T x) const
        {
            return !m_high && m_lo < x;
        }

        inline operator bool() const { return m_lo || m_high; }
        inline operator uint8_t() const { return static_cast<uint8_t>(m_lo); }
        inline operator uint16_t() const { return static_cast<uint16_t>(m_lo); }
        inline operator uint32_t() const { return static_cast<uint32_t>(m_lo); }
        inline operator uint64_t() const { return static_cast<uint64_t>(m_lo); }
        inline operator uint128_t() const { return m_lo; }
};

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator==(T number, const uint256_t& x) { return x.operator==(number); }

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator!=(T number, const uint256_t& x) { return x.operator!=(number); }

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator<(T number, const uint256_t& x) { return x.operator>(number); }

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator<=(T number, const uint256_t& x) { return x.operator>=(number); }

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator>(T number, const uint256_t& x) { return x.operator<(number); }

template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, T>::type >
inline bool operator>=(T number, const uint256_t& x) { return x.operator<=(number); }

std::ostream& operator<<(std::ostream& os, const uint256_t& x);

} // end namespace

namespace std {
    template<> struct is_arithmetic<sdsl::uint256_t> : ::std::true_type {};
    template<> struct is_integral<sdsl::uint256_t> : ::std::true_type {};
    template<> struct is_unsigned<sdsl::uint256_t> : ::std::true_type {};
    template<> struct make_unsigned<sdsl::uint256_t> { typedef sdsl::uint256_t type; };
} // end namespace

#endif
