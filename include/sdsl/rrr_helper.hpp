/* sdsl - succinct data structures library
    Copyright (C) 2011 Simon Gog

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
/*! \file rrr_helper.hpp
   \brief rrr_helper.hpp contains the sdsl::binomial class,
          a class which contains informations about the binomial coefficients
   \author Simon Gog, Stefan Arnold
*/
#ifndef SDSL_RRR_HELPER
#define SDSL_RRR_HELPER

#include <algorithm> // for next permutation

namespace sdsl
{

// Helper class
template<uint8_t n>
class binomial
{
    private:

        static class impl
        {
            public:
                static const int MAX_SIZE=32;
                uint8_t m_space_for_bt[n+1];
                uint8_t m_space_for_bt_pair[256*(n==15)];
                uint64_t m_C[MAX_SIZE];
                int_vector<32> m_nr_to_bin;
                int_vector<32> m_bin_to_nr;

                impl() {
                    m_nr_to_bin.resize(1<<n);
                    m_bin_to_nr.resize(1<<n);
                    for (int i=0, cnt=0, class_cnt=0; i<=n; ++i) {
                        m_C[i] = cnt;
                        class_cnt = 0;
                        std::vector<bool> b(n,0);
                        for (int j=0; j<i; ++j) b[n-j-1] = 1;
                        do {
                            uint32_t x=0;
                            for (int k=0; k<n; ++k)
                                x |= ((uint32_t)b[n-k-1])<<(n-1-k);
                            m_nr_to_bin[cnt] = x;
                            m_bin_to_nr[x] = class_cnt;
                            ++cnt;
                            ++class_cnt;
                        } while (next_permutation(b.begin(), b.end()));
                        if (class_cnt == 1)
                            m_space_for_bt[i] = 0;
                        else
                            m_space_for_bt[i] = bit_magic::l1BP(class_cnt)+1;
                        //			cout << cnt << " " << class_cnt << " " << (int)m_space_for_bt[i] << endl;
                    }
                    if (n == 15) {
                        for (int x=0; x<256; ++x) {
                            m_space_for_bt_pair[x] = m_space_for_bt[x>>4] + m_space_for_bt[x&0x0F];
                        }
                    }
                }
        } iii;

    public:

        static inline uint8_t space_for_bt(uint32_t i) {
            return iii.m_space_for_bt[i];
        }

        static inline uint32_t nr_to_bin(uint8_t k, uint32_t nr) {
            return iii.m_nr_to_bin[iii.m_C[k]+nr];
        }

        static inline uint32_t bin_to_nr(uint32_t bin) {
            return iii.m_bin_to_nr[bin];
        }

        static inline uint8_t space_for_bt_pair(uint8_t x) {
            return iii.m_space_for_bt_pair[x];
        }
};
template<uint8_t n>
typename binomial<n>::impl binomial<n>::iii;


// Second helper class.
// Idea based on the paper
// Gonzalo Navarro and Eliana Providel: Fast, Small, Simple Rank/Select on Bitmaps, SEA 2012
// TODO: Further improve by precalc decoding of blocks with small k
template<uint8_t n>
class binomial2
{
    private:

        static class impl
        {
            public:
                static const int MAX_SIZE=64;
                uint8_t m_space_for_bt[n+1];
                uint64_t m_coefficients[MAX_SIZE][MAX_SIZE]; // m_coefficient[n][k] stores /n
                //                            \k/
                // TODO: this table is the same for all possible n
                // the different binomial2 classes should share it
                impl() {
                    for (int k=0; k<MAX_SIZE; ++k) {
                        m_coefficients[0][k] = 0;
                    }
                    for (int nn=0; nn<MAX_SIZE; ++nn) {
                        m_coefficients[nn][0] = 1;
                    }
                    for (int nn=1; nn<MAX_SIZE; ++nn) {
                        for (int k=1; k<MAX_SIZE; ++k) {
                            m_coefficients[nn][k] = m_coefficients[nn-1][k-1] + m_coefficients[nn-1][k];
                            // check overflow
                            if (m_coefficients[nn][k] < m_coefficients[nn-1][k-1] or
                                m_coefficients[nn][k] < m_coefficients[nn-1][k-1]) {
                                std::cout << "ERROR: for nn="<< nn <<" k="<< k << std::endl;
                            }
                        }
                    }
                    for (uint8_t k=0; k<=n; ++k) {
                        if (m_coefficients[n][k] == 1) {
                            m_space_for_bt[k] = 0;
                        } else {
                            m_space_for_bt[k] = bit_magic::l1BP(m_coefficients[n][k])+1;
                        }
//					std::cout<<"m_space_for_bt["<<(int)k<<"]="<<(int)m_space_for_bt[k]<<" m_coefficients[n][k]="<<m_coefficients[n][k]<<std::endl;
                    }
                }
        } iii;

    public:

        static inline uint8_t space_for_bt(uint64_t i) {
            return iii.m_space_for_bt[i];
        }

        static inline uint64_t nr_to_bin(uint8_t k, uint64_t nr) {
            if (k == n) {
                return bit_magic::Li1Mask[n];
            } else if (k == 0) {
                return 0;
            } else if (k == 1) {
                return 1ULL<<(n-1-nr);
            }
            uint64_t bin = 0;
            uint64_t mask = 1;
            uint8_t nn = n;
            while (k > 1) {
//		while ( k  ){
                if (nr >= iii.m_coefficients[nn-1][k]) {
                    nr -= iii.m_coefficients[nn-1][k];
                    --k;
                    bin |= mask;
                }
                --nn;
                mask <<= 1;
            }
            // now: k == 1
            return bin | (1ULL<<(n-1-nr));
//		return bin;
        }

        static inline uint64_t bin_to_nr(uint64_t bin) {
            if (bin == 0 or bin == bit_magic::Li1Mask[n]) {  // handle special cases
                return 0;
            }
            uint64_t nr = 0;
            uint8_t  k  = bit_magic::b1Cnt(bin); // get number of ones
            uint8_t  nn = n; // size of the block
            while (bin) {
                if (bin&1) {
                    nr += iii.m_coefficients[nn-1][k];
                    --k; // go to the case (n-1, k-1)
                }// else go to the case (n-1, k)
                bin >>= 1;
                --nn;
            }
            return nr;
        }

        static inline uint8_t space_for_bt_pair(uint8_t x) {
            return 0;//iii.m_space_for_bt_pair[x];
        }
};
template<uint8_t n>
typename binomial2<n>::impl binomial2<n>::iii;

typedef unsigned int uint128_t __attribute__((mode(TI)));

template<uint8_t n>
class binomial3
{
    private:

        static class impl
        {
            public:
                static const int MAX_SIZE=128;
                uint8_t m_space_for_bt[n+1];
                uint128_t m_coefficients[MAX_SIZE][MAX_SIZE]; // m_coefficient[n][k] stores /n
                //                            \k/
                // TODO: this table is the same for all possible n
                // the different binomial2 classes should share it
                impl() {
                    for (int k=0; k<MAX_SIZE; ++k) {
                        m_coefficients[0][k] = 0;
                    }
                    for (int nn=0; nn<MAX_SIZE; ++nn) {
                        m_coefficients[nn][0] = 1;
                    }
                    for (int nn=1; nn<MAX_SIZE; ++nn) {
                        for (int k=1; k<MAX_SIZE; ++k) {
                            m_coefficients[nn][k] = m_coefficients[nn-1][k-1] + m_coefficients[nn-1][k];
                        }
                    }
                    for (uint8_t k=0; k<=n; ++k) {
                        uint64_t high = (m_coefficients[n][k] >> 64);
                        uint64_t low = m_coefficients[n][k];
                        if (high == 0 && low == 1) {
                            m_space_for_bt[k] = 0;
                        } else {
                            if (high == 0) m_space_for_bt[k] = bit_magic::l1BP(low)+1;
                            else m_space_for_bt[k] = bit_magic::l1BP(high) + 64 + 1;
                        }
                    }
                }
        } iii;

    public:

        static inline uint8_t space_for_bt(uint64_t i) {
            return iii.m_space_for_bt[i];
        }

        static inline const uint128_t nr_to_bin(uint8_t k,uint128_t& nr) {
            if (k == n) {
                return ((uint128_t)bit_magic::Li1Mask[n-64]<<64) + bit_magic::Li1Mask[64];
            } else if (k == 0) {
                return 0;
            } else if (k == 1) {
                return ((uint128_t)1ULL<<(n-nr-1));
            }

            uint128_t bin = 0;
            uint128_t mask = 1;
            uint8_t nn = n;

            while (k > 1) {
                /*
                print_m128a((__m128i*)&mask,"mask");
                print_m128a((__m128i*)&nr,"nr");
                print_m128a((__m128i*)&bin,"bin");*/
                if (nr >= iii.m_coefficients[nn-1][k]) {
                    nr -= iii.m_coefficients[nn-1][k];
                    --k;
                    bin |= mask;
                }
                --nn;
                mask = (mask << 1ULL);
            }
            /* k == 1 */
            return bin | ((uint128_t)1<<(n-nr-1));
        };

        static inline uint128_t bin_to_nr(uint128_t& bin) {
            if (bin == 0 or bin == (((uint128_t)bit_magic::Li1Mask[n-64]<<64) +
                                    bit_magic::Li1Mask[64])) {   // handle special cases
                return 0;
            }
            uint128_t nr = 0;
            uint64_t high = (bin>>64);
            uint64_t low = bin;
            uint8_t  k  = bit_magic::b1Cnt(high) + bit_magic::b1Cnt(low); // get number of ones
            uint8_t  nn = n; // size of the block
            while (bin) {
                if (bin&1) {
                    nr += iii.m_coefficients[nn-1][k];
                    --k; // go to the case (n-1, k-1)
                }// else go to the case (n-1, k)
                bin >>= 1;
                --nn;
            }
            return nr;
        }

        static inline uint8_t space_for_bt_pair(uint8_t x) {
            return 0;//iii.m_space_for_bt_pair[x];
        }
};

template<uint8_t n>
typename binomial3<n>::impl binomial3<n>::iii;


} // end namespace
#endif
