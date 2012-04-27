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
   \author Simon Gog, Matthias Petri, Stefan Arnold
*/
#ifndef SDSL_RRR_HELPER
#define SDSL_RRR_HELPER

#include <algorithm> // for next permutation
#include <iostream>
#include "bitmagic.hpp"

namespace sdsl
{

// Helper class
template<uint8_t n>
class binomial
{
	public:
		typedef uint32_t number_type;
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

		static inline uint8_t popcount(number_type x) {
			return bit_magic::b1Cnt(x);
		}

		static inline uint8_t select(number_type x, uint32_t i) {
           	return bit_magic::i1BP(x, i);
		}

		template<class bit_vector_type>
		static inline number_type decode_btnr(const bit_vector_type &bv, 
											  typename bit_vector_type::size_type btnrp,
											  uint8_t btnrlen){
			return bv.get_int(btnrp, btnrlen);
		}

		template<class bit_vector_type>
		static inline uint8_t get_bt(const bit_vector_type &bv, 
									 typename bit_vector_type::size_type pos,
									 uint8_t block_size){
			return bit_magic::b1Cnt( bv.get_int(pos, block_size) );
		}

		template<class bit_vector_type>
		static void set_bt(bit_vector_type &bv, 
						   typename bit_vector_type::size_type pos,
						   number_type bt,
						   uint8_t space_for_bt){
			bv.set_int(pos, bt, space_for_bt);
		}

		static inline number_type Li1Mask(uint8_t off){
			return bit_magic::Li1Mask[off];
		}
};
template<uint8_t n>
typename binomial<n>::impl binomial<n>::iii;


// Second helper class.
// Idea based on the paper
// Gonzalo Navarro and Eliana Providel: Fast, Small, Simple Rank/Select on Bitmaps, SEA 2012
// TODO: Further improve by precalc decoding of blocks with small k
template<uint8_t n>
class binomial64
{
	public:
		typedef uint64_t number_type;
    private:

        static class impl
        {
            public:
                static const int MAX_SIZE=64;
                uint8_t m_space_for_bt[MAX_SIZE+1];
                uint64_t m_coefficients[MAX_SIZE+1][MAX_SIZE+1]; // m_coefficient[n][k] stores /n
                //                            \k/
                // TODO: this table is the same for all possible n
                // the different binomial2 classes should share it
                impl() {
                    for (int k=0; k<=MAX_SIZE; ++k) {
                        m_coefficients[0][k] = 0;
                    }
                    for (int nn=0; nn<=MAX_SIZE; ++nn) {
                        m_coefficients[nn][0] = 1;
                    }
                    for (int nn=1; nn<=MAX_SIZE; ++nn) {
                        for (int k=1; k<=MAX_SIZE; ++k) {
                            m_coefficients[nn][k] = m_coefficients[nn-1][k-1] + m_coefficients[nn-1][k];
                            // check overflow
                            if (m_coefficients[nn][k] < m_coefficients[nn-1][k-1] or
                                m_coefficients[nn][k] < m_coefficients[nn-1][k-1]) {
                                std::cout << "ERROR: for nn="<< nn <<" k="<< k << std::endl;
                            }
                        }
                    }
                    for (uint8_t k=0; k<=MAX_SIZE; ++k) {
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
            } else if (k == 1) { // optimization if only on bit is set
                return 1ULL<<(n-1-nr); // (n-1-nr) is always < 64 for n<=64
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

        static inline uint8_t space_for_bt_pair(uint8_t x) { return 0; }

		static inline uint8_t popcount(number_type x) {
			return bit_magic::b1Cnt(x);
		}

		static inline uint8_t select(number_type x, uint32_t i) {
           	return bit_magic::i1BP(x, i);
		}

		template<class bit_vector_type>
		static inline number_type decode_btnr(const bit_vector_type &bv, 
											  typename bit_vector_type::size_type btnrp,
											  uint8_t btnrlen){
			return bv.get_int(btnrp, btnrlen);
		}

		template<class bit_vector_type>
		static inline uint8_t get_bt(const bit_vector_type &bv, 
									 typename bit_vector_type::size_type pos,
									 uint8_t block_size){
			return bit_magic::b1Cnt( bv.get_int(pos, block_size) );
		}

		template<class bit_vector_type>
		static void set_bt(bit_vector_type &bv, 
						   typename bit_vector_type::size_type pos,
						   number_type bt,
						   uint8_t space_for_bt){
			bv.set_int(pos, bt, space_for_bt);
		}

		static inline number_type Li1Mask(uint8_t off){
			return bit_magic::Li1Mask[off];
		}
};


template<uint8_t n>
typename binomial64<n>::impl binomial64<n>::iii;

typedef unsigned int uint128_t __attribute__((mode(TI)));

std::ostream& operator<<(std::ostream &os, const uint128_t &x){
	uint64_t X[2] = {(uint64_t)(x >> 64), (uint64_t)x};
	for( int j=0; j < 2; ++j ){
		for( int i=0; i < 16; ++i ){
			os << std::hex << ((X[j]>>60)&0xFULL);
			X[j] <<= 4;
		}	
	}
	return os;
};

template<uint8_t n>
class binomial128
{
	public:
		typedef uint128_t number_type;
    private:

        static class impl
        {
            public:
                static const int MAX_SIZE=128;
                uint8_t m_space_for_bt[n+1];
                uint128_t m_coefficients[MAX_SIZE+1][MAX_SIZE+1]; // m_coefficient[n][k] stores /n

				uint128_t m_L1Mask[MAX_SIZE+1];
                //                            \k/
                // TODO: this table is the same for all possible n
                // the different binomial2 classes should share it
                impl() {
                    for (int k=0; k <= MAX_SIZE; ++k) {
                        m_coefficients[0][k] = 0;
                    }
                    for (int nn=0; nn <= MAX_SIZE; ++nn) {
                        m_coefficients[nn][0] = 1;
                    }
                    for (int nn=1; nn <= MAX_SIZE; ++nn) {
                        for (int k=1; k <= MAX_SIZE; ++k) {
                            m_coefficients[nn][k] = m_coefficients[nn-1][k-1] + m_coefficients[nn-1][k];
                        }
                    }
                    for (uint8_t k=0; k<=n; ++k) {
                        if (m_coefficients[n][k] == (uint128_t)1) {
                            m_space_for_bt[k] = 0;
                        } else {
							m_space_for_bt[k] = l1BP( m_coefficients[n][k] )+1;
                        }
                    }
					m_L1Mask[0] = 0;
					uint128_t mask = 1;
					for (int i=1; i<=128; ++i){
						m_L1Mask[i] = mask;
						mask <<= 1;
						mask |= (uint128_t)1;
					}
                }
        } iii;

    public:

        static inline uint8_t space_for_bt(uint64_t i) {
            return iii.m_space_for_bt[i];
        }

		// TODO: is a speed-up possible if we only decode as many bit as we need?
        static inline const uint128_t nr_to_bin(uint8_t k,uint128_t& nr) {
            if (k == n) {
                return iii.m_L1Mask[n]; 
            } else if (k == 0) {
                return 0;
            } else if (k == 1) {
                return ((uint128_t)1ULL<<(n-nr-1));
            }

            uint128_t bin = 0;
            uint128_t mask = 1;
            uint8_t nn = n;
            while (k > 1) {
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
            if (bin == 0 or bin == iii.m_L1Mask[n] ){ // handle special case
			   	//if ( bin == 0 or (((uint128_t)bit_magic::Li1Mask[n-64]<<64) + bit_magic::Li1Mask[64])) {   
                return 0;
            }
            uint128_t nr = 0;
            uint8_t  k  = popcount( bin ); //bit_magic::b1Cnt(high) + bit_magic::b1Cnt(low); // get number of ones
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

        static inline uint8_t space_for_bt_pair(uint8_t x) { return 0; }

		static inline uint8_t popcount(number_type x) {
			return bit_magic::b1Cnt(x >> 64) + bit_magic::b1Cnt(x);
		}

		static inline uint8_t l1BP(number_type x) {
			if ( (x >> 64) ){
				return bit_magic::l1BP( x >> 64 ) + 64;
			}else{
				return bit_magic::l1BP( x );
			}
		}

		static inline uint8_t select(number_type x, uint32_t i) {
			uint64_t low = x;
			uint64_t poplow = bit_magic::b1Cnt(low);
			if ( poplow >= i ){
            	return bit_magic::i1BP(low, i);
			}else{
				uint64_t high = x>>64;
            	return 64 + bit_magic::i1BP( high, i-poplow);
			}
		}

		template<class bit_vector_type>
		static inline number_type decode_btnr(const bit_vector_type &bv, 
											  typename bit_vector_type::size_type btnrp,
											  uint8_t btnrlen){
			if ( btnrlen <= 64 ){
				return bv.get_int(btnrp, btnrlen);
			}else{
                return ((((uint128_t) bv.get_int(btnrp+64, btnrlen-64))<<64) + bv.get_int(btnrp, 64));
			}
		}

		template<class bit_vector_type>
		static inline uint8_t get_bt(const bit_vector_type &bv, 
									 typename bit_vector_type::size_type pos,
									 uint8_t block_size){
			if ( block_size <= 64 ){
				return bit_magic::b1Cnt( bv.get_int(pos, block_size) );
			}else{
                return bit_magic::b1Cnt( bv.get_int(pos+64, block_size-64) )+ 
					   bit_magic::b1Cnt( bv.get_int(pos, 64) );
			}
		}

		template<class bit_vector_type>
		static void set_bt(bit_vector_type &bv, 
						   typename bit_vector_type::size_type pos,
						   number_type bt,
						   uint8_t space_for_bt){
			if ( space_for_bt <= 64 ){
				bv.set_int(pos, bt, space_for_bt);
			}else{
                bv.set_int(pos, (uint64_t)bt, 64);
			    bv.set_int(pos+64, bt>>64, space_for_bt-64 );
			}
		}

		static inline number_type Li1Mask(uint8_t off){
			return iii.m_L1Mask[off];
		}
};

template<uint8_t n>
typename binomial128<n>::impl binomial128<n>::iii;

class uint256_t{
	public:
		friend std::ostream& operator << (std::ostream &, const uint256_t&);
	private:
		uint64_t m_lo;
		uint64_t m_mid;
		uint128_t m_high;

	public:
		uint256_t(uint64_t lo=0, uint64_t mid=0, 
				  uint128_t high=0):m_lo(lo), 
									m_mid(mid),
									m_high(high){}

		uint256_t(const uint256_t &x):m_lo(x.m_lo), m_mid(x.m_mid), m_high(x.m_high){}

		uint16_t popcount(){
			return ((uint16_t)bit_magic::b1Cnt(m_lo)) + bit_magic::b1Cnt(m_mid)
				 + bit_magic::b1Cnt(m_high>>64) + bit_magic::b1Cnt(m_high);
		}

		uint16_t l1BP() {
			if ( m_high == 0 ){
				if ( m_mid ){
					return bit_magic::l1BP( m_mid ) + 64;
				}else{
					return bit_magic::l1BP( m_lo );
				}
			}else{
				uint64_t hh = (m_high >> 64);
				if ( hh ){
					return bit_magic::l1BP( hh ) + 192;
				}else{
					return bit_magic::l1BP( m_high ) + 128;
				}
			}
		}

		uint16_t select(uint32_t i) {
			uint16_t x = 0;
			if ( (x=bit_magic::b1Cnt(m_lo)) >= i ){
				return bit_magic::i1BP(m_lo, i);
			}
			i -= x;
			if ( (x=bit_magic::b1Cnt(m_mid)) >= i ){
				return bit_magic::i1BP( m_mid, i ) + 64;
			}
			i -= x;
			uint64_t hh = m_high >> 64;
			uint64_t lh = m_high;
			if ( (x=bit_magic::b1Cnt(lh)) >= i ){
				return bit_magic::i1BP( lh, i) + 128;
			}
			i -= x;
			return bit_magic::i1BP( hh, i) + 192;
		}

		uint256_t& operator+=(const uint256_t &x){
			uint128_t lo = (uint128_t)m_lo + x.m_lo;
			uint128_t mid = (uint128_t)m_mid + x.m_mid + (lo >> 64);
			m_lo = lo; m_mid = mid;
			m_high += x.m_high + (mid >> 64);
			return *this;
//			return uint256_t(lo, mid, m_high + x.m_high + (mid >> 64));
		}

		uint256_t operator+(uint256_t x){
			uint128_t lo = (uint128_t)m_lo + x.m_lo;
			uint128_t mid = (uint128_t)m_mid + x.m_mid + (lo >> 64);
			return uint256_t(lo, mid, m_high + x.m_high + (mid >> 64));
		}

		uint256_t operator-(uint256_t x){
//			add two's complement of x 
			uint128_t lo = (uint128_t)m_lo + (~x.m_lo) + 1;
			uint128_t mid = (uint128_t)m_mid + (~x.m_mid) + (lo >> 64);
			return uint256_t(lo, mid, m_high + (~x.m_high) + (mid >> 64));
		}

		uint256_t& operator-=(const uint256_t &x){
//			add two's complement of x 
			uint128_t lo = (uint128_t)m_lo + (~x.m_lo) + 1;
			uint128_t mid = (uint128_t)m_mid + (~x.m_mid) + (lo >> 64);
			m_lo = lo; 
			m_mid = mid;
			m_high += (~x.m_high) + (mid >> 64);
			return *this;
		}
		

		uint256_t operator|(const uint256_t &x){
			return uint256_t(m_lo|x.m_lo, m_mid|x.m_mid, m_high|x.m_high);
		}

		uint256_t& operator|=(const uint256_t &x){
			m_lo |= x.m_lo; m_mid |= x.m_mid; m_high |= x.m_high;
			return *this;
		}

		uint256_t operator&(const uint256_t &x){
			return uint256_t(m_lo&x.m_lo, m_mid&x.m_mid, m_high&x.m_high);
		}
/* // is not needed since we can convert uint256_t to uint64_t
		uint64_t operator&(uint64_t x){
			return m_lo & x;	
		}
*/		

		uint256_t operator<<(int x){
			if ( x < 128){
				uint128_t high = m_high << x;
				uint128_t low  = (((uint128_t)m_mid<<64) | m_lo);
				high |= (low >> (128-x));
				low = low << x;
				return uint256_t(low, low>>64, high );
			}else{ // x >= 128
				uint128_t high = (((uint128_t)m_mid<<64) | m_lo) << (x-128);
				return uint256_t(0, 0, high);
			}
		}

		uint256_t operator>>(int x){
			if ( x < 128 ){
				uint128_t low  = (((uint128_t)m_mid<<64) | m_lo) >> x;
				low |= ((m_high << (127-x))<<1);
				return uint256_t(low, low>>64, m_high>>x );
			}else{ // x >= 128
				uint128_t low = (m_high >> (x-128));
				return uint256_t(low, low>>64, 0);
			}
		}

		uint256_t& operator=(const uint64_t &x){
			m_high = 0;
			m_mid = 0;
			m_lo = x;
			return *this;
		}

		bool operator==(const uint256_t &x) const{
			return (m_lo == x.m_lo) and (m_mid == x.m_mid) and (m_high == x.m_high);
		}

		bool operator!=(const uint256_t &x) const{
			return !(*this == x);
		}

		bool operator>=(const uint256_t &x) const{
			if ( m_high != x.m_high ){
				return m_high > x.m_high;
			}
			if( m_mid != x.m_mid ){
				return m_mid > x.m_mid;
			}else{
				return m_lo >= x.m_lo;
			}
		}

		bool operator<=(const uint256_t &x) const{
			if ( m_high != x.m_high ){
				return m_high < x.m_high;
			}
			if( m_mid != x.m_mid ){
				return m_mid < x.m_mid;
			}else{
				return m_lo <= x.m_lo;
			}
		}

		bool operator>(const uint256_t &x) const{
			if ( m_high != x.m_high ){
				return m_high > x.m_high;
			}
			if( m_mid != x.m_mid ){
				return m_mid > x.m_mid;
			}else{
				return m_lo > x.m_lo;
			}
		}

		bool operator>(const uint64_t &x) const{
			if ( m_high > 0 or m_mid > 0 ){
				return true;
			}
			return m_lo > x;
		}

		bool operator<(const uint256_t &x) const{
			if ( m_high != x.m_high ){
				return m_high < x.m_high;
			}
			if( m_mid != x.m_mid ){
				return m_mid < x.m_mid;
			}else{
				return m_lo < x.m_lo;
			}
		}

		operator uint64_t(){
			return m_lo;
		}
};

std::ostream& operator<<(std::ostream &os, const uint256_t &x){
	uint64_t X[4] = {(uint64_t)(x.m_high >> 64), (uint64_t)x.m_high, x.m_mid, x.m_lo};
	for( int j=0; j < 4; ++j ){
		for( int i=0; i < 16; ++i ){
			os << std::hex << ((X[j]>>60)&0xFULL) << std::dec;
			X[j] <<= 4;
		}	
	}
	return os;
};


template<uint16_t n>
class binomial256
{
	public:
		typedef uint256_t number_type;
    private:

        static class impl
        {
            public:
                static const int MAX_SIZE=255;
                uint8_t m_space_for_bt[n+1];
                uint256_t m_coefficients[MAX_SIZE+1][MAX_SIZE+1]; 

				uint256_t m_L1Mask[MAX_SIZE+1];
                impl() {
                    for (int k=0; k <= MAX_SIZE; ++k) {
                        m_coefficients[0][k] = uint256_t(0);
                    }
                    for (int nn=0; nn <= MAX_SIZE; ++nn) {
                        m_coefficients[nn][0] = uint256_t(1);
                    }
                    for (int nn=1; nn <= MAX_SIZE; ++nn) {
                        for (int k=1; k <= MAX_SIZE; ++k) {
                            m_coefficients[nn][k] = (m_coefficients[nn-1][k-1]) + (m_coefficients[nn-1][k]);
                        }
                    }
                    for (int k=0; k<=n; ++k) {
                        if (m_coefficients[n][k] == uint256_t(1) ) {
                            m_space_for_bt[k] = 0;
                        } else {
							m_space_for_bt[k] = l1BP( m_coefficients[n][k] )+1;
                        }
                    }
					m_L1Mask[0] = uint256_t(0);
					uint256_t mask = uint256_t(1);
					for (int i=1; i<=256; ++i){
						m_L1Mask[i] = mask;
						mask =  (mask << 1);
						mask |= uint256_t(1);
					}
                }
        } iii;

    public:

        static inline uint8_t space_for_bt(uint64_t i) {
            return iii.m_space_for_bt[i];
        }

		// TODO: is a speed-up possible if we only decode as many bit as we need?
        static inline uint256_t nr_to_bin(uint16_t k, uint256_t nr) {
            if (k == n) {
                return iii.m_L1Mask[n]; 
            } else if (k == 0) {
                return uint256_t(0);
            }
//		   	else if (k == 1) {
//              return (((uint256_t)1ULL) << (int)((uint64_t)n-(uint64_t)nr-1) );
//        }

            uint256_t bin = uint256_t(0);
            uint256_t mask = uint256_t(1);
            uint16_t nn = n;
//			while (k > 1) {
			while (k > 0) {
                if ( nr >= iii.m_coefficients[nn-1][k] ) {
                    nr -= iii.m_coefficients[nn-1][k];
                    --k;
                    bin |= mask;
                }
                --nn;
                mask = (mask << 1);
            }
            /* k == 1 */
			return bin;
//            return bin | (((uint256_t(1)) << (int)(((uint64_t)n)-((uint64_t)nr)-1));
        };

        static inline uint256_t bin_to_nr(uint256_t bin) {
            if (bin == uint256_t(0) or bin == iii.m_L1Mask[n] ){ // handle special case
                return uint256_t(0);
            }
            uint256_t nr = 0;
            uint16_t  k  = popcount(bin);
            uint16_t  nn = n; // size of the block
            while ( (const uint256_t)bin > 0ULL ) {
                if ( 1ULL & bin ) {
                    nr += (iii.m_coefficients[nn-1][k]);
                    --k; // go to the case (n-1, k-1)
                }// else go to the case (n-1, k)
                bin = (bin >> 1);
                --nn;
            }
            return nr;
        }

        static inline uint8_t space_for_bt_pair(uint8_t x) { return 0; }

		static inline uint16_t popcount(number_type x) {
			return x.popcount();
		}

		static inline uint16_t l1BP(number_type x) {
			return x.l1BP();
		}

		static inline uint16_t select(number_type x, uint32_t i) {
			return x.select(i);
		}

		template<class bit_vector_type>
		static inline number_type decode_btnr(const bit_vector_type &bv, 
											  typename bit_vector_type::size_type btnrp,
											  uint16_t btnrlen){
			if ( btnrlen <= 64 ){
				return number_type( bv.get_int(btnrp, btnrlen) );
			}else if( btnrlen <= 128 ){
                return number_type( bv.get_int(btnrp, 64),
									bv.get_int(btnrp+64, btnrlen-64) );
			}else if( btnrlen <= 192 ){
                return number_type( bv.get_int(btnrp, 64),
									bv.get_int(btnrp + 64, 64),
									(uint128_t)bv.get_int(btnrp + 128, btnrlen-128) );
			}else{ // > 192
                return number_type( bv.get_int(btnrp, 64),
									bv.get_int(btnrp+64, 64),
									(((uint128_t)bv.get_int(btnrp+192, btnrlen-192))<<64) | bv.get_int(btnrp+128, 64) );
			}
		}

		template<class bit_vector_type>
		static inline uint16_t get_bt(const bit_vector_type &bv, 
									 typename bit_vector_type::size_type pos,
									 uint16_t block_size){
			if ( block_size <= 64 ){
				return bit_magic::b1Cnt( bv.get_int(pos, block_size) );
			}else if( block_size <= 128 ){
                return bit_magic::b1Cnt( bv.get_int(pos+64, block_size-64) )+ 
					   bit_magic::b1Cnt( bv.get_int(pos, 64) );
			}else if( block_size <= 192 ){
		          return bit_magic::b1Cnt( bv.get_int(pos+128, block_size-128) )+ 
						bit_magic::b1Cnt( bv.get_int(pos+64, 64) )+ 
					    bit_magic::b1Cnt( bv.get_int(pos, 64) );
			}else{ // > 192
			      return bit_magic::b1Cnt( bv.get_int(pos+192, block_size-192) )+ 
						bit_magic::b1Cnt( bv.get_int(pos+128, 64) )+ 
						bit_magic::b1Cnt( bv.get_int(pos+64, 64) )+ 
					    bit_magic::b1Cnt( bv.get_int(pos, 64) );
			}
		}

		template<class bit_vector_type>
		static void set_bt(bit_vector_type &bv, 
						   typename bit_vector_type::size_type pos,
						   number_type bt,
						   uint16_t space_for_bt){
			if ( space_for_bt <= 64 ){
				bv.set_int(pos, bt, space_for_bt);
			}else if ( space_for_bt <= 128 ){
                bv.set_int(pos, bt, 64);
			    bv.set_int(pos+64, bt>>64, space_for_bt-64 );
			}else if ( space_for_bt <= 192 ){
                bv.set_int(pos, bt, 64);
			    bv.set_int(pos+64, (uint64_t)(bt>>64), 64 );
			    bv.set_int(pos+128, (uint64_t)(bt>>128), space_for_bt-128 );
			}else{ // > 192
                bv.set_int(pos, bt, 64);
			    bv.set_int(pos+64, bt>>64, 64 );
			    bv.set_int(pos+128, bt>>128, 64 );
			    bv.set_int(pos+192, bt>>192, space_for_bt-192 );
			}
		}

		static inline number_type Li1Mask(uint16_t off){
			return iii.m_L1Mask[off];
		}
};

template<uint16_t n>
typename binomial256<n>::impl binomial256<n>::iii;



} // end namespace
#endif
