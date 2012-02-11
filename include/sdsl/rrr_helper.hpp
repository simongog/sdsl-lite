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

namespace sdsl{

// Helper class
template<uint8_t n>	
class binomial{
	private:

	static class impl{	
		public:
		static const int MAX_SIZE=32;	
		uint8_t m_space_for_bt[MAX_SIZE];
		uint8_t m_space_for_bt_pair[256*(n==15)];
		uint64_t m_C[MAX_SIZE];
		int_vector<32> m_nr_to_bin;
		int_vector<32> m_bin_to_nr;

		impl(){
			m_nr_to_bin.resize( 1<<n );
			m_bin_to_nr.resize( 1<<n );
			for(int i=0, cnt=0, class_cnt=0; i<=n; ++i)
			{
				m_C[i] = cnt;
				class_cnt = 0;
				std::vector<bool> b(n,0); 
				for(int j=0; j<i; ++j) b[n-j-1] = 1;
				do{
					uint32_t x=0;
					for(int k=0; k<n; ++k)
						x |= ((uint32_t)b[n-k-1])<<(n-1-k);
					m_nr_to_bin[cnt] = x;
					m_bin_to_nr[x] = class_cnt;
					++cnt;
					++class_cnt;
				}
				while( next_permutation(b.begin(), b.end()) );
				if( class_cnt == 1 )
					m_space_for_bt[i] = 0;
				else
					m_space_for_bt[i] = bit_magic::l1BP(class_cnt)+1;
	//			cout << cnt << " " << class_cnt << " " << (int)m_space_for_bt[i] << endl;
			}	
			if( n == 15 ){
				for(int x=0; x<256; ++x){
					m_space_for_bt_pair[x] = m_space_for_bt[x>>4] + m_space_for_bt[x&0x0F];
				}
			}
		} 
	} iii;

	public:
	
	static inline uint8_t space_for_bt(uint32_t i){
		return iii.m_space_for_bt[i];
	}

	static inline uint32_t nr_to_bin(uint8_t k, uint32_t nr){
		return iii.m_nr_to_bin[iii.m_C[k]+nr];	
	}

	static inline uint32_t bin_to_nr(uint32_t bin){
		return iii.m_bin_to_nr[bin];	
	}

	static inline uint8_t space_for_bt_pair(uint8_t x){
		return iii.m_space_for_bt_pair[x];
	}
};
template<uint8_t n>
typename binomial<n>::impl binomial<n>::iii;

}
#endif
