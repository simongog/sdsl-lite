#include "sdsl/algorithms_for_suffix_array_construction.hpp"

namespace sdsl{

namespace algorithm{	
// \param c Char array pointing to the text
// \param n Length of the text
bool shift_text(char *c, uint64_t n, bool shift){
	uint64_t cnt_c[256] = {0};
	for(uint64_t i=0; i<n; ++i){
		++cnt_c[(unsigned char)c[i]];	
	}
	int cnt_cc=0;
	int min_cnt=1;
	for(int i=0; i<256; ++i){
		cnt_cc += cnt_c[i]>0;
		if( i > 0 and cnt_c[i] < cnt_c[min_cnt] )
			min_cnt = i;
	}
	if(cnt_cc == 256){
		std::cerr << "# ERROR: Alphabet size>255!!!" << std::endl;
		std::cerr << "min_cnt = "<< min_cnt << " C[min_cnt]="<< cnt_c[min_cnt] << std::endl;
		for(uint64_t i=0; i<n; ++i)
			if( !c[i] )
				c[i] = min_cnt;
		return false;
	}
	else{
		if(shift){
			cnt_cc=0;
			for(int i=0; i<256; ++i)
				if(cnt_c[i])
					cnt_c[i] = ++cnt_cc;
			for(uint64_t i=0; i<n; ++i)
				c[i] = cnt_c[(unsigned char)c[i]];
		}
	}
	return true;
}

} // end namespace algorithm

} // end namespace sdsl

