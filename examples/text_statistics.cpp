#include <sdsl/suffixtrees.hpp>
#include <iostream>

using namespace std;
using namespace sdsl;

int main(int argc, char* argv[]){
	if ( argc < 2 ){
		cout << "Usage: "<<argv[0]<<" file"<<endl;
		return 1;
	}
	cst_sct3<> cst;
	typedef cst_sct3<>::size_type size_type;
	typedef cst_sct3<>::char_type char_type;
	construct_cst(argv[1], cst);

	long double red = 0;
	long double lcp_mean = 0;
	if( cst.csa.size() ){
		char_type bwti_1 = cst.csa.bwt[0];
		for(size_type i=1; i<cst.csa.size(); ++i){
			char_type bwti = cst.csa.bwt[i];
			if ( bwti_1 == bwti  ){
				red += 1.0;
			}else{
				bwti_1 = bwti;
			}
			lcp_mean += cst.lcp[i];
		}
		lcp_mean /= cst.csa.size();
		cout << "reducible ratio in percent: " << (red/cst.csa.size())*100 << endl;
		cout << "lcp mean: " << lcp_mean << endl;
	}
}
