#include <sdsl/wavelet_trees.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[]){
//	wt_int<> wt;
//	wt<> wt;
//	wt_rlmn<> wt;
//	wt_rlg<> wt;
//	wt_rlg8<> wt;
	wt_huff<> wt;
	construct(wt, argv[1], 1);
	for (wt_int<>::size_type i=0; i<wt.size(); ++i){
		cout<<(char)wt[i];
	}
	cout<<endl;
}
