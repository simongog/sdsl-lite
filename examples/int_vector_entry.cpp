//! Program for debugging. Load vector v of type int_vector<> from
// file argv[1]. Then the program outputs for each integer i,
// which is read from stdin, the value v[i].
#include <sdsl/int_vector.hpp>
#include <iostream>

using namespace std;
using namespace sdsl;

int main(int argc, char* argv[]){
	int_vector<> v;
	util::load_from_file(v, argv[1]);
	cout<<"loaded vector of size "<<v.size()<<" and "<<(int)v.width()<<"-bit integers from file";
	cout<<" "<<argv[1]<<endl;
	if ( v.size() == 0 )
		return 0;
	if (argc>3){
		size_t a=atoi(argv[2]);
		size_t b=atoi(argv[3]);
		if ( b >= v.size() )
			b = v.size()-1;
		if ( a > b )
			a = b;
		for(size_t i=a; i<=b; ++i)
			cout<<"v["<<i<<"]="<<v[i]<<endl;
	}
	size_t i;
	while (cin>>i){
		cout<<v[i]<<endl;
	}
}
