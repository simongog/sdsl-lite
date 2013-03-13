#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[]){
	if ( argc < 2 ){
		cout << "Usage: " << argv[0] << " int32_file nr" << endl;
		cout << "     Reads a file containing 32-bit integers and outputs the first nr of them." << endl;
		return 1;
	}
	int_vector<> v;

	util::load_vector_from_file(v, argv[1], 4);
	
	cout << "v.size() = " << v.size() << endl;
	cout << "v.width() = "  << (int)v.width() << endl;
	cout << "v =";
	size_t max_out = 10;
	if ( argc > 2 ){
		max_out = atoll(argv[2]);
	}
	for (size_t i=0; i < max_out and i < v.size(); ++i)
		cout << " " << v[i];
	if (v.size() > max_out)
		cout << "....";
	cout << endl;
}
