#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sdsl/int_vector.hpp>

using namespace sdsl;
using namespace std;

template<typename T>
void write_inc_file(const char *name, size_t len){
	ofstream out(name);
	for(size_t i=0; i<len; ++i){
		T x = i;
		out.write((char*)&x,sizeof(x));
	}
	out.close();
}

template<typename T>
void write_const_file(const char *name, size_t len){
	int_vector<8> v(len, 42);
	util::store_to_plain_array<T>(v, name);
}

int main(int argc, char* argv[]){
	if ( argc < 2 ){
		cout << "Usage: "<<argv[0]<<" len" << endl;
		cout << "Generates 8 files. The first 4 contain an increasing sequence\n";
		cout << "[0..len-1], where each number is represented with a fixed\n";
		cout << "x=8, 16, 32, and 64-bit integer. The result files are named\n";
		cout << "v.xbit.\n";
		cout << "The remaining files contain a sequence of length len and all";
		cout << "elements are set to 42. The result riles are named v42.xbit.";
		return 1;
	}
	size_t len = atoll(argv[1]);
	write_inc_file<uint8_t>("v.8bit",   len);
	write_inc_file<uint16_t>("v.16bit", len);
	write_inc_file<uint32_t>("v.32bit", len);
	write_inc_file<uint64_t>("v.64bit", len);

	write_const_file<uint8_t>("v42.8bit",   len);
	write_const_file<uint16_t>("v42.16bit", len);
	write_const_file<uint32_t>("v42.32bit", len);
	write_const_file<uint64_t>("v42.64bit", len);


}
