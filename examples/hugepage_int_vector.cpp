#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vector_interleaved.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main(){	
	{
		bit_vector v(1ULL<<37);
		v.map_it();
		cout << "Hello World" << endl;
		v.unmap_it();
		util::clear(v);
	}
	{
		bit_vector_interleaved<> v(1ULL<<36);
		v.map_it();
		cout << "Hello World" << endl;
		v.unmap_it();

	}
}
