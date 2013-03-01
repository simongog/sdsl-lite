#include <sdsl/rrr_vector.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main(){
	bit_vector b(100,0);
	b[5]=1;
	rrr_vector<15> v(b);
	cout << v[0] << endl;
	cout << v[5] << endl;
}
