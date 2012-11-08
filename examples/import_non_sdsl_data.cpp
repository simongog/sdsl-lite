#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[]){
	if ( argc < 2 ){
		cout << "Usage: " << argv[0] << " text_file nr" << endl;
		cout << "     Reads a text file and outputs the first nr symbols." << endl;
		return 1;
	}
	int_vector<8> text;

	util::load_vector_from_file(text, argv[1], 1);
	
	cout << "text.size() = " << text.size() << endl;
	cout << "text =";
	size_t max_out = 10;
	if ( argc > 2 ){
		max_out = atoll(argv[2]);
	}
	for (size_t i=0; i < max_out and i < text.size(); ++i)
		cout << " " << text[i];
	if (text.size() > max_out)
		cout << "....";
	cout << endl;
}
