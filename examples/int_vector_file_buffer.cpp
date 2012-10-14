#include <sdsl/int_vector.hpp>
#include <string>
#include <iostream>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[]){
	size_t size = 10000000;
	// create a int_vector (each int of fixed size 8 bits) of length 10000000
	int_vector<8> v(size);

	// initialize vector with random bits, seed of the random process = 42
	util::set_random_bits(v, 42);

	// store int_vector to temporary file
	string tmp_file = "file_buffer_example.int_vector";
	util::store_to_file(v, tmp_file.c_str());

	// open file buffer
	int_vector_file_buffer<8> v_buf(tmp_file.c_str());
		
	// stream vector data from disk and compare it with in-memory data
	for (size_t i = 0, r_sum = 0, r = v_buf.load_next_block(); i < r_sum;){
		for(; i < r_sum +r; ++i){
			if( v[i] != v_buf[i-r_sum] ){
				std::cerr << "ERROR: v["<< i << "] != v_buf[" << i-r_sum << "]" << std::endl;
				return 1;
			}
		}
		r_sum += r; r = v_buf.load_next_block();
	}

	// remove temporary file
	std::remove(tmp_file.c_str());
}
