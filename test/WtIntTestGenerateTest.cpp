#include "sdsl/int_vector.hpp"
#include "sdsl/config.hpp" // for CMAKE_SOURCE_DIR

using namespace std;
using namespace sdsl;

int main(){
	std::string test_cases_dir = std::string(SDSL_XSTR(CMAKE_SOURCE_DIR)) + "/test/test_cases/crafted";
	util::store_to_file(int_vector<>(), (test_cases_dir + "/empty.iv").c_str());
	util::store_to_file(int_vector<>(1023,0,1), (test_cases_dir + "/1023-0bits.iv").c_str());
	util::store_to_file(int_vector<>(100023,0,1), (test_cases_dir + "/100023-0bits.iv").c_str());
	util::store_to_file(int_vector<>(64,0,2), (test_cases_dir + "/64-0bitpairs.iv").c_str());
	int_vector<> iv(1000000,0,18);
	util::store_to_file(iv,  (test_cases_dir + "/1000000-random-bits.iv").c_str());
}
