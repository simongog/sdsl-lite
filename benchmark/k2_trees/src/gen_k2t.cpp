#include <iostream>
#include <fstream>
#include <tuple>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/k2_tree.hpp>
#include <stdexcept>
#include <string>
#include <vector>

using namespace sdsl;

typedef K2_TYPE::idx_type idx_type;
typedef K2_TYPE::size_type size_type;

template<class k2_type>
void load_edges_from_file(
		const std::string &file_path,
		std::vector<std::tuple<typename k2_type::idx_type,
							   typename k2_type::idx_type>>&vec,
		typename k2_type::idx_type &max)
{
	std::ifstream infile(file_path);
	std::string line;
	typedef std::tuple<typename k2_type::idx_type,
					   typename k2_type::idx_type> t_tuple;
	max = 0;

	while(std::getline(infile, line)) {
		std::istringstream iss(line);
		idx_type v, u;
		if(!(iss >> v >> u)) {
			throw std::invalid_argument("Not expected line at construct");
		}
		if(v > max)
			max = v;
		if(u > max)
			max = u;
		vec.push_back(t_tuple {v, u});
	}
	max++;
}

int main(int argc, char* argv[])
{
	if(argc < 3) {
        std::cout<<"Usage: input_file output_file temp_dir" << std::endl;
	}

	std::vector<std::tuple<typename K2_TYPE::idx_type,
						   typename K2_TYPE::idx_type>> v;
	K2_TYPE::idx_type max;
    load_edges_from_file<K2_TYPE>(argv[1], v, max);
	K2_TYPE k2(v, max);

	std::ofstream fs;
	fs.open(argv[2]);
	k2.serialize(fs);
	fs.close();
}
