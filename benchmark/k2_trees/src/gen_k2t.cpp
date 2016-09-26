#include <iostream>
#include <fstream>
#include <tuple>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/k2_tree.hpp>
#include <stdexcept>
#include <string>
#include <vector>

using namespace sdsl;

void generate_bit_vector_buffers(const std::string& idx_file,
								 const std::string& output_x_file,
								 const std::string& output_y_file)
{
	std::ifstream infile(idx_file);
	std::string line;
	uint64_t cnt = 0;

	for (int i = 0; std::getline(infile, line); ++i)
		    cnt++;

	infile.clear();
	infile.seekg(0, std::ios::beg);
	// Set size of vector to the amount of lines in the input file.
	int_vector<>xv(cnt), yv(cnt);
	cnt = 0;

	while(std::getline(infile, line)) {
		sdsl::k2_tree_ns::idx_type x, y;
		std::istringstream iss(line);
		if(!(iss >> x >> y))
			throw std::invalid_argument("Not expected line at construct");
		xv[cnt] = x;
		yv[cnt++] = y;
	}

	store_to_file(xv, output_x_file);
	store_to_file(yv, output_y_file);
}

inline bool exists(const std::string& name) {
	std::ifstream f(name.c_str());
	return f.good();
}

int main(int argc, char* argv[])
{
	if(argc < 4) {
        std::cout<<"Usage: input_file output_file_prefix output_k2_file" << std::endl;
	}

	std::string out_x(argv[2]);
	out_x.append(".x");
	std::string out_y(argv[2]);
	out_y.append(+ ".y");

	if(!exists(out_x) || !exists(out_y))
		generate_bit_vector_buffers(argv[1], out_x, out_y);

	K2_TYPE k2(argv[2]);
	std::ofstream fs;
	fs.open(argv[3]);
	k2.serialize(fs);
}
