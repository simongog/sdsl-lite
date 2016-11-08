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

void generate_query_file(const std::string& idx_file, std::string& output_file, uint32_t amount_of_queries = 100000, std::string option="random"){
    std::ifstream infile(idx_file);
    uint number_of_nodes = 0;
    read_member(number_of_nodes, infile);

    std::ofstream outfile (output_file.c_str() , std::ofstream::binary);

    srand(0);//seed for random for reproduceability

    write_member(amount_of_queries, outfile);
    if(option == "random"){
        for(uint i=0;i<amount_of_queries;i++) {
            write_member(rand()%number_of_nodes, outfile);
            //cout<<a[i]<<endl;
        }
    } else {
        for(uint i=0;i<amount_of_queries;i++) {
            write_member(i, outfile);
            //cout<<a[i]<<endl;
        }
    }

    outfile.close();
}

inline bool exists(const std::string& name) {
	std::ifstream f(name.c_str());
	return f.good();
}

int main(int argc, char* argv[])
{
	if(argc < 3) {
        std::cout<<"Usage: input_file output_file_prefix" << std::endl;
	}

    std::string output_file(argv[2]);
    output_file.append(".queries");

	if(!exists(output_file))
		generate_query_file(argv[1], output_file);
}
