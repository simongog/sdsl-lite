#include <iostream>
#include <fstream>
#include <sdsl/suffix_arrays.hpp>
#include <string>
#include <vector>

using namespace sdsl;

//routine to save a vector in different formats, see lower implementations
template<class INT_VECTOR, uint8_t num_byte>
void saveVector(const INT_VECTOR &v, const char *dest);

//main function to generate MTF of BWT of an integer vector.
// CSA_WT: used wavelet - tree - based suffix array implementation
// INT_VECTOR: used integer vector for extracting BWT
// num_byte: value indicating how result has to be opened / saved
// srcfile: file from which to generate
// destfile: file where to save result
// tmpdir: directory used for temporary results
// conf_bwt_key: key what is able to fetch bwt after suffix array construction
template<class CSA_WT, class INT_VECTOR, uint8_t num_byte>
void gen_bwt_mtf(const char *srcfile, const char *destfile, const char *tmpdir,
                 const char *conf_bwt_key) {
	//utility for CSA generation	
	cache_config cc(false, tmpdir, "gen_bwt_mtf_");
	INT_VECTOR bwt;

	//create suffix array
	CSA_WT wt;
	construct(wt, srcfile, cc, num_byte);

	//compute alphabet table from suffix array
	std::vector<uint64_t> alph_tbl( wt.sigma );
	for (uint64_t i = 0; i < wt.sigma; i++) {
		alph_tbl.push_back( wt.comp2char[i] );
	}

	//fetch bwt
	load_from_file(bwt, cache_file_name(conf_bwt_key, cc));

	//create mtf
	for (uint64_t i = 0; i < bwt.size(); i++) {
		uint64_t c = bwt[i];
		//find c in alphabet table and move it to front
		uint64_t j = 0;
		do {
			uint64_t tmp = alph_tbl[j];
			alph_tbl[j++] = c;
			c = tmp;
		} while (c != alph_tbl.front());
		//and write it's index to mtf transform of bwt
		bwt[i] = j-1;
	}
	
	//save everything
	saveVector<INT_VECTOR, num_byte>( bwt, destfile );
	
	//and free resources
	util::delete_all_files(cc.file_map);
}

//functions for saving an integer vector in different formats
//generic version (raw output)
template<class INT_VECTOR, uint8_t num_byte>
void saveVector(const INT_VECTOR &v, const char *dest) {
	std::ofstream out(dest);
	out.write((char *)v.data(), num_byte * v.size());
}
//serialization of integer vector
template<>
void saveVector<int_vector<>, 0>(const int_vector<> &v, const char *dest) {
	store_to_file(v, dest);
}
//decimal digits
template<>
void saveVector<int_vector<>, 'd'>(const int_vector<> &v, const char *dest) {
	std::ofstream out(dest);
	if (v.size())	out << v[0];
	for (uint64_t i = 1; i < v.size(); i++) {
		out << " " << v[i];
	}
}

//main function
int main(int argc, char* argv[]) {
	if (argc != 5) {
        	std::cout<<"Usage: input_file output_file temp_dir num_byte" << std::endl;
		return 1;
    	}
	std::cout << "Calculate MTF Transform of BWT of " << argv[1] 
	          << " and store it to " << argv[2] << std::endl;

	typedef csa_wt<> csa_wt_byte;
	typedef csa_wt<wt_int<>, 64, 64, sa_order_sa_sampling<>, int_vector<>, int_alphabet<>> csa_wt_int;

	switch (argv[4][0]) {
	case 'd':	//decimal digits
		gen_bwt_mtf<csa_wt_int, int_vector<>, 'd'>(argv[1], argv[2], argv[3], conf::KEY_BWT_INT);
		return 0;
	case '0':	//serialized integer vector
		gen_bwt_mtf<csa_wt_int, int_vector<>, 0>(argv[1], argv[2], argv[3], conf::KEY_BWT_INT);
		return 0;
	case '1':	//byte integer vector
		gen_bwt_mtf<csa_wt_byte, int_vector<8>, 1>(argv[1], argv[2], argv[3], conf::KEY_BWT);
		return 0;
	case '2':	//2 byte integer vector
		gen_bwt_mtf<csa_wt_int, int_vector<>, 2>(argv[1], argv[2], argv[3], conf::KEY_BWT_INT);
		return 0;
	case '4':	//4 byte integer vector
		gen_bwt_mtf<csa_wt_int, int_vector<>, 4>(argv[1], argv[2], argv[3], conf::KEY_BWT_INT);
		return 0;
	case '8':	//8 byte integer vector
		gen_bwt_mtf<csa_wt_int, int_vector<>, 8>(argv[1], argv[2], argv[3], conf::KEY_BWT_INT);
		return 0;
	default:
		std::cout << "Illegal num_byte, allowed are 'd', 0, 1, 2, 4, 8" << std::endl;
		return 1;
	}
}
