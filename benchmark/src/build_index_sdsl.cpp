#include "csa_typedefs.hpp"
#include <sdsl/csa_wt.hpp>
#include <sdsl/csa_sada.hpp>
#include <sdsl/cst_construct.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <iostream>
#include <string>

using namespace sdsl;
using namespace std;

int main(int argc, char **argv) {
	if ( argc < 2 ){
		cout << "Usage ./" << argv[0] << " input_file [tmp_dir]" << endl;
		return 0;
	}
	CSA_TYPE csa;
	if ( argc < 3 ){
		construct_csa(argv[1], csa);
	}else{
		tMSS file_map;
		construct_csa(argv[1], csa, file_map, false, argv[2], util::basename(argv[1]));
	}
	util::store_to_file( csa, (string(argv[1]) + "." + string(SUF)).c_str() );
}
