#include "csa_typedefs.hpp"
#include <sdsl/csa_wt.hpp>
#include <sdsl/csa_sada.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/construct.hpp>
#include <iostream>
#include <string>

using namespace sdsl;
using namespace std;

int main(int argc, char **argv) {
	if ( argc < 2 ){
		cout << "Usage ./" << argv[0] << " input_file tmp_dir" << endl;
		return 0;
	}
	CSA_TYPE csa;
	if ( argc < 3 ){
		construct(csa, argv[1], 1);
	}else{
		// config: do not delete files, tmp_dir=argv[2], id=basename(argv[1])
		cache_config cconfig(false, argv[2], util::basename(argv[1]));	
		construct(csa, argv[1], cconfig, 1);
	}
	util::store_to_file( csa, (string(argv[1]) + "." + string(SUF)).c_str() );
}
