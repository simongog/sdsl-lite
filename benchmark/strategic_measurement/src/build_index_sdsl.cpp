#include <sdsl/suffix_arrays.hpp>
#include <iostream>
#include <string>

using namespace sdsl;
using namespace std;

int main(int argc, char** argv)
{
    if (argc < 4) {
        cout << "Usage ./" << argv[0] << " input_file tmp_dir output_file" << endl;
        return 0;
    }
    CSA_TYPE csa;
    if (argc < 3) {
        construct(csa, argv[1], 1);
    } else {
        // config: do not delete files, tmp_dir=argv[2], id=basename(argv[1])
        cache_config cconfig(false, argv[2], util::basename(argv[1]));
        construct(csa, argv[1], cconfig, 1);
    }
    store_to_file(csa, argv[3]);
}
