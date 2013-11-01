#include <iostream>
#include <fstream>
#include <sdsl/suffix_arrays.hpp>

using namespace sdsl;
int main(int argc, char* argv[])
{

    if (argc > 3) {
        cache_config cc(false, argv[3],"gen_bwt_");
        std::cout << cache_file_name(conf::KEY_BWT, cc) << std::endl;
        {
            csa_wt<> wt;
            construct(wt, argv[1], cc, 1);
        }
        int_vector<8> bwt;
        load_from_file(bwt, cache_file_name(conf::KEY_BWT, cc));
        std::ofstream out(argv[2]);
        out.write((char*)bwt.data(), bwt.size());
        util::delete_all_files(cc.file_map);
    } else {
        std::cout<<"Usage: input_file output_file temp_dir" << std::endl;
    }
}
