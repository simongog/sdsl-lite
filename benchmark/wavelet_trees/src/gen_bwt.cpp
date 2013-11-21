#include <iostream>
#include <fstream>
#include <sdsl/suffix_arrays.hpp>
#include <string>

using namespace sdsl;
int main(int argc, char* argv[])
{
    if (argc == 6) {
        uint8_t num_byte = argv[5][0]=='d' ? 'd' : argv[5][0]-'0';
        if (strcmp(argv[4], "BWT")==0) {
            std::cout<<"Calculate BWT of " << argv[1] << " and store it to " << argv[2] << std::endl;
            cache_config cc(false, argv[3], "gen_bwt_");
            if (1 == num_byte) {
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
                {
                    csa_wt<wt_int<>, 64, 64, sa_order_sa_sampling<>, int_vector<>, int_alphabet<>> wt;
                    construct(wt, argv[1], cc, num_byte);
                }
                int_vector<> bwt;
                load_from_file(bwt, cache_file_name(conf::KEY_BWT_INT, cc));
                std::ofstream out(argv[2]);
                if ('d' == num_byte) {
                    if (bwt.size()) {
                        out << bwt[0];
                    }
                    for (uint64_t i=1; i<bwt.size(); ++i) {
                        out << " " << bwt[i];
                    }
                } else if (0 == num_byte) {
                    store_to_file(bwt, argv[2]);
                } else {
                    out.write((char*)bwt.data(), num_byte*bwt.size());
                }
                util::delete_all_files(cc.file_map);
            }
        }
    } else {
        std::cout<<"Usage: input_file output_file temp_dir create_bwt num_byte" << std::endl;
    }
}
