#include <sdsl/suffix_trees.hpp>
#include <iostream>

using namespace sdsl;

typedef CST_TYPE cst_type;

int main(int argc, char **argv) {
    if(argc < 6) {
        std::cout << "Usage: " << argv[0] << " input_file tmp_dir cst_file time_file struct_file" << std::endl;
        std::cout << " Generates a CST and stores it in cst_file. The time and space required for the" << std::endl;
        std::cout << " construction of the CST is stored in time_file. A visualization of the" << std::endl;
        std::cout << " data structure is stored in struct_file. Temporary files are stored in tmp_dir." << std::endl;
        return 1;
    }

    const char* input_file = argv[1];
    const char* tmp_dir = argv[2];
    const char* cst_file = argv[3];
    const char* time_file = argv[4];
    const char* struct_file = argv[5];

    memory_monitor::start();

    cst_type cst;
    cache_config config(true, tmp_dir);
    construct(cst, input_file, config, cst_type::alphabet_category::WIDTH == 8 ? 1 : 0);

    memory_monitor::stop();

    store_to_file(cst, cst_file);

    std::ofstream out_time(time_file);
    memory_monitor::write_memory_log<HTML_FORMAT>(out_time);
    out_time.close();

    std::ofstream out_struct(struct_file);
    write_structure<HTML_FORMAT>(cst, out_struct);
    out_struct.close();

    return 0;
}
