#include <sdsl/suffix_trees.hpp>
#include <sdsl/io.hpp>
#include <iostream>
#include <chrono>

using timer = std::chrono::high_resolution_clock;
using namespace std::chrono;
using namespace sdsl;

int main(int argc, char** argv)
{
    if (argc < 2) {
      std::cout << "Usage: " << argv[0] << " file" << std::endl;
      std::cout << " Creates a CST for a byte file and visualizes the space used by the data structure." << std::endl;
      return 1;
    }
    cst_sct3<> cst;
    auto start = timer::now();
    std::cout << "constructing cst..." << std::endl;
    construct(cst, argv[1], 1);
    std::cout << "construction cst time in seconds: " << duration_cast<seconds>(timer::now()-start).count() << std::endl;

    std::ofstream ofs("cst-space-usage.html");
    std::cout << "writing storage visualization to cst-space-usage.html" << std::endl;
    sdsl::write_structure<HTML_FORMAT>(cst,ofs);
}
