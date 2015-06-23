#include <sdsl/suffix_arrays.hpp>
#include <iostream>
#include <fstream>

using namespace sdsl;
using namespace std;

int main()
{
    csa_wt<> csa;
    construct_im(csa, "This is how it works", 1);
    // write only one object to std::cout
    write_structure<HTML_FORMAT>(csa, cout);
    wt_int<> wt;
    construct_im(wt, int_vector<>(1000,3));
    // write multiple objects into one file
    {
        ofstream out("write_structure.html");
        write_structure<HTML_FORMAT>(out, csa, wt);
    }
    // write one object into a file
    write_structure<HTML_FORMAT>(csa, "csa_structure.html");
}
