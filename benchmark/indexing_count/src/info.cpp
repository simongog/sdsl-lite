/*
 * This program outputs the structure of an index.
 */
#include <sdsl/suffix_arrays.hpp>
#include <string>

#include <cstdlib>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 2)  {
        cout << "./" << argv[0] << " index_file " << endl;
        return 1;
    }
    CSA_TYPE csa;
    load_from_file(csa, argv[1]);
    write_structure<HTML_FORMAT>(csa, cout);
}
