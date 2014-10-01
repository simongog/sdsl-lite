/* This example shows how the representation of the alphabet dependent
 * part of a CST can be altered by using policy classes.
 */
#include <sdsl/suffix_arrays.hpp>
#include <iostream>
#include <string>

using namespace sdsl;
using namespace std;

template<class csa_t>
void csa_info(csa_t& csa, const char* file, bool json)
{
    cout << "file          : " << file << endl;
    construct(csa, file, 1);
    cout << "csa of type   : " << util::demangle(typeid(csa).name()) << endl;
    cout << "size in bytes : " << size_in_bytes(csa) << endl;
    if (json) {
        cout << "---------------" << endl;
        cout << "json output: " << endl;
        write_structure<JSON_FORMAT>(csa, cout);
        cout << endl;
    }
    cout << "---------------" << endl;
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " file [JSON]" << endl;
        cout << " (1) Constructs CSAs 2 csa_sada and 2 csa_wt, with" << endl;
        cout << "     alphabet strategies." << endl;
        cout << " (2) Outputs type and size. If JSON is specified," << endl;
        cout << "     also the structure in JSON-format." << endl;
        return 1;
    }
    bool json = false;
    if (argc > 2) {
        json = true;
    }
    csa_sada<enc_vector<>, 32, 32, sa_order_sa_sampling<>, isa_sampling<>, byte_alphabet> csa1;
    csa_sada<enc_vector<>, 32, 32, sa_order_sa_sampling<>, isa_sampling<>, succinct_byte_alphabet<> > csa2;
    csa_wt<wt_huff<>, 32, 32, sa_order_sa_sampling<>, isa_sampling<>, byte_alphabet> csa3;
    csa_wt<wt_huff<>, 32, 32, sa_order_sa_sampling<>, isa_sampling<>, succinct_byte_alphabet<> > csa4;
    csa_info(csa1, argv[1], json);
    csa_info(csa2, argv[1], json);
    csa_info(csa3, argv[1], json);
    csa_info(csa4, argv[1], json);
}
