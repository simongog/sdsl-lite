#include <iostream>
#include <fstream>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <string>

using namespace std;
using namespace sdsl;

int main(int argc, char* argv[])
{
    if (argc < 4) {
        cout << "Usage: " << argv[0] << " length density file" << endl;
        cout << " generates a bit_vector of`length` bits, populates it with " << endl;
        cout << " `density`\% set bits and saves it to `file`.\n" << endl;
        cout << "`length` format:\n";
        cout << "  X, where X is the length in number of bytes\n";
        cout << "  XkB, where X is the length in number of kilobytes\n";
        cout << "  XMB, where X is the length in number of megabytes\n";
        cout << "  XGB, where X is the length in number of gigabytes\n";
        return 1;
    }
    uint64_t length = 0;
    length = atoll(argv[1])*8; // length in bits
    string length_str(argv[1]);
    size_t Bpos = length_str.find_first_of("B");
    if (Bpos != string::npos and Bpos > 1) {
        char order = length_str.substr(Bpos-1,1)[0];
        if (order == 'k' or order == 'K')
            length <<= 10;
        else if (order == 'm' or order == 'M')
            length <<= 20;
        else if (order == 'g' or order == 'G')
            length <<= 30;
    }
    const uint64_t density = atoi(argv[2]);
    cout << length << endl;
    bit_vector v(length);
    srand(17);
    for (uint64_t i=0; i<v.size(); ++i) {
        if ((uint64_t)(rand()%100) < density) {
            v[i] = 1;
        }
    }
    if (!store_to_file(v, argv[3]))
        return 1;
}
