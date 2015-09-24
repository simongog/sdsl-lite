#include "sdsl/int_vector.hpp"
#include <string>
#include <cstdlib>

using namespace std;
using namespace sdsl;

int main(int argc, char* argv[])
{
    if (argc < 4) {
        cout << "Usage: " << argv[0] << " FILE X Y" << endl;
        cout << "  Loads an int_vector<>(SIZE, DEFAULT_VALUE, WIDTH)" << endl;
        cout << "  and replaces all values X with Y and stores the " << endl;
        cout << "  result in the same file." << endl;
        return 1;
    }
    uint64_t x = atoi(argv[2]);
    uint64_t y = atoi(argv[3]);
    int_vector_buffer<> v(argv[1]);
    if (v.good()) {
        for (size_t i=0; i<v.size(); ++i) {
            uint64_t val = v[i];
            if (val == x) {
                v[i] = y;
            }
        }
    } else {
        cerr << "Could not open int_vector<> file " << argv[1] << endl;
        return 1;
    }
}
