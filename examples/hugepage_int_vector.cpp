#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vector_il.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main()
{
    {
        bit_vector v(1ULL<<30);
        bool mapped = mm::map_hp();
        cout << v[0] << endl;
        if (mapped) {
            mm::unmap_hp();
        }
    }
    {
        bit_vector_il<> v(1ULL<<30);
        bool mapped = mm::map_hp();
        cout << v[0] << endl;
        if (mapped) {
            mm::unmap_hp();
        }
    }
}
