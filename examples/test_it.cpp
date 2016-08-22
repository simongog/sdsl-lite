#include <vector>
#include <stxxl/vector>
#include <bitset>
#include "sdsl/k2_tree_helper.hpp"

using namespace std;

int main()
{
 /*   for (uint i = 0; i < 512; i++){
        uint a = 0;
        a = a | (i << (32-9) >> (32-3) << 12)
              | (i << (32-6) >> (32-3) << 6)
              | (i << (32-3) >> (32-3));
        std::cout << "0x" << hex << a << ", \t";
    }*/

    using namespace sdsl;
    using namespace k2_treap_ns;


    uint q = 22576;
    uint p = 51760;

    std::cout << "q: " << bitset<32>(q) << std::endl;
    std::cout << "p: " << bitset<32>(p) << std::endl;

    uint64_t z = access_shortcut_helper<8>::corresponding_subtree(q, p, 16, 2);

    std::cout << bitset<32>(z) << std::endl;

    q = 51760;
    p = 22576;
    z = access_shortcut_helper<8>::corresponding_subtree(q, p, 16, 2);

    std::cout << bitset<32>(z) << std::endl;
    return 0;
}
