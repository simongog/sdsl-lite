#include <sdsl/vectors.hpp>
#include <iostream>
#include <algorithm>

using namespace sdsl;
using namespace std;

int main()
{
    string file = "int_vector_buffer_iterator.sdsl";
    {
        int_vector<> iv = {1,2,3,4,5,9,8,7,6};
        store_to_file(iv, file);
    }
    int_vector_buffer<> ivb(file);
    cout << std::accumulate(ivb.begin(), ivb.end(), 0) << endl;
    cout << std::count_if(ivb.begin(), ivb.end(), [](uint64_t x) {
        return 0 == x%2;
    }) << endl;
}


