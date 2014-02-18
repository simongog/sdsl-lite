#include <sdsl/vectors.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main()
{
    std::string  tmp_file = "tmp_file.sdsl";
    {
        // generate int_vector and store to a file
        int_vector<> v = {1,3,5,7,2,3,4,9,8,7,10,1};
        store_to_file(v, tmp_file);
    }
    {
        // use int_vector_buffer to open the file
        int_vector_buffer<> ivb(tmp_file);
        // output elements and assign 42 to each of them
for (auto x : ivb) {
            cout << x << endl;
            x = 42;
        }
        // int_vector_buffer is destroy at the end of the
        // scope and all value are written to disk
    }
    {
        // read vector from file and output it
        int_vector<> v;
        load_from_file(v, tmp_file);
        cout << v << endl;
    }
    // delete temporary file
    sdsl::remove(tmp_file);
}
