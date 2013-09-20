#include <sdsl/rmq_support.hpp>
#include <iostream>

using namespace std;
using namespace sdsl;

int main()
{
    //                0 1 2 3 4 5 6 7 8 9 0
    int_vector<> v = {5,3,8,9,1,2,5,3,9,0,7};
    rmq_succinct_sct<> rmq(&v); // <- pointer to b
    util::clear(v);
    cout << "v.size() = " << v.size() << endl;
    cout << rmq(0, 10) << ", " << rmq(2, 7) << endl;
}
