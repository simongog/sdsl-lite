#include <iostream>
#include <sdsl/rmq_support.hpp>

using namespace sdsl;
using namespace std;


int main()
{
    rmq_succinct_sct<> rmq;
    {
        // initilaize int_vector object A
        int_vector<> A = {6,7,3,2,4,2,1,5,8,2};
        // output a
        cout << "A = " << A << endl;

        // generate range-minimum structure by passing
        // pointer to ordered array
        rmq = rmq_succinct_sct<>(&A);

        // calculate index of minimal element in a[0..a.size()-1]
        auto min_idx = rmq(0, A.size()-1);

        cout << "A[" << min_idx << "]=" << A[min_idx];
        cout << " is the samllest element in A[0.." << A.size()-1 << "]" << endl;

        // calculate index of minimal element in a[2..5]
        min_idx = rmq(2,5);

        cout << "A[" << min_idx << "]=" << A[min_idx];
        cout << " is the smallest element in A[2..5]" << endl;
    }
    // no A does not exist any more, queries can still be answered
    auto min_idx = rmq(1,5);
    cout << "A[" << min_idx << "] is the smallest in A[1..5]" << endl;
    return 0;
}
