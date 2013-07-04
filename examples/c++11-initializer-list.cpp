#include <sdsl/suffix_trees.hpp>
#include <iostream>
#include <algorithm>

using namespace sdsl;
using namespace std;

int main()
{
    int_vector<8> b = {'c','d','a','b','d'};
    int_vector<> a = {1,8,3,15,5};

    for (auto c : b) {
        cout << c;
    } cout << endl;
    for (auto x : a) {
        cout << x << ",";
    }
    cout << endl;
    sort(a.begin(), a.end());
    for (auto x: a) {
        cout << x << ",";
    }
}
