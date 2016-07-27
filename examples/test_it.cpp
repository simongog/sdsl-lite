#include <vector>
#include <stxxl/vector>

using namespace std;

int main()
{

    typedef typename stxxl::VECTOR_GENERATOR<pair<uint,uint>>::result stxxl_pair_vector;
    vector<stxxl_pair_vector> buffers;
    //buffers.reserve(4);

    stxxl_pair_vector one;
    one.push_back(make_pair<uint,uint>(2,2));
    buffers.emplace_back(one);

    return 0;
}
