#include <vector>
#include <stxxl/vector>

using namespace std;

int main()
{

    typedef typename stxxl::VECTOR_GENERATOR<pair<uint,uint>>::result stxxl_pair_vector;
    vector<stxxl_pair_vector> buffers;
    buffers.resize(4);

    buffers[1].push_back(make_pair<uint,uint>(2,2));

    return 0;
}
