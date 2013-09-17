#include "sdsl/rmq_support.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <stack>

using namespace std;
using namespace sdsl;

namespace
{

typedef vector<uint64_t> tVUI;

string  test_file;
string  temp_file;

template<class T>
class RMQTest : public ::testing::Test { };


using testing::Types;

typedef Types<sdsl::rmq_succinct_sct<>,
        sdsl::rmq_succinct_sada<>
        > Implementations;

TYPED_TEST_CASE(RMQTest, Implementations);


TYPED_TEST(RMQTest, ConstructAndStore)
{
    int_vector<> v;
    load_from_file(v, test_file);
    TypeParam rmq(&v);
    bool success = sdsl::store_to_file(rmq, temp_file);
    ASSERT_EQ(success, true);
}

// helper class for next test
class state
{
    public:
        uint64_t l, r; // left and right border of interval
        uint64_t idx;  // index of the min value
        uint64_t min;  // min value in the interval
        state(uint64_t fl=0, uint64_t fr=0, uint64_t fidx = 0, uint64_t fmin=0) :
            l(fl), r(fr), idx(fidx), min(fmin) {}
};

//! Test range minimum queries
TYPED_TEST(RMQTest, RmqLoadAndQuery)
{
    int_vector<> v;
    ASSERT_TRUE(load_from_file(v, test_file));
    TypeParam rmq;
    ASSERT_TRUE(load_from_file(rmq, temp_file));
    ASSERT_EQ(v.size(), rmq.size());
    if (rmq.size() > 0) {
        stack<state> s;
        uint64_t idx = rmq(0, rmq.size()-1);
        ASSERT_TRUE(idx >= (uint64_t)0); ASSERT_TRUE(idx < rmq.size());
        s.push(state(0, rmq.size()-1, idx,  v[idx]));
        while (!s.empty()) {
            state st = s.top(); s.pop();
            if (st.l < st.idx) {
                idx = rmq(st.l, st.idx-1);
                ASSERT_TRUE(idx >= st.l); ASSERT_TRUE(idx <= st.idx-1);
                ASSERT_TRUE(v[idx] >= v[st.idx])
                        << "v["<<idx<<"]="<< v[idx]
                        << " < " << "v["<<st.idx<<"]="
                        << v[st.idx] << endl
                        << "[" << st.l << "," << st.r << "]" << endl;
                s.push(state(st.l, st.idx-1, idx, v[idx]));
            }
            if (st.idx < st.r) {
                idx = rmq(st.idx+1, st.r);
                ASSERT_TRUE(idx >= st.idx+1); ASSERT_TRUE(idx <= st.r);
                ASSERT_TRUE(v[idx] >= v[st.idx])
                        << "v["<<idx<<"]="<< v[idx]
                        << " < " << "v["<<st.idx<<"]="
                        << v[st.idx] << endl
                        << "[" << st.l << "," << st.r << "]" << endl;
                s.push(state(st.idx+1, st.r, idx, v[idx]));
            }
        }
    }
}


TYPED_TEST(RMQTest, DeleteTest)
{
    sdsl::remove(temp_file);
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 3) {
        // LCOV_EXCL_START
        cout << "Usage: " << argv[0] << " test_file temp_file" << endl;
        cout << " (1) Generates a RMQ out of the serialized int_vector in test_file." << endl;
        cout << "     in test_file. Stores RMQ in temp_file." << endl;
        cout << " (2) Performs tests." << endl;
        cout << " (3) Deletes temp_file." << endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    test_file = argv[1];
    temp_file = argv[2];

    return RUN_ALL_TESTS();
}
