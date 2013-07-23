#include <sdsl/suffix_arrays.hpp>
#include "gtest/gtest.h"
#include <iostream>
#include <vector>
#include <string>

namespace
{

using namespace sdsl;
using namespace std;

typedef sdsl::int_vector<>::size_type size_type;

string test_file;
string test_file_rev;

template<class T>
class SearchBidirectionalTest : public ::testing::Test { };

// Compare bidirectional search and backward search
template<class csa_type>
void test_bidirectional()
{
    bool debug = false;

    csa_type csa1;
    csa_type csa1_rev;
    construct(csa1, test_file.c_str(), 1);
    construct(csa1_rev, test_file_rev.c_str(), 1);

    bool patterntest = true;
    bool error = false;

    std::mt19937_64 rng(13);
    std::uniform_int_distribution<uint64_t> distribution(0, csa1.size()-1);


    for (size_type h = 0; h<1000; ++h) {
        //search for existing pattern forward and backward using bidirectional_search:
        size_type x = 4; // number of characters that are added to the pattern in each step
        size_type steps = 10; // maximal number of patternparts that are searched for
        if (debug) cout << "forward and backward search: alternating every " << x << " characters" << endl;
        size_type start = distribution(rng); //inclusive
        size_type end = start;  //exclusive
        bool forward = false;

        size_type l_rev = 0;
        size_type r_rev = csa1_rev.size()-1;
        size_type l = 0;
        size_type r = csa1.size()-1;
        size_type intervalsize = csa1.size();
        size_type i,pos;

        // alternating forward and backward search using bidirectional_search: alternating every x characters
        for (size_type rep = 0; rep<steps and start>=x and end<=csa1.size()-x ; ++rep) {
            string newpat = "";
            if (forward) {
                //forward
                if (debug) cout << "FORWARD" << endl;
                i = 0;
                pos = end;
                if (patterntest) {
                    if (debug) cout << "new pattern part: ";
                    for (size_type j=0; j<x; ++j) {
                        newpat.push_back(csa1.text[pos+j]);
                        if (debug) cout << newpat[j];
                    }
                    if (debug) cout << endl;
                    intervalsize = bidirectional_search_forward(csa1, csa1_rev, l, r, l_rev, r_rev, newpat.begin(), newpat.end(), l, r, l_rev, r_rev);
                    i = newpat.size();
                } else {
                    while (i < x and intervalsize) {
                        intervalsize = bidirectional_search(csa1_rev, l_rev, r_rev, l, r, csa1.text[pos+i], l_rev, r_rev, l, r);
                        ++i;
                    }
                }
                end += i;
            } else {
                //backward
                if (debug) cout << "BACKWARD" << endl;
                i = 0;
                pos = start-1;
                if (patterntest) {
                    if (debug) cout << "new pattern part: ";
                    for (size_type j=0; j<x; ++j) {
                        newpat.push_back(csa1.text[pos-x+1+j]);
                        if (debug) cout << newpat[j];
                    }
                    if (debug) cout << endl;
                    intervalsize = bidirectional_search_backward(csa1, csa1_rev, l, r, l_rev, r_rev, newpat.begin(), newpat.end(), l, r, l_rev, r_rev);
                    i = newpat.size();
                } else {
                    while (i < x and intervalsize) {
                        intervalsize = bidirectional_search(csa1, l, r, l_rev, r_rev, csa1.text[pos-i], l, r, l_rev, r_rev);
                        ++i;
                    }
                }
                start -= i;
            }

            //output
            if (debug) {
                cout << "pattern: ";
                for (size_type j=start; j<end; ++j) {
                    cout << csa1.text[j];
                }
                cout << " (at text[" << start << ".." << end-1 << "])" << endl;
                if (intervalsize) {
                    cout << "interval of pattern in csa1 is [" << l   << ".." << r   << "]" << endl;
                    cout << "interval of reverse pattern in csa1_rev is [" << l_rev << ".." << r_rev << "]" << endl;
                } else {
                    cout << "Pattern not found in input." << endl;
                }
            }
            ASSERT_EQ(true, (bool)intervalsize); // make sure pattern was found in input (it has to be because we took part of the input as pattern)

            {
                //check using backward_search
                string pat = "";
                if (debug) cout << "pattern: ";
                for (size_type j=0; j<end-start; ++j) {
                    pat.push_back(csa1.text[start+j]);
                    if (debug) cout << pat[j];
                }
                if (debug) cout << endl;
                size_type b_l,b_r;
                size_type res = backward_search(csa1, 0, csa1.size()-1, pat.begin(), pat.end(), b_l, b_r);
                if (debug) {
                    if (intervalsize == res and l == b_l and r == b_r) {
                        cout << "correct: Interval of pattern is correct." << endl;
                    } else {
                        error = true;
                        cout << "error: Intervals of pattern differ." << endl;
                        if (intervalsize != res) cout << "intervalsize bidirectional = " << intervalsize << ", intervalsize backward = " << res << endl;
                        if (l != b_l) cout << "l bidirectional = " << l << ", b_l backward = " << b_l << endl;
                        if (r != b_r) cout << "r bidirectional = " << r << ", b_r backward = " << b_r << endl;
                    }
                }
                ASSERT_EQ(res, intervalsize);
                ASSERT_EQ(b_l, l);
                ASSERT_EQ(b_r, r);
            }

            //change direction
            forward = !forward;
        }
    }
    if (debug) {
        if (error) cout << "error" << endl;
        else cout << "all tests correct" << endl;
    }
}

//! Test Bidirectional Search
TEST(SearchBidirectionalTest, Bidirectional1)
{
    typedef sdsl::csa_wt<wt<>, 32, 32, sa_order_sa_sampling<>, int_vector<>, byte_alphabet> csa_type;
    test_bidirectional< csa_type >();
}
//! Test Bidirectional Search
TEST(SearchBidirectionalTest, Bidirectional2)
{
    typedef sdsl::csa_wt<wt<>, 32, 32, sa_order_sa_sampling<>, int_vector<>, succinct_byte_alphabet<> > csa_type;
    test_bidirectional< csa_type >();
}
//! Test Bidirectional Search
TEST(SearchBidirectionalTest, Bidirectional3)
{
    typedef sdsl::csa_wt<wt_hutu<>, 32, 32, sa_order_sa_sampling<>, int_vector<>, byte_alphabet> csa_type;
    test_bidirectional< csa_type >();
}
//! Test Bidirectional Search
TEST(SearchBidirectionalTest, Bidirectional4)
{
    typedef sdsl::csa_wt<wt_hutu<>, 32, 32, sa_order_sa_sampling<>, int_vector<>, succinct_byte_alphabet<> > csa_type;
    test_bidirectional< csa_type >();
}
//! Test Bidirectional Search
TEST(SearchBidirectionalTest, Bidirectional5)
{
    typedef sdsl::csa_wt<wt_hutu<bit_vector_il<> >, 32, 32, sa_order_sa_sampling<>, int_vector<>, byte_alphabet> csa_type;
    test_bidirectional< csa_type >();
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    if (argc!=2) {
        cout << "Usage: " << argv[0] << " test_file" << endl;
        cout << " (1) Reverses test_file; stores it in test_file_rev." << endl;
        cout << " (2) Performs tests." << endl;
        cout << " (3) Deletes test_file_reverse." << endl;
        return 1;
    }
    test_file = argv[1];
    test_file_rev = test_file + "_rev";

    {
        //reverse input
        int_vector<8> text;
        load_vector_from_file(text, test_file, 1);
        size_type n = text.size();
        int_vector<8> text_rev(n);
        for (size_type i=0; i<n; i++) {
            text_rev[n-1-i] = text[i];
        }
        char* text2 = (char*)text_rev.data();
        ofstream of(test_file_rev.c_str(),ofstream::binary);
        of.write(text2,n);
        of.close();
    }
    int result = RUN_ALL_TESTS();
    sdsl::remove(test_file_rev);
    return result;
}
