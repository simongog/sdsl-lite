#include <sdsl/int_vector.hpp>
#include <sdsl/select_support.hpp>
#include <future>
#include <iostream>

using namespace sdsl;

size_t
do_work(size_t i,size_t n)
{
    size_t count = 0;
    /* build */
    //for(size_t k=0;k<i;k++) {
    for (size_t l=0; l<n; l++) {
        size_t steps = 500;
        size_t m = 765432;
        bit_vector bv(765432);
        for (size_t j=0; j<m; j+=steps) {
            bv[j] = 1;
        }
        select_support_mcl<1> select_one(&bv);
        /* query */

        for (size_t j=1; j<(m/steps)-5; j++) {
            count += select_one.select(j);
        }
    }
    //}
    return count;
}


int main(int argc, char** argv)
{
    std::vector< std::future<size_t> > v;
    srand(getpid());

    size_t n = 6234;
    srand(4711);
    for (size_t i=0; i<10; i++) {
        v.push_back(std::async(std::launch::async,[&,i] {
            return do_work(2+i,1+rand()%n);
        }));
    }

    size_t count = 0;
    for (auto& f : v) count += f.get();
    std::cout << "count = " << count << std::endl;

}
