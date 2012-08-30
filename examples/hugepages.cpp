#include <sdsl/int_vector.hpp>
#include <sdsl/suffixtrees.hpp>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/testutils.hpp>
#include <sdsl/bit_vector_interleaved.hpp>
#include <iostream>
#include <string>

using namespace sdsl;
using namespace std;

template<class tCsa>
void do_something(const tCsa &csa){
	stop_watch sw;
	uint64_t sum=0;
	sw.start();
	for(size_t i=0; i<csa.size() and i<10000000;++i){
		sum+=csa.psi(i);
	}
	sw.stop();
	cout << sw.get_real_time() << endl;
	cout<<"sum="<<sum<<endl;
}

int main(int argc, char** argv){
//	util::verbose = true;
	{
		bit_vector b(32,1);
		cout << b << endl;
	}
	{
		bit_vector b(32,0);
		bool mapped = mm::map_hp();
		if ( mapped ) mm::unmap_hp();
		cout << b << endl;
	}
	bool mapped;
//	bool mapped = mm::map_hp();
//	if( mapped ) mm::unmap_hp();

	csa_wt<wt_huff<bit_vector_interleaved<> > > csa;
	construct_csa(string(argv[1]), csa);
	do_something(csa); // before it is mapped
	mapped = mm::map_hp();
	do_something(csa); // while it is mapped
	if ( mapped ) mm::unmap_hp();
	do_something(csa); // after it is unmapped 
	util::clear(csa);
}
