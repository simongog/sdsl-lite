#include <sdsl/construct.hpp>
#include <sdsl/suffixtrees.hpp>
#include <sdsl/algorithms.hpp>

using namespace sdsl;
using namespace std;

typedef csa_wt<wt_int<rrr_vector<63> >, 32, 16, sa_order_sa_sampling<>, int_vector<>, int_alphabet_strategy<> > tCsa;
typedef cst_sada<tCsa, lcp_dac<> > tCst;
typedef int_vector<>::size_type size_type;
typedef tCsa::char_type char_type;

void print_suffix(const tCsa &csa, size_type i, size_type extract_len, std::ostream &out){
	size_type len=0;
	char_type *buf = new char_type[extract_len+2];
	algorithm::extract(csa, csa[i], csa[i]+extract_len, buf, len);
	for (size_type j=0; j<len; ++j){ if(sizeof(char_type)>1) out<<" "; out<<buf[j]; }
	if (csa[i]+extract_len <= csa.size()){ out<<"..."; }
	delete [] buf;
}

void print_csa(const tCsa &csa, size_type sp, size_type ep){
	cout<<" i   BWT  PSI   LF    SA    TEXT\n";
	for(int_vector_size_type i=sp; i<=ep; ++i){
		cout<<setw(3)<<i<<" "<<setw(4)<<csa.bwt[i]
			            <<" "<<setw(4)<<csa.psi[i]
			            <<" "<<setw(4)<<csa.psi(i)
						<<" "<<setw(4)<<csa[i] << "  ";
		print_suffix(csa, i, 10, cout);
		cout<<endl;
	}
}


void print_cst(const tCst &cst, size_type sp, size_type ep){
	cout<<" i   BWT  PSI   LF  LCP   SA    TEXT\n";
	for (size_type i=sp; i<=ep; ++i){
		cout<<setw(3)<<i<<" "<<setw(4)<<cst.csa.bwt[i]
			            <<" "<<setw(4)<<cst.csa.psi[i]
			            <<" "<<setw(4)<<cst.csa.psi(i)
			            <<" "<<setw(4)<<cst.lcp[i]
						<<" "<<setw(4)<<cst.csa[i] << "  ";
		print_suffix(cst.csa, i, 10, cout);
		cout<<endl;
	}
}

int main(int argc, char *argv[]){
	if ( argc < 2 ){
		cout << "Usage: "<<argv[0]<<" file"<<endl;
		cout <<"        Format of the file sequence of 64-bit integers."<<endl;
	}
	tCsa csa;
	construct(csa, argv[1], 8);
	cout << "csa.size() = " << csa.size() << endl;
	tCst cst;
	construct(cst, argv[1], 8);
	cout << "cst.size() = " << csa.size() << endl;
	size_type sp,ep;
	while ( cin >> sp >> ep ){
		print_cst(cst, sp, ep);
	}
}
