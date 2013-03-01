#include <sdsl/construct.hpp>
#include <sdsl/suffixtrees.hpp>
#include <sdsl/algorithms.hpp>
#include <string>

using namespace sdsl;
using namespace std;

typedef csa_wt<wt_int<rrr_vector<63> >, 32, 128, sa_order_sa_sampling<>, int_vector<>, int_alphabet_strategy<> > tCsa;
//typedef csa_wt<wt_int<>, 32, 128, sa_order_sa_sampling<>, int_vector<>, int_alphabet_strategy<> > tCsa;
//typedef csa_wt<wt_int<>, 32, 16, sa_order_sa_sampling<>, int_vector<>, int_alphabet_strategy<> > tCsa;
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
	cout<<"         i        BWT        PSI         LF         SA   TEXT\n";
	for(int_vector_size_type i=sp; i<=ep; ++i){
		cout<<setw(10)<<i<<" "<<setw(10)<<csa.bwt[i]
			             <<" "<<setw(10)<<csa.psi[i]
			             <<" "<<setw(10)<<csa.psi(i)
						 <<" "<<setw(10)<<csa[i] << "  ";
		print_suffix(csa, i, 10, cout);
		cout<<endl;
	}
}

/*
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
*/
int main(int argc, char *argv[]){
	if ( argc < 2 ){
		cout << "Usage: "<<argv[0]<<" file"<<endl;
		cout <<"        Format of the file sequence of 64-bit integers."<<endl;
	}
/*	
	tCsa csa;
	{
		string csa_file = string(argv[1])+"."+util::class_to_hash(csa)+".csa";
		ifstream in(csa_file.c_str());
		if ( !in ){
			cache_config config(false, "./tmp/", util::basename(argv[1]));
			construct(csa, argv[1], config, 8);
			util::store_to_file(csa, csa_file.c_str());
		}else{
			util::load_from_file(csa, csa_file.c_str());
		}
	}
	cout << "csa.size() = " << csa.size() << endl;
	cout << "csa.sigma = " << csa.sigma << endl;
	cout << "csa.rank_bwt(csa.size(), 1) = " << csa.rank_bwt(csa.size(), 1) << endl;
*/	
/*	
	cout << "Test permutation property of LF" << endl;
	bit_vector visited(csa.size(),0);
	bool valid = true;
	for (size_type i=0; i<csa.size() and valid; ++i){
		size_type lf = csa.psi(i);
		if ( !visited[lf] ){ 
			visited[lf] = 1; 
		}else{
			cout<<"csa.psi("<<i<<")="<<lf<<" but was already mapped"<<endl;
			valid = false;
		}
	}
	cout << "Test permutation property of Psi" << endl;
	util::assign(visited, bit_vector(csa.size(),0));
	valid = true;
	for (size_type i=0; i<csa.size() and valid; ++i){
		size_type lf = csa.psi[i];
		if ( !visited[lf] ){ 
			visited[lf] = 1; 
		}else{
			cout<<"csa.psi["<<i<<"]="<<lf<<" but was already mapped"<<endl;
			valid = false;
		}
	}
	cout<<"done"<<endl;
*/	
	
	tCst cst;
	{
		string cst_file = string(argv[1])+"."+util::class_to_hash(cst)+".cst";
		ifstream in(cst_file.c_str());
		if ( !in ){
			cache_config config(false, "./tmp/", util::basename(argv[1]));
			construct(cst, argv[1], config, 8);
			util::store_to_file(cst, cst_file.c_str());
		}else{
			util::load_from_file(cst, cst_file.c_str());
		}
	}
	
//	cout << "cst.size() = " << csa.size() << endl;
//	util::write_structure<JSON_FORMAT>(cst, cout);
	cout << "\nk     Hk      context/n\n";
	for(size_t i=0; i<10; ++i){
		size_type context;
		double x = Hk(cst, i, context);
		cout << i<<" "<<x<<" "<<((double)context)/cst.csa.size()<<endl;
	}
	cout << endl;
/*	
	size_type sp,ep;
	while ( cin >> sp >> ep ){
		print_csa(csa, sp, ep);
	}
*/	
}
