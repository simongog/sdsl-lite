#include <iostream>
#include <vector>
#include <sdsl/vectors.hpp>
#include <fstream>
#include <sdsl/bits.hpp>
#include <sdsl/coder_elias_delta.hpp>
#include <sdsl/suffix_arrays.hpp>

using namespace std;
using namespace sdsl;
using namespace sdsl::coder;

int main(int argc, char* argv[])
{	
	int numBytes = 1;
	
	if (argc < 2) {
		cout<< "Bitte eine input Datein angeben\n";
		return 0;
	}
	if (argc == 3)
		numBytes = atoi(argv[2]);
	
	construct_config::byte_algo_sa = SE_SAIS;
	csa_sada<> csa_enc;

	construct(csa_enc, argv[1], numBytes);		
	
    // COUNT TEST
    //string pattern = "size_in_bytes";
    vector<string> patterns;
    patterns.push_back("aa");
    patterns.push_back("a");
    patterns.push_back("e");
    patterns.push_back("r");
    patterns.push_back("abracadabrabarbara");
    patterns.push_back("abracadabrabarbaraa");
    patterns.push_back("abracadabrabarbarae");
    patterns.push_back("abracadabrabarbaraaaa");
    patterns.push_back("");
    for (int i = 0; i < patterns.size(); i++)
    {
		string pattern = patterns[i];
		
		cout << "TEST PATTERN: " << pattern << endl;
		cout << "CSA-ENC-COUNT: " << count(csa_enc, pattern.begin(), pattern.end()) << endl;
		
		unsigned long int res_l, res_r, res_size;
		
		res_size = backward_search(csa_enc, 
			0, csa_enc.size() - 1, 
			pattern.begin(), pattern.end(),
			res_l, res_r);
		cout << "BWS: " << res_l << " " << res_r << " = " << res_size << endl;
		res_size = forward_search(csa_enc, 
			0, csa_enc.size() - 1, 
			pattern.begin(), pattern.end(),
			res_l, res_r);
		cout << "FWS: " << res_l << " " << res_r << " = " << res_size << endl;
	}

	return 0;
}
