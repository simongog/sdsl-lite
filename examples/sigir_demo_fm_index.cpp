#include <sdsl/csa_wt.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <string>
#include <iostream>
#include <algorithm>

using namespace sdsl;
using namespace std;

int main(int argc, char **argv){
	if( argc <  2 ){
		cout << "Usage " << argv[0] << " text_file [max_locations] [extract_length]" << endl;
		cout << "      This program constructs a very compact FM-index" << endl;
		cout << "      which supports count, locate, and extract queries." << endl;
		cout << "      text_file      Original text file." << endl;
		cout << "      max_locations  Maximal number of location to report." <<endl;
		cout << "      extract_length Maximal length of the reported context." << endl;
		return 1;
	}
	size_t max_locations = 5;
	size_t extract_length = 10;
	if ( argc >= 3 ){ max_locations = atoi(argv[2]); }
	if ( argc >= 4 ){ extract_length = atoi(argv[3]); }
	string index_suffix = ".fm9";
	string index_file   = string(argv[1]) + index_suffix;
	csa_wt<wt_huff<rrr_vector<255> >, 512, 1024> fm_index;

	if( !util::load_from_file(fm_index, index_file.c_str()) ){
		// generate index
		construct_csa(argv[1], fm_index);
		// save it
		util::store_to_file(fm_index, index_file.c_str());
	}
	cout << "File is indexed. Input search terms and press Ctrl-D to exit." << endl;
	string query;
	while ( getline(cin, query) ){
		size_t occs = algorithm::count(fm_index, (const unsigned char*)query.c_str(), query.size());
		cout << "# of occcurences: " << occs << endl;
		if ( occs > 0 ){
			cout << "Location and context of first occurences: " << endl;
			int_vector<> locations;
			algorithm::locate(fm_index, (const unsigned char*)query.c_str(), query.size(), locations);
			sort(locations.begin(), locations.end());
			for(size_t i = 0; i < min(occs, max_locations); ++i){
				cout << locations[i] << ": ";
				string extract = algorithm::extract(fm_index, locations[i], 
						                            min(locations[i] + query.size() - 1 + extract_length, fm_index.size()-1) );
				cout << "\e[1;31m"; 
				cout << extract.substr(0, query.size());
				cout << "\e[0m";
				string context = extract.substr(query.size());
				cout << context.substr(0, context.find_first_of('\n')) << endl;
			}
		}
	}
}

