#include "sdsl/construct.hpp"

namespace sdsl{

	bool contains_no_zero_symbol(const int_vector<> &text, const char* file){
		for (int_vector_size_type i=0; i < text.size(); ++i){
			if ( (uint64_t)0 == text[i] ){
				throw std::logic_error((std::string("Error: File \"")+std::string(file)+std::string("\" contains zero symbol.")).c_str());
				return false;
			}
		}
		return true;
	}

	void append_zero_symbol(int_vector<> &text){
		text.resize(text.size()+1);
		text[text.size()-1] = 0;
	}
}
