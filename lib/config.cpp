#include "sdsl/config.hpp"
#include "sdsl/util.hpp"

namespace sdsl{
	cache_config::cache_config(bool f_delete_files, std::string f_dir, std::string f_id, tMSS f_file_map) : delete_files(f_delete_files), dir(f_dir), id(f_id), file_map(f_file_map) { 
		if ( "" == id ){
			id = util::to_string(util::get_pid())+"_"+util::to_string(util::get_id());
		}
	}

	template<>
	const char* key_text_trait<0>::KEY_TEXT = constants::KEY_TEXT_INT;
	template<>
	const char* key_text_trait<8>::KEY_TEXT = constants::KEY_TEXT;

	template<>
	const char* key_bwt_trait<0>::KEY_BWT = constants::KEY_BWT_INT;
	template<>
	const char* key_bwt_trait<8>::KEY_BWT = constants::KEY_BWT;

}// end namespace sdsl
