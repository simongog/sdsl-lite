#include "sdsl/config.hpp"
#include "sdsl/util.hpp"

namespace sdsl{
	cache_config::cache_config(bool f_delete_files, std::string f_dir, std::string f_id, tMSS f_file_map) : delete_files(f_delete_files), dir(f_dir), id(f_id), file_map(f_file_map) { 
		if ( "" == id ){
			id = util::to_string(util::get_pid())+"_"+util::to_string(util::get_id());
		}
	}
}// end namespace sdsl
