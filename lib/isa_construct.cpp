/* sdsl - succinct data structures library
    Copyright (C) 2010 Simon Gog 
*/	
#include "sdsl/isa_construct.hpp"

namespace sdsl{

	bool construct_isa( tMSS &file_map, const std::string &dir, const std::string &id){
		typedef int_vector<>::size_type size_type;
		// TODO: bei csa_uncompressed ist die Berechnung des isa nicht notwendig; die Speicherung schon
		if( file_map.find("isa") == file_map.end() ){ // if isa is not already on disk => calculate it
			write_R_output("cst", "construct ISA", "begin", 1, 0);
			int_vector<> isa;
			if( !util::load_from_file(isa, file_map["sa"].c_str() ) ){
				throw std::ios_base::failure("cst_construct: Cannot load SA from file system!");
			}
			{
				int_vector_file_buffer<> sa_buf(file_map["sa"].c_str());

				for(size_type i=0, r_sum=0, r = sa_buf.load_next_block(); r_sum < isa.size(); ){
					for(; i<r_sum+r; ++i){
						isa[ sa_buf[i-r_sum] ] = i;
					}
					r_sum += r; r = sa_buf.load_next_block();
				}
			}
			if( !util::store_to_file(isa, (dir+"isa_"+id).c_str() ) ){
				throw std::ios_base::failure( "cst_construct: Cannot store ISA to file system!" );
				return false;
			}else{
				file_map["isa"] = dir+"isa_"+id;
			}
			write_R_output("cst", "construct ISA", "end", 1, 0);
/*			
			int_vector<32> isa2(isa.size());
			isa.resize(0);
			write_R_output("cst", "construct ISA2", "begin", 1, 0);
			{
				int_vector_file_buffer<> sa_buf(file_map["sa"].c_str());

				for(size_type i=0, r_sum=0, r = sa_buf.load_next_block(); r_sum < isa2.size(); ){
					for(; i<r_sum+r; ++i){
						isa2[ sa_buf[i-r_sum] ] = i;
					}
					r_sum += r; r = sa_buf.load_next_block();
				}
			}
			if( !util::store_to_file(isa2, (dir+"isa2_"+id).c_str() ) ){
				throw std::ios_base::failure( "cst_construct: Cannot store ISA to file system!" );
				return false;
			}else{
				file_map["isa2"] = dir+"isa2_"+id;
			}
			write_R_output("cst", "construct ISA2", "end", 1, 0);
*/			
		}
		return true;
	}	

	// TODO methode die den isa auf der bwt berechnet

}// end namespace

