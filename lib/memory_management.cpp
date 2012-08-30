#include "sdsl/memory_management.hpp"

#include <cstdlib> // for malloc and free
#include <sys/mman.h>

#define HUGE_LEN 1073741824 
#define HUGE_PROTECTION (PROT_READ | PROT_WRITE)
#define HUGE_FLAGS (MAP_HUGETLB | MAP_ANONYMOUS | MAP_PRIVATE)

//! Namespace for the succinct data structure library
namespace sdsl
{
	static int nifty_counter;
	std::map<uint64_t, mm_item_base*> mm::m_items;
	uint64_t mm::m_total_memory;
	uint64_t *mm::m_data;

	mm_initializer::mm_initializer(){
		if ( 0 == nifty_counter++ ){
			// initialize static members object here
			// mm::m_items.clear();
			mm::m_total_memory = 0;
			mm::m_data = NULL;
		}
	}
	mm_initializer::~mm_initializer(){
		if ( 0 == --nifty_counter ){
			// clean up
		}
	}

	bool mm::map_hp(){
		m_total_memory = 0; // memory of all int_vectors
		for(tMVecItem::const_iterator it=m_items.begin(); it!=m_items.end(); ++it){
			m_total_memory += it->second->size();
		}				
		if(util::verbose){
			std::cout<<"m_total_memory"<<m_total_memory<<std::endl;
		}
		size_t hpgs= (m_total_memory+HUGE_LEN-1)/HUGE_LEN; // number of huge pages required to store the int_vectors 
		m_data = (uint64_t*)mmap(NULL, hpgs*HUGE_LEN, HUGE_PROTECTION, HUGE_FLAGS, 0, 0);
		if (m_data == MAP_FAILED) {
			std::cout << "mmap was not successful" << std::endl;
			return false;
		}else{
			if( util::verbose ){
				std::cerr<<"map " << m_total_memory << " bytes" << std::endl; 
			}
		}
		// map int_vectors
		uint64_t *addr = m_data;
		bool success = true;
		for(tMVecItem::const_iterator it=m_items.begin(); it!=m_items.end(); ++it){
//			std::cerr<<"addr = "<< addr << std::endl;
			success = success && it->second->map_hp( addr );
		}
		return success;
	}

	bool mm::unmap_hp(){
		size_t hpgs= (m_total_memory+HUGE_LEN-1)/HUGE_LEN; // number of huge pages 
		if (  util::verbose ){
			std::cerr<<"unmap "<< m_total_memory << " bytes" <<std::endl;
			std::cerr.flush();
		}
		bool success = true;
		for(tMVecItem::const_iterator it=m_items.begin(); it!=m_items.end(); ++it){
			success = success && it->second->unmap_hp();
		}
//		uint64_t* tmp_data = (uint64_t*)malloc(m_total_memory); // allocate memory for int_vectors 
//		memcpy(tmp_data, m_data, len); // copy data from the mmapped region
		int ret = munmap((void*)m_data, hpgs*HUGE_LEN ); 
		if ( ret == -1 ){
			perror("Unmap failed");
			return false;
		}
		return success;
	}

} // end namespace
