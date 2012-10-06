/*!\file memory_management.hpp
   \brief memory_management.hpp contains two function for allocating and deallocating memory
   \author Your name
*/
#ifndef INCLUDED_MEMORY_MANAGEMENT 
#define INCLUDED_MEMORY_MANAGEMENT

#include "uintx_t.hpp"
#include "util.hpp"
#include <map>
#include <iostream>
using std::cout;
using std::endl;

namespace sdsl{
	
class mm_item_base{
	public:
		mm_item_base(){};
		virtual bool map_hp(uint64_t*&){return false;};// const = 0;
		virtual bool unmap_hp(){return false;};// const = 0;
		virtual ~mm_item_base(){};
		virtual uint64_t size(){return 0;};
};

template<class int_vector_type>
class mm_item : public mm_item_base{
	private:
		int_vector_type *m_v;
	public:
		explicit mm_item(int_vector_type *v):m_v(v){}
		~mm_item(){ }

		//! Map content of int_vector to a hugepage starting at address addr
		/*! 
		 *  Details: The content of the corresponding int_vector of mm_item 
		 *           is copied into a hugepage starting at address addr.
		 *           The string position in the hugepage is then increased
		 *           by the number of bytes used by the int_vector.
		 *           So addr is the new starting address for the next
		 *           mm_item which have to be mapped.
		 */
		bool map_hp(uint64_t*& addr){
			uint64_t len = size();
			if ( m_v->m_data != NULL ){
				memcpy((char*)addr, m_v->m_data, len); // copy old data
				free(m_v->m_data);
				m_v->m_data = addr;
				addr += (len/8);		
			}
			return true;	
		}

		//! 
		bool unmap_hp(){
			uint64_t len = size();
			if (  util::verbose ){
				std::cerr<<"unmap int_vector of size "<< len <<std::endl;
				std::cerr<<"m_data="<<m_v->m_data<<std::endl;
			}
			uint64_t* tmp_data = (uint64_t*)malloc(len); // allocate memory for m_data
			memcpy(tmp_data, m_v->m_data, len); // copy data from the mmapped region
			m_v->m_data = tmp_data;
			return true;
		}

		uint64_t size(){
			return ((m_v->bit_size()+63)>>6)<<3;
		}
};

class mm_initializer; // forward declaration of initialization helper

// memory management class
class mm{
	friend class mm_initializer;
	typedef std::map<uint64_t, mm_item_base*> tMVecItem;
	static tMVecItem m_items; 
	static uint64_t m_total_memory;
	static uint64_t *m_data;
	public:
		mm();

		template<class int_vector_type>
		static void add(int_vector_type *v){
			if( mm::m_items.find((uint64_t)v) == mm::m_items.end() ){
				mm_item_base* item = new mm_item<int_vector_type>(v); 
				if(false and util::verbose) { 
					cout << "mm::add: add vector " << v << endl;
					cout.flush();
				}
				mm::m_items[(uint64_t)v] = item;
			}else{
				if(false and util::verbose) cout << "mm::add: mm_item is already in the set" << endl;
			}
		}

		template<class int_vector_type>
		static void remove(int_vector_type *v){
			if( mm::m_items.find((uint64_t)v) != mm::m_items.end() ){
				if( false and util::verbose ){ cout << "mm:remove: remove vector " << v << endl; };
				mm_item_base* item = m_items[(uint64_t)v];
				mm::m_items.erase((uint64_t)v);
				delete item;
			}else{
				if( false and util::verbose ){ cout << "mm:remove: mm_item is not in the set" << endl; };
			}
		}

		static bool map_hp();
		static bool unmap_hp();
};

static class mm_initializer{
	public:
		mm_initializer ();
		~mm_initializer ();
} initializer;

} // end namespace

#endif
