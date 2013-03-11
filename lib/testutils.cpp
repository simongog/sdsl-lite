/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
#include "sdsl/testutils.hpp"
#include <cxxabi.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <unistd.h> // for file_size, also contains clock_gettime 


namespace sdsl
{

timeval stop_watch::m_first_t = {0,0};
rusage stop_watch::m_first_r = {{0,0},{0,0}};

void stop_watch::start()
{
    gettimeofday(&m_timeOfDay1, 0);
    getrusage(RUSAGE_SELF, &m_ruse1);
}

void stop_watch::stop()
{
    getrusage(RUSAGE_SELF, &m_ruse2);
    gettimeofday(&m_timeOfDay2, 0);
}

double stop_watch::get_user_time()
{
    timeval t1, t2;
    t1 = m_ruse1.ru_utime;
    t2 = m_ruse2.ru_utime;
    return ((double)(t2.tv_sec*1000000 + t2.tv_usec - (t1.tv_sec*1000000 + t1.tv_usec)))/1000.0;
}

double stop_watch::get_sys_time()
{
    timeval t1, t2;
    t1 = m_ruse1.ru_stime;
    t2 = m_ruse2.ru_stime;
    return ((double)(t2.tv_sec*1000000 + t2.tv_usec - (t1.tv_sec*1000000 + t1.tv_usec)))/1000.0;
}

double stop_watch::get_real_time()
{
    double result = ((double)((m_timeOfDay2.tv_sec*1000000 + m_timeOfDay2.tv_usec)-(m_timeOfDay1.tv_sec*1000000 + m_timeOfDay1.tv_usec)))/1000.0;
    if (result < get_sys_time() + get_user_time())
        return get_sys_time()+get_user_time();
    return result;
}

uint64_t stop_watch::get_abs_real_time()
{
    uint64_t result = (((m_timeOfDay2.tv_sec*1000000 + m_timeOfDay2.tv_usec - (m_first_t.tv_sec*1000000 + m_first_t.tv_usec))))/1000;
    return result;
}

uint64_t stop_watch::get_abs_user_time()
{
    timeval t1, t2;
    t1 = m_first_r.ru_utime;
    t2 = m_ruse2.ru_utime;
    return (t2.tv_sec*1000000 + t2.tv_usec - (t1.tv_sec*1000000 + t1.tv_usec))/1000;
}


uint64_t stop_watch::get_abs_sys_time()
{
    timeval t1, t2;
    t1 = m_first_r.ru_stime;
    t2 = m_ruse2.ru_stime;
    return (t2.tv_sec*1000000 + t2.tv_usec - (t1.tv_sec*1000000 + t1.tv_usec))/1000;
}

uint64_t stop_watch::get_abs_page_faults()
{
    return m_ruse2.ru_majflt - m_first_r.ru_majflt; // does not work on my platform
}

std::string clock::get_time_string()
{
    time_t rawtime;
    struct tm* timeinfo;
    char buffer[1024];
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer, 1024, "%Y-%m-%d-%H%M%S", timeinfo);
    return buffer;
}




void file::write_text(const char* file_name, const char* c, uint64_t len)
{
    std::ofstream out(file_name);
    if (out) {
        out.write(c, len);
        out.close();
    }
}

uint64_t file::read_text(const char* file_name, char*& c, bool trunc, uint64_t lim)
{
    if (c != NULL) {
        delete [] c;
        c = NULL;
    }
    uint64_t n = util::file_size(file_name) + 1; // add one for the 0 byte
    if (trunc and lim+1 < n) {
        n = lim+1;
    }
//std::cerr<<"file has size "<< n <<std::endl;
    std::ifstream in;
    in.open(file_name);
    if (in) {
        c = new char[n];
        c[n-1] = '\0';
        char* cp = c;
        in.read(cp, n-1);
        return n; // added 0 byte
    }
    return 0;
    /*
    	std::ifstream in;
    	in.open(file_name);
    	if( in ){
    		const uint64_t BLOCK_SIZE = (1<<20);
    		uint64_t n=0, read = 0;
    		char buf[BLOCK_SIZE], *cp;
    		do{
    			in.read(buf, BLOCK_SIZE);
    			read = in.gcount();
    			n+=read;
    		}while( BLOCK_SIZE == read );
    		if(n==0)
    			return 0;
    		c = new char[n+2+BLOCK_SIZE]; //TODO checken warum das nicht gut ist
    //		*(c+(n+1)) = 0;
    		in.close();
    		in.open(file_name);
    		if(!in){
    			delete [] c;
    			c = NULL;
    			return 0;
    		}
    		cp=c;
    		do{
    			in.read(cp, BLOCK_SIZE);
    			read = in.gcount();
    			cp+= read;
    		}while( BLOCK_SIZE == read );
    		*(c+n) = '\0';
    		return n;
    	}
    	else
    		return 0;
    */
}


std::vector<std::string> paths_from_config_file(const std::string &file, const char *prefix){
		std::ifstream config_in(file.c_str());
		if ( config_in ){ // opened file successfully
			std::vector<std::string> result;
			const size_t name_max_size = 1024;
			char * name = new char [name_max_size];
			while ( config_in.getline( name, name_max_size ) ) {
				if ( strlen(name) > 0 and '#' != name[0] ){ // check empty line and comment
					std::string path = std::string(name);
					if ( prefix != NULL ){
						path = std::string(prefix) + "/" + path;
					}
					result.push_back( path );
				}
			}
			delete [] name;
			return result;
		} else {
			std::cerr << "WARNING: Could not open config file: `";
			std::cerr << file << "`" << std::endl;
			return std::vector<std::string>();
		}
}

} // end of namespace sdsl
