#include "sdsl/io.hpp"
#include "sdsl/sfstream.hpp"
#include "sdsl/util.hpp"
#include <vector>

namespace sdsl
{


bool store_to_file(const char* v, const std::string& file)
{
    osfstream out(file, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out) {
        if (util::verbose) {
            std::cerr<<"ERROR: store_to_file(const char *v, const std::string&)"<<std::endl;
            return false;
        }
    }
    uint64_t n = strlen((const char*)v);
    out.write(v, n);
    out.close();
    return true;
}


template<>
size_t write_member<std::string>(const std::string& t, std::ostream& out, structure_tree_node* v, std::string name)
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(t));
    size_t written_bytes = 0;
    written_bytes += write_member(t.size(), out, child, "length");
    out.write(t.c_str(), t.size());
    written_bytes += t.size();
    structure_tree::add_size(v, written_bytes);
    return written_bytes;
}

template<>
void read_member<std::string>(std::string& t, std::istream& in)
{
    std::string::size_type size;
    read_member(size, in);
    char* buf = new char[size];
    in.read(buf, size);
    std::string temp(buf, size);
    delete [] buf;
    t.swap(temp);
}

template<>
bool load_from_file(void*& v, const std::string& file_name)
{
    return true;
}

bool load_from_file(char*& v, const std::string& file_name)
{
    if (v != NULL) {
        delete [] v;
        v = NULL;
    }
    isfstream in(file_name, std::ios::binary | std::ios::in);
    if (in) {
        const uint64_t SDSL_BLOCK_SIZE = (1<<20);
        uint64_t n=0, read = 0;
        char buf[SDSL_BLOCK_SIZE], *cp;
        do {
            in.read(buf, SDSL_BLOCK_SIZE);
            read = in.gcount();
            n+=read;
        } while (SDSL_BLOCK_SIZE == read);
        if (n==0)
            return false;
        v = new char[n+1];
        in.close();
        in.open(file_name.c_str());
        if (!in) {
            delete [] v;
            v = NULL;
            return false;
        }
        cp=v;
        do {
            in.read(cp, SDSL_BLOCK_SIZE);
            read = in.gcount();
            cp+= read;
        } while (SDSL_BLOCK_SIZE == read);
        *(v+n) = '\0';
        return true;
    } else
        return false;
}

std::string cache_file_name(const std::string& key, const cache_config& config)
{
    return config.dir+"/"+key+"_"+config.id+".sdsl";
}

void register_cache_file(const std::string& key, cache_config& config)
{
    std::string file_name = cache_file_name(key, config);
    isfstream in(file_name);
    if (in) {  // if file exists, register it.
        config.file_map[key] = file_name;
    }
}


bool cache_file_exists(const std::string& key, const cache_config& config)
{
    std::string file_name = cache_file_name(key, config);
    isfstream in(file_name);
    if (in) {
        in.close();
        return true;
    }
    return false;
}

std::string tmp_file(const cache_config& config, std::string name_part)
{
    return config.dir+"/"+ util::to_string(util::pid()) + "_" + util::to_string(util::id()) + name_part + ".sdsl";
}

}// end namespace sdsl

