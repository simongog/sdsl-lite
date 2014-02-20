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

bool store_to_checked_file(const char* v, const std::string& file)
{
    std::string checkfile = file+"_check";
    osfstream out(checkfile, std::ios::binary | std::ios::trunc | std::ios::out);
    if (!out) {
        if (util::verbose) {
            std::cerr<<"ERROR: store_to_checked_file(const char *v, const std::string&)"<<std::endl;
            return false;
        }
    }
    add_hash(v, out);
    out.close();
    return store_to_file(v, file);
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

uint64_t _parse_number(std::string::const_iterator& c, const std::string::const_iterator& end)
{
    std::string::const_iterator s = c;
    while (c != end and isdigit(*c)) ++c;
    if (c > s) {
        return stoull(std::string(s,c));
    } else {
        return 0;
    }
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

std::string tmp_file(const std::string& filename, std::string name_part)
{
    return util::dirname(filename) + "/" + util::to_string(util::pid()) + "_" +
           util::to_string(util::id()) + name_part + ".sdsl";
}

}// end namespace sdsl

