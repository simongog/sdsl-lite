#include "sdsl/ram_fs.hpp"
#include <cstdio>
#include <iostream>

static int nifty_counter = 0;

sdsl::ram_fs::mss_type sdsl::ram_fs::m_map;

sdsl::ram_fs_initializer::ram_fs_initializer()
{
    if (0 == nifty_counter++) {
        ram_fs::m_map = ram_fs::mss_type();
    }
}

sdsl::ram_fs_initializer::~ram_fs_initializer()
{
    if (0 == --nifty_counter) {
        // clean up
    }
}

namespace sdsl
{

ram_fs::ram_fs() {};

void
ram_fs::store(const std::string& name, const std::string& data)
{
    m_map[name] = data;
}

const std::string&
ram_fs::content(const std::string& name)
{
    std::cout<<"content of `"<<name<<"`"<<std::endl;
    return m_map[name];
}

int
ram_fs::remove(const std::string& name)
{
    m_map.erase(name);
    return 0;
}

int
ram_fs::rename(const std::string old_filename, const std::string new_filename)
{
    // TODO: this is expensive with map
    m_map[new_filename] = m_map[old_filename];
    remove(old_filename);
    return 0;
}

bool is_ram_file(const std::string& file)
{
    if (file.size() > 0) {
        if (file[0]=='@') {
            return true;
        }
    }
    return false;
}

std::string ram_file_name(const std::string& file)
{
    if (is_ram_file(file)) {
        return file;
    } else {
        return "@" + file;
    }
}

std::string disk_file_name(const std::string& file)
{
    if (!is_ram_file(file)) {
        return file;
    } else {
        return file.substr(1);
    }
}

int remove(const std::string& file)
{
    if (is_ram_file(file)) {
        return ram_fs::remove(file);
    } else {
        return std::remove(file.c_str());
    }
}

int rename(const std::string& old_filename, const std::string& new_filename)
{
    if (is_ram_file(old_filename)) {
        // TODO: check if new_file is also RAM file name
        return ram_fs::rename(old_filename, new_filename);
    } else {
        return std::rename(old_filename.c_str(), new_filename.c_str());
    }
}

} // end namespace sdsl
