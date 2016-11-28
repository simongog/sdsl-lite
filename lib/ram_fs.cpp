#include "sdsl/ram_fs.hpp"
#include "sdsl/util.hpp"
#include "sdsl/memory_management.hpp"
#include <cstdio>
#include <iostream>
#include <algorithm>

namespace sdsl
{

ram_fs::ram_fs() {
    m_fd_map[-1] = "";
}

void
ram_fs::store(const std::string& name, content_type data)
{
    auto& r= ram_fs::the_ramfs();
    std::lock_guard<std::recursive_mutex> lock(r.m_rlock);
    if (!exists(name)) {
        std::string cname = name;
        r.m_map.insert(std::make_pair(std::move(cname), std::move(data)));
    } else {
        r.m_map[name] = std::move(data);
    }
}

bool
ram_fs::exists(const std::string& name)
{
    auto& r= ram_fs::the_ramfs();
    std::lock_guard<std::recursive_mutex> lock(r.m_rlock);
    return r.m_map.find(name) != r.m_map.end();
}

ram_fs::content_type&
ram_fs::content(const std::string& name)
{
    auto& r= ram_fs::the_ramfs();
    std::lock_guard<std::recursive_mutex> lock(r.m_rlock);
    return r.m_map[name];
}

size_t
ram_fs::file_size(const std::string& name)
{
    auto& r= ram_fs::the_ramfs();
    std::lock_guard<std::recursive_mutex> lock(r.m_rlock);
    if (exists(name)) {
        return r.m_map[name].size();
    } else {
        return 0;
    }
}

int
ram_fs::remove(const std::string& name)
{
    auto& r= ram_fs::the_ramfs();
    std::lock_guard<std::recursive_mutex> lock(r.m_rlock);
    r.m_map.erase(name);
    return 0;
}

int
ram_fs::rename(const std::string old_filename, const std::string new_filename)
{
    auto& r= ram_fs::the_ramfs();
    std::lock_guard<std::recursive_mutex> lock(r.m_rlock);
    r.m_map[new_filename] = std::move(r.m_map[old_filename]);
    remove(old_filename);
    return 0;
}

ram_fs::content_type&
ram_fs::content(const int fd)
{
    auto& r= ram_fs::the_ramfs();
    std::lock_guard<std::recursive_mutex> lock(r.m_rlock);
    auto name = r.m_fd_map[fd];
    return r.m_map[name];
}

int
ram_fs::truncate(const int fd,size_t new_size)
{
    auto& r= ram_fs::the_ramfs();
    std::lock_guard<std::recursive_mutex> lock(r.m_rlock);
    if(r.m_fd_map.count(fd) == 0) return -1;
    auto name = r.m_fd_map[fd];
    std::cout<<"truncate("<<fd<<","<<name<<","<<r.m_map[name].size()<<")"<<std::endl;
    r.m_map[name].reserve(new_size);
    r.m_map[name].resize(new_size,0);
    std::cout<<"truncate "<<new_size<<", "<<r.m_map[name].size()<<std::endl;
    return 0;
}

size_t
ram_fs::file_size(const int fd)
{
    auto& r= ram_fs::the_ramfs();
    std::lock_guard<std::recursive_mutex> lock(r.m_rlock);
    if(r.m_fd_map.count(fd) == 0) return 0;
    auto name = r.m_fd_map[fd];
    return r.m_map[name].size();
}

int
ram_fs::open(const std::string& name)
{
    auto& r= ram_fs::the_ramfs();
    std::lock_guard<std::recursive_mutex> lock(r.m_rlock);
std::cout<<"name="<<name<<std::endl;
    if(!exists(name)) {
std::cout<<"exists not"<<name<<std::endl;
        store(name,content_type{});
    }
std::cout<<"exists"<<name<<std::endl;
    int fd = -2;
std::cout<<"r.m_fd_map.size()="<<r.m_fd_map.size()<<std::endl;
    auto largest_fd = r.m_fd_map.rbegin()->first;
std::cout<<"largest_fd="<<largest_fd<<std::endl;
    if( largest_fd < 0 ) {
        auto smallest_fd = r.m_fd_map.begin()->first;
        fd = smallest_fd - 1;
    } else {
        r.m_fd_map.erase(largest_fd);
        fd = - largest_fd;
    }
std::cout<<"fd="<<fd<<std::endl;
    r.m_fd_map[fd] = name;
    return fd;
}

int
ram_fs::close(const int fd)
{
    auto& r= ram_fs::the_ramfs();
    std::lock_guard<std::recursive_mutex> lock(r.m_rlock);
    if( fd >= -1 ) return -1;
    if(r.m_fd_map.count(fd) == 0) {
        return -1;
    } else {
        r.m_fd_map.erase(fd);
        r.m_fd_map[-fd] = "";
    }
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

bool is_ram_file(const int fd)
{
    return fd < -1;
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
        if (!is_ram_file(new_filename)) {  // error, if new file is not also RAM-file
            return -1;
        }
        return ram_fs::rename(old_filename, new_filename);
    } else {
        return std::rename(old_filename.c_str(), new_filename.c_str());
    }
}

} // end namespace sdsl
