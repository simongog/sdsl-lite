#include "sdsl/ram_filebuf.hpp"
#include <iostream>

namespace sdsl
{


ram_filebuf::~ram_filebuf() {}

ram_filebuf::ram_filebuf() {}

ram_filebuf::ram_filebuf(std::vector<char>& ram_file) : m_ram_file(&ram_file)
{
    char* begin = m_ram_file->data();
    char* end   = begin + m_ram_file->size();
    setg(begin, begin, end); // set get pointers eback(), eptr(), egptr()
}

std::streambuf*
ram_filebuf::open(const std::string name, std::ios_base::openmode mode)
{
    // open ram_file
    if ((mode & std::ios_base::in) and !(mode & std::ios_base::trunc)) {
        // file must exist, initial position at the start
        if (!ram_fs::exists(name)) {
            m_ram_file = nullptr;
        } else {
            m_ram_file = &ram_fs::content(name);
        }
    } else { // existence of file not required
        if (!ram_fs::exists(name)) {
            // create empty file, if it does not yet exist
            ram_fs::store(name, ram_fs::content_type());// TODO: create method in ram_fs?? or store w 1 arg?
        }
        m_ram_file = &ram_fs::content(name);
        if ((mode & std::ios_base::out) and !(mode & std::ios_base::app)) {
            m_ram_file->clear();
        }
    }

    if (m_ram_file and (mode & std::ios_base::trunc)) {
        m_ram_file->clear();
    }
    if (m_ram_file) {
        if (mode & std::ios_base::ate) {
            // TODO: move put pointer to the end of the file
        } else {

        }
        setg(m_ram_file->data(), m_ram_file->data(), m_ram_file->data()+m_ram_file->size());
        setp(m_ram_file->data(), m_ram_file->data()+m_ram_file->size());
    }
// ATTENTION: if m_ram_file->size() == 0, then data might be nullptr !!!
    return m_ram_file ? this : nullptr;
}

bool
ram_filebuf::is_open()
{
    return m_ram_file!=nullptr;
}

ram_filebuf*
ram_filebuf::close()
{
    if (!this->is_open())
        return nullptr;
    return this;
}

ram_filebuf::pos_type
ram_filebuf::seekpos(pos_type sp, std::ios_base::openmode mode)
{
    if (sp >= 0 and sp < m_ram_file->size()) {
        setg(eback(), eback()+sp, egptr());
        setp(pbase(), epptr());
        pbump(pbase()+sp-pptr()); // pptr should be pbase() anyway after the setp call?
    } else {
        if (mode & std::ios_base::out) {
            // extend buffer
            m_ram_file->resize(sp, 0);
            setg(m_ram_file->data(), m_ram_file->data()+sp, m_ram_file->data()+m_ram_file->size());
            setp(m_ram_file->data(), m_ram_file->data()+m_ram_file->size());
            pbump(sp);
        } else {
            return pos_type(off_type(-1));
        }
    }
    return 0;
}

ram_filebuf::pos_type
ram_filebuf::pubseekoff(off_type off, std::ios_base::seekdir way,
                        std::ios_base::openmode which)
{
    if (std::ios_base::beg == way) {
        if (seekpos(off, which) == pos_type(-1)) {
            return pos_type(-1);
        }
    } else if (std::ios_base::cur == way) {
        if (seekpos(gptr()-eback()+off, which) == pos_type(-1)) {
            return pos_type(-1);
        }
    } else if (std::ios_base::end == way) {
        if (seekpos(egptr()-eback()+off, which) == pos_type(-1)) {
            return pos_type(-1);
        }
    }
    return gptr()-eback();
}


ram_filebuf::pos_type
ram_filebuf::pubseekpos(pos_type sp, std::ios_base::openmode which)
{
    if (seekpos(sp, which) == pos_type(-1)) {
        return pos_type(-1);
    } else {
        return gptr()-eback();
    }
}

/*
std::streamsize
ram_filebuf::xsputn(const char_type* s, std::streamsize n){
    std::cerr<<"call xsputn("<<s<<","<<n<<")"<<std::endl;
    if ( m_ram_file ){
        if (  n < epptr()-pptr() ){
            memcpy(pptr(), s, n);
            pbump(n);
            return n;
        } else {
            m_ram_file->resize( (pptr()+n)-epptr() );
            memcpy(pptr(), s, n);
            setp(m_ram_file->data(), m_ram_file->data()+m_ram_file->size());
            pbump(epptr()-pbase());
            return n;
        }
    } else {
        return 0;
    }
}
*/

int
ram_filebuf::sync()
{
    return 0; // we are always in sync, since buffer is sink
}

ram_filebuf::int_type
ram_filebuf::overflow(int_type c)
{
    if (m_ram_file) {
        m_ram_file->push_back(c);
        setp(m_ram_file->data(), m_ram_file->data()+m_ram_file->size());
        pbump(epptr()-pbase());
    }
    return traits_type::to_int_type(c);
}


}

