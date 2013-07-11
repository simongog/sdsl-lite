#include "sdsl/sfstream.hpp"
#include "sdsl/ram_fs.hpp"
#include <iostream>

#define DEBUG_STREAM 0

namespace sdsl
{

//  IMPLEMENTATION OF OSFSTREAM

osfstream::osfstream() : std::ostream(m_streambuf), m_streambuf(NULL), m_file(""), m_use_ram(true), m_closed(true)
{
    this->init(m_streambuf);
}

osfstream::osfstream(const std::string& file, std::ios_base::openmode mode) : std::ostream(m_streambuf),
    m_streambuf(NULL), m_file(""), m_use_ram(true), m_closed(true)
{
    this->init(m_streambuf);
    open(file, mode);
}

std::streambuf*
osfstream::open(const std::string& file, std::ios_base::openmode mode)
{
    if (DEBUG_STREAM)std::cerr<<"OSFSTREAM: try to open " << file <<std::endl;
    if (NULL != m_streambuf) {
        delete m_streambuf;
        m_streambuf = NULL;
    }
    m_file = file;
    m_use_ram = is_ram_file(file);
    if (m_use_ram) {
        m_streambuf = new std::stringbuf(ram_fs::content(m_file), mode | std::ios_base::out);
        if (DEBUG_STREAM)std::cerr<<"m_stringbuf="<<(size_t)m_streambuf<<std::endl;
    } else {
        std::filebuf* f_buf = new std::filebuf();
        if (!f_buf->open(m_file.c_str(), mode | std::fstream::out)) {
            this->setstate(std::ios_base::failbit);
            if (NULL != f_buf) {
                delete f_buf;
            }
            f_buf = NULL;
            if (DEBUG_STREAM)std::cerr<<"ERROR OSFSTREAM: could not open `"<<m_file<<"`"<<std::endl;
        } else {
            this->clear();
        }
        m_streambuf = f_buf;
    }
    this->rdbuf(m_streambuf);
    m_closed = false;
    return m_streambuf;
}

bool
osfstream::is_open()
{
    if (NULL == m_streambuf or m_closed) {
        return false;
    } else {
        if (m_use_ram) {
            return true;
        } else {
            return ((std::filebuf*)m_streambuf)->is_open();
        }
    }
}

void
osfstream::close()
{
    if (NULL == m_streambuf) {
        this->setstate(std::ios::failbit);
        return;
    }
    if (m_use_ram and !m_closed) {
        ram_fs::store(m_file, ((std::stringbuf*)m_streambuf)->str());
        m_closed = true;
    } else {
        if (((std::filebuf*)m_streambuf)->close()) {
            this->clear();
            m_closed = true;
        } else {
            if (DEBUG_STREAM)std::cerr<<"STREAM: set failbit"<<std::endl;
            this->setstate(std::ios::failbit);
        }
    }
}

osfstream::~osfstream()
{
    if (is_open()) {
        close();
    }
    if (m_streambuf != NULL)
        delete m_streambuf;

}

osfstream::operator voidptr()const
{
    // TODO: && is->valid() ...
    return m_streambuf;
}

//  IMPLEMENTATION OF ISFSTREAM

isfstream::isfstream() : std::istream(m_streambuf), m_streambuf(NULL), m_file(""), m_use_ram(true), m_closed(true)
{
    this->init(m_streambuf);
}

isfstream::isfstream(const std::string& file, std::ios_base::openmode mode) : std::istream(m_streambuf),
    m_streambuf(NULL), m_file(""), m_use_ram(true), m_closed(true)
{
    this->init(m_streambuf);
    open(file, mode);
}

std::streambuf*
isfstream::open(const std::string& file, std::ios_base::openmode mode)
{
    if (DEBUG_STREAM)std::cout<<"ISFSTREAM: try to open "<<file<<std::endl;
    if (NULL != m_streambuf) {
        delete m_streambuf;
        m_streambuf = NULL;
    }
    m_file = file;
    m_use_ram = is_ram_file(file);
    if (m_use_ram) {
        if (ram_fs::exists(m_file)) {
            m_streambuf = new std::stringbuf(ram_fs::content(m_file), mode | std::ios_base::in);
            if (NULL == m_streambuf) {
                this->setstate(std::ios_base::failbit);
                m_closed = true;
                if (DEBUG_STREAM)std::cout<<"ERROR: opening failed "<<m_file<<std::endl;
            } else {
                this->clear();
                this->rdbuf(m_streambuf);
                m_closed = false;
                if (DEBUG_STREAM)std::cout<<"opened "<<m_file<<std::endl;
            }
        } else {
            m_closed = true;
            this->setstate(std::ios_base::failbit);
        }
    } else {
        std::filebuf* f_buf = new std::filebuf();
        if (!f_buf->open(file.c_str(), mode | std::ios_base::in)) {
            this->setstate(std::ios_base::failbit);
            if (NULL != f_buf) {
                delete f_buf;
            }
        } else {
            this->clear();
            if (NULL != m_streambuf) {
                delete m_streambuf;
                m_streambuf = NULL;
            }
            m_streambuf = f_buf;
            this->rdbuf(m_streambuf);
            m_closed = false;
        }
    }
    return m_streambuf;
}

bool
isfstream::is_open()
{
    if (NULL == m_streambuf or m_closed) {
        return false;
    } else {
        if (m_use_ram) {
            return true;
        } else {
            return ((std::filebuf*)m_streambuf)->is_open();
        }
    }
}

void
isfstream::close()
{
    if (NULL == m_streambuf) {
        this->setstate(std::ios::failbit);
        return;
    }
    if (!m_use_ram) {
        if (((std::filebuf*)m_streambuf)->close()) {
            this->clear();
            m_closed = true;
        } else {
            this->setstate(std::ios::failbit);
        }
    }
}

isfstream::~isfstream()
{
    if (m_streambuf != NULL)
        delete m_streambuf;
}

isfstream::operator voidptr()const
{
    // TODO: && is->valid() ...
    return m_streambuf;
}

}// end namespace sdsl
