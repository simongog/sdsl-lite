#include "sdsl/wt.hpp"

namespace sdsl
{

unsigned char& unsigned_char_map::operator[](unsigned char i)
{
    return *(m_map+i);
}

unsigned char unsigned_char_map::operator[](unsigned char i)const
{
    return m_map[i];
}

void unsigned_char_map::clear()
{
    for (uint16_t i=0; i<256; ++i)
        m_map[i] = 0;
}

uint16_t unsigned_char_map::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    uint16_t written_bytes = 0;
    for (uint16_t i=0; i<256; ++i) {
        out.write((char*)&m_map[i], sizeof(m_map[i]));
        written_bytes += sizeof(m_map[256]);
    }
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void unsigned_char_map::load(std::istream& in)
{
    for (uint16_t i=0; i<256; ++i) {
        in.read((char*) &m_map[i], sizeof(m_map[i]));
    }
}

void unsigned_char_map::swap(unsigned_char_map& map)
{
    if (this != &map) {
        for (uint16_t i=0; i<256; ++i) {
            std::swap(m_map[i], map.m_map[i]);
        }
    }
}

} // end namespace
