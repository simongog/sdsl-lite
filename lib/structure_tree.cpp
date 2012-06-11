#include "sdsl/structure_tree.hpp"
#include <algorithm> // for std::swap
#include <string>
#include <iostream>

namespace sdsl
{

structure_tree_node::structure_tree_node():m_parent(NULL), parent(m_parent), children(m_children), key_values(m_key_values) {}

// User defined operator.
structure_tree_node::structure_tree_node(structure_tree_node* v):m_parent(v), parent(m_parent), children(m_children), key_values(m_key_values) {}

// User defined operator.
structure_tree_node::structure_tree_node(structure_tree_node* v, const string& name, const string& class_name):m_parent(v), parent(m_parent), children(m_children), key_values(m_key_values)
{
    add_key_value("name", name);
    add_key_value("class_name", class_name);
}

// copy constructor
structure_tree_node::structure_tree_node(const structure_tree_node& v):parent(m_parent), children(m_children), key_values(m_key_values)
{
    if (this != &v) { // if it is not the same object
        m_parent     = v.m_parent;
        m_children   = v.m_children;
        m_key_values = v.m_key_values;
    }
}

structure_tree_node::~structure_tree_node()
{
    // recursively delete children
    for (size_t i=0; i < m_children.size(); ++i) {
        delete m_children[i];
    }
}

// swap operator
void structure_tree_node::swap(structure_tree_node& v)
{
    std::swap(m_parent, v.m_parent);
    m_children.swap(v.m_children);
    m_key_values.swap(v.m_key_values);
}

void structure_tree_node::add_key_value(const string& key, const string& value)
{
    m_key_values[key] = value;
}

structure_tree_node* structure_tree::add_child(structure_tree_node* v, const string& name, const string& class_name)
{
    if (NULL == v) {
        return NULL;
    } else {
        structure_tree_node* v_new = new structure_tree_node(v, name, class_name);
        v->m_children.push_back(v_new);
        return v_new;
    }
}

structure_tree_node* structure_tree::parent(const structure_tree_node* v)
{
    if (NULL == v) {
        return NULL;
    } else {
        return v->parent;
    }
}


template<>
void write_structure_tree<JSON_FORMAT>(const structure_tree_node* v, std::ostream& out)
{
    if (NULL == v or (v->children.size()==0 and v->key_values.size()==0)) {
        return;
    }
    typedef structure_tree_node::tKeyValue::const_iterator const_iterator;
    size_t written_elements = 0;
    out << "{"; // begin json element
    for (const_iterator it = v->key_values.begin(); it != v->key_values.end(); ++it) {
        if (written_elements++ > 0) {
            out << ",";
        }
        out << "\"" << it->first << "\":\"" << it->second << "\"";
    }
    if (v->children.size() > 0) {
        if (written_elements++ > 0) {
            out << ",";
        }
        out << "\"children\":["; // open children
        size_t written_child_elements = 0;
        for (size_t i = 0; i < v->children.size(); ++i) {
            if (written_child_elements++ > 0) {
                out << ",";
            }
            write_structure_tree<JSON_FORMAT>(v->children[i], out);
        }
        out << "]";// close children
    }
    out << "}"; // end json element
}

template<>
void write_structure_tree<R_FORMAT>(const structure_tree_node* v, std::ostream& out)
{
    if (NULL == v or (v->children.size()==0 and v->key_values.size()==0)) {
        return;
    }
    typedef structure_tree_node::tKeyValue::const_iterator const_iterator;
    size_t written_elements = 0;
    out << "list("; // begin R list
    for (const_iterator it = v->key_values.begin(); it != v->key_values.end(); ++it) {
        if (written_elements++ > 0) {
            out << ",";
        }
        out << it->first << " = " << it->second;
    }
    if (v->children.size() > 0) {
        if (written_elements++ > 0) {
            out << ",";
        }
        out << "list("; // open children
        size_t written_child_elements = 0;
        for (size_t i = 0; i < v->children.size(); ++i) {
            if (written_child_elements++ > 0) {
                out << ",";
            }
            write_structure_tree<JSON_FORMAT>(v->children[i], out);
        }
        out << ")";// close children
    }
    out << ")"; // end R list
}


} // end namespace
