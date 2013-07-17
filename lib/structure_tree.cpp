#include "sdsl/structure_tree.hpp"
#include "sdsl/util.hpp"
#include <algorithm> // for std::swap

namespace sdsl
{

structure_tree_node::structure_tree_node():m_parent(NULL), parent(m_parent), children(m_children), key_values(m_key_values) {}

// User defined operator.
structure_tree_node::structure_tree_node(structure_tree_node* v):m_parent(v), parent(m_parent), children(m_children), key_values(m_key_values) {}

// User defined operator.
structure_tree_node::structure_tree_node(structure_tree_node* v, const std::string& name, const std::string& class_name):m_parent(v), parent(m_parent), children(m_children), key_values(m_key_values)
{
    add_key_value(STRUCTURE_TREE_NAME_KEY, name);
    add_key_value(STRUCTURE_TREE_CLASS_NAME_KEY, class_name);
}

// copy method
void structure_tree_node::copy(const structure_tree_node& v)
{
    if (this != &v) { // if it is not the same object
        m_parent     = v.m_parent;
        m_children   = v.m_children;
        m_key_values = v.m_key_values;
        for (vector<structure_tree_node*>::size_type i=0; i < m_children.size(); ++i) { // make a deep copy
            m_children[i] = new structure_tree_node(*(v.m_children[i]));
        }
    }
}

// copy constructor
structure_tree_node::structure_tree_node(const structure_tree_node& v):parent(m_parent), children(m_children), key_values(m_key_values)
{
    copy(v);
}

structure_tree_node& structure_tree_node::operator=(const structure_tree_node& v)
{
    copy(v);
    return *this;
}

structure_tree_node::~structure_tree_node()
{
    delete_children();
}

void structure_tree_node::delete_children()
{
    // recursively delete children
    for (size_t i=0; i < m_children.size(); ++i) {
        delete m_children[i];
    }
    m_children.clear();
}

// swap operator
void structure_tree_node::swap(structure_tree_node& v)
{
    std::swap(m_parent, v.m_parent);
    m_children.swap(v.m_children);
    m_key_values.swap(v.m_key_values);
}

void structure_tree_node::add_key_value(const std::string& key, const std::string& value)
{
    m_key_values[key] = value;
}

structure_tree_node* structure_tree::add_child(structure_tree_node* v, const std::string& name, const std::string& class_name)
{
    if (nullptr == v) {
        return nullptr;
    } else {
        structure_tree_node* v_new = new structure_tree_node(v, name, class_name);
        v->m_children.push_back(v_new);
        return v_new;
    }
}

structure_tree_node* structure_tree::parent(const structure_tree_node* v)
{
    if (nullptr == v) {
        return nullptr;
    } else {
        return v->parent;
    }
}

void structure_tree_node::add_size(uint64_t value)
{
    m_key_values[STRUCTURE_TREE_SIZE_KEY] = std::to_string(value);
}

void structure_tree::add_size(structure_tree_node* v, uint64_t value)
{
    if (nullptr != v) {
        v->add_size(value);
    }
}

bool structure_tree_node::merge(const structure_tree_node& v)
{
    if (equal_structure(v)) {
        rec_merge(v);
        return true;
    }
    return false;
}

bool structure_tree_node::equal_key_value_pairs(const key_val_t& kv1, const key_val_t& kv2)const
{
    if (kv1.size() != kv2.size())
        return false;
    for (auto it1 = kv1.begin(), it2 = kv2.begin(); it1 != kv1.end(); ++it1, ++it2) {
        if (it1->first != it2->first) {
            return false;
        }
        if (it1->first != "size") {
            if (it1->second != it2->second)
                return false;
        }
    }
    return true;
}

bool structure_tree_node::equal_structure(const structure_tree_node& v)const
{
    if (equal_key_value_pairs(m_key_values, v.m_key_values)) {
        if (m_children.size() != v.m_children.size())
            return false;
        // NOTE: The children are not sorted.
        for (size_t i=0; i < m_children.size(); ++i) {
            if (!m_children[i]->equal_structure(*(v.m_children[i])))
                return false;
        }
        return true;
    } else {
        std::cout<<"structure not the same: key_values_pairs differ"<<std::endl;
        std::cout<<"m_key_values.size()="<<m_key_values.size()<<std::endl;
        std::cout<<"v.m_key_values.size()="<<v.m_key_values.size()<<std::endl;
        for (auto x : v.m_key_values) {
            std::cout << "(" << x.first << "," << x.second << ")" << std::endl;
        }
        std::cout<<"======="<<std::endl;
        for (auto x : v.m_key_values) {
            std::cout << "(" << x.first << "," << x.second << ")" << std::endl;
        }
        std::cout<<"-------"<<std::endl;
    }
    return false;
}

void structure_tree_node::rec_merge(const structure_tree_node& v)
{
    for (size_t i=0; i < m_children.size(); ++i) {
        m_children[i]->rec_merge(*(v.m_children[i]));
    }
    if (m_key_values.find(STRUCTURE_TREE_SIZE_KEY) != m_key_values.end()) {
        uint64_t size, vsize;
        std::stringstream(m_key_values[STRUCTURE_TREE_SIZE_KEY]) >> size;
        std::stringstream(v.m_key_values.find(STRUCTURE_TREE_SIZE_KEY)->second) >> vsize;
        m_key_values[STRUCTURE_TREE_SIZE_KEY] = std::to_string(size + vsize);
    }
}

bool structure_tree::merge_children(structure_tree_node* v)
{
    if (nullptr != v and v->m_children.size() > 0) {
        // copy the first child
        structure_tree_node* merged_child = new structure_tree_node(*(v->m_children[0]));
        bool success = true;
        for (vector<structure_tree_node*>::size_type i = 1; i < v->m_children.size() and success; ++i) {
            success = merged_child->merge(*(v->m_children[i]));
        }
        if (success) {
            uint64_t m_1 = v->m_children.size()-1;
            v->delete_children();
            merged_child->m_parent = v; // set parent of merged child to v
            v->m_children.push_back(merged_child);
            v->key_values[STRUCTURE_TREE_NAME_KEY] = "[0.." + std::to_string(m_1) + "]";
        }
    }
    return true;
}


template<>
void write_structure_tree<JSON_FORMAT>(const structure_tree_node* v, std::ostream& out)
{
    if (NULL == v or (v->children.size()==0 and v->key_values.size()==0)) {
        return;
    }
    size_t written_elements = 0;
    out << "{"; // begin json element
    for (auto x : v->key_values) {
        if (written_elements++ > 0) {
            out << ",";
        }
        out << "\"" << x.first << "\":\"" << x.second << "\"";
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
    size_t written_elements = 0;
    out << "list("; // begin R list
    for (auto it = v->key_values.begin(); it != v->key_values.end(); ++it) {
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
