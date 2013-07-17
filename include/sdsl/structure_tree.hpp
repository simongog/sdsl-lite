/*!\file structure_tree.hpp
   \brief structure_tree.hpp contains a helper class which can represent the memory structure of a class.
   \author Simon Gog
*/
#ifndef INCLUDED_SDSL_STRUCTURE_TREE
#define INCLUDED_SDSL_STRUCTURE_TREE

#include "uintx_t.hpp"
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <memory>

using std::map;

//! Namespace for the succinct data structure library
namespace sdsl
{

class structure_tree_node
{
    private:
        using map_type = std::map<size_t,std::unique_ptr<structure_tree_node>>;
        map_type            m_children;
    public:
        const map_type& children = m_children;
        size_t              size = 0;
        std::string         name;
        std::string         type;
        size_t              hash;
    private:
        static size_t hash_object(const std::string& name,const std::string& type) {
            static std::hash<std::string> hash_func;
            return hash_func(name + type);
        }
    public:
        structure_tree_node(const std::string& n, const std::string& t) : name(n) , type(t) {
            hash = structure_tree_node::hash_object(n,t);
        }
        structure_tree_node* add_child(const std::string& n, const std::string& t) {
            auto hash = structure_tree_node::hash_object(n,t);
            auto child_itr = m_children.find(hash);
            if (child_itr == m_children.end()) {
                // add new child as we don't have one of this type yet
                structure_tree_node* new_node = new structure_tree_node(n,t);
                m_children[hash] = std::unique_ptr<structure_tree_node>(new_node);
                return new_node;
            } else {
                // child of same type and name exists
                return (*child_itr).second.get();
            }
        }
        void add_size(size_t s) { size += s; }
        void append_node(structure_tree_node* other) {
            // make sure they are the same kind of node
            if (this->hash != other->hash) {
                throw std::logic_error("structure_tree_node: trying to append non-equal node.");
            } else {
                add_size(other->size);
                for (const auto& child : m_children) {
                    const auto& child_hash = child.second.get()->hash;
                    auto other_child = other->m_children.find(child_hash);
                    if (other_child == other->m_children.end()) {
                        throw std::logic_error("structure_tree_node: child of equal node not existant.");
                    } else {
                        child.second.get()->append_node((*other_child).second.get());
                    }
                }
            }
        }
};

class structure_tree
{
    public:
        static structure_tree_node* add_child(structure_tree_node* v, const std::string& name, const std::string& type) {
            if (v) return v->add_child(name,type);
            return NULL;
        };
        static void add_size(structure_tree_node* v, uint64_t value) {
            if (v) v->add_size(value);
        };
};


enum format_type {JSON_FORMAT, R_FORMAT};
template<format_type F>
void write_structure_tree(const structure_tree_node* v, std::ostream& out, size_t level = 0);


}
#endif
