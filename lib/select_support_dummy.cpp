#include "sdsl/select_support_dummy.hpp"

namespace sdsl
{

select_support_dummy::select_support_dummy(const bit_vector* v):select_support(v) { }


select_support_dummy::select_support_dummy(const select_support_dummy& ss):select_support(ss.m_v) { }


select_support_dummy& select_support_dummy::operator=(const select_support_dummy& ss) {
    return *this;
}

void select_support_dummy::swap(select_support_dummy& ss) { }

select_support_dummy::~select_support_dummy() { }

const select_support_dummy::size_type select_support_dummy::select(size_type i)const {
    throw std::logic_error("select_dummy does not implement select");
    return 0;
}


const select_support_dummy::size_type select_support_dummy::operator()(size_type i)const {
    return select(i);
}

void select_support_dummy::set_vector(const bit_vector* v) { }

select_support_dummy::size_type select_support_dummy::serialize(std::ostream& out, structure_tree_node* v, std::string name)const {
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    structure_tree::add_size(child, 0);
    return 0;
}

void select_support_dummy::load(std::istream& in, const bit_vector* v) { }

}
