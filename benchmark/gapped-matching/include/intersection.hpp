#pragma once

#include "sdsl/int_vector.hpp"

template<class t_itr,class t_itr2>
intersection_result
intersect(t_itr fbegin,t_itr fend,t_itr2 sbegin,t_itr2 send,int64_t offset = 0)
{
    auto n = std::distance(fbegin,fend);
    auto m = std::distance(sbegin,send);
    intersection_result res(std::min(n,m));

    // std::vector<uint64_t> ss(m);
    // auto tmp = sbegin;
    // size_t k=0;
    // while(tmp != send) {
    //     ss[k++] = *tmp;
    //     ++tmp;
    // }

    size_t i=0;
    if (n < m) {
        size_t value_offset = (size_t) std::min(int64_t(0),offset);

        while (fbegin != fend) {
            int64_t cur = *fbegin;
            if (cur + offset >= 0 && sbegin.skip(cur+offset)) {
                // std::cout << "1 FOUND " << cur+offset << "(" << cur+value_offset << ")" << std::endl;
                res[i++] = cur+value_offset;
            }
            if (sbegin == send) break;
            ++fbegin;
        }
    } else {
        size_t value_offset = (size_t) std::max(int64_t(0),offset);
        while (sbegin != send) {
            int64_t cur = *sbegin;
            if (cur-offset >= 0 && fbegin.skip(cur-offset)) {
                // std::cout << "1 FOUND " << cur-offset << "(" << cur-value_offset << ")" << std::endl;
                res[i++] = cur-value_offset;
            }
            if (fbegin == fend) break;
            ++sbegin;
        }
    }
    res.resize(i);
    return res;
}

template<class t_list1,class t_list2>
intersection_result
intersect(const t_list1& first,const t_list2& second,int64_t offset = 0)
{
    return intersect(first.begin(),first.end(),second.begin(),second.end(),offset);
}


template<class t_list>
intersection_result
intersect(std::vector<t_list> lists)
{
    // sort by size
    std::sort(lists.begin(),lists.end());

    // perform SvS intersection
    auto res = intersect(lists[0],lists[1]);
    for (size_t i=2; i<lists.size(); i++) {
        res = intersect(res,lists[i]);
        if (res.size()==0) break;
    }
    return res;
}

template<class t_list>
intersection_result
pos_intersect(std::vector<t_list> lists)
{
    // sort by size
    std::sort(lists.begin(),lists.end());

    // perform SvS intersection
    auto offset = lists[1].offset() - lists[0].offset();
    auto res = intersect(lists[0],lists[1],offset);
    res.offset = std::min(lists[0].offset(),lists[1].offset());
    for (size_t i=2; i<lists.size(); i++) {
        auto offset = lists[i].offset() - res.offset;
        auto new_offset = std::min(res.offset,lists[i].offset());
        res = intersect(res,lists[i],offset);
        res.offset = new_offset;
        if (res.size()==0) break;
    }
    return res;
}


template<class t_list>
intersection_result
pos_intersect(std::vector<t_list> lists,uint64_t thres)
{
    // sort by size
    std::sort(lists.begin(),lists.end());

    // perform SvS intersection
    auto offset = lists[1].offset() - lists[0].offset();
    auto res = intersect(lists[0],lists[1],offset);
    res.offset = std::min(lists[0].offset(),lists[1].offset());
    for (size_t i=2; i<lists.size(); i++) {
        auto offset = lists[i].offset() - res.offset;
        auto new_offset = std::min(res.offset,lists[i].offset());
        res = intersect(res,lists[i],offset);
        res.offset = new_offset;
        if (res.size()<=thres) break;
    }
    return res;
}
