#include <chrono>
#include "sdsl/memory_management.hpp"

using namespace std::chrono;

namespace sdsl
{

template<>
void write_mem_log<CSV>(std::ostream& out,const memory_monitor& m)
{
    // write header
    out << "timestamp;memory_usage;event" << std::endl;

    auto first_ts = m.mem_events[0].timestamp;
    auto cur_event = m.events[0];
for (const auto& event : m.mem_events) {
        out << duration_cast<milliseconds>(event.timestamp-first_ts).count() << ";" << event.usage << ";";
    }
}

template<>
void write_mem_log<JSON>(std::ostream& out,const memory_monitor& m)
{
    out << "[";

    auto first_ts = m.mem_events[0].timestamp;
    for (size_t i=0; i<m.mem_events.size()-1; i++) {
        const auto& event = m.mem_events[i];
        out << "[" << duration_cast<milliseconds>(event.timestamp-first_ts).count() << "," << event.usage << "], ";
    }
    const auto& event = m.mem_events[m.mem_events.size()-1];
    out << "[" << duration_cast<milliseconds>(event.timestamp-first_ts).count() << "," << event.usage << "] ";
    out << "]" << std::endl;

    out << "[";

    for (size_t i=0; i<m.events.size()-1; i++) {
        const auto& ev = m.events[i];
        out << "[" << duration_cast<milliseconds>(std::get<0>(ev)-first_ts).count() << ","<< std::get<1>(ev) <<
            ",\"" << std::get<2>(ev) << "\"], ";
    }
    const auto& ev = m.events[m.events.size()-1];
    out << "[" << duration_cast<milliseconds>(std::get<0>(ev)-first_ts).count() << ","<< std::get<1>(ev) <<
        ",\"" << std::get<2>(ev) << "\"]";
    out << "]";
}

#define ALIGNMENT             sizeof(uint64_t)
#define ALIGNSPLIT(size)      (((size)) & ~0x7)
#define ALIGN(size)           (((size) + (ALIGNMENT-1)) & ~0x7)
#define MM_BLOCK_OVERHEAD     (sizeof(size_t)+sizeof(size_t))
#define MIN_BLOCKSIZE         (ALIGN(sizeof(mm_block_t)+sizeof(mm_block_foot_t)))
#define UNMASK_SIZE(size)     ((size)&~1)
#define ISFREE(size)          ((size)&1)
#define SETFREE(size)         ((size)|1)
#define SPLIT_THRESHOLD       (MIN_BLOCKSIZE)


/* from a memory location get the corresponding block header */
using namespace sdsl;

mm_block_t*
block_cur(void* ptr)
{
    mm_block_t* bptr = (mm_block_t*)((uint8_t*)ptr - sizeof(size_t));
    return bptr;
}

/* given a block retrieve the previous block if any. nullptr otherwise */
mm_block_t*
block_prev(mm_block_t* cur_bptr,mm_block_t* first)
{
    /* start of the heap? */
    if (cur_bptr == first) return nullptr;
    mm_block_foot_t* prev_foot = (mm_block_foot_t*)((uint8_t*)cur_bptr - sizeof(mm_block_foot_t));
    mm_block_t* prev_bptr = (mm_block_t*)((uint8_t*)cur_bptr - UNMASK_SIZE(prev_foot->size));
    return prev_bptr;
}

/* given a block retrieve the next block if any. nullptr otherwise */
mm_block_t*
block_next(mm_block_t* cur_bptr,uint8_t* top)
{
    /* end of the heap? */
    if ((uint8_t*)((uint8_t*)cur_bptr+UNMASK_SIZE(cur_bptr->size)) >= top) return nullptr;

    mm_block_t* next_bptr = (mm_block_t*)((uint8_t*)cur_bptr + UNMASK_SIZE(cur_bptr->size));
    return next_bptr;
}

/* calculate the size of a memory block */
size_t
block_size(void* ptr)
{
    mm_block_t* bptr = block_cur(ptr);
    return UNMASK_SIZE(bptr->size);
}

bool
block_isfree(mm_block_t* ptr)
{
    ;
    return ((ptr->size)&1ULL);
}

/* is the next block free */
bool
block_nextfree(mm_block_t* ptr,uint8_t* top)
{
    mm_block_t* next = block_next(ptr,top);
    if (next && block_isfree(next)) return true;
    return false;
}

/* is the prev block free */
bool
block_prevfree(mm_block_t* ptr,mm_block_t* begin)
{
    mm_block_t* prev = block_prev(ptr,begin);
    if (prev && block_isfree(prev)) return 1;
    return 0;
}

/* update the footer with a new size */
void
foot_update(mm_block_t* ptr,size_t size)
{
    mm_block_foot_t* fptr = (mm_block_foot_t*)((uint8_t*)ptr+
                            UNMASK_SIZE(size)-sizeof(mm_block_foot_t));
    fptr->size = size;
}

/* update the block with a new size */
void
block_update(mm_block_t* ptr,size_t size)
{
    ptr->size = size;
    foot_update(ptr,size);
}

/* return the pointer to the "data" */
void*
block_data(mm_block_t* ptr)
{
    return (void*)((uint8_t*)ptr+sizeof(size_t));
}

/* return size of the data that can be stored in the block */
size_t
block_getdatasize(mm_block_t* ptr)
{
    size_t blocksize = UNMASK_SIZE(ptr->size);
    return blocksize - sizeof(size_t) - sizeof(mm_block_foot_t);
}

/* mark the block as free */
void
block_markfree(mm_block_t* ptr)
{
    block_update(ptr,SETFREE(ptr->size));
}

/* mark the block as used */
void
block_markused(mm_block_t* ptr)
{
    block_update(ptr,UNMASK_SIZE(ptr->size));
}


void
hugepage_allocator::coalesce_block(mm_block_t* block)
{
    mm_block_t* newblock = block;
    if (block_nextfree(block,m_top)) {
        mm_block_t* next = block_next(block,m_top);
        /* remove the "next" block from the free list */
        remove_from_free_set(next);
        /* add the size of our block */
        block_update(block,UNMASK_SIZE(block->size)+UNMASK_SIZE(next->size));
    }
    if (block_prevfree(block,m_first_block)) {
        mm_block_t* prev = block_prev(block,m_first_block);
        /* we remove the old prev block and readd it to the correct
           size list if necessary */
        remove_from_free_set(prev);
        newblock = prev;
        block_update(prev,UNMASK_SIZE(prev->size)+UNMASK_SIZE(block->size));
    }
    if (newblock) {
        block_markfree(newblock);
        insert_into_free_set(newblock);
    }
}


void
hugepage_allocator::split_block(mm_block_t* bptr,size_t size)
{
    size_t blocksize = UNMASK_SIZE(bptr->size);
    /* only split if we get at least a small block
       out of it */
    int64_t newblocksize = ALIGNSPLIT(blocksize - ALIGN(size+MM_BLOCK_OVERHEAD));
    if (newblocksize >= (int64_t)SPLIT_THRESHOLD) {
        /* update blocksize of old block */
        block_update(bptr,blocksize-newblocksize);
        mm_block_t* newblock = (mm_block_t*)((char*)bptr+(blocksize-newblocksize));
        block_update(newblock,newblocksize);
        block_markfree(newblock);
        insert_into_free_set(newblock);
    }
}


uint8_t*
hugepage_allocator::hsbrk(size_t size)
{
    ptrdiff_t left = (ptrdiff_t) m_total_size - (m_top - m_base);
    if (left < (ptrdiff_t) size) {  // enough space left?
        throw std::system_error(ENOMEM,std::system_category(),
                                "hugepage_allocator: not enough hugepage memory available");
    }
    uint8_t* new_mem = m_top;
    m_top += size;
    return new_mem;
}

mm_block_t*
hugepage_allocator::new_block(size_t size)
{
    size = ALIGN(size+MM_BLOCK_OVERHEAD);
    if (size < MIN_BLOCKSIZE) size = MIN_BLOCKSIZE;
    mm_block_t* ptr = (mm_block_t*) hsbrk(size);
    block_update(ptr,size);
    return ptr;
}

mm_block_t*
hugepage_allocator::last_block()
{
    mm_block_t* last = nullptr;
    if (m_top != m_base) {
        mm_block_foot_t* fptr = (mm_block_foot_t*)(m_top - sizeof(size_t));
        last = (mm_block_t*)(((uint8_t*)fptr) - UNMASK_SIZE(fptr->size) + sizeof(size_t));
    }
    return last;
}

void
block_print(int id,mm_block_t* bptr)
{
    fprintf(stdout, "%d addr=%p size=%lu (%lu) free=%d\n",id,((void*)bptr),
            UNMASK_SIZE(bptr->size),bptr->size,block_isfree(bptr));
    fflush(stdout);
}

void
hugepage_allocator::print_heap()
{
    mm_block_t* bptr = m_first_block;
    size_t id = 0;
    while (bptr) {
        block_print(id,bptr);
        id++;
        bptr = block_next(bptr,m_top);
    }
}

void
hugepage_allocator::remove_from_free_set(mm_block_t* block)
{
    auto eq_range = m_free_large.equal_range(block->size);
    // find the block amoung the blocks with equal size
    auto itr = eq_range.first;
    auto last = eq_range.second;
    auto found = m_free_large.end();
    while (itr != last) {
        if (itr->second == block) {
            found = itr;
        }
        ++itr;
    }
    if (found == m_free_large.end()) {
        found = last;
    }
    m_free_large.erase(found);
}

void
hugepage_allocator::insert_into_free_set(mm_block_t* block)
{
    m_free_large.emplace(block->size,block);
}

mm_block_t*
hugepage_allocator::find_free_block(size_t size_in_bytes)
{
    mm_block_t* bptr = nullptr;
    auto free_block = m_free_large.lower_bound(size_in_bytes);
    if (free_block != m_free_large.end()) {
        bptr = free_block->second;
        m_free_large.erase(free_block);
    }
    return bptr;
}

void*
hugepage_allocator::mm_alloc(size_t size_in_bytes)
{
    mm_block_t* bptr = nullptr;
    if ((bptr=find_free_block(size_in_bytes + MM_BLOCK_OVERHEAD)) != nullptr) {
        block_markused(bptr);
        /* split if we have a block too large for us? */
        split_block(bptr,size_in_bytes);
    } else {
        // check if last block is free
        bptr = last_block();
        if (bptr && block_isfree(bptr)) {
            // extent last block as it is free
            size_t blockdatasize = block_getdatasize(bptr);
            size_t needed = ALIGN(size_in_bytes - blockdatasize);
            hsbrk(needed);
            remove_from_free_set(bptr);
            block_update(bptr,UNMASK_SIZE(bptr->size)+needed);
            insert_into_free_set(bptr);
        } else {
            bptr = new_block(size_in_bytes);
        }
    }
    return block_data(bptr);
}

void
hugepage_allocator::mm_free(void* ptr)
{
    if (ptr) {
        mm_block_t* bptr = block_cur(ptr);
        block_markfree(bptr);
        /* coalesce if needed. otherwise just add */
        coalesce_block(bptr);
    }
}

void*
hugepage_allocator::mm_realloc(void* ptr, size_t size)
{
    /* handle special cases first */
    if (ptr==NULL) return mm_alloc(size);
    if (size==0) {
        mm_free(ptr);
        return NULL;
    }
    mm_block_t* bptr = block_cur(ptr);

    bool need_malloc = 0;
    size_t blockdatasize = block_getdatasize(bptr);
    /* we do nothing if the size is equal to the block */
    if (size == blockdatasize)
        return ptr; /* do nothing if size fits already */
    if (size < blockdatasize) {
        /* we shrink */
        /* do we shrink enough to perform a split? */
        split_block(bptr,size);
    } else {
        /* we expand */
        /* if the next block is free we could use it! */
        mm_block_t* next = block_next(bptr,m_top);
        if (!next) {
            // we are the last block so we just expand
            blockdatasize = block_getdatasize(bptr);
            size_t needed = ALIGN(size - blockdatasize);
            hsbrk(needed);
            block_update(bptr,UNMASK_SIZE(bptr->size)+needed);
            return block_data(bptr);
        } else {
            // we are not the last block
            if (next && block_isfree(next)) {
                /* do we have enough space if we use the next block */
                if (blockdatasize + UNMASK_SIZE(next->size) >= size) {
                    /* the next block is enough! */
                    /* remove the "next" block from the free list */
                    remove_from_free_set(next);
                    /* add the size of our block */
                    block_update(bptr,UNMASK_SIZE(bptr->size)+UNMASK_SIZE(next->size));
                } else {
                    /* the next block is not enough. we allocate a new one instead */
                    need_malloc = true;
                }
            } else {
                /* try combing the previous block if free */
                mm_block_t* prev = block_prev(bptr,m_first_block);
                if (prev && block_isfree(prev)) {
                    if (blockdatasize + UNMASK_SIZE(prev->size) >= size) {
                        remove_from_free_set(prev);
                        size_t newsize = UNMASK_SIZE(prev->size)+UNMASK_SIZE(bptr->size);
                        block_update(prev,newsize);
                        /* move the data into the previous block */
                        ptr = memmove(block_data(prev),ptr,blockdatasize);
                    } else {
                        /* not enough in the prev block */
                        need_malloc = true;
                    }
                } else {
                    /* prev block not free. get more memory */
                    need_malloc = true;
                }
            }
        }
    }
    if (need_malloc) {
        void* newptr = mm_alloc(size);
        memcpy(newptr,ptr,size);
        mm_free(ptr);
        ptr = newptr;
    }
    return ptr;
}


}
