#include <chrono>
#include <algorithm>
#include "sdsl/memory_management.hpp"

using namespace std::chrono;

namespace sdsl
{

void output_event_json(std::ostream& out,const memory_monitor::mm_event& ev,const memory_monitor& m)
{
    out << "\t\t" << "\"name\" : " << "\"" << ev.name << "\",\n";
    out << "\t\t" << "\"usage\" : [" << "\n";
    for (size_t j=0; j<ev.allocations.size(); j++)  {
        out << "\t\t\t[" << duration_cast<milliseconds>(ev.allocations[j].timestamp-m.start_log).count()
            << "," << ev.allocations[j].usage << "]";
        if (j+1<ev.allocations.size()) {
            out << ",\n";
        } else {
            out << "\n";
        }
    }
    out << "\t\t" << "]\n";
}

template<>
void write_mem_log<JSON_FORMAT>(std::ostream& out,const memory_monitor& m)
{
    auto events = m.completed_events;
    std::sort(events.begin(),events.end());

    // output
    out << "[\n";
    for (size_t i=0; i<events.size(); i++) {
        out << "\t{\n";
        output_event_json(out,events[i],m);
        if (i<events.size()-1) {
            out << "\t},\n";
        } else {
            out << "\t}\n";
        }
    }
    out << "]\n";
}

std::string create_mem_html_header()
{
    std::stringstream jsonheader;
    jsonheader
            << "<html>\n"
            << "<head>\n"
            << "<meta charset=\"utf-8\">\n"
            << "<style>\n"
            << "    body { font: 11px sans-serif; }\n"
            << "    .rule { height: 90%; position: absolute; border-right: 1px dotted #000; text-align: right; }\n"
            << "</style>\n"
            << "<title>sdsl memory usage visualization</title>\n"
            << "<script src=\"http://d3js.org/d3.v3.js\"></script>\n"
            << "</head>\n"
            << "<body marginwidth=\"0\" marginheight=\"0\">\n"
            << "<button><a id=\"download\">Save as SVG</a></button>\n"
            << "<div class=\"chart\"><div id=\"visualization\"></div></div><script>\n";
    return jsonheader.str();
}

std::string create_mem_js_body(const std::string& jsonObject)
{
    std::stringstream jsonbody;
    jsonbody
            << "var events = " << jsonObject << ";\n"
            << "var w = window,d = document,e = d.documentElement,g = d.getElementsByTagName('body')[0],\n"
            << "  xw = w.innerWidth || e.clientWidth || g.clientWidth,\n"
            << "  yh = w.innerHeight || e.clientHeight || g.clientHeight;\n\n"
            << "var margin = {top: 20,right: 80,bottom: 120,left: 120},\n"
            << "  width = xw - margin.left - margin.right,height = yh - margin.top - margin.bottom;\n"
            << "var x = d3.scale.linear().range([0, width]);\n"
            << "var y = d3.scale.linear().range([height, 0]);\n"
            << "var xAxis = d3.svg.axis().scale(x).orient(\"bottom\");\n"
            << "var yAxis = d3.svg.axis().scale(y).orient(\"left\").ticks(5);\n"
            << "var color = d3.scale.category10();\n"
            << "var x_max = d3.max(events, function (d) { return d3.max(d.usage, function (u) { return u[0] / 1000;})})\n"
            << "var y_max = d3.max(events, function (d) { return d3.max(d.usage, function (u) { return 1.1 * u[1] / (1024 * 1024);})})\n"
            << "var peak = d3.max(events, function (d) { return d3.max(d.usage, function (u) { return u[1]; })})\n"
            << "var data = []\nevents.forEach(function (d) { data = data.concat(d.usage); });\n"
            << "var peakelem = data.filter(function (a) { return a[1] == peak; });\n"
            << "var peakelem = peakelem.splice(0,1);\n"
            << "x.domain([0, x_max]);\n y.domain([0, y_max]);\n"
            << "var svg = d3.select(\"#visualization\").append(\"svg\")\n"
            << "  .attr(\"width\", width + margin.left + margin.right)\n"
            << "  .attr(\"height\", height + margin.top + margin.bottom)\n"
            << "  .attr(\"xmlns\", \"http://www.w3.org/2000/svg\")\n"
            << "  .append(\"g\").attr(\"transform\",\"translate(\" + margin.left + \",\" + margin.top + \")\");\n\n"
            << "  svg.append(\"g\").attr(\"class\", \"xaxis\").attr(\"transform\", \"translate(0,\" + height + \")\")\n"
            << "  .call(xAxis).append(\"text\").attr(\"text-anchor\", \"end\")\n"
            << "  .attr(\"shape-rendering\", \"crispEdges\").attr(\"x\", width / 2 + 50).attr(\"y\", 70).attr(\"shape-rendering\", \"crispEdges\")\n"
            << "  .attr(\"font-family\", \"sans-serif\").attr(\"font-size\", \"20px\").text(\"Time (seconds)\");\n\n"
            << "svg.append(\"g\").attr(\"class\", \"yaxis\").call(yAxis).append(\"text\").attr(\"transform\", \"rotate(-90)\").attr(\"x\", -height / 2 + 50)\n"
            << "  .attr(\"y\", -80).attr(\"shape-rendering\", \"crispEdges\").attr(\"font-family\", \"sans-serif\").attr(\"font-size\", \"20px\").style(\"text-anchor\", \"end\")\n"
            << "  .text(\"Memory Usage (MiB)\");\n\n"
            << "svg.selectAll(\".tick text\").style(\"font-size\", \"20px\");\n"
            << "svg.selectAll(\".xaxis .tick text\").attr(\"dy\", 23);\nsvg.selectAll(\".yaxis .tick text\").attr(\"dx\", -10);\n"
            << "svg.selectAll(\"line\").attr(\"fill\", \"none\").attr(\"stroke\", \"black\")\nsvg.selectAll(\"path\").attr(\"fill\", \"none\").attr(\"stroke\", \"black\")\n\n"
            << "svg.selectAll(\"line.horizontalGrid\").data(y.ticks(5)).enter().append(\"line\")\n"
            << "  .attr({\"class\": \"horizontalGrid\",\"x1\": 0,\"x2\": width,\"y1\": function (d) { return y(d);},\n"
            << "     \"y2\": function (d) { return y(d); }, \"fill\": \"none\", \"shape-rendering\": \"crispEdges\",\n"
            << "     \"stroke\": \"lightgrey\",\"stroke-dasharray\": \"10,10\",\"stroke-width\": \"1.5px\"});\n\n"
            << "var area = d3.svg.area().x(function (d) { return x(d[0] / 1000);}).y0(height).y1(function (d) { return y(d[1] / (1024 * 1024))});\n\n"
            << "var ev = svg.selectAll(\".event\").data(events).enter().append(\"svg:path\").attr(\"class\", \"area\")\n"
            << "  .attr(\"fill\", function (d) { return d3.rgb(color(d.name)); })\n"
            << "  .attr(\"d\", function (d) { return area(d.usage) })\n"
            << "  .style(\"stroke\", function (d) { return d3.rgb(color(d.name)).darker(2);}).style(\"stroke-width\", \"2px\")\n\n"
            << "svg.selectAll(\".dot\").data(peakelem).enter().append(\"circle\").attr(\"r\", 3).attr(\"fill\", \"red\")\n"
            << "  .attr(\"cx\", function (d) {return x(d[0] / 1000)})\n"
            << "  .attr(\"cy\", function (d) {return y(d[1] / (1024 * 1024))})\n"
            << "  .attr(\"fill\", \"red\").attr(\"stroke-width\", 2).attr(\"stroke\", \"#cc0000\")\n\n"
            << "svg.selectAll(\".dot\").data(peakelem).enter().append(\"svg:text\")\n"
            << "  .attr(\"x\", function (d) {return x(d[0] / 1000)}).attr(\"y\", function (d) {return y(d[1] / (1024 * 1024) * 1.025)})\n"
            << "  .text(function (d) {return \"Peak Usage: \" + Math.round(d[1] / (1024 * 1024)) + \" MB\"})\n"
            << "  .attr(\"font-size\", 12).attr(\"fill\", \"red\");\n\n"
            << "svg.selectAll(\".dot\").data(peakelem).enter().append(\"circle\")\n"
            << "  .attr(\"r\", 5).attr(\"fill\", \"red\")\n"
            << "  .attr(\"cx\", function (d) {return x(d[0] / 1000)})\n"
            << "  .attr(\"cy\", function (d) {return y(d[1] / (1024 * 1024))})\n"
            << "  .attr(\"fill\", \"none\").attr(\"stroke-width\", 2).attr(\"stroke\", \"#cc0000\").each(pulsepeak());\n\n"
            << "function pulsepeak() { return function (d, i, j) {\n"
            << "  d3.select(this).attr(\"r\", 5).style(\"stroke-opacity\", 1.0).transition()\n"
            << "    .ease(\"linear\").duration(1000).attr(\"r\", 10).style(\"stroke-opacity\", 0.0).each(\"end\", pulsepeak());};}\n\n"
            << "var vertical = d3.select(\".chart\").append(\"div\").attr(\"class\", \"remove\")\n"
            << "  .style(\"position\", \"absolute\").style(\"z-index\", \"19\").style(\"width\", \"1px\")\n"
            << "  .style(\"height\", height - margin).style(\"top\", \"30px\").style(\"bottom\", \"50px\")\n"
            << "  .style(\"left\", \"0px\").style(\"opacity\", \"0.4\").style(\"background\", \"black\");\n\n"
            << "var tooltip = d3.select(\".chart\").append(\"div\").attr(\"class\", \"remove\")\n"
            << "  .style(\"position\", \"absolute\").style(\"z-index\", \"20\").style(\"visibility\", \"hidden\").style(\"top\", \"10px\");\n\n"
            << "var circle = svg.append(\"circle\").attr(\"cx\", 100).attr(\"cy\", 350).attr(\"r\", 3).attr(\"fill\", \"black\").style(\"opacity\", \"0\")\n\n"
            << "d3.select(\"svg\").on(\"mousemove\", function () {\n"
            << "  mousex = d3.mouse(this);\n"
            << "  if (mousex[0] < margin.left + 3 || mousex[0] >= xw - margin.right) {\n"
            << "    vertical.style(\"opacity\", \"0\"); tooltip.style(\"opacity\", \"0\"); circle.style(\"opacity\", \"0\")\n"
            << "  } else {\n"
            << "    var xvalue = x.invert(mousex[0] - margin.left); var pos = findPosition(xvalue)\n"
            << "    vertical.style(\"opacity\", \"0.4\"); tooltip.style(\"opacity\", \"1\"); circle.style(\"opacity\", \"1\")\n"
            << "    circle.attr(\"cx\", pos.x).attr(\"cy\", pos.y); vertical.style(\"left\", mousex[0] + \"px\");tooltip.style(\"left\", mousex[0] + 15 + \"px\")\n"
            << "    tooltip.html(\"<p>\" + xvalue.toFixed(2) + \" Seconds <br>\" + Math.round(pos.mem) + \" MiB <br> \" + pos.name + "
            << "  \"<br> Phase Time: \" + pos.ptime + \" Seconds </p>\").style(\"visibility\", \"visible\");\n"
            << "  }\n})"
            << ".on(\"mouseover\", function () {\n"
            << "  mousex = d3.mouse(this);\n  if (mousex[0] < margin.left + 3 || mousex[0] > xw - margin.right) {\n"
            << "    vertical.style(\"opacity\", \"0\")\n  } else {\n    vertical.style(\"opacity\", \"0.4\");vertical.style(\"left\", mousex[0] + 7 + \"px\")\n}})\n"
            << "d3.select(\"#download\").on(\"click\", function () {\n"
            << "d3.select(this).attr(\"href\", 'data:application/octet-stream;base64,' + btoa(d3.select(\"#visualization\").html())).attr(\"download\", \"viz.svg\")})\n\n"
            << "function findPosition(e){correctArea=d3.selectAll(\".area\").filter(function(t){if(t.usage[0][0]<=e*1e3&&t.usage[t.usage.length-1][0]>=e*1e3){return true}"
            << "return false});if(correctArea.empty()){return 0}var t=new Array;correctArea[0].forEach(function(n){t.push(findYValueinArea(n,e))});"
            << "max_elem=d3.max(t,function(e){return e.mem});var n=t.filter(function(e){return e.mem==max_elem});return n[0]}"
            << "function findYValueinArea(e,t){len=e.getTotalLength();var n=0;var r=len;for(var i=0;i<=len;i+=50){var s=e.getPointAtLength(i);"
            << "var o=x.invert(s.x);var u=y.invert(s.y);if(u>0&&o>t){n=Math.max(0,i-50);r=i;break}}var a=e.getPointAtLength(0);"
            << "var f=1;while(n<r){var l=(r+n)/2;a=e.getPointAtLength(l);target_x=x.invert(a.x);if((l==n||l==r)&&Math.abs(target_x-t)>.01){break}if(target_x>t)r=l;"
            << "else if(target_x<t)n=l;else{break}if(f>50){break}f++}var c=new function(){this.mem=y.invert(a.y);this.name=e.__data__.name;"
            << "this.min=d3.min(e.__data__.usage,function(e){return e[0]/1e3});this.max=d3.max(e.__data__.usage,function(e){return e[0]/1e3});"
            << "this.ptime=Math.round(this.max-this.min);this.x=a.x;this.y=a.y};return c}\n</script></body></html>";
    return jsonbody.str();
}


template<>
void write_mem_log<HTML_FORMAT>(std::ostream& out,const memory_monitor& m)
{
    std::stringstream json_data;
    write_mem_log<JSON_FORMAT>(json_data,m);

    out << create_mem_html_header();
    out << create_mem_js_body(json_data.str());
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
    //std::cout << "coalesce_block()" << std::endl;
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
    //std::cout << "split_block("<< (void*)bptr << ")" << std::endl;
    size_t blocksize = UNMASK_SIZE(bptr->size);
    //std::cout << "cur_block_size = " << blocksize << std::endl;
    /* only split if we get at least a small block
       out of it */
    int64_t newblocksize = ALIGNSPLIT(blocksize - ALIGN(size+MM_BLOCK_OVERHEAD));
    //std::cout << "new_block_size = " << newblocksize << std::endl;
    if (newblocksize >= (int64_t)SPLIT_THRESHOLD) {
        /* update blocksize of old block */
        //std::cout << "block_update = " << blocksize-newblocksize << std::endl;
        block_update(bptr,blocksize-newblocksize);
        mm_block_t* newblock = (mm_block_t*)((char*)bptr+(blocksize-newblocksize));
        //std::cout << "new block ptr = " << (void*)newblock << std::endl;
        block_update(newblock,newblocksize);
        coalesce_block(newblock);
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
    //std::cout << "new_block(" << size << ")" << std::endl;
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
    //std::cout << "m_top = " << (void*)m_top << std::endl;
    //std::cout << "m_base = " << (void*)m_base << std::endl;
    if (m_top != m_base) {
        mm_block_foot_t* fptr = (mm_block_foot_t*)(m_top - sizeof(size_t));
        //std::cout << "foot of last = " << (void*)fptr << std::endl;
        //std::cout << "size of last = " << UNMASK_SIZE(fptr->size) << std::endl;
        last = (mm_block_t*)(((uint8_t*)fptr) - UNMASK_SIZE(fptr->size) + sizeof(size_t));
        //std::cout << "last = " << (void*)last << std::endl;
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
    //std::cout << "remove_from_free_set()" << std::endl;
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
    //std::cout << "insert_into_free_set("<< (void*)block << "," << UNMASK_SIZE(block->size) << ")" << std::endl;
    //std::cout << "insert_into_free_set("<< (void*)block << "," << block->size << ")" << std::endl;
    m_free_large.insert({block->size,block});
}

mm_block_t*
hugepage_allocator::find_free_block(size_t size_in_bytes)
{
    //std::cout << "find_free_block(" << size_in_bytes << ")" << std::endl;

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
    //std::cout << "ALLOC(" << size_in_bytes << ")" << std::endl;
    mm_block_t* bptr = nullptr;
    if ((bptr=find_free_block(size_in_bytes + MM_BLOCK_OVERHEAD)) != nullptr) {
        //std::cout << "found free block = " << (void*)bptr << std::endl;
        block_markused(bptr);
        /* split if we have a block too large for us? */
        split_block(bptr,size_in_bytes);
    } else {
        //std::cout << "no free block found that is big enough!" << std::endl;
        // check if last block is free
        //std::cout << "check last block" << std::endl;
        bptr = last_block();
        if (bptr && block_isfree(bptr)) {
            //std::cout << "last block is free. -> extend!" << std::endl;
            // extent last block as it is free
            size_t blockdatasize = block_getdatasize(bptr);
            size_t needed = ALIGN(size_in_bytes - blockdatasize);
            hsbrk(needed);
            remove_from_free_set(bptr);
            block_update(bptr,blockdatasize+needed+sizeof(size_t)+sizeof(mm_block_foot_t));
            //insert_into_free_set(bptr);
            block_markused(bptr);
        } else {
            bptr = new_block(size_in_bytes);
        }
    }
    //print_heap();
    //void* ptr = block_data(bptr);
    //std::cout << "return ptr = " << ptr << std::endl;
    return block_data(bptr);
}

void
hugepage_allocator::mm_free(void* ptr)
{
    //print_heap();
    //std::cout << "FREE(" << ptr << ")" << std::endl;
    if (ptr) {
        mm_block_t* bptr = block_cur(ptr);
        block_markfree(bptr);
        /* coalesce if needed. otherwise just add */
        coalesce_block(bptr);
    }
    //print_heap();
}

void*
hugepage_allocator::mm_realloc(void* ptr, size_t size)
{
    //print_heap();
    //std::cout << "REALLOC(" << ptr << "," << size << ")" << std::endl;
    /* handle special cases first */
    if (nullptr==ptr) return mm_alloc(size);
    if (size==0) {
        mm_free(ptr);
        return nullptr;
    }
    mm_block_t* bptr = block_cur(ptr);

    bool need_malloc = 0;
    size_t blockdatasize = block_getdatasize(bptr);
    /* we do nothing if the size is equal to the block */
    if (size == blockdatasize) {
        //std::cout << "return ptr = " << ptr << std::endl;
        return ptr; /* do nothing if size fits already */
    }
    if (size < blockdatasize) {
        /* we shrink */
        /* do we shrink enough to perform a split? */
        //std::cout << "shrink!" << std::endl;
        split_block(bptr,size);
    } else {
        //std::cout << "expand!" << std::endl;
        /* we expand */
        /* if the next block is free we could use it! */
        mm_block_t* next = block_next(bptr,m_top);
        if (!next) {
            //std::cout << "no next! -> expand!" << std::endl;
            // we are the last block so we just expand
            blockdatasize = block_getdatasize(bptr);
            size_t needed = ALIGN(size - blockdatasize);
            hsbrk(needed);
            block_update(bptr,UNMASK_SIZE(bptr->size)+needed);
            return block_data(bptr);
        } else {
            // we are not the last block
            //std::cout << "try combine next" << std::endl;
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
                //std::cout << "try combine prev" << std::endl;
                mm_block_t* prev = block_prev(bptr,m_first_block);
                if (prev && block_isfree(prev)) {
                    if (blockdatasize + UNMASK_SIZE(prev->size) >= size) {
                        remove_from_free_set(prev);
                        size_t newsize = UNMASK_SIZE(prev->size)+UNMASK_SIZE(bptr->size);
                        block_update(prev,newsize);
                        block_markused(prev);
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
        //std::cout << "need_alloc in REALLOC!" << std::endl;
        void* newptr = mm_alloc(size);
        memcpy(newptr,ptr,size);
        mm_free(ptr);
        ptr = newptr;
    }
    //print_heap();
    //std::cout << "return ptr = " << ptr << std::endl;
    return ptr;
}

uint64_t extract_number(std::string& line)
{
    std::string num_str;
    for (size_t i=line.size()-1; i+1>=1; i--) {
        if (isdigit(line[i])) {
            num_str.insert(num_str.begin(),line[i]);
        } else {
            if (num_str.size() > 0) {
                break;
            }
        }
    }
    return std::strtoull(num_str.c_str(),nullptr,10);
}

uint64_t extract_multiplier(std::string& line)
{
    uint64_t num = 1;
    if (line[line.size()-2] == 'k' || line[line.size()-2] == 'K') {
        num = 1024;
    }
    if (line[line.size()-2] == 'm' || line[line.size()-2] == 'M') {
        num = 1024*1024;
    }
    if (line[line.size()-2] == 'g' || line[line.size()-2] == 'G') {
        num = 1024*1024*1024;
    }
    return num;
}

size_t
hugepage_allocator::determine_available_hugepage_memory()
{
    size_t size_in_bytes = 0;
    size_t page_size_in_bytes = 0;
    size_t num_free_pages = 0;
    const std::string meminfo_file = "/proc/meminfo";
    const std::string ps_str = "Hugepagesize:";
    const std::string pf_str = "HugePages_Free:";
    std::ifstream mifs(meminfo_file);
    if (mifs.is_open()) {
        // find size of one page
        std::string line;
        while (std::getline(mifs, line)) {
            auto ps = std::mismatch(ps_str.begin(),ps_str.end(), line.begin());
            if (ps.first == ps_str.end()) {
                page_size_in_bytes = extract_number(line) * extract_multiplier(line);
            }
            auto pf = std::mismatch(pf_str.begin(),pf_str.end(), line.begin());
            if (pf.first == pf_str.end()) {
                num_free_pages = extract_number(line);
            }
        }
        size_in_bytes = page_size_in_bytes*num_free_pages;
    } else {
        throw std::system_error(ENOMEM,std::system_category(),
                                "hugepage_allocator could not automatically determine available hugepages");
    }
    return size_in_bytes;
}


}
