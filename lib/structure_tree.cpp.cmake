#include "sdsl/structure_tree.hpp"
#include "sdsl/util.hpp"

namespace sdsl
{

void output_tab(std::ostream& out,size_t level)
{
    for (size_t i=0; i<level; i++) out << "\t";
}

template<>
void write_structure_tree<JSON_FORMAT>(const structure_tree_node* v, std::ostream& out, size_t level)
{
    if (v) {
        output_tab(out,level); out << "{" << std::endl;
        output_tab(out,level+1); out << "\"class_name\":" << "\"" << v->type << "\"," << std::endl;
        output_tab(out,level+1); out << "\"name\":" << "\"" << v->name << "\"," << std::endl;
        output_tab(out,level+1); out << "\"size\":" << "\"" << v->size << "\"";

        if (v->children.size()) {
            out << "," << std::endl; // terminate the size tag from before if there are children
            output_tab(out,level+1); out << "\"children\":[" << std::endl;
            size_t written_child_elements = 0;
            for (const auto& child : v->children) {
                if (written_child_elements++ > 0) {
                    out << "," << std::endl;
                }
                write_structure_tree<JSON_FORMAT>(child.second.get(), out,level+2);
            }
            out << std::endl;
            output_tab(out,level+1); out << "]" << std::endl;
        } else {
            out << std::endl;
        }
        output_tab(out,level); out << "}";
    }
}

std::string create_html_header(const char* file_name)
{
    std::stringstream jsonheader;
    jsonheader
            << "<html>\n"
            << "   <head>\n"
            << "    <meta http-equiv=\"Content-Type\" content=\"text/html;charset=utf-8\">\n"
            << "    <title>" << file_name << "</title>\n"
            << "    <script src=\"file://@PROJECT_SOURCE_DIR@/external/d3/d3.min.js\"></script>"
            << "    <script src=\"http://d3js.org/d3.v2.js\"></script>\n"
            << "    <style type=\"text/css\">\n"
            << "      path { stroke: #000; stroke-width: 0.8; cursor: pointer; }\n"
            << "      text { font: 11px sans-serif; cursor: pointer; }\n"
            << "      body { width: 900; margin: 0 auto; }\n"
            << "      h1 { text-align: center; margin: .5em 0; }\n"
            << "      #breadcrumbs { display: none; }\n"
            << "      svg { font: 10px sans-serif; }\n"
            << "     </style>\n"
            << "  </head>\n"
            << "<body marginwidth=\"0\" marginheight=\"0\">\n"
            << "<button><a id=\"download\">Save as SVG</a></button>\n"
            << "  <div id=\"chart\"></div>" << std::endl;
    return jsonheader.str();
}

std::string create_js_body(const std::string& jsonsize)
{
    std::stringstream jsonbody;
    jsonbody
            << "<script type=\"text/javascript\">" << std::endl
            << ""
            "var w = 800,\n"
            "  h = w,\n"
            "  r = w / 2,\n"
            "  x = d3.scale.linear().range([0, 2 * Math.PI]),\n"
            "  y = d3.scale.pow().exponent(1.3).domain([0, 1]).range([0, r]),\n"
            "  p = 5,\n"
            "  color = d3.scale.category20c(),\n"
            "  duration = 1000;\n"
            "\n"
            "var vis = d3.select(\"#chart\").append(\"svg:svg\")\n"
            "  .attr(\"width\", w + p * 2)\n"
            "  .attr(\"height\", h + p * 2)\n"
            "  .append(\"g\")\n"
            "  .attr(\"transform\", \"translate(\" + (r + p) + \",\" + (r + p) + \")\");\n"
            "\n"
            "vis.append(\"p\")\n"
            "  .attr(\"id\", \"intro\")\n"
            "  .text(\"Click to zoom!\");\n"
            "\n"
            "var partition = d3.layout.partition()\n"
            "  .sort(null)\n"
            "  .size([2 * Math.PI, r * r])\n"
            "  .value(function(d) { return d.size; });\n"
            "\n"
            "var arc = d3.svg.arc()\n"
            "  .startAngle(function(d) { return Math.max(0, Math.min(2 * Math.PI, x(d.x))); })\n"
            "  .endAngle(function(d) { return Math.max(0, Math.min(2 * Math.PI, x(d.x + d.dx))); })\n"
            "  .innerRadius(function(d) { return Math.max(0, d.y ? y(d.y) : d.y); })\n"
            "  .outerRadius(function(d) { return Math.max(0, y(d.y + d.dy)); });\n"
            "\n"
            "      " << std::endl
            << "var spaceJSON = " << jsonsize << ";" << std::endl << std::endl
            << "\n"
            "\n"
            "  var nodes = partition.nodes(spaceJSON);\n"
            "\n"
            "  var path = vis.selectAll(\"path\").data(nodes);\n"
            "  path.enter().append(\"path\")\n"
            "    .attr(\"id\", function(d, i) { return \"path-\" + i; })\n"
            "    .attr(\"d\", arc)\n"
            "    .attr(\"fill-rule\", \"evenodd\")\n"
            "    .style(\"fill\", colour)\n"
            "    .on(\"click\", click);\n"
            "\n"
            "  path.append(\"title\").text(function(d) { return 'class name: ' + d.class_name + '\\nmember_name: ' + d.name + '\\n size: ' + sizeMB(d) });\n"
            "\n"
            "  var text = vis.selectAll(\"text\").data(nodes);\n"
            "  var textEnter = text.enter().append(\"text\")\n"
            "    .style(\"opacity\", 1)\n"
            "    .style(\"fill\", function(d) {\n"
            "    return brightness(d3.rgb(colour(d))) < 125 ? \"#eee\" : \"#000\";\n"
            "    })\n"
            "    .attr(\"text-anchor\", function(d) {\n"
            "    return x(d.x + d.dx / 2) > Math.PI ? \"end\" : \"start\";\n"
            "    })\n"
            "    .attr(\"dy\", \".2em\")\n"
            "    .attr(\"transform\", function(d) {\n"
            "    var multiline = (d.name || \"\").split(\" \").length > 1,\n"
            "      angle = x(d.x + d.dx / 2) * 180 / Math.PI - 90,\n"
            "      rotate = angle + (multiline ? -.5 : 0);\n"
            "    return \"rotate(\" + rotate + \")translate(\" + (y(d.y) + p) + \")rotate(\" + (angle > 90 ? -180 : 0) + \")\";\n"
            "    })\n"
            "    .on(\"click\", click);\n"
            "\n"
            "  textEnter.append(\"title\").text(function(d) { return 'class name: ' + d.class_name + '\\nmember_name: ' + d.name + '\\n size: ' + sizeMB(d) });\n"
            "\n"
            "  textEnter.append(\"tspan\")\n"
            "    .attr(\"x\", 0)\n"
            "    .text(function(d) { return d.dx < 0.05 ? \"\" : d.depth ? d.name.split(\" \")[0] : \"\"; });\n"
            "  textEnter.append(\"tspan\")\n"
            "    .attr(\"x\", 0)\n"
            "    .attr(\"dy\", \"1em\")\n"
            "    .text(function(d) { return d.dx < 0.05 ? \"\" : d.depth ? d.name.split(\" \")[1] || \"\" : \"\"; });\n"
            "\n"
            "  function click(d) {\n"
            "  path.transition()\n"
            "    .duration(duration)\n"
            "    .attrTween(\"d\", arcTween(d));\n"
            "\n"
            "  // Somewhat of a hack as we rely on arcTween updating the scales.\n"
            "  text\n"
            "    .style(\"visibility\", function(e) {\n"
            "    return isParentOf(d, e) ? null : d3.select(this).style(\"visibility\");\n"
            "    })\n"
            "    .transition().duration(duration)\n"
            "    .attrTween(\"text-anchor\", function(d) {\n"
            "    return function() {\n"
            "      return x(d.x + d.dx / 2) > Math.PI ? \"end\" : \"start\";\n"
            "    };\n"
            "    })\n"
            "    .attrTween(\"transform\", function(d) {\n"
            "    var multiline = (d.name || \"\").split(\" \").length > 1;\n"
            "    return function() {\n"
            "      var angle = x(d.x + d.dx / 2) * 180 / Math.PI - 90,\n"
            "        rotate = angle + (multiline ? -.5 : 0);\n"
            "      return \"rotate(\" + rotate + \")translate(\" + (y(d.y) + p) + \")rotate(\" + (angle > 90 ? -180 : 0) + \")\";\n"
            "    };\n"
            "    })\n"
            "    .style(\"opacity\", function(e) { return isParentOf(d, e) ? 1 : 1e-6; })\n"
            "    .each(\"end\", function(e) {\n"
            "    d3.select(this).style(\"visibility\", isParentOf(d, e) ? null : \"hidden\");\n"
            "    });\n"
            "  }\n"
            "\n"
            "\n"
            "function sizeMB(d) {\n"
            "//  if (d.children) {\n"
            "//  var sum = calcSum(d);\n"
            "//  return (sum / (1024*1024)).toFixed(2) + 'MB';\n"
            "//  } else {\n"
            "  return (d.size / (1024*1024)).toFixed(2) + 'MB';\n"
            "//  }\n"
            "}\n"
            "\n"
            "function calcSum(d) {\n"
            "  if(d.children) {\n"
            "  var sum = 0;\n"
            "  function recurse(d) {\n"
            "    if(d.children) d.children.forEach( function(child) { recurse(child); } );\n"
            "    else sum += d.size;\n"
            "  }\n"
            "  recurse(d,sum);\n"
            "  console.log(sum);\n"
            "  console.log(d.children);\n"
            "  return sum;\n"
            "  } else {\n"
            "  console.log(d.size);\n"
            "  return d.size;\n"
            "  }\n"
            "}\n"
            "\n"
            "function isParentOf(p, c) {\n"
            "  if (p === c) return true;\n"
            "  if (p.children) {\n"
            "  return p.children.some(function(d) {\n"
            "    return isParentOf(d, c);\n"
            "  });\n"
            "  }\n"
            "  return false;\n"
            "}\n"
            "\n"
            "function colour(d) {\n"
            "  return color(d.name);\n"
            "}\n"
            "\n"
            "// Interpolate the scales!\n"
            "function arcTween(d) {\n"
            "  var my = maxY(d),\n"
            "    xd = d3.interpolate(x.domain(), [d.x, d.x + d.dx]),\n"
            "    yd = d3.interpolate(y.domain(), [d.y, my]),\n"
            "    yr = d3.interpolate(y.range(), [d.y ? 20 : 0, r]);\n"
            "  return function(d) {\n"
            "  return function(t) { x.domain(xd(t)); y.domain(yd(t)).range(yr(t)); return arc(d); };\n"
            "  };\n"
            "}\n"
            "\n"
            "// Interpolate the scales!\n"
            "function arcTween2(d) {\n"
            "  var xd = d3.interpolate(x.domain(), [d.x, d.x + d.dx]),\n"
            "    yd = d3.interpolate(y.domain(), [d.y, 1]),\n"
            "    yr = d3.interpolate(y.range(), [d.y ? 20 : 0, radius]);\n"
            "  return function(d, i) {\n"
            "  return i\n"
            "    ? function(t) { return arc(d); }\n"
            "    : function(t) { x.domain(xd(t)); y.domain(yd(t)).range(yr(t)); return arc(d); };\n"
            "  };\n"
            "}\n"
            "\n"
            "function maxY(d) {\n"
            "  return d.children ? Math.max.apply(Math, d.children.map(maxY)) : d.y + d.dy;\n"
            "}\n"
            "\n"
            "// http://www.w3.org/WAI/ER/WD-AERT/#color-contrast\n"
            "function brightness(rgb) {\n"
            "  return rgb.r * .299 + rgb.g * .587 + rgb.b * .114;\n"
            "}\n"
            "d3.select(\"#download\").on(\"click\", function () {\n"
            "d3.select(this).attr(\"href\", 'data:application/octet-stream;base64,' + btoa(d3.select(\"#chart\").html())).attr(\"download\", \"memorysun.svg\")})\n\n"
            "click(nodes[0]);\n"
            "    "   << std::endl
            << "</script>" << std::endl
            << "</body>" << std::endl
            << "</html>" << std::endl;
    return jsonbody.str();
}

template<>
void write_structure_tree<HTML_FORMAT>(const structure_tree_node* v, std::ostream& out, SDSL_UNUSED size_t level)
{
    std::stringstream json_data;
    write_structure_tree<JSON_FORMAT>(v, json_data);

    out << create_html_header("sdsl data structure visualization");
    out << create_js_body(json_data.str());
}

} // namespace end
