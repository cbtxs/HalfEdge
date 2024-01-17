#ifndef _FIGURE_
#define _FIGURE_

#include <algorithm>
#include <cairomm/cairomm.h>
#include <iostream>
#include <functional>

class Figure
{
public:
  struct Color
  {
    double r,g,b;
  }; 
public:
  Figure(std::string name = "out", std::array<double, 4> param = {0, 0, 400, 400})
  {
    double x = param[0], y = param[1], dx = param[2], dy = param[3];
    scal_ = (dx>1920?dx:1920) / dx;
    dx = dx * scal_; dy = dy * scal_;

    surface_ = Cairo::SvgSurface::create(name+".svg", dx, dy);
    cr_ = Cairo::Context::create(surface_);
    cr_->set_line_join(Cairo::LINE_JOIN_BEVEL);
    cr_->scale(scal_, -scal_);
    cr_->translate(0, -param[3]);
    cr_->translate(-x, -y);
  }

  Color get_color(double v,double vmin,double vmax)
  {
    Color c = {1.0,1.0,1.0}; // white
    double dv;

    if (v < vmin)
      v = vmin;
    if (v > vmax)
      v = vmax;
    dv = vmax - vmin;

    if (v<(vmin + 0.25*dv)) 
    {
      c.r = 0;
      c.g = 4*(v - vmin)/dv;
    } 
    else if (v<(vmin + 0.5*dv)) 
    {
      c.r = 0;
      c.b = 1 + 4*(vmin + 0.25*dv - v)/dv;
    } 
    else if (v<(vmin + 0.75*dv)) 
    {
      c.r = 4*(v - vmin - 0.5 * dv)/dv;
      c.b = 0;
    } 
    else 
    {
      c.g = 1 + 4*(vmin + 0.75*dv - v)/dv;
      c.b = 0;
    }
    return c;
  }

  template<typename Mesh>
  void draw_mesh(Mesh & m, bool showindex = false)
  {
    draw_mesh<Mesh, double>(m, showindex);
  }

  template<typename Mesh, typename Data>
  void draw_mesh(Mesh & m, bool showindex = false, std::string name="");

  template<typename Mesh>
    void draw_halfedge(Mesh & m, bool showindex = false);

  template<typename Mesh>
  void draw_node(Mesh & m, bool showindex = false);

  template<typename Mesh>
  void draw_edge(Mesh & m, bool showindex = false);

private:
  void _show_label(double x, double y, std::string label, double ax,  double ay, 
      double size, std::array<double, 4> pcolor, std::array<double, 4> textcolor);

private:
  Cairo::RefPtr<Cairo::SvgSurface> surface_;
  Cairo::RefPtr<Cairo::Context> cr_;
  double scal_;
};

void Figure::_show_label(double x, double y, std::string label, double ax,  double ay, double size, 
    std::array<double, 4> pcolor, std::array<double, 4> textcolor)
{
  cr_->set_source_rgba(pcolor[0], pcolor[1], pcolor[2], pcolor[3]);  /**< RGBA */

  cr_->arc(ax, ay, size, 0, 2*M_PI); /**< 绘制点 */
  cr_->close_path();
  cr_->fill_preserve();
  cr_->set_line_width(0.0);
  cr_->stroke();

  cr_->set_source_rgba(textcolor[0], textcolor[1], textcolor[2], textcolor[3]);  /**< RGBA */

  cr_->select_font_face("Sans", Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL);
  cr_->set_font_matrix(Cairo::Matrix(size*4, 0, 0, -size*4, 0, 0));

  cr_->move_to(x, y); /**< 绘制编号 */
  cr_->show_text(label);
}

template<typename Mesh, typename Data>
void Figure::draw_mesh(Mesh & m, bool showindex, std::string dname)
{
  uint32_t n = 0;
  uint32_t NC = m.number_of_cells();
  auto & cell = *(m.get_cell());
  std::vector<double> ch(NC);
  std::function<Color(typename Mesh::Cell &)> color_fun;
  if(dname.size()<1)
    color_fun = [](typename Mesh::Cell &)->Color { return Color{128.0/256.0, 230.0/256.0, 115.0/256.0};};
  else
  {
    auto & data = *(m.template get_cell_data<Data>(dname));
    double max = -1e10, min = 1e10;
    for(auto d : data)
    {
      max = std::max(max, (double)d);
      min = std::min(min, (double)d);
    }
    color_fun = [max, min, &data, this](typename Mesh::Cell & c)->Color
    {
      return this->get_color((double)data[c.index()], min, max);
    };
  }

  for(auto & c : cell)
  {
    auto * start = c.halfedge();
    Color color = color_fun(c);
    cr_->set_source_rgb(color.r, color.g, color.b);
    cr_->move_to(start->node()->coordinate().x, start->node()->coordinate().y);

    for(auto * h = start->next(); h != start; h=h->next())
      cr_->line_to(h->node()->coordinate().x, h->node()->coordinate().y);

    ch[n] = std::sqrt(c.area()); 
    cr_->set_line_width(0.008*ch[n++]); /**< 获取多边形的尺寸并调整 linewidth 为尺寸的 0.008 倍*/

    cr_->close_path();
    cr_->fill_preserve();

    cr_->set_source_rgb(0, 0, 0);
    cr_->stroke();  // Stroke to draw the edges
  }
  if(showindex)
  {
    n = 0;
    std::array<double, 4> pcolor = {0, 0, 1, 0.5};
    std::array<double, 4> textcolor = {0, 0, 0, 1};
    for(auto & c : cell)
    {
      auto p = c.barycenter();
      _show_label(p.x, p.y, std::to_string(n), p.x, p.y, 0.024*ch[n], pcolor, textcolor);
      n++;
    }
  }
}

template<typename Mesh>
void Figure::draw_halfedge(Mesh & m, bool showindex)
{
  uint32_t n = 0;
  auto & halfedge = *(m.get_halfedge());
  for(auto & h : halfedge)
  {
    double l = std::sqrt(std::min(h.cell()->area(), h.opposite()->cell()->area()));
    double l0 = std::sqrt(std::max(h.cell()->area(), h.opposite()->cell()->area()));

    auto normal = h.normal()*0.8;
    auto tangen = h.tangential()*0.8;
    auto p0 = h.previous()->node()->coordinate()+normal.normalize()*(l*0.025+l0*0.005) + tangen*(1.0/8.0);
    auto p1 = p0 + tangen;
    if(h.edge()->halfedge()==&h)
      cr_->set_source_rgb(1, 0, 0);
    else
      cr_->set_source_rgb(0, 0, 0);

    cr_->set_line_width(0.01*l); /**< 获取多边形的尺寸并调整 linewidth 为尺寸的 0.008 倍*/

    /** [p0, p1] 半边的杆 */
    cr_->move_to(p0.x, p0.y);
    cr_->line_to(p1.x, p1.y);
    cr_->stroke();  

    normal = normal.normalize()*l;
    tangen = tangen.normalize()*l;

    /** [p4, p2, p3]  半边的头 */
    auto p2 = p1 - tangen*0.05 - normal*0.005; 
    auto p4 = p1 + tangen*0.024 - normal*0.005;
    cr_->move_to(p4.x, p4.y);
    cr_->line_to(p2.x, p2.y);

    auto p3 = p2 + normal*0.03;

    cr_->line_to(p3.x, p3.y);
    cr_->close_path();
    cr_->fill_preserve();

    cr_->set_line_width(0.0);
    cr_->stroke();  // Stroke to draw the edges
  }
  if(showindex)
  {
    std::array<double, 4> pcolor = {0, 0, 1, 0.5};
    std::array<double, 4> textcolor = {0, 0, 0, 1};
    for(auto & h : halfedge)
    {
      double l = std::sqrt(std::min(h.cell()->area(), h.opposite()->cell()->area()));
      double l0 = std::sqrt(std::max(h.cell()->area(), h.opposite()->cell()->area()));

      auto p = h.barycenter();
      auto normal = h.normal();
      auto tangen = h.tangential();
      double size = 0.018*l;
      p += normal.normalize()*(l*0.025+l0*0.005);

      double x = p.x, y = p.y;
      normal = normal.normalize();

      /** 下面四行是获得字符的大小，确定字符放什么位置 */
      cr_->set_font_size(size*4);
      cr_->select_font_face("Sans", Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL);
      Cairo::TextExtents extents;
      cr_->get_text_extents(std::to_string(n), extents);
      if(tangen.x > 0)
      {
        if(tangen.y > 0)
          x -= extents.width*1.1;
      }
      else
      {
        y -= extents.height*1.1;
        if(tangen.y > 0)
          x -= extents.width*1.1;
      }
      _show_label(x, y, std::to_string(n), p.x, p.y, size, pcolor, textcolor);
      n++;
    }
  }
}

template<typename Mesh>
void Figure::draw_node(Mesh & m, bool showindex)
{
  auto & nodes = *(m.get_node());
  if(showindex)
  {
    uint32_t n = 0;
    std::array<double, 4> pcolor = {1, 0, 0, 1};
    std::array<double, 4> textcolor = {0, 0, 0, 1};
    for(auto & node : nodes)
    {
      auto p = node.coordinate();
      double size = 1e100;
      uint32_t N = node.get_top();
      for(uint32_t i = 0; i < N; i++)
        size = std::min(size, node.node2cell[i]->area());
      size = 0.018*std::sqrt(size);
      _show_label(p.x+size*0.5, p.y+size*0.5, std::to_string(n), p.x, p.y, size, pcolor, textcolor);
      n++;
    }
  }
}

#endif /* _FIGURE_ */ 
