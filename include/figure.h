#ifndef _FIGURE_
#define _FIGURE_

#include <cairomm/cairomm.h>
#include <iostream>

class Figure
{
public:
  Figure(std::string name = "out", std::array<double, 4> param = {0, 0, 400, 400})
  {
    double x = param[0], y = param[1], dx = param[2], dy = param[3];
    double scal = (dx>600?dx:600) / dx;
    dx = dx * scal; dy = dy * scal;

    surface_ = Cairo::SvgSurface::create(name+".svg", dx, dy);
    cr_ = Cairo::Context::create(surface_);
    cr_->set_line_join(Cairo::LINE_JOIN_BEVEL);
    cr_->scale(scal, -scal);
    cr_->translate(0, -param[3]);
    cr_->translate(-x, -y);
  }

  template<typename Mesh>
  void draw_mesh(Mesh & m, bool showindex = false);

  template<typename Mesh>
  void draw_halfedge(Mesh & m, bool showindex = false);

  template<typename Mesh>
  void draw_node(const Mesh & m, bool showindex = false);

  template<typename Mesh>
  void draw_edge(const Mesh & m, bool showindex = false);

private:
  void _show_label(double x, double y, std::string label, double ax,  double ay, double size, double a);

private:
  Cairo::RefPtr<Cairo::SvgSurface> surface_;
  Cairo::RefPtr<Cairo::Context> cr_;
};

void Figure::_show_label(double x, double y, std::string label, double ax,  double ay, double size, double a)
{
  cr_->set_source_rgba(0, 0, 1, a);  /**< RGBA */

  cr_->arc(ax, ay, size, 0, 2*M_PI); /**< 绘制点 */
  cr_->close_path();
  cr_->fill_preserve();
  cr_->set_line_width(0.0);
  cr_->stroke();

  cr_->set_source_rgb(0, 0, 0); /**< 设置文本的颜色 */

  cr_->select_font_face("Sans", Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL);
  cr_->set_font_matrix(Cairo::Matrix(size*4, 0, 0, -size*4, 0, 0));

  cr_->move_to(x, y); /**< 绘制编号 */
  cr_->show_text(label);
}

template<typename Mesh>
void Figure::draw_mesh(Mesh & m, bool showindex)
{
  uint32_t n = 0;
  uint32_t NC = m.number_of_cells();
  auto & cell = *(m.get_cell());
  std::vector<double> ch(NC);
  for(auto & c : cell)
  {
    auto * start = c.halfedge();
    cr_->set_source_rgb(128.0/256.0, 230.0/256.0, 115.0/256.0);
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
    for(auto & c : cell)
    {
      auto p = c.barycentary();
      _show_label(p.x, p.y, std::to_string(n), p.x, p.y, 0.02*ch[n], 0.5);
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
    auto normal = h.normal()*0.8;
    auto tangen = h.tangential()*0.8;
    auto p0 = h.previous()->node()->coordinate()+normal*0.02 + tangen*0.1;
    auto p1 = p0 + tangen;
    if(h.edge()->halfedge()==&h)
      cr_->set_source_rgb(1, 0, 0);
    else
      cr_->set_source_rgb(0, 0, 0);

    double l = std::sqrt(std::min(h.cell()->area(), h.opposite()->cell()->area()));
    cr_->set_line_width(0.01*l); /**< 获取多边形的尺寸并调整 linewidth 为尺寸的 0.008 倍*/

    cr_->move_to(p0.x, p0.y);
    cr_->line_to(p1.x, p1.y);
    cr_->stroke();  

    normal = normal.normalize()*l;
    tangen = tangen.normalize()*l;

    auto p2 = p1 - tangen*0.05 - normal*0.005; 
    auto p4 = p1 + tangen*0.03 - normal*0.005;
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
    for(auto & h : halfedge)
    {
      auto p = h.barycentary();
      auto normal = h.normal();
      auto tangen = h.tangential();
      double size = 0.02*std::sqrt(std::min(h.cell()->area(), h.opposite()->cell()->area()));
      p += normal*0.016;

      double s = std::log(n)/std::log(10)+1;
      double x = p.x, y = p.y;
      normal = normal.normalize();
      if(tangen.y > 0)
      {
        std::cout << n << " " << s << std::endl;
        x += normal.x*size*s*2;
        y += normal.y*size*s*2;
      }
      _show_label(x, y, std::to_string(n), p.x, p.y, size, 0.5);
      n++;
    }
  }
}

#endif /* _FIGURE_ */ 
