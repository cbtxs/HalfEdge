#include <cairomm/cairomm.h>
#include <cmath>

class Figure
{
public:
  Figure(std::string name)
  {
    surface_ = Cairo::SvgSurface::create(name+".svg", 400, 400);
    cr_ = Cairo::Context::create(surface_);
  }

  template<typename Mesh>
  void draw_mesh(const Mesh & m, bool showindex = false);

  template<typename Mesh>
  void draw_halfedge(const Mesh & m, bool showindex = false);

  template<typename Mesh>
  void draw_node(const Mesh & m, bool showindex = false);

  template<typename Mesh>
  void draw_edge(const Mesh & m, bool showindex = false);

private:
  Cairo::RefPtr<Cairo::SvgSurface> surface_;
  Cairo::RefPtr<Cairo::Context> cr_;
};

template<typename Mesh>
void Figure::draw_mesh(const Mesh & m, bool showindex)
{
  auto & node = *(m.get_node());
  auto & edge = *(m.get_edge());
  auto & cell = *(m.get_cell());
  auto & halfedge = *(m.get_halfedge());

  for(auto & c : cell)
  {
    auto * start = c.halfedge();
    cr_->set_source_rgb(0, 0, 1);
    cr_->move_to(start->node()->coordinate().x*10, start->node()->coordinate().y*10);

    for(auto * h = start->next(); h != start; h=h->next())
      cr_->line_to(h->node()->coordinate().x*10, h->node()->coordinate().y*10);

    cr_->close_path();
    cr_->fill_preserve();

    cr_->set_source_rgb(1, 0, 0);
    cr_->stroke();  // Stroke to draw the edges
  }
  if(showindex)
  {
    uint32_t n = 0;
    for(auto & c : cell)
    {
      auto p = c.barycentery()*10;
      // 设置点的颜色和不透明度
      cr_->set_source_rgba(0, 0, 1, 0.5);  // RGBA，这里的透明度设置为0.5
      cr_->set_line_width(5.0);

      // 绘制点
      cr_->arc(p.x, p.y, 5.0, 0, 2 * M_PI);
      cr_->fill_preserve();  // 填充并保留路径
      cr_->stroke();

      // 设置文本的颜色
      cr_->set_source_rgb(0, 0, 0);

      // 设置文本的字体和大小
      cr_->select_font_face("Sans", Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL);
      cr_->set_font_size(12);

      // 绘制编号
      char text[10];
      snprintf(text, sizeof(text), "%d", n++);
      cr_->move_to(p.x + 10, p.y + 5);
      cr_->show_text(text);
    }
  }
}

class Point {
public:
    double x, y;

    Point(double x, double y) : x(x), y(y) {}
};

class Arrow {
public:
    Point start, end;

    Arrow(const Point& start, const Point& end) : start(start), end(end) {}
};

class Polygon {
public:
    std::vector<Point> points;
    std::vector<Arrow> arrows;

    void addPoint(double x, double y) {
        points.emplace_back(x, y);
    }

    void calculateArrows() {
        arrows.clear();
        for (std::size_t i = 0; i < points.size(); ++i) 
        {
          arrows.emplace_back(points[i], points[(i + 1) % points.size()]);
          arrows.emplace_back(points[(i + 1) % points.size()], points[i]);
        }
    }
};
void draw_half_arrow(const Cairo::RefPtr<Cairo::Context>& cr, const Arrow& arrow) {
    double angle = std::atan2(arrow.end.y - arrow.start.y, arrow.end.x - arrow.start.x);
    //double length = std::hypot(arrow.end.y - arrow.start.y, arrow.end.x - arrow.start.x);
    double arrow_length = 10.0;  // Adjust as needed
    double half_arrow_length = arrow_length / 2.0;

    // Calculate the midpoint of the arrow
    double mid_x = (arrow.start.x + arrow.end.x) / 2.0;
    double mid_y = (arrow.start.y + arrow.end.y) / 2.0;

    // Draw the arrow line
    cr->move_to(arrow.start.x, arrow.start.y);
    cr->line_to(mid_x, mid_y);
    cr->stroke();

    // Draw the half-arrowhead
    cr->move_to(mid_x - half_arrow_length * std::cos(angle - M_PI / 6),
                mid_y - half_arrow_length * std::sin(angle - M_PI / 6));
    cr->line_to(mid_x, mid_y);
    cr->line_to(mid_x - half_arrow_length * std::cos(angle + M_PI / 6),
                mid_y - half_arrow_length * std::sin(angle + M_PI / 6));
    cr->stroke();
}


void draw_arrow(const Cairo::RefPtr<Cairo::Context>& cr, const Arrow& arrow) {
    double angle = std::atan2(arrow.end.y - arrow.start.y, arrow.end.x - arrow.start.x);
    //double length = std::hypot(arrow.end.y - arrow.start.y, arrow.end.x - arrow.start.x);
    double arrow_length = 10.0;  // Adjust as needed

    // Draw the arrow line
    cr->move_to(arrow.start.x, arrow.start.y);
    cr->line_to(arrow.end.x, arrow.end.y);
    cr->stroke();

    // Draw the arrowhead
    cr->move_to(arrow.end.x, arrow.end.y);
    cr->line_to(arrow.end.x - arrow_length * std::cos(angle - M_PI / 6),
                arrow.end.y - arrow_length * std::sin(angle - M_PI / 6));
    cr->line_to(arrow.end.x - arrow_length * std::cos(angle + M_PI / 6),
                arrow.end.y - arrow_length * std::sin(angle + M_PI / 6));
    cr->close_path();
    cr->fill();
}

void draw_polygon(const Cairo::RefPtr<Cairo::Context>& cr, const Polygon& polygon) {
    // Set blue color for the fill
    cr->set_source_rgb(0, 0, 1);
    cr->move_to(polygon.points[0].x, polygon.points[0].y);

    for (const auto& point : polygon.points) {
        cr->line_to(point.x, point.y);
    }

    cr->close_path();
    cr->fill_preserve();

    // Set red color for the edges
    cr->set_source_rgb(1, 0, 0);
    cr->stroke();  // Stroke to draw the edges

    // Draw two arrows on each edge
    cr->set_source_rgb(0, 0, 0);  // Set green color for the arrows

    for (const auto& arrow : polygon.arrows) {
        draw_arrow(cr, arrow);
    }
}

void save_svg(const char *filename, const Polygon& polygon) 
{
  auto surface = Cairo::SvgSurface::create(filename, 400, 400);
  auto cr = Cairo::Context::create(surface);

  draw_polygon(cr, polygon);
}

int main() {
    Polygon polygon;
    polygon.addPoint(200, 50);
    polygon.addPoint(350, 150);
    polygon.addPoint(350, 300);
    polygon.addPoint(200, 400);
    polygon.addPoint(50, 300);

    polygon.calculateArrows();

    save_svg("output.svg", polygon);

    return 0;
}

