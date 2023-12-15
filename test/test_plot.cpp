#include <cairomm/cairomm.h>
#include <cmath>

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
    double length = std::hypot(arrow.end.y - arrow.start.y, arrow.end.x - arrow.start.x);
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
    double length = std::hypot(arrow.end.y - arrow.start.y, arrow.end.x - arrow.start.x);
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

void save_svg(const char *filename, const Polygon& polygon) {
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

