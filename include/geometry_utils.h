#ifndef GEOMETRY_UTILS2D_H 
#define GEOMETRY_UTILS2D_H

#include <vector>
#include <stdint.h>

#include "geometry.h"

namespace HEM
{

/**
 * @brief The GeometryUtils class
 * This class contains utility functions for geometry operations
 */
class GeometryUtils2D
{
public:
  /**
   * @brief Constructor 
   */
  GeometryUtils2D(double tol = 1e-6): tol_(tol) {}

  /**
   * @brief Determines if two points are equal within a tolerance, if the
   * distance between the points is less than the tolerance.
   * @param p1 The first point
   * @param p2 The second point
   */
  bool points_equal(const Point2d& p1, const Point2d& p2) const
  {
    double dist = (p1 - p2).length();
    return dist < tol_;
  }

  /**
   * @brief Determines if two points are equal within a tolerance, if the
   * distance between the points is less than the tolerance.
   * @param p1 The first point
   * @param p2 The second point
   */
  double dist_point_to_segment(const Point2d & p0, const Point2d & p1, const Point2d & p) const
  {
    auto v = p1 - p0;
    auto w = p  - p0;

    double c1 = w.dot(v);
    if ( c1 <= 0 )
      return w.length();

    double c2 = v.dot(v);
    if ( c2 <= c1 )
      return (p-p1).length();

    double b = c1 / c2;
    Point2d pb = p0 + v*b;
    return (p-pb).length();
  }

  /**
   * @brief Determines if a point is on a line within a tolerance, if the
   * distance between the point and the line is less than the tolerance.
   * @param p0 The first point on the line
   * @param p1 The second point on the line
   * @param p The point to check
   */
  bool point_on_segment(const Point2d& p0, const Point2d & p1, const Point2d& p) const
  {
    double dist = dist_point_to_segment(p0, p1, p); 
    return std::abs(dist) < tol_;
  }

  /**
   * @brief Project a point to a line
   * @param p0 The first point on the line
   * @param p1 The second point on the line
   * @param p The point to project
   */
  void project_point_to_line(const Point2d& p0, const Point2d & p1, Point2d& p) const
  {
    auto v = p1 - p0;
    p = p0 + v*(v.dot(p-p0)/v.dot(v));
  }

  /**
   * @brief Determines if a point is on a polygon
   * @param polygon The vertices of the polygon
   * @param p The point to check
   */
  bool point_in_polygon(const std::vector<Point2d*>& polygon, const Point2d& p) const
  {
    int n = polygon.size();
    int count = 0;
    for (int i = 0; i < n; i++)
    {
      Point2d & p0 = *polygon[i];
      Point2d & p1 = *polygon[(i + 1) % n];

      if((p.x>p0.x) != (p.x>p1.x))
      {
        double ft = (p0.x-p.x)/(p0.x-p1.x);
        if(p.y < p0.y-ft*(p0.y-p1.y))
          count += 1; 
      }
    }
    return count % 2 == 1;
  }
  
  /**
   * @brief relative of a point and a polygon
   * @param polygon The vertices of the polygon
   * @param p The point to check
   * @param index The index of the vertex or edge
   * @return 0 if the point is the index-th polygon vertex,
   *         1 if the point is on the index-th polygon edge,
   *         2 if the point is inside the polygon
   *         3 if the point is outside the polygon
   */
  uint8_t relative_position_of_point_and_polygon(
      const std::vector<Point2d*>& polygon, 
      const Point2d& p,
      uint32_t & index) const
  {
    index = 0;
    int n = polygon.size();
    for (int i = 0; i < n; i++)
    {
      Point2d & p0 = *polygon[i];
      if (points_equal(p0, p))
      {
        index = i;
        return 0;
      }
    }
    for (int i = 0; i < n; i++)
    {
      Point2d & p0 = *polygon[i];
      Point2d & p1 = *polygon[(i + 1) % n];
      if (point_on_segment(p0, p1, p))
      { 
        index = i;
        return 1;
      }
    }
    if (point_in_polygon(polygon, p))
    {
      return 2;
    }
    return 3;
  }


  /**
   * @brief Computes the intersection point of two line segments
   * @param p0 The first point of the first line segment
   * @param p1 The second point of the first line segment
   * @param q0 The first point of the second line segment
   * @param q1 The second point of the second line segment
   * @param p The intersection point
   */
  bool intersection_of_two_segments(const Point2d& p0, 
                                          const Point2d& p1, 
                                          const Point2d& q0, 
                                          const Point2d& q1,
                                                Point2d& p) const
  {
    auto v0 = p1-p0;
    auto v1 = q0-q1;
    auto v2 = q0-p0;
    double den = v0.cross(v1);
    if (den == 0)
    {
      return false;
    }
    double t = (v2.cross(v1))/den;
    double s = (v2.cross(v0))/den;
    if (t >= 0 && t <= 1 && s >= 0 && s <= 1)
    {
      p = p0 + v0*t;
      return true;
    }
    return false;
  }

  /**
   * @brief releative position of two line segments
   * @param p0 The first point of the first line segment
   * @param p1 The second point of the first line segment
   * @param q0 The first point of the second line segment
   * @param q1 The second point of the second line segment
   * @param p The intersection point
   * @return 0 if the two segments are overlapping, 
   *         1 if the two segments intersect at the vertex q0, 
   *         2 if the two segments intersect at the vertex q1,
   *         3 if the two segments intersect at a point,
   *         4 if the two segments are not intersecting
   */
  uint8_t relative_position_of_two_segments(const Point2d& p0, 
                                                  const Point2d& p1, 
                                                  const Point2d& q0, 
                                                  const Point2d& q1,
                                                        Point2d& p) const
  {
    bool flag0 = point_on_segment(p0, p1, q0);
    bool flag1 = point_on_segment(p0, p1, q1);
    if (flag0 && flag1)
    {
      p = q0;
      return 0; /**< overlapping */
    }
    else if(flag0)
    {
      p = q0;
      return 1; /**< intersecting at the vertex q0 */
    }
    else if(flag1)
    {
      p = q1;
      return 2; /**< intersecting at the vertex q1 */
    }
    bool flag2 = point_on_segment(q0, q1, p0);
    bool flag3 = point_on_segment(q0, q1, p1);
    if(flag2 || flag3)
    {
      p = flag2 ? p0 : p1;
      return 3; /**< intersecting at a point */
    }

    auto v0 = p1-p0;
    auto v1 = q0-q1;
    auto v2 = q0-p0;
    double den = v0.cross(v1);
    if (den == 0)
    {
      return 4; /**< parallel, not intersecting */
    }
    double t = (v2.cross(v1))/den;
    double s = (v0.cross(v2))/den;
    if (t >= 0 && t <= 1 && s >= 0 && s <= 1)
    {
      p = p0 + v0*t;
      return 3; /**< intersecting at a point */
    }
    return 4; /**< not intersecting */
  }

  uint8_t relative_position_of_two_segments(const Point2d& p0, 
                                                  const Point2d& p1, 
                                                  const Point2d& q0, 
                                                  const Point2d& q1) const
  {
    Point2d p;
    return relative_position_of_two_segments(p0, p1, q0, q1, p);
  }

  /**
   * @brief Determines the quadrant of a vector
   * @param v The vector
   * @return The quadrant of the vector
   */
  uint8_t quadrant_of_vector(const Vector2d & v)
  {
    uint8_t a = v.x<0;
    uint8_t b = v.y<0;
    return a+3*b-2*a*b;
  }

private:
  double tol_;

};



}

#endif // GEOMETRY_UTILS2D_H
