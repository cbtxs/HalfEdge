#ifndef GEOMETRY_UTILS2D_H 
#define GEOMETRY_UTILS2D_H

#include <utility>
#include <vector>
#include <iostream>
#include <algorithm>

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
  bool dist_point_to_line(const Point2d& p0, const Point2d& p1, const Point2d& p) const
  {
    double dist = (p1 - p0).cross(p - p0) / (p1 - p0).length();
    return std::abs(dist);
  }

  /**
   * @brief Determines if a point is on a line within a tolerance, if the
   * distance between the point and the line is less than the tolerance.
   * @param p0 The first point on the line
   * @param p1 The second point on the line
   * @param p The point to check
   */
  bool point_on_line(const Point2d& p0, const Point2d & p1, const Point2d& p) const
  {
    double dist = dist_point_to_line(p0, p1, p); 
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
    p = p0 + (p1 - p0).normalize() * (p - p0).dot(p1 - p0);
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
   * @return 0 if the point is outside the polygon,
   *         1 if the point is in the polygon,
   *         2 if the point is on a polygon edge,
   *         3 if the point is a polygon vertex
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
        return 3;
      }
    }
    for (int i = 0; i < n; i++)
    {
      Point2d & p0 = *polygon[i];
      Point2d & p1 = *polygon[(i + 1) % n];
      if (point_on_line(p0, p1, p))
      { 
        index = i;
        return 2;
      }
    }
    if (point_in_polygon(polygon, p))
    {
      return 1;
    }
    return 0;
  }


  /**
   * @brief Computes the intersection point of two line segments
   * @param p0 The first point of the first line segment
   * @param p1 The second point of the first line segment
   * @param q0 The first point of the second line segment
   * @param q1 The second point of the second line segment
   * @param p The intersection point
   */
  bool intersection_point_of_two_segments(const Point2d& p0, 
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
   *         1 if the two segments are intersecting at a point, 
   *         2 if the two segments are not intersecting
   */
  uint8_t relative_position_of_two_segments(const Point2d& p0, 
                                                  const Point2d& p1, 
                                                  const Point2d& q0, 
                                                  const Point2d& q1,
                                                        Point2d& p) const
  {
    bool flag0 = point_on_line(q0, q1, p0);
    bool flag1 = point_on_line(q0, q1, p1);
    if (flag0 && flag1)
    {
      return 0;
    }
    auto v0 = p1-p0;
    auto v1 = q0-q1;
    auto v2 = q0-p0;
    double den = v0.cross(v1);
    if (den == 0)
    {
      return 2;
    }
    double t = (v2.cross(v1))/den;
    double s = (v2.cross(v0))/den;
    if (t >= 0 && t <= 1 && s >= 0 && s <= 1)
    {
      p = p0 + v0*t;
      return 1;
    }
    return 2;
  }

  /**
   * @brief intersection points of a line segment and a polygon
   * @param polygon The vertices of the polygon
   * @param p0 The first point of the line segment
   * @param p1 The second point of the line segment
   * @param points The intersection points
   * @param edges The edges of the polygon that intersect with the line segment
   * @return The number of intersection points
   */
  uint8_t intersection_points_of_segments_and_polygon(const std::vector<Point2d*>& polygon, 
                                             const Point2d& p0, 
                                             const Point2d& p1, 
                                             std::vector<std::pair<Point2d, uint32_t>>& points) const
  {
    uint32_t n = polygon.size(); 
    for (uint32_t i = 0; i < n; i++)
    {
      Point2d & q0 = *polygon[i];
      Point2d & q1 = *polygon[(i + 1) % n];
      Point2d p;
      uint8_t flag = relative_position_of_two_segments(p0, p1, q0, q1, p);
      if (flag == 1)
      {
        points.emplace_back(p, i);
      }
      else if (flag == 0)
      {
        points.emplace_back(p0, i);
        points.emplace_back(p1, (i + 1) % n);
      }
    }
    /** 根据交点在 p0, p1 之间的位置排序 */
    std::sort(points.begin(), points.end(), 
      [&](const std::pair<Point2d, uint32_t>& a, const std::pair<Point2d, uint32_t>& b)
      {
        return (a.first - p0).length() < (b.first - p0).length();
      }
    );
    return points.size();
  }





private:
  double tol_;

};



}

#endif // GEOMETRY_UTILS2D_H
