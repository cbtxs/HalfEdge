#ifndef _GEOMETRY_
#define _GEOMETRY_

#include <cmath>

namespace HEM 
{

class Vector2d {
public:
    // 构造函数
    Vector2d(double x_ = 0.0, double y_ = 0.0) : x(x_), y(y_) {}

    // 复制构造函数
    Vector2d(const Vector2d& other) : x(other.x), y(other.y) {}

    Vector2d rotcw() { return Vector2d(y, -x); }

    // 向量相加
    Vector2d operator+(const Vector2d& other) const {
        return Vector2d(x + other.x, y + other.y);
    }

    // 向量相减
    Vector2d operator-(const Vector2d& other) const {
        return Vector2d(x - other.x, y - other.y);
    }

    // 向量与标量相乘
    Vector2d operator*(double scalar) const {
        return Vector2d(x * scalar, y * scalar);
    }

    // 向量与标量相除
    Vector2d operator/(double scalar) const 
    {
            return Vector2d(x / scalar, y / scalar);
    }

    // 复合赋值运算符 +=
    Vector2d& operator+=(const Vector2d& other) {
        x += other.x;
        y += other.y;
        return *this;
    }

    // 复合赋值运算符 -=
    Vector2d& operator-=(const Vector2d& other) {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    // 重载赋值运算符 =
    Vector2d& operator=(const Vector2d& other) {
        if (this != &other) {
            x = other.x;
            y = other.y;
        }
        return *this;
    }

    // 向量模长
    double length() const 
    {
        return std::sqrt(x * x + y * y);
    }

    // 向量点积
    double dot(const Vector2d& other) const {
        return x * other.x + y * other.y;
    }

    // 向量叉积
    double cross(const Vector2d& other) const {
        return x * other.y - y * other.x;
    }

    // 归一化向量
    Vector2d normalize() const {
      double mag = length();
      return *this / mag;
    }

public:
    double x;
    double y;
};

using Vector = Vector2d;
using Point  = Vector2d;

/** 计算两个线段 [p0, p1], [q0, q1] 的交点，交点为 (1-t)*p0 + t*p1 */
bool intersection_point_of_two_segments(const Point & p0, const Point & p1, 
    const Point & q0, const Point & q1, double & t, double eps)
{
  Vector v0 = p1-p0;
  Vector v1 = q0-q1;
  Vector v2 = q0-p0;
  double l = v0.length();
  double v = v0.cross(v1);
  if(std::abs(v) < 1e-14)
  {
    return false;
  }
  t = v2.cross(v1)/v;
  if(t*l<l+eps && t*l>-eps)
  {
    l = v1.length();
    double s = v0.cross(v2)/v;
    if(s*l<l+eps && s*l>-eps)
      return true;
  }
  return false;
}

/** 判断点 p 是否在 [p0, p1] 所在直线的左边 */
inline bool is_on_the_left(const Point & p0, const Point & p1, const Point & p)
{
  Vector v0 = p1-p0;
  Vector v1 = p-p0;
  return v0.cross(v1)>0;
}

}
#endif /* _GEOMETRY_ */ 
