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

    /**
     * @brief 计算一个旋转矩阵 [[cost, sint], [-sint, cost]] 作用到自身的结果
     */
    Vector2d rotate(const double & sint, const double & cost)
    {
      double nx = x*cost + y*sint;
      double ny = -x*sint + y*cost;
      return Vector2d(nx, ny);
    }

public:
    double x;
    double y;
};

using Vector = Vector2d;
using Point  = Vector2d;

}
#endif /* _GEOMETRY_ */ 
