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

class Vector3d {
public:
    // 构造函数
    Vector3d(double x_ = 0.0, double y_ = 0.0, double z_ = 0.0) : x(x_), y(y_), z(z_) {}

    // 复制构造函数
    Vector3d(const Vector3d& other) : x(other.x), y(other.y), z(other.z) {}

    Vector3d rotcwxy() { return Vector3d(y, -x, z); }
    Vector3d rotcwxz() { return Vector3d(-z, y, x); }
    Vector3d rotcwyz() { return Vector3d(x, -z, y); }

    // 向量相加
    Vector3d operator+(const Vector3d& other) const {
        return Vector3d(x + other.x, y + other.y, z + other.z);
    }

    // 向量相减
    Vector3d operator-(const Vector3d& other) const {
        return Vector3d(x - other.x, y - other.y, z - other.z);
    }

    // 向量与标量相乘
    Vector3d operator*(double scalar) const {
        return Vector3d(x * scalar, y * scalar, z * scalar);
    }

    // 向量与标量相除
    Vector3d operator/(double scalar) const {
        return Vector3d(x / scalar, y / scalar, z / scalar);
    }

    // 复合赋值运算符 +=
    Vector3d& operator+=(const Vector3d& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    // 复合赋值运算符 -=
    Vector3d& operator-=(const Vector3d& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    // 重载赋值运算符 =
    Vector3d& operator=(const Vector3d& other) {
        if (this != &other) {
            x = other.x;
            y = other.y;
            z = other.z;
        }
        return *this;
    }

    // 向量模长
    double length() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    // 向量点积
    double dot(const Vector3d& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    // 向量叉积
    Vector3d cross(const Vector3d& other) const {
        return Vector3d(y * other.z - z * other.y,
                        z * other.x - x * other.z,
                        x * other.y - y * other.x);
    }

    // 归一化向量
    Vector3d normalize() const {
        double mag = length();
        return *this / mag;
    }

    /**
     * @brief 计算一个绕任意轴旋转的结果
     * @param axis 旋转轴的单位向量
     * @param angle 旋转角度（弧度）
     */
    Vector3d rotate(const Vector3d& axis, double angle) const {
        double cosAngle = std::cos(angle);
        double sinAngle = std::sin(angle);
        double oneMinusCos = 1.0 - cosAngle;

        double xRot = (cosAngle + axis.x * axis.x * oneMinusCos) * x +
                      (axis.x * axis.y * oneMinusCos - axis.z * sinAngle) * y +
                      (axis.x * axis.z * oneMinusCos + axis.y * sinAngle) * z;

        double yRot = (axis.y * axis.x * oneMinusCos + axis.z * sinAngle) * x +
                      (cosAngle + axis.y * axis.y * oneMinusCos) * y +
                      (axis.y * axis.z * oneMinusCos - axis.x * sinAngle) * z;

        double zRot = (axis.z * axis.x * oneMinusCos - axis.y * sinAngle) * x +
                      (axis.z * axis.y * oneMinusCos + axis.x * sinAngle) * y +
                      (cosAngle + axis.z * axis.z * oneMinusCos) * z;

        return Vector3d(xRot, yRot, zRot);
    }

public:
    double x;
    double y;
    double z;
};

using Point2d  = Vector2d;
using Point3d  = Vector3d;

}
#endif /* _GEOMETRY_ */ 
