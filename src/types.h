#ifndef TYPES_H
#define TYPES_H

#include "header.h"
#include "cmath"

template<typename T>
class Vector2 {
public:
  Vector2() { set( 0,0 ); }
  Vector2(const T& x, const T& y) { set(x,y); }

  void set(const T& x, const T& y) {
    m_x = x;
    m_y = y;
  }
  T x() { return m_x; }
  T y() { return m_y; }
  T mod2() { return (m_x*m_x+m_y*m_y); }
  T mod() { return sqrt(mod2()); }
protected:
  T m_x;
  T m_y;
};

template<typename T>
class Vector3 {
public:
  Vector3() { set(0,0,0); }
  Vector3(const T& x, const T& y, const T& z) { set(x,y,z); }
  void set(const T& x, const T& y, const T& z) {
    m_x = x;
    m_y = y;
    m_z = z;
  }
  T mod2() { return (m_x*m_x+m_y*m_y+m_z*m_z); }
  T mod() { return sqrt(mod2()); }

protected:
  T m_x;
  T m_y;
  T m_z;
};

typedef Vector2<double> Vector2d;
typedef Vector2<int>    Vector2i;

typedef Vector3<double> Vector3d;
typedef Vector3<int>    Vector3i;

#endif // TYPES_H
