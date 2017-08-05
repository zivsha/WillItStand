#ifndef vec3_h__
#define vec3_h__

#include <ostream>
#include "mat4.h"

class vec3 {
    double vec[3];
public:
    vec3(double x = 0.0, double y = 0.0, double z = 0.0);
    vec3(const vec3& v);
    const double& operator[](int i) const;
    double& operator[](int i);

    vec3 operator-() const;
    vec3& operator+=(const vec3& v);
    vec3& operator-=(const vec3& v);
    vec3& operator*=(double a);
    vec3& operator*=(vec3& v);
    vec3& operator/=(double a);

    double x() const;
    double y() const;
    double z() const;

    vec3& normalize();
    bool Equals(const vec3& other);

    friend std::ostream& operator<<(std::ostream &os, const vec3& v);

};

vec3 operator-(const vec3& v1, const vec3& v2);
vec3 operator*(const vec3& p, const mat4& m);
vec3 operator+(const vec3& v1, const vec3& v2);

#endif // vec3_h__
