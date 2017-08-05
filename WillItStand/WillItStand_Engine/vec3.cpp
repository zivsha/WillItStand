#include "vec3.h"
#include <ostream>
#include "mat4.h"

std::ostream& operator<<(std::ostream &os, const vec3& v)
{
    return os << v.x() << " " << v.y() << " " << v.z();
}

vec3 operator*(const vec3& p, const mat4& m)
{

    double x = p.x() * m.at(0, 0) + p.y() * m.at(0, 1) + p.z() * m.at(0, 2);
    double y = p.x() * m.at(1, 0) + p.y() * m.at(1, 1) + p.z() * m.at(1, 2);
    double z = p.x() * m.at(2, 0) + p.y() * m.at(2, 1) + p.z() * m.at(2, 2);

    return vec3(x, y, z);
}

vec3 operator-(const vec3& v1, const vec3& v2)
{
    return vec3(v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z());
}

vec3 operator+(const vec3& v1, const vec3& v2)
{
    return vec3(v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z());
}
vec3::vec3(double x, double y, double z)
{
    vec[0] = x;
    vec[1] = y;
    vec[2] = z;
}
vec3::vec3(const vec3& v)
{
    vec[0] = v[0];
    vec[1] = v[1];
    vec[2] = v[2];
}
const double& vec3::operator[](int i) const
{
    if (i<0 || i >2)
        throw new std::exception();
    return vec[i];
}
double& vec3::operator[](int i)
{
    if (i<0 || i >2)
        throw new std::exception("Invalid index");
    return vec[i];
}

vec3 vec3::operator-() const
{
    return vec3(-vec[0], -vec[1], -vec[2]);
}
vec3& vec3::operator+=(const vec3& v)
{
    for (int i = 0; i < 3; i++)
        vec[i] += v[i];
    return *this;
}
vec3& vec3::operator-=(const vec3& v)
{
    return (*this) += -v;
}
vec3& vec3::operator*=(double a)
{
    for (int i = 0; i < 3; i++)
        vec[i] *= a;
    return *this;
}
vec3& vec3::operator*=(vec3& v)
{
    for (int i = 0; i < 3; i++)
        vec[i] *= v[i];
    return *this;
}
vec3& vec3::operator/=(double a)
{
    for (int i = 0; i < 3; i++)
        vec[i] /= a;
    return *this;
}

double vec3::x() const
{
    return vec[0];
}
double vec3::y() const
{
    return vec[1];
}
double vec3::z() const
{
    return vec[2];
}

vec3& vec3::normalize()
{
    double w = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    if (w == 0)
    {
        return *this;
    }
    for (int i = 0; i < 3; i++)
    {
        vec[i] /= w;
    }
    return *this;
}
bool vec3::Equals(const vec3& other)
{
    return     abs(vec[0] - other.x()) <= std::numeric_limits<double>::epsilon()
            && abs(vec[1] - other.y()) <= std::numeric_limits<double>::epsilon()
            && abs(vec[2] - other.z()) <= std::numeric_limits<double>::epsilon();
}