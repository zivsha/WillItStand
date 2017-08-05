#pragma once
#include "vec3.h"
#include "mat4.h"
#include <vector>
#include "Polygon3D.h"
#include <utility>

using std::vector;
using std::pair;
class GeometricTools
{
public:
    static bool IsPointInPolygon(int numVertices, std::vector<double>& polygonXCoords, std::vector<double>& polygonYCoords, double pointX, double pointY);
    static vec3 CalculateCentroid(const vector<vec3>& points);
    static vec3 CalculateProjectionPoint(vec3 point, vec3 normal, double d);
    static vector<Polygon3D> GetIntersectingPolygons(const vector<Polygon3D>& originalPolygons, const vector<pair<vec3, double> >& facetsNormals);
    static double dot(const vec3& a, const vec3& b);
    static vec3 cross(const vec3& a, const vec3& b);
    static vec3 normalize(const vec3& a);
    static mat4 CreateFromAxisAngle(const vec3& axis, double angle);
    static vec3 CalcCenterOfMassTriangle(const vec3 & p0, const vec3 & p1, const vec3 & p2);
    static vec3 CalcCenterOfMassTriangle(const vec3 & p0, const vec3 & p1, const vec3 & p2, const vec3& n);
    static void CalcMassAndCenterOfMassTriangle(const vec3 & v0, const vec3 & v1, const vec3 & v2, vec3& c,  double& m);
    static void CalcMassAndCenterOfMassTriangle(const vec3 & v0, const vec3 & v1, const vec3 & v2, const vec3& n, vec3& c, double& m);
    static vec3 calcCOM_StackOverflow(const vector<Polygon3D>& polygons);
    static vec3 CalcCOM_3(const vector<Polygon3D>& polygons);
};

