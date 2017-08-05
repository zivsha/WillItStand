#pragma once
#include <vector>
#include "vec3.h"

using std::vector;

class Polygon3D {

    vec3 normal;
    double normalOffset;
    vector<vec3> vertices;
public:
    Polygon3D(vector<vec3> vertices, vec3 normal, double normalOffset);
    int numVertices();
    vector<vec3> Vertices();
    vec3 Normal();
    double NormalOffset();

};
