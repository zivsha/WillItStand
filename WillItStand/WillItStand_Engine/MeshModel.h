#pragma once

#include <string>
#include <vector>
#include "vec3.h"
#include "Polygon3D.h"

using namespace std;

class MeshModel
{
private:
    vector<vec3> vertex_positions;
    vector<vec3> normal_positions;
    vector<vec3> vertices;
    vector<vec3> face_normal_positions;
    vector<Polygon3D> polygons;
public:
    MeshModel() {}
    void Parse(string fileName);
   // vec3 GetCenteroid();
    vector<vec3> GetVertices();
    vector<Polygon3D> GetPolygons();
    vector<vec3> GetVerticesSet();
};
