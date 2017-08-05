#include "Polygon3D.h"


Polygon3D::Polygon3D(vector<vec3> vertices, vec3 normal, double normalOffset) : vertices(vertices), normal(normal), normalOffset(normalOffset){}

int Polygon3D::numVertices() 
{
    return vertices.size(); 
}
vector<vec3> Polygon3D::Vertices()
{
    return vertices;
}
vec3 Polygon3D::Normal() 
{ 
    return normal; 
}
double Polygon3D::NormalOffset()
{ 
    return normalOffset;
}
