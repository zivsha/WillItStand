#include "MeshModel.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "GeometricTools.h"
#include <shlwapi.h>

using namespace std;

struct FaceIdcs
{
	int v[4];
	int vn[4];

	FaceIdcs()
	{
		for (int i=0; i<4; i++)
			v[i] = vn[i] = 0;
	}

	FaceIdcs(std::istream & aStream)
	{
		for (int i=0; i<4; i++)
			v[i] = vn[i] = 0;

		char c;
		for(int i = 0; i < 3; i++)
		{
			aStream >> std::ws >> v[i] >> std::ws;

			if (aStream.peek() != '/')
				continue;

			aStream >> c >> std::ws;
			
            if (aStream.peek() == '/')
			{
				aStream >> c >> std::ws >> vn[i];
				continue;
			}

			if (aStream.peek() != '/')
				continue;

			aStream >> c >> vn[i];
		}
	}
};

vec3 vec3fFromStream(std::istream & aStream)
{
	double x, y, z;
	aStream >> x >> std::ws >> y >> std::ws >> z;
    return vec3(x, y, z);
}

void MeshModel::Parse(string fileName)
{
    char* ext = PathFindExtensionA(fileName.c_str());
    if (_strcmpi(ext, ".obj") != 0)
    {
        string excep("Invalid input file: ");
        excep += fileName + "\nFile must be a valid OBJ file.";
        throw exception(excep.c_str());
    }
	ifstream ifile(fileName.c_str());
	vector<FaceIdcs> faces;
    //vector<vec3> vertices;
    vector<vec3> normals;
	
	while (!ifile.eof())
	{
		string curLine;
		getline(ifile, curLine);
		istringstream issLine(curLine);
		string lineType;
		issLine >> std::ws >> lineType;

		if (lineType == "v")
            vertices.push_back(vec3fFromStream(issLine));
		else if (lineType == "vn")
			normals.push_back(vec3fFromStream(issLine));
		else if (lineType == "f")
			faces.push_back(issLine);
		else if (lineType == "#" || lineType == "") {}
			// comment / empty line
		else
			cout<< "Found unknown line Type \"" << lineType << "\"" << endl;
	}

    /*
        vector<FaceIdcs> faces - Holds the indices of the vertices for every face;
        vector<vec3> vertices - Holds all vertices
        vector<vec3> normals;
    */
	for (vector<FaceIdcs>::iterator it = faces.begin(); it != faces.end(); ++it)
	{
		for (int i = 0; i < 3; i++)
		{
			vertex_positions.push_back(vertices[it->v[i]-1]);
		}
	}

	if(normals.size() != 0)
	{
		for (vector<FaceIdcs>::iterator it = faces.begin(); it != faces.end(); ++it)
		{
			for (int i = 0; i < 3; i++)
			{
				normal_positions.push_back(normals[it->vn[i]-1]);
			}
		}
	}

    for (vector<vec3>::iterator it = vertex_positions.begin(); it != vertex_positions.end(); ++it)
	{
        vec3 v1 = *it;
        vec3 v2 = *(++it);
        vec3 v3 = *(++it);

        vec3 n = GeometricTools::normalize(GeometricTools::cross(v1 - v2, v1 - v3));
		face_normal_positions.push_back(n);
		face_normal_positions.push_back(n);
		face_normal_positions.push_back(n);
        //-(Ax+By+Cz)
        double d = -(n.x()*v1.x() + n.y()*v1.y() + n.z()*v1.z());
        polygons.push_back(Polygon3D({ v1, v2, v3 }, n, d));
	}
}
vector<vec3> MeshModel::GetVertices(){
    return vertex_positions;
}
vector<Polygon3D> MeshModel::GetPolygons(){
    return polygons;
}

vector<vec3> MeshModel::GetVerticesSet()
{
    return vertices;
}
