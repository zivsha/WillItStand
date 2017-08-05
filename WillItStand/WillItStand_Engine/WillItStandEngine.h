#pragma once
#include "MeshModel.h"
#include "Qhull.h"
#include <vector>
#include "polygon3D.h"

typedef struct MyStruct
{
    int a, b, c;
}ThreeInts;


class WillItStand
{
private:
    bool m_isReady;
    bool m_exportCH;
    bool m_exportCenterOfMass;
    bool m_exportResults;
    string m_inputFileName;
    string m_outCHFileName;
    string m_outCenterOfMassFileName;
    string m_outResultsFileName;
    orgQhull::Qhull m_qhull;
    vec3 m_centroid;
    vector<Polygon3D> m_intersectedPolygons;
    vector<Polygon3D> m_result;
    vector<Polygon3D> m_centerOfMassMesh;
    MeshModel m_mesh;

    void ProcessMesh();
    void CalcConvexHull(const char* args);
    void CalcCentroid();
    void CalcSilhouetteFacets();
    void CalcBalancingFacets(vector<Polygon3D> model);
    void CalcCenterOfMass();
    void PrintUsage();
    void OutputResult(ostream& output);
    vector<Polygon3D> CreateMeshFromCH();
    void CreatMeshFromCOM();
public:
    WillItStand();
    WillItStand(string filename, bool outputResult, bool outputCH, bool outputCenterOfMass, string outputFolder);
    orgQhull::Qhull* GetConvexHull();
    vec3 GetCenterOfMass();
    vector<Polygon3D> GetBalancingFacets();
    vector<Polygon3D> GetFacetsIntersectingCHull();
    vector<Polygon3D> GetCenterOfMassAsMesh();
    string GetFileName();

    void ParseInput(int argc, char** argv);
    bool Run(ostream& output);

    void ExportPolygonsToOBJFile(const vector<Polygon3D>& polygons, string fileName);
    void ExportQhullToOBJFile(orgQhull::Qhull* qhull, string fileName);


};
