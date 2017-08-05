#include <fstream>
#include <iostream>

#include "WillItStandEngine.h"
#include "mat4.h"
#include "GeometricTools.h"
#include "RboxPoints.h"
#include "QhullError.h"
#include "Qhull.h"
#include "QhullQh.h"
#include "QhullFacet.h"
#include "QhullFacetList.h"
#include "QhullLinkedList.h"
#include "QhullVertex.h"
#include "QhullSet.h"
#include "QhullVertexSet.h"
#include <map>
#include <fstream>

using namespace std;
using namespace orgQhull;

/************************************************************************/
/* Constructors                                                          */
/************************************************************************/
WillItStand::WillItStand()
{
    m_exportCH = false;
    m_exportCenterOfMass = false;
    m_exportResults = false;
    m_isReady = false;
}
WillItStand::WillItStand(string filename, bool outputResult, bool outputCH, bool outputCenterOfMass, string outputFolder) :
m_inputFileName(filename), m_exportResults(outputResult), m_exportCH(outputCH), m_exportCenterOfMass(outputCenterOfMass)
{
    string fileNameWithExt = filename.substr(filename.find_last_of('\\') + 1, filename.size() - filename.find_last_of('\\'));
    string fileNameWithoutExt = fileNameWithExt.substr(0, fileNameWithExt.find_last_of('.'));

    m_outCHFileName =           outputFolder + "\\" + fileNameWithoutExt + "_ConvexHull.obj";
    m_outCenterOfMassFileName = outputFolder + "\\" + fileNameWithoutExt + "_CenterOfMass.obj";
    m_outResultsFileName =      outputFolder + "\\" + fileNameWithoutExt + "_BalancingFaces.obj";
}

/************************************************************************/
/* Public Getters                                                       */
/************************************************************************/
string WillItStand::GetFileName()
{
    return m_inputFileName;
}
Qhull* WillItStand::GetConvexHull()
{
    if (!m_isReady)
    {
        throw std::exception("WillItStand did not run yet or did not run properly, please use Run(); method");
    }
    return &m_qhull; //Very bad, I know, it's just to get it out
}
vec3 WillItStand::GetCenterOfMass()
{
    if (!m_isReady)
    {
        throw std::exception("WillItStand did not run yet or did not run properly, please use Run(); method");
    }
    return m_centroid;
}
vector<Polygon3D> WillItStand::GetBalancingFacets()
{
    if (!m_isReady)
    {
        throw std::exception("WillItStand did not run yet or did not run properly, please use Run(); method");
    }
    return m_result;
}
vector<Polygon3D> WillItStand::GetFacetsIntersectingCHull()
{
    if (!m_isReady)
    {
        throw std::exception("WillItStand did not run yet or did not run properly, please use Run(); method");
    }
    return m_intersectedPolygons;
}
vector<Polygon3D> WillItStand::GetCenterOfMassAsMesh()
{
    if (!m_isReady)
    {
        throw std::exception("WillItStand did not run yet or did not run properly, please use Run(); method");
    }
    return m_centerOfMassMesh;
}

/************************************************************************/
/* Main Functions                                                       */
/************************************************************************/
void WillItStand::ParseInput(int argc, char** argv)
{
    if (argc < 3)
    {
        PrintUsage();
        throw exception("Invalid input: not enough parameters");
    }
    bool inValidInput = false;
    for (int i = 1; i < argc; i++)
    {
        if (!string("-f").compare(argv[i]))
        {
            if (i + 1 == argc)
            {
                throw exception("Invalid input: no file name specified");
            }

            string fileName = argv[++i];
            if (FILE *file = fopen(fileName.c_str(), "r")) {
                fclose(file);
                m_inputFileName = fileName;
                inValidInput = true;
            }
            else {
                throw exception((string("Can't open file: ") + fileName).c_str());
            }
        }
        else  if (!string("--export_ch").compare(argv[i]))
        {
            if (i + 1 == argc)
            {
                throw exception("Invalid input: no file name specified");
            }
            m_outCHFileName = argv[++i];
            m_exportCH = true;
        }
        else  if (!string("--export_com").compare(argv[i]))
        {
            if (i + 1 == argc)
            {
                throw exception("Invalid input: no file name specified");
            }
            m_outCenterOfMassFileName = argv[++i];
            m_exportCenterOfMass = true;
        }
        else  if (!string("--export_result").compare(argv[i]))
        {
            if (i + 1 == argc)
            {
                throw exception("Invalid input: no file name specified");
            }
            m_outResultsFileName = argv[++i];
            m_exportResults = true;
        }
        else
        {
            throw exception((string("Invalid input. No such option exists: ") + argv[i]).c_str());
        }
    }

    if (!inValidInput)
    {
        PrintUsage();
        throw exception("Invalid input: Failed to get input file");

    }
}
bool WillItStand::Run(ostream& output)
{
    m_isReady = false;
    if (m_inputFileName.empty())
    {
        throw std::exception("Please run ParseInput(..); prior to Run();");
    }
    try
    {
        output << "Processing Mesh...\n";
        ProcessMesh();
        output << "Done Processing Mesh.\n";

        output << "Calculating Convex Hull...\n";
        CalcConvexHull("Qt");
        output << "Done Calculating Convex Hull.\n";

        CalcCenterOfMass();
        CreatMeshFromCOM();
        CalcBalancingFacets(CreateMeshFromCH());
        m_isReady = true;

        OutputResult(output);

        if (m_exportCH)
        {
            ExportQhullToOBJFile(&m_qhull, m_outCHFileName);
            output << "Exported Convex Hull Obj to: " << m_outCHFileName << "\n";
        }
        if (m_exportCenterOfMass)
        {
            ExportPolygonsToOBJFile(m_centerOfMassMesh, m_outCenterOfMassFileName);
            output << "Exported Center Of Mass as Obj to: " << m_outCenterOfMassFileName << "\n";

        }
        if (m_exportResults)
        {
            ExportPolygonsToOBJFile(m_result, m_outResultsFileName);
            output << "Exported Result Obj to: " << m_outResultsFileName << "\n";
        }
    }
    catch (std::exception& e)
    {
        output << "Caught an exception while running WillItStand:\n" << e.what() << "\n";
    }
    return m_isReady;
}

/************************************************************************/
/* Private Auxiliary Methods                                         */
/************************************************************************/

void WillItStand::OutputResult(ostream& output)
{
    output << "\n\n============================================================\n";
    output << "                        Results\n";
    output << "============================================================\n\n";
    output << "Input file: " << m_inputFileName.substr(m_inputFileName.find_last_of('\\') + 1, m_inputFileName.size() - m_inputFileName.find_last_of('\\')) << "\n";

    output << "Center of Mass = (" << m_centroid << ")\n";

    output << "Object can stand on the following: \n\n";
    int i = 0;
    for each(Polygon3D poly in m_result)
    {
        output << i++ << ")\n";
        output << "\tNormal: " << poly.Normal() << "\n";
        output << "\tVertices:\n";
        int j = 0;
        for each (vec3 p in poly.Vertices())
        {
            output << "\t\t" << j++ << ") " << p << "\n";
        }
        output << "--------------------------------------\n";
    }
}
void WillItStand::PrintUsage()
{
    cout << "Will It Stand - By Ziv Shahaf & Kujan Lauz:\n";
    cout << "Usage:\n";
    cout << "         -f <inputFileFullPath>\n";
    cout << "Optional:\n";
    cout << "         --export_ch           <outputFileFullPath>\n";
    cout << "         --export_com          <outputFileFullPath>\n";
    cout << "         --export_result       <outputFileFullPath>\n";

}

void WillItStand::ProcessMesh()
{
    m_mesh.Parse(m_inputFileName);
}
void WillItStand::CalcConvexHull(const char* args)
{
    //This runs the qhull algorithm
    m_qhull.runQhull3D(m_mesh.GetVertices(), args);

    //Output qhull result for sanity:
    m_qhull.outputQhull();
    if (m_qhull.hasQhullMessage())
    {
        cout << "\nResults:\n" << m_qhull.qhullMessage();
        m_qhull.clearQhullMessage();
    }
}
void WillItStand::CalcCenterOfMass()
{
    m_centroid = GeometricTools::CalcCOM_3(m_mesh.GetPolygons());
   // m_centroid = GeometricTools::calcCOM_StackOverflow(m_mesh.GetPolygons());
    /*
    vector<vec3> com;
    double mass = 0;
    for each (Polygon3D poly in m_mesh.GetPolygons())
    {
        double m = 0;
        vec3 c;
        vec3 v0 = poly.Vertices()[0];
        vec3 v1 = poly.Vertices()[1];
        vec3 v2 = poly.Vertices()[2];
        GeometricTools::CalcMassAndCenterOfMassTriangle(v0, v1, v2, c, m);
        com.push_back(c);
        mass += m;
    }

    for each (vec3 v in com)
    {
        m_centroid[0] += v.x();
        m_centroid[1] += v.y();
        m_centroid[2] += v.z();
    }
    m_centroid /= 24.0;
    m_centroid /= mass;
    */
}
void WillItStand::CreatMeshFromCOM()
{
    vec3 v1 = m_centroid + vec3(-0.05, -0.05, 0.05);
    vec3 v2 = m_centroid + vec3(-0.05, 0.05, 0.05);
    vec3 v3 = m_centroid + vec3(0.05, 0.05, 0.05);
    vec3 v4 = m_centroid + vec3(0.05, -0.05, 0.05);
    vec3 v5 = m_centroid + vec3(-0.05, -0.05, -0.05);
    vec3 v6 = m_centroid + vec3(-0.05, 0.05, -0.05);
    vec3 v7 = m_centroid + vec3(0.05, 0.05, -0.05);
    vec3 v8 = m_centroid + vec3(0.05, -0.05, -0.05);

    Polygon3D p1({ v4, v3, v2, v1 }, vec3(), 0);
    Polygon3D p2({ v1, v2, v6, v5 }, vec3(), 0);
    Polygon3D p3({ v7, v6, v2, v3 }, vec3(), 0);
    Polygon3D p4({ v4, v8, v7, v3 }, vec3(), 0);
    Polygon3D p5({ v5, v8, v4, v1 }, vec3(), 0);
    Polygon3D p6({ v5, v6, v7, v8 }, vec3(), 0);

    m_centerOfMassMesh.push_back(p1);
    m_centerOfMassMesh.push_back(p2);
    m_centerOfMassMesh.push_back(p3);
    m_centerOfMassMesh.push_back(p4);
    m_centerOfMassMesh.push_back(p5);
    m_centerOfMassMesh.push_back(p6);
}
void WillItStand::CalcBalancingFacets(vector<Polygon3D> model)
{
    //Test every intersected polygon with projected centroid
    for each (Polygon3D poly in model)
    {
        vec3 projectedPoint = GeometricTools::CalculateProjectionPoint(m_centroid, poly.Normal(), poly.NormalOffset());

        std::vector<double> polygonXCoords;
        std::vector<double> polygonYCoords;

        //Apply transformation for all points to be on XY Plane (according to current plane)
        //The rotation-axis is not the z-axis, it is the cross product of the normal of the triangle and the z-axis:
        vec3 rot_axis = GeometricTools::normalize(GeometricTools::cross(vec3(0, 0, 1), poly.Normal()));
        //    the rot_angle that you already calculated :
        double rot_angle = acos(GeometricTools::dot(vec3(0, 0, 1), poly.Normal()));
        mat4 tranformMatrix = GeometricTools::CreateFromAxisAngle(rot_axis, rot_angle);

        for (int i = 0; i < poly.numVertices(); i++)
        {
            polygonXCoords.push_back((poly.Vertices()[i] * tranformMatrix).x());
            polygonYCoords.push_back((poly.Vertices()[i] * tranformMatrix).y());
        }

        double x, y;
        x = (projectedPoint *tranformMatrix).x();
        y = (projectedPoint *tranformMatrix).y();
		if (GeometricTools::IsPointInPolygon(poly.numVertices(), polygonXCoords, polygonYCoords, x, y))
		{

            //double dist = 0;
            //double cx = 0;
            //double cy = 0;
            //for (int i = 0; i < poly.numVertices(); i++)
            //{
            //	cx += polygonXCoords[i];
            //	cy += polygonYCoords[i];
            //}
            //cx = cx / poly.numVertices();
            //cy = cy / poly.numVertices();
            //double distance = std::sqrt((cx - x)*(cx - x) + (cy - y)* (cy - y));
            //double longsestDist = 0;
            //for (int i = 0, j = 1; i < poly.numVertices(); i++, j = (j + 1) % poly.numVertices())
            //{
            //	double tempDist = std::sqrt((polygonXCoords[i] - polygonXCoords[j])*(polygonXCoords[i] - polygonXCoords[j]) + (polygonYCoords[i] - polygonYCoords[j])* (polygonYCoords[i] - polygonYCoords[j]));
            //	if (tempDist > longsestDist)
            //		longsestDist = tempDist;
            //}
            //if (distance < distanceThershold * longsestDist)
            //{
            	m_result.push_back(poly);
            //}
		}

    }
}
vector<Polygon3D> WillItStand::CreateMeshFromCH()
{
    if (!m_qhull.initialized())
    {
        throw exception("CreateMeshFromCH cannot run before Qhull.Run is invoked");
    }
    vector<Polygon3D> chull;
    QhullFacetList facets = m_qhull.facetList();
    for (QhullFacetList::iterator it = facets.begin(); it != facets.end(); ++it)
    {
        if (!(*it).isGood()) continue;
        QhullFacet f = *it;
        QhullVertexSet vSet = f.vertices();
        vector<vec3> vertices;
        for (QhullVertexSet::iterator vIt = vSet.begin(); vIt != vSet.end(); ++vIt)
        {

            QhullVertex v = *vIt;
            QhullPoint p = v.point();

            double * coords = p.coordinates();
            vertices.push_back(vec3(coords[0], coords[1], coords[2]));
        }

        if (f.hyperplane().isDefined())
        {
            auto coord = f.hyperplane().coordinates();
            vec3 normal(coord[0], coord[1], coord[2]);
            double offset = f.hyperplane().offset();
            chull.push_back(Polygon3D(vertices, normal, offset));
        }
    }
    return chull;
}

void WillItStand::ExportQhullToOBJFile(orgQhull::Qhull* qhull, string fileName)
{
    map<int, vec3> vertices_by_ids;
    vector<ThreeInts> vertices_ids__per_face;
    int numVertices = qhull->vertexList().toStdVector().size();

    ofstream outFile;
    outFile.open(fileName);


    QhullFacetList facets = qhull->facetList();
    for (QhullFacetList::iterator it = facets.begin(); it != facets.end(); ++it)
    {
        if (!(*it).isGood()) continue;
        QhullFacet f = *it;
        vector<int> Ids;
        QhullVertexSet vSet = f.vertices();
        vec3 aPoint;
        for (QhullVertexSet::iterator vIt = vSet.begin(); vIt != vSet.end(); ++vIt)
        {
            QhullVertex v = *vIt;
            QhullPoint p = v.point();
            double * coords = p.coordinates();
            Ids.push_back(v.id());
            vertices_by_ids[v.id()] = vec3(coords[0], coords[1], coords[2]);
            aPoint = vec3(coords[0], coords[1], coords[2]);
        }

        vec3 normal = GeometricTools::cross(vertices_by_ids[Ids[0]] - vertices_by_ids[Ids[1]], vertices_by_ids[Ids[0]] - vertices_by_ids[Ids[2]]);
        vec3 diffVec = m_centroid - aPoint;
        if (GeometricTools::dot(diffVec, normal) > 0)//To Make the Obj file normals consistent
            vertices_ids__per_face.push_back({ Ids[2], Ids[1], Ids[0] }); //We are sure there are only there because we used "Qt" argument for qhul.Run()
        else
            vertices_ids__per_face.push_back({ Ids[0], Ids[1], Ids[2] }); //We are sure there are only there because we used "Qt" argument for qhul.Run()
        
    }
    int idx = 0;

    outFile << "#Vertices" << endl;

    for each (auto v in vertices_by_ids)
    {
        int vId = v.first;
        if (vId == idx)
        {
            vec3 vec = v.second;
            outFile << "v " << vec << endl;
        }
        else
        {
            if (idx < vId)
            {
                for (int j = idx; j < vId; j++)
                {
                    outFile << "v " << "000000000000000000 000000000000000000 000000000000000000" << endl; //no one will use this...
                }
                idx = vId;
            }
            vec3 vec = v.second;
            outFile << "v " << vec << endl;
        }
        idx++;
    }

    outFile << "#Faces" << endl;
    for each (auto f in vertices_ids__per_face)
    {
        outFile << "f " << f.a + 1 << " " << f.b + 1 << " " << f.c + 1 << endl;
    }
    outFile.close();
}
void WillItStand::ExportPolygonsToOBJFile(const vector<Polygon3D>& polygons, string fileName)
{
    ofstream outFile;
    outFile.open(fileName);

    //Go over m_mesh's polygons and points and save to file
    outFile << "#vertices:\n";
    int i = 1;
    for each (Polygon3D p in polygons)
    {
        for each (vec3 v in p.Vertices())
        {
            outFile << "v " << v.x() << " " << v.y() << " " << v.z() << endl;
        }

        outFile << "f "; // << i++  << " " << i++ << " " << i++ << endl;
        vec3 normal = GeometricTools::cross(p.Vertices()[0] - p.Vertices()[1], p.Vertices()[0] - p.Vertices()[2]);
        vec3 diffVec = m_centroid - p.Vertices()[0];
        vector<int> ids;
        for (unsigned int j = 0; j < p.Vertices().size(); j++)
        {
            ids.push_back(i++);
        }

        if (GeometricTools::dot(diffVec, normal) > 0)//To Make the Obj file normals consistent
        {//2 1 0
            for (int j = p.Vertices().size() - 1; j >= 0; j--)
            {
                outFile << ids[j] << " ";
            }
        }
        else
        {           // 0 1 2
            for (unsigned int j = 0; j < p.Vertices().size(); j++)
            {
                outFile << ids[j] << " ";
            }
        }




        outFile << endl;
    }
    outFile.close();
}


/************************************************************************/
/*                       Not in used methods                            */
/************************************************************************/

void WillItStand::CalcCentroid()
{
    m_centroid = GeometricTools::CalculateCentroid(m_mesh.GetVertices());
}
void WillItStand::CalcSilhouetteFacets()
{
    //Get all convex hull's normals and offsets
    vector<std::pair<vec3, double> > facetsNormals;
    for each (QhullFacet facet in m_qhull.facetList().toStdVector())
    {
        if (facet.hyperplane().isDefined())
        {
            auto coord = facet.hyperplane().coordinates();
            vec3 normal(coord[0], coord[1], coord[2]);
            double offset = facet.hyperplane().offset();
            facetsNormals.push_back(pair<vec3, double>(normal, offset));
        }
    }

    //Get Intersecting polygons with convex hull
    m_intersectedPolygons = GeometricTools::GetIntersectingPolygons(m_mesh.GetPolygons(), facetsNormals);
}

