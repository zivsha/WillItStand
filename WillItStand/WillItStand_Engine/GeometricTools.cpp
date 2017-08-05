#include "GeometricTools.h"
#include <limits>
#include <iostream>
#include <ctime>

using namespace std;

bool GeometricTools::IsPointInPolygon(int numVertices, std::vector<double>& polygonXCoords, std::vector<double>& polygonYCoords, double pointX, double pointY)
{
    bool isPointIn = false;
    int i, j;
    for (i = 0, j = numVertices - 1; i < numVertices; j = i++)
    {
        if (((polygonYCoords[i]>pointY) != (polygonYCoords[j] > pointY)) &&
            (pointX < (polygonXCoords[j] - polygonXCoords[i]) * (pointY - polygonYCoords[i]) / (polygonYCoords[j] - polygonYCoords[i]) + polygonXCoords[i]))
            isPointIn = !isPointIn;
    }
    return isPointIn;
}
vec3 GeometricTools::CalculateProjectionPoint(vec3 point, vec3 normal, double d){
    ;
    double distancePointFromPlane = GeometricTools::dot(point, normal) + d;
    double x = point.x() - distancePointFromPlane * normal.x();
    double y = point.y() - distancePointFromPlane * normal.y();
    double z = point.z() - distancePointFromPlane * normal.z();
    //cout << "Projection Point is: " << x << "," << y << "," << z << endl;
    vec3 projectionPoint(x, y, z);
    return projectionPoint;
}
vector<Polygon3D> GeometricTools::GetIntersectingPolygons(const vector<Polygon3D>& originalPolygons, const vector<pair<vec3, double> >& facetsNormals)
{
    vector<Polygon3D> retVector;
    for each (Polygon3D poly in originalPolygons)
    {
        for each (auto normalOffsetPair in facetsNormals)
        {
           // if (normalOffsetPair.first.Equals(poly.Normal()) || normalOffsetPair.first.Equals(-poly.Normal()))
            if (normalOffsetPair.first.Equals(poly.Normal()) || normalOffsetPair.first.Equals(-poly.Normal()))
            {
                bool polygonIsIntersected = true;
                //check if every point in polygon applies to the plane equation
                for each (vec3 v in poly.Vertices())
                {
                    //Ax + By + Cz + d = 0 ?
                    if (abs(normalOffsetPair.first.x() * v.x() + normalOffsetPair.first.y() * v.y() + normalOffsetPair.first.z() * v.z() + normalOffsetPair.second) > numeric_limits<double>::epsilon())
                    {
                        polygonIsIntersected = false;
                        break;
                    }
                }
                //if so, add it to result
                if (polygonIsIntersected)
                {
                    retVector.push_back(poly);
                }
            }
        }
    }
    return retVector;
}
double GeometricTools::dot(const vec3& a, const vec3& b)
{
    return a.x()*b.x() + a.y()*b.y() + a.z()*b.z();
}
vec3 GeometricTools::cross(const vec3& a, const vec3& b)
{
    return vec3(a.y()*b.z() - a.z()*b.y(), a.z()*b.x() - a.x()*b.z(), a.x()*b.y() - a.y()*b.x());
}
vec3 GeometricTools::normalize(const vec3& a)
{
    double length = sqrt(a.x()*a.x() + a.y()* a.y() + a.z() * a.z());
    if (abs(length) <= std::numeric_limits<double>::epsilon())
    {
        return vec3(0, 0, 0);
    }
    return vec3(a.x() / length, a.y() / length, a.z() / length);
}
mat4 GeometricTools::CreateFromAxisAngle(const vec3& axis, double angle)
{
    double c = cos(angle);
    double s = sin(angle);
    double t = 1.0f - c;

    double x = axis.x();
    double y = axis.y();
    double z = axis.z();

    mat4 C = mat4(
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 0
        ).scale(c);

    C.putAt(3, 3, 1);

    mat4 T = mat4(
        x * x, x * y, x * z, 0,
        x * y, y * y, y * z, 0,
        x * z, y * z, z * z, 0,
        0, 0, 0, 0
        ).scale(t);

    mat4 S = mat4(
        0, z, -y, 0,
        -z, 0, x, 0,
        y, -x, 0, 0,
        0, 0, 0, 0
        ).scale(s);

    C += T;
    C += S;
    return C;
}
vec3 GeometricTools::CalcCenterOfMassTriangle(const vec3 & v0, const vec3 & v1, const vec3 & v2, const vec3& n)
{
    vec3 c;
    vec3 f1 = v0 + v1 + v2;
    vec3 f2;
    f2[0] = f1[0] * f1[0] - (v0[0] * v1[0] + v1[0] * v2[0] + v2[0] * v0[0]);
    f2[1] = f1[1] * f1[1] - (v0[1] * v1[1] + v1[1] * v2[1] + v2[1] * v0[1]);
    f2[2] = f1[2] * f1[2] - (v0[2] * v1[2] + v1[2] * v2[2] + v2[2] * v0[2]);

    // center of mass
    c[0] = f2[0] * n[0];
    c[1] = f2[1] * n[1];
    c[2] = f2[2] * n[2];
    return c;
}
vec3 GeometricTools::CalcCenterOfMassTriangle(const vec3 & v0, const vec3 & v1, const vec3 & v2) {
    // get edges and cross product of edges
    vec3 e1 = v1 - v0;
    vec3 e2 = v0 - v2;
    vec3 e3 = v2 - v1;
    vec3 n = -cross(e1,e2);
    return CalcCenterOfMassTriangle(v0, v1, v2, n);
}


void GeometricTools::CalcMassAndCenterOfMassTriangle(const vec3 & v0, const vec3 & v1, const vec3 & v2, const vec3& n, /*[out]*/ vec3& c, /*[out]*/ double& m)
{
    vec3 f1 = v0 + v1 + v2;
    vec3 f2;
    f2[0] = f1[0] * f1[0] - (v0[0] * v1[0] + v1[0] * v2[0] + v2[0] * v0[0]);
    f2[1] = f1[1] * f1[1] - (v0[1] * v1[1] + v1[1] * v2[1] + v2[1] * v0[1]);
    f2[2] = f1[2] * f1[2] - (v0[2] * v1[2] + v1[2] * v2[2] + v2[2] * v0[2]);

    // mass
    m = f1[0] * n[0];

    // center of mass
    c[0] = f2[0] * n[0];
    c[1] = f2[1] * n[1];
    c[2] = f2[2] * n[2];
}

void GeometricTools::CalcMassAndCenterOfMassTriangle(const vec3 & v0, const vec3 & v1, const vec3 & v2, /*[out]*/ vec3& c, /*[out]*/ double& m) {
    // get edges and cross product of edges
    vec3 e1 = v1 - v0;
    vec3 e2 = v0 - v2;
    vec3 e3 = v2 - v1;
    vec3 n = -cross(e1, e2);
    return CalcMassAndCenterOfMassTriangle(v0, v1, v2, n, c, m);
}
vec3 GeometricTools::CalculateCentroid(const vector<vec3>& points){
    std::clock_t start;
    double duration;

    start = std::clock();
    cout << "Calculating Centroid...";
    double centerX = 0;
    double centerY = 0;
    double centerZ = 0;
    int counter = 0;
    int samplesAmount = 10000;
    for (unsigned int i = 0; i < points.size(); i++)
    {
        double x0 = points[i].x();
        double y0 = points[i].y();
        double z0 = points[i].z();
        for (unsigned int j = i + 1; j < points.size(); j++)
        {
            double x1 = points[j].x();
            double y1 = points[j].y();
            double z1 = points[j].z();
            double dirX = x1 - x0;
            double dirY = y1 - y0;
            double dirZ = z1 - z0;
            for (int k = 0; k < samplesAmount; k++)
            {
                centerX += x0 + dirX * k / samplesAmount;
                centerY += y0 + dirY * k / samplesAmount;
                centerZ += z0 + dirZ * k / samplesAmount;
                counter++;
            }
        }
        duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        if (duration > 5)
        {
            start = std::clock();
            cout << ".";
        }
    }
    cout << endl;
    centerX = centerX / counter;
    centerY = centerY / counter;
    centerZ = centerZ / counter;
    cout << "Centroid = (" << centerX << ", " << centerY << ", " << centerZ << ")" << endl;
    vec3 centroid(centerX, centerY, centerZ);
    return centroid;
}

class Triangle // 3 vertices of each triangle
{
public:
    double x1, y1, z1;
    double x2, y2, z2;
    double x3, y3, z3;
};


vec3 GeometricTools::calcCOM_StackOverflow(const vector<Polygon3D>& polygons)
{
    double totalVolume = 0, currentVolume;
    double xCenter = 0, yCenter = 0, zCenter = 0;
    for each (Polygon3D poly in polygons)
    {
        Triangle triangle;
        triangle.x1 = poly.Vertices()[0].x();
        triangle.y1 = poly.Vertices()[0].y();
        triangle.z1 = poly.Vertices()[0].z();

        triangle.x2 = poly.Vertices()[1].x();
        triangle.y2 = poly.Vertices()[1].y();
        triangle.z2 = poly.Vertices()[1].z();

        triangle.x3 = poly.Vertices()[2].x();
        triangle.y3 = poly.Vertices()[2].y();
        triangle.z3 = poly.Vertices()[2].z();

        totalVolume += currentVolume = (triangle.x1*triangle.y2*triangle.z3 - triangle.x1*triangle.y3*triangle.z2 - triangle.x2*triangle.y1*triangle.z3 + triangle.x2*triangle.y3*triangle.z1 + triangle.x3*triangle.y1*triangle.z2 - triangle.x3*triangle.y2*triangle.z1) / 6;
        xCenter += ((triangle.x1 + triangle.x2 + triangle.x3) / 4) * currentVolume;
        yCenter += ((triangle.y1 + triangle.y2 + triangle.y3) / 4) * currentVolume;
        zCenter += ((triangle.z1 + triangle.z2 + triangle.z3) / 4) * currentVolume;
    }

    cout << endl << "Total Volume = " << totalVolume << endl;
    cout << endl << "X center = " << xCenter / totalVolume << endl;
    cout << endl << "Y center = " << yCenter / totalVolume << endl;
    cout << endl << "Z center = " << zCenter / totalVolume << endl;

    return vec3(xCenter / totalVolume, yCenter / totalVolume, zCenter / totalVolume);
}

//#define MAX_VERTS 100     /* maximum number of polyhedral vertices */
//#define MAX_FACES 100     /* maximum number of polyhedral faces */
//#define MAX_POLYGON_SZ 10 /* maximum number of verts per polygonal face */

#define X 0
#define Y 1
#define Z 2

/*
============================================================================
macros
============================================================================
*/

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

/*
============================================================================
globals
============================================================================
*/

static int A;   /* alpha */
static int B;   /* beta */
static int C;   /* gamma */

/* projection integrals */
static double P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;

/* face integrals */
static double Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca;

/* volume integrals */
static double T0;
vec3 T1, T2, TP;


void compProjectionIntegrals(Polygon3D& face)
{
    double a0, a1, da;
    double b0, b1, db;
    double a0_2, a0_3, a0_4, b0_2, b0_3, b0_4;
    double a1_2, a1_3, b1_2, b1_3;
    double C1, Ca, Caa, Caaa, Cb, Cbb, Cbbb;
    double Cab, Kab, Caab, Kaab, Cabb, Kabb;

    P1 = Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0.0;
    for (unsigned i = 0; i < face.Vertices().size(); i++)
     {
        a0 = face.Vertices()[i][A];
        b0 = face.Vertices()[i][B];
        a1 = face.Vertices()[(i + 1) % face.Vertices().size()][A];
        b1 = face.Vertices()[(i + 1) % face.Vertices().size()][B];
        da = a1 - a0;
        db = b1 - b0;
        a0_2 = a0 * a0; a0_3 = a0_2 * a0; a0_4 = a0_3 * a0;
        b0_2 = b0 * b0; b0_3 = b0_2 * b0; b0_4 = b0_3 * b0;
        a1_2 = a1 * a1; a1_3 = a1_2 * a1;
        b1_2 = b1 * b1; b1_3 = b1_2 * b1;

        C1 = a1 + a0;
        Ca = a1*C1 + a0_2; Caa = a1*Ca + a0_3; Caaa = a1*Caa + a0_4;
        Cb = b1*(b1 + b0) + b0_2; Cbb = b1*Cb + b0_3; Cbbb = b1*Cbb + b0_4;
        Cab = 3 * a1_2 + 2 * a1*a0 + a0_2; Kab = a1_2 + 2 * a1*a0 + 3 * a0_2;
        Caab = a0*Cab + 4 * a1_3; Kaab = a1*Kab + 4 * a0_3;
        Cabb = 4 * b1_3 + 3 * b1_2*b0 + 2 * b1*b0_2 + b0_3;
        Kabb = b1_3 + 2 * b1_2*b0 + 3 * b1*b0_2 + 4 * b0_3;

        P1 += db*C1;
        Pa += db*Ca;
        Paa += db*Caa;
        Paaa += db*Caaa;
        Pb += da*Cb;
        Pbb += da*Cbb;
        Pbbb += da*Cbbb;
        Pab += db*(b1*Cab + b0*Kab);
        Paab += db*(b1*Caab + b0*Kaab);
        Pabb += da*(a1*Cabb + a0*Kabb);
    }

    P1 /= 2.0;
    Pa /= 6.0;
    Paa /= 12.0;
    Paaa /= 20.0;
    Pb /= -6.0;
    Pbb /= -12.0;
    Pbbb /= -20.0;
    Pab /= 24.0;
    Paab /= 60.0;
    Pabb /= -60.0;
}

void compFaceIntegrals(Polygon3D& face)
{
    double  w;
    vec3 n;
    double k1, k2, k3, k4;

    compProjectionIntegrals(face);

    w = face.NormalOffset();
    n = face.Normal();
    k1 = 1 / n[C]; k2 = k1 * k1; k3 = k2 * k1; k4 = k3 * k1;

    Fa = k1 * Pa;
    Fb = k1 * Pb;
    Fc = -k2 * (n[A] * Pa + n[B] * Pb + w*P1);

    Faa = k1 * Paa;
    Fbb = k1 * Pbb;
    Fcc = k3 * (SQR(n[A])*Paa + 2 * n[A] * n[B] * Pab + SQR(n[B])*Pbb
        + w*(2 * (n[A] * Pa + n[B] * Pb) + w*P1));

    Faaa = k1 * Paaa;
    Fbbb = k1 * Pbbb;
    Fccc = -k4 * (CUBE(n[A])*Paaa + 3 * SQR(n[A])*n[B] * Paab
        + 3 * n[A] * SQR(n[B])*Pabb + CUBE(n[B])*Pbbb
        + 3 * w*(SQR(n[A])*Paa + 2 * n[A] * n[B] * Pab + SQR(n[B])*Pbb)
        + w*w*(3 * (n[A] * Pa + n[B] * Pb) + w*P1));

    Faab = k1 * Paab;
    Fbbc = -k2 * (n[A] * Pabb + n[B] * Pbbb + w*Pbb);
    Fcca = k3 * (SQR(n[A])*Paaa + 2 * n[A] * n[B] * Paab + SQR(n[B])*Pabb
        + w*(2 * (n[A] * Paa + n[B] * Pab) + w*Pa));
}

void compVolumeIntegrals(const vector<Polygon3D>& p)
{
    double nx, ny, nz;

    T0 = T1[X] = T1[Y] = T1[Z]
        = T2[X] = T2[Y] = T2[Z]
        = TP[X] = TP[Y] = TP[Z] = 0;

    
    /*for (i = 0; i < p->numFaces; i++) {*/
    for each (Polygon3D f in p)
    {
        //f = &p->faces[i];

        nx = fabs(f.Normal().x());
        ny = fabs(f.Normal().y());
        nz = fabs(f.Normal().z());
        if (nx > ny && nx > nz) C = X;
        else C = (ny > nz) ? Y : Z;
        A = (C + 1) % 3;
        B = (A + 1) % 3;

        compFaceIntegrals(f);

        T0 += f.Normal().x() * ((A == X) ? Fa : ((B == X) ? Fb : Fc));

        T1[A] += f.Normal()[A] * Faa;
        T1[B] += f.Normal()[B] * Fbb;
        T1[C] += f.Normal()[C] * Fcc;
        T2[A] += f.Normal()[A] * Faaa;
        T2[B] += f.Normal()[B] * Fbbb;
        T2[C] += f.Normal()[C] * Fccc;
        TP[A] += f.Normal()[A] * Faab;
        TP[B] += f.Normal()[B] * Fbbc;
        TP[C] += f.Normal()[C] * Fcca;
    }

    T1[X] /= 2; T1[Y] /= 2; T1[Z] /= 2;
    T2[X] /= 3; T2[Y] /= 3; T2[Z] /= 3;
    TP[X] /= 2; TP[Y] /= 2; TP[Z] /= 2;
}


vec3 GeometricTools::CalcCOM_3(const vector<Polygon3D>& polygons)
{


    /*******************************************************
    *                                                      *
    *  volInt.c                                            *
    *                                                      *
    *  This code computes volume integrals needed for      *
    *  determining mass properties of polyhedral bodies.   *
    *                                                      *
    *  For more information, see the accompanying README   *
    *  file, and the paper                                 *
    *                                                      *
    *  Brian Mirtich, "Fast and Accurate Computation of    *
    *  Polyhedral Mass Properties," journal of graphics    *
    *  tools, volume 1, number 1, 1996.                    *
    *                                                      *
    *  This source code is public domain, and may be used  *
    *  in any way, shape or form, free of charge.          *
    *                                                      *
    *  Copyright 1995 by Brian Mirtich                     *
    *                                                      *
    *  mirtich@cs.berkeley.edu                             *
    *  http://www.cs.berkeley.edu/~jfc/mirtich/            *
    *                                                      *
    *******************************************************/

    /*
    Revision history

    26 Jan 1996	Program creation.

    3 Aug 1996	Corrected bug arising when polyhedron density
    is not 1.0.  Changes confined to function main().
    Thanks to Zoran Popovic for catching this one.

    27 May 1997     Corrected sign error in translation of inertia
    product terms to center of mass frame.  Changes
    confined to function main().  Thanks to
    Chris Hecker.
    */


    double density, mass;
    vec3 r;            /* center of mass */
    vector<vec3> J;
    J.push_back(vec3());
    J.push_back(vec3());
    J.push_back(vec3());


    compVolumeIntegrals(polygons);

    density = 1.0;  /* assume unit density */

    mass = density * T0;

    /* compute center of mass */
    r[X] = T1[X] / T0;
    r[Y] = T1[Y] / T0;
    r[Z] = T1[Z] / T0;

    /* compute inertia tensor */
    J[X][X] = density * (T2[Y] + T2[Z]);
    J[Y][Y] = density * (T2[Z] + T2[X]);
    J[Z][Z] = density * (T2[X] + T2[Y]);
    J[X][Y] = J[Y][X] = -density * TP[X];
    J[Y][Z] = J[Z][Y] = -density * TP[Y];
    J[Z][X] = J[X][Z] = -density * TP[Z];

    /* translate inertia tensor to center of mass */
    J[X][X] -= mass * (r[Y] * r[Y] + r[Z] * r[Z]);
    J[Y][Y] -= mass * (r[Z] * r[Z] + r[X] * r[X]);
    J[Z][Z] -= mass * (r[X] * r[X] + r[Y] * r[Y]);
    J[X][Y] = J[Y][X] += mass * r[X] * r[Y];
    J[Y][Z] = J[Z][Y] += mass * r[Y] * r[Z];
    J[Z][X] = J[X][Z] += mass * r[Z] * r[X];

   //printf("center of mass:  (%+12.6f,%+12.6f,%+12.6f)\n\n", r[X], r[Y], r[Z]);
   //
   //printf("inertia tensor with origin at c.o.m. :\n");
   //printf("%+15.6f  %+15.6f  %+15.6f\n", J[X][X], J[X][Y], J[X][Z]);
   //printf("%+15.6f  %+15.6f  %+15.6f\n", J[Y][X], J[Y][Y], J[Y][Z]);
   //printf("%+15.6f  %+15.6f  %+15.6f\n\n", J[Z][X], J[Z][Y], J[Z][Z]);
    return r;
}