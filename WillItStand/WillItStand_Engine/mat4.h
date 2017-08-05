#ifndef matr4_h__
#define matr4_h__

#include <iostream>
#include <assert.h>

class mat4
{
    double mat[4][4];

public:
    mat4();
    mat4(double onDiag);
    mat4(const mat4& m);
    mat4(double a00, double a01, double a02, double a03,
        double a10, double a11, double a12, double a13,
        double a20, double a21, double a22, double a23,
        double a30, double a31, double a32, double a33);
    double at(int row, int col) const;
    void putAt(int row, int col, double val);
    mat4 operator-() const;
    mat4& operator+=(const mat4& m);
    mat4& operator-=(const mat4& m);
    mat4& operator*=(const mat4& m);
    double sumRowCol(int m1Row, int m2Col, const mat4& m1, const mat4& m2);
    mat4 transpose() const;
    double determinant() const;
    mat4 inverse() const;
    mat4& scale(const double scalar);

    //friend bool operator==(const mat4& m1, const mat4& m2);
    //friend std::ostream& operator<<(const std::ostream& os, const mat4& m);

};

#endif // matr4_h__
