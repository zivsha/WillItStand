#include "mat4.h"
mat4::mat4()
{
    for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
        mat[i][j] = 0.0;
}

mat4::mat4(double onDiag)
{
    for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
        mat[i][j] = (i == j) ? onDiag : 0.0;
}

mat4::mat4(const mat4& m)
{
    for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
        mat[i][j] = m.at(i, j);
}
mat4::mat4(double a00, double a01, double a02, double a03,
    double a10, double a11, double a12, double a13,
    double a20, double a21, double a22, double a23,
    double a30, double a31, double a32, double a33)
{
    mat[0][0] = a00;
    mat[0][1] = a01;
    mat[0][2] = a02;
    mat[0][3] = a03;
    mat[1][0] = a10;
    mat[1][1] = a11;
    mat[1][2] = a12;
    mat[1][3] = a13;
    mat[2][0] = a20;
    mat[2][1] = a21;
    mat[2][2] = a22;
    mat[2][3] = a23;
    mat[3][0] = a30;
    mat[3][1] = a31;
    mat[3][2] = a32;
    mat[3][3] = a33;
}

double mat4::at(int row, int col) const
{
    if (row<0 || row>3 || col<0 || col>3)
    {
        throw new std::exception();
    }
    return mat[row][col];
}

void mat4::putAt(int row, int col, double val)
{
    if (row<0 || row>3 || col<0 || col>3)
    {
        throw new std::exception();
    }
    mat[row][col] = val;
}

mat4 mat4::operator-() const
{
    mat4 m;
    for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
        m.putAt(i, j, -(mat[i][j]));
    return m;
}

mat4& mat4::operator+=(const mat4& m)
{
    for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
        mat[i][j] += m.at(i, j);
    return *this;
}

mat4& mat4::operator-=(const mat4& m)
{
    return *this += -m;
}

mat4& mat4::operator*=(const mat4& m)
{
    double sum = 0;
    mat4 tempMat;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            tempMat.putAt(i, j, sumRowCol(i, j, *this, m));
        }
    }
    *this = tempMat;
    return *this;
}

double mat4::sumRowCol(int m1Row, int m2Col, const mat4& m1, const mat4& m2)
{
    double sum = 0;
    for (int i = 0; i < 4; i++)
    {
        sum += (m1.at(m1Row, i) * m2.at(i, m2Col));
    }
    return sum;
}

mat4 mat4::transpose() const
{
    mat4 retMat;
    for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
        retMat.putAt(j, i, mat[i][j]);
    return retMat;
}

double mat4::determinant() const
{
    return (mat[3][0] * mat[2][1] * mat[1][2] * mat[0][3] - mat[2][0] * mat[3][1] * mat[1][2] * mat[0][3] -
        mat[3][0] * mat[1][0] * mat[2][2] * mat[0][3] + mat[1][0] * mat[3][1] * mat[2][2] * mat[0][3] +
        mat[2][0] * mat[1][0] * mat[3][2] * mat[0][3] - mat[1][0] * mat[2][1] * mat[3][2] * mat[0][3] -
        mat[3][0] * mat[2][1] * mat[0][2] * mat[1][3] + mat[2][0] * mat[3][1] * mat[0][2] * mat[1][3] +
        mat[3][0] * mat[0][1] * mat[2][2] * mat[1][3] - mat[0][0] * mat[3][1] * mat[2][2] * mat[1][3] -
        mat[2][0] * mat[0][1] * mat[3][2] * mat[1][3] + mat[0][0] * mat[2][1] * mat[3][2] * mat[1][3] +
        mat[3][0] * mat[1][0] * mat[0][2] * mat[2][3] - mat[1][0] * mat[3][1] * mat[0][2] * mat[2][3] -
        mat[3][0] * mat[0][1] * mat[1][2] * mat[2][3] + mat[0][0] * mat[3][1] * mat[1][2] * mat[2][3] +
        mat[1][0] * mat[0][1] * mat[3][2] * mat[2][3] - mat[0][0] * mat[1][0] * mat[3][2] * mat[2][3] -
        mat[2][0] * mat[1][0] * mat[0][2] * mat[3][3] + mat[1][0] * mat[2][1] * mat[0][2] * mat[3][3] +
        mat[2][0] * mat[0][1] * mat[1][2] * mat[3][3] - mat[0][0] * mat[2][1] * mat[1][2] * mat[3][3] -
        mat[1][0] * mat[0][1] * mat[2][2] * mat[3][3] + mat[0][0] * mat[1][0] * mat[2][2] * mat[3][3]);
}

mat4 mat4::inverse() const {
    mat4 m;
    m.putAt(0, 0, (mat[1][2] * mat[2][3] * mat[3][1] - mat[1][3] * mat[2][2] * mat[3][1] + mat[1][3] * mat[2][1] * mat[3][2] - mat[1][1] * mat[2][3] * mat[3][2] - mat[1][2] * mat[2][1] * mat[3][3] + mat[1][1] * mat[2][2] * mat[3][3]));
    m.putAt(0, 1, (mat[0][3] * mat[2][2] * mat[3][1] - mat[0][2] * mat[2][3] * mat[3][1] - mat[0][3] * mat[2][1] * mat[3][2] + mat[0][1] * mat[2][3] * mat[3][2] + mat[0][2] * mat[2][1] * mat[3][3] - mat[0][1] * mat[2][2] * mat[3][3]));
    m.putAt(0, 2, (mat[0][2] * mat[1][3] * mat[3][1] - mat[0][3] * mat[1][2] * mat[3][1] + mat[0][3] * mat[1][1] * mat[3][2] - mat[0][1] * mat[1][3] * mat[3][2] - mat[0][2] * mat[1][1] * mat[3][3] + mat[0][1] * mat[1][2] * mat[3][3]));
    m.putAt(0, 3, (mat[0][3] * mat[1][2] * mat[2][1] - mat[0][2] * mat[1][3] * mat[2][1] - mat[0][3] * mat[1][1] * mat[2][2] + mat[0][1] * mat[1][3] * mat[2][2] + mat[0][2] * mat[1][1] * mat[2][3] - mat[0][1] * mat[1][2] * mat[2][3]));
    m.putAt(1, 0, (mat[1][3] * mat[2][2] * mat[3][0] - mat[1][2] * mat[2][3] * mat[3][0] - mat[1][3] * mat[2][0] * mat[3][2] + mat[1][0] * mat[2][3] * mat[3][2] + mat[1][2] * mat[2][0] * mat[3][3] - mat[1][0] * mat[2][2] * mat[3][3]));
    m.putAt(1, 1, (mat[0][2] * mat[2][3] * mat[3][0] - mat[0][3] * mat[2][2] * mat[3][0] + mat[0][3] * mat[2][0] * mat[3][2] - mat[0][0] * mat[2][3] * mat[3][2] - mat[0][2] * mat[2][0] * mat[3][3] + mat[0][0] * mat[2][2] * mat[3][3]));
    m.putAt(1, 2, (mat[0][3] * mat[1][2] * mat[3][0] - mat[0][2] * mat[1][3] * mat[3][0] - mat[0][3] * mat[1][0] * mat[3][2] + mat[0][0] * mat[1][3] * mat[3][2] + mat[0][2] * mat[1][0] * mat[3][3] - mat[0][0] * mat[1][2] * mat[3][3]));
    m.putAt(1, 3, (mat[0][2] * mat[1][3] * mat[2][0] - mat[0][3] * mat[1][2] * mat[2][0] + mat[0][3] * mat[1][0] * mat[2][2] - mat[0][0] * mat[1][3] * mat[2][2] - mat[0][2] * mat[1][0] * mat[2][3] + mat[0][0] * mat[1][2] * mat[2][3]));
    m.putAt(2, 0, (mat[1][1] * mat[2][3] * mat[3][0] - mat[1][3] * mat[2][1] * mat[3][0] + mat[1][3] * mat[2][0] * mat[3][1] - mat[1][0] * mat[2][3] * mat[3][1] - mat[1][1] * mat[2][0] * mat[3][3] + mat[1][0] * mat[2][1] * mat[3][3]));
    m.putAt(2, 1, (mat[0][3] * mat[2][1] * mat[3][0] - mat[0][1] * mat[2][3] * mat[3][0] - mat[0][3] * mat[2][0] * mat[3][1] + mat[0][0] * mat[2][3] * mat[3][1] + mat[0][1] * mat[2][0] * mat[3][3] - mat[0][0] * mat[2][1] * mat[3][3]));
    m.putAt(2, 2, (mat[0][1] * mat[1][3] * mat[3][0] - mat[0][3] * mat[1][1] * mat[3][0] + mat[0][3] * mat[1][0] * mat[3][1] - mat[0][0] * mat[1][3] * mat[3][1] - mat[0][1] * mat[1][0] * mat[3][3] + mat[0][0] * mat[1][1] * mat[3][3]));
    m.putAt(2, 3, (mat[0][3] * mat[1][1] * mat[2][0] - mat[0][1] * mat[1][3] * mat[2][0] - mat[0][3] * mat[1][0] * mat[2][1] + mat[0][0] * mat[1][3] * mat[2][1] + mat[0][1] * mat[1][0] * mat[2][3] - mat[0][0] * mat[1][1] * mat[2][3]));
    m.putAt(3, 0, (mat[1][2] * mat[2][1] * mat[3][0] - mat[1][1] * mat[2][2] * mat[3][0] - mat[1][2] * mat[2][0] * mat[3][1] + mat[1][0] * mat[2][2] * mat[3][1] + mat[1][1] * mat[2][0] * mat[3][2] - mat[1][0] * mat[2][1] * mat[3][2]));
    m.putAt(3, 1, (mat[0][1] * mat[2][2] * mat[3][0] - mat[0][2] * mat[2][1] * mat[3][0] + mat[0][2] * mat[2][0] * mat[3][1] - mat[0][0] * mat[2][2] * mat[3][1] - mat[0][1] * mat[2][0] * mat[3][2] + mat[0][0] * mat[2][1] * mat[3][2]));
    m.putAt(3, 2, (mat[0][2] * mat[1][1] * mat[3][0] - mat[0][1] * mat[1][2] * mat[3][0] - mat[0][2] * mat[1][0] * mat[3][1] + mat[0][0] * mat[1][2] * mat[3][1] + mat[0][1] * mat[1][0] * mat[3][2] - mat[0][0] * mat[1][1] * mat[3][2]));
    m.putAt(3, 3, (mat[0][1] * mat[1][2] * mat[2][0] - mat[0][2] * mat[1][1] * mat[2][0] + mat[0][2] * mat[1][0] * mat[2][1] - mat[0][0] * mat[1][2] * mat[2][1] - mat[0][1] * mat[1][0] * mat[2][2] + mat[0][0] * mat[1][1] * mat[2][2]));
    double det = this->determinant();
    if (det == 0)
    {
        assert(0);
    }
    else
    {
        return m.scale(1 / det);
    }
    return m.scale(1 / det);
}
mat4& mat4::scale(const double scalar)
{
    for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
        mat[i][j] *= scalar;
    return *this;
}
