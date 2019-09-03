/*
 * Copyright 2016 Computer Graphics Group, RWTH Aachen University
 * Author: Max Lyon <lyon@cs.rwth-aachen.de>
 *
 * This file is part of HexEx.
 *
 * HexEx is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * HexEx is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HexEx.  If not, see <http://www.gnu.org/licenses/>.
 */


#pragma once

#include <array>
#include <ostream>

#include <OpenVolumeMesh/Geometry/VectorT.hh>

namespace HexEx{

template <typename T>
class Matrix4x4T;

template <typename T>
bool operator==(const Matrix4x4T<T>& _mat, const Matrix4x4T<T>& _other);

template <typename T>
class Matrix4x4T
{
    using Matrix4x4 = Matrix4x4T<T>;

    friend bool operator==<>(const Matrix4x4T<T>& _mat, const Matrix4x4T<T>& _other);

public:
    Matrix4x4T()
    {
        for (unsigned int i = 0; i < 4; ++i)
            for (unsigned int j = 0; j < 4; ++j)
                operator()(i,j) = (i == j) ? 1 : 0;
    }

    Matrix4x4T(T* _entries)
    {
        for (unsigned int i = 0; i < 16; ++i)
            entries[i] = *_entries++;
    }

    void clear()
    {
        for (unsigned int i = 0; i < 16; ++i)
            entries[i] = 0;
    }

    T& operator()(unsigned int row, unsigned int col)             { return  entries[row+col*4]; }
    const T& operator()(unsigned int row, unsigned int col) const { return  entries[row+col*4]; }

    void operator+=(const Matrix4x4& _other)
    {
        for (unsigned int i = 0; i < 16; ++i)
            entries[i] += _other.entries[i];
    }

    template <typename Vec>
    Vec operator*(const Vec& _vec)
    {
        auto res = Vec(0,0,0,0);

        for (unsigned int i = 0; i < 4; ++i)
            for (unsigned int j = 0; j < 4; ++j)
                res[i] += operator()(i, j) * _vec[j];

        return res;
    }

    Matrix4x4 operator*(const Matrix4x4& _mat)
    {
        auto res = Matrix4x4();
        res.clear();

        for (unsigned int i = 0; i < 4; ++i)
            for (unsigned int j = 0; j < 4; ++j)
                for (unsigned int k = 0; k < 4; ++k)
                    res(i,j) += operator()(i,k) * _mat(k,j);

        return res;
    }

    bool operator==(const Matrix4x4& _other)
    {
      return entries == _other.entries;
    }

    template <typename Vec>
    Vec transform_point(const Vec& _vec)
    {
        auto res = Vec(0,0,0);
        auto vec = OpenVolumeMesh::Geometry::Vec4d(_vec[0],_vec[1],_vec[2],1);
        for (unsigned int i = 0; i < 3; ++i)
            for (unsigned int j = 0; j < 4; ++j)
                res[i] += operator()(i,j) * vec[j];

        return res;
    }

    template <typename Vec>
    Vec transform_vector(const Vec& _vec)
    {
        auto res = Vec(0,0,0);
        for (unsigned int i = 0; i < 3; ++i)
            for (unsigned int j = 0; j < 3; ++j)
                res[i] += operator()(i,j) * _vec[j];

        return res;
    }

    /*
     * Compute inverse of 4x4 transformation matrix.
     * Taken from OpenFlipper.
     * Originally taken from Mesa3.1
     * Code contributed by Jacques Leroy jle@star.be */
    bool invert()
    {
        auto& M = *this;
#define SWAP_ROWS(a, b) { T *_tmp = a; (a)=(b); (b)=_tmp; }

        T wtmp[4][8];
        T m0, m1, m2, m3, s;
        T *r0, *r1, *r2, *r3;

        r0 = wtmp[0], r1 = wtmp[1], r2 = wtmp[2], r3 = wtmp[3];

        r0[0] = M(0,0); r0[1] = M(0,1);
        r0[2] = M(0,2); r0[3] = M(0,3);
        r0[4] = 1.0, r0[5] = r0[6] = r0[7] = 0.0;

        r1[0] = M(1,0); r1[1] = M(1,1);
        r1[2] = M(1,2); r1[3] = M(1,3);
        r1[5] = 1.0, r1[4] = r1[6] = r1[7] = 0.0;

        r2[0] = M(2,0); r2[1] = M(2,1);
        r2[2] = M(2,2); r2[3] = M(2,3);
        r2[6] = 1.0, r2[4] = r2[5] = r2[7] = 0.0;

        r3[0] = M(3,0); r3[1] = M(3,1);
        r3[2] = M(3,2); r3[3] = M(3,3);
        r3[7] = 1.0, r3[4] = r3[5] = r3[6] = 0.0;


        /* choose pivot - or die */
        if (fabs(r3[0])>fabs(r2[0])) SWAP_ROWS(r3, r2);
        if (fabs(r2[0])>fabs(r1[0])) SWAP_ROWS(r2, r1);
        if (fabs(r1[0])>fabs(r0[0])) SWAP_ROWS(r1, r0);
        if (0.0 == r0[0])  return false;


        /* eliminate first variable     */
        m1 = r1[0]/r0[0]; m2 = r2[0]/r0[0]; m3 = r3[0]/r0[0];
        s = r0[1]; r1[1] -= m1 * s; r2[1] -= m2 * s; r3[1] -= m3 * s;
        s = r0[2]; r1[2] -= m1 * s; r2[2] -= m2 * s; r3[2] -= m3 * s;
        s = r0[3]; r1[3] -= m1 * s; r2[3] -= m2 * s; r3[3] -= m3 * s;
        s = r0[4];
        if (s != 0.0) { r1[4] -= m1 * s; r2[4] -= m2 * s; r3[4] -= m3 * s; }
        s = r0[5];
        if (s != 0.0) { r1[5] -= m1 * s; r2[5] -= m2 * s; r3[5] -= m3 * s; }
        s = r0[6];
        if (s != 0.0) { r1[6] -= m1 * s; r2[6] -= m2 * s; r3[6] -= m3 * s; }
        s = r0[7];
        if (s != 0.0) { r1[7] -= m1 * s; r2[7] -= m2 * s; r3[7] -= m3 * s; }


        /* choose pivot - or die */
        if (fabs(r3[1])>fabs(r2[1])) SWAP_ROWS(r3, r2);
        if (fabs(r2[1])>fabs(r1[1])) SWAP_ROWS(r2, r1);
        if (0.0 == r1[1])  return false;


        /* eliminate second variable */
        m2 = r2[1]/r1[1]; m3 = r3[1]/r1[1];
        r2[2] -= m2 * r1[2]; r3[2] -= m3 * r1[2];
        r2[3] -= m2 * r1[3]; r3[3] -= m3 * r1[3];
        s = r1[4]; if (0.0 != s) { r2[4] -= m2 * s; r3[4] -= m3 * s; }
        s = r1[5]; if (0.0 != s) { r2[5] -= m2 * s; r3[5] -= m3 * s; }
        s = r1[6]; if (0.0 != s) { r2[6] -= m2 * s; r3[6] -= m3 * s; }
        s = r1[7]; if (0.0 != s) { r2[7] -= m2 * s; r3[7] -= m3 * s; }


        /* choose pivot - or die */
        if (fabs(r3[2])>fabs(r2[2])) SWAP_ROWS(r3, r2);
        if (0.0 == r2[2])  return false;

        /* eliminate third variable */
        m3 = r3[2]/r2[2];
        r3[3] -= m3 * r2[3];
        r3[4] -= m3 * r2[4];
        r3[5] -= m3 * r2[5];
        r3[6] -= m3 * r2[6];
        r3[7] -= m3 * r2[7];

        /* last check */
        if (0.0 == r3[3]) return false;

        s = 1.0/r3[3];              /* now back substitute row 3 */
        r3[4] *= s; r3[5] *= s; r3[6] *= s; r3[7] *= s;

        m2 = r2[3];                 /* now back substitute row 2 */
        s  = 1.0/r2[2];
        r2[4] = s * (r2[4] - r3[4] * m2); r2[5] = s * (r2[5] - r3[5] * m2);
        r2[6] = s * (r2[6] - r3[6] * m2); r2[7] = s * (r2[7] - r3[7] * m2);
        m1 = r1[3];
        r1[4] -= r3[4] * m1; r1[5] -= r3[5] * m1;
        r1[6] -= r3[6] * m1; r1[7] -= r3[7] * m1;
        m0 = r0[3];
        r0[4] -= r3[4] * m0; r0[5] -= r3[5] * m0;
        r0[6] -= r3[6] * m0; r0[7] -= r3[7] * m0;

        m1 = r1[2];                 /* now back substitute row 1 */
        s  = 1.0/r1[1];
        r1[4] = s * (r1[4] - r2[4] * m1); r1[5] = s * (r1[5] - r2[5] * m1);
        r1[6] = s * (r1[6] - r2[6] * m1); r1[7] = s * (r1[7] - r2[7] * m1);
        m0 = r0[2];
        r0[4] -= r2[4] * m0; r0[5] -= r2[5] * m0;
        r0[6] -= r2[6] * m0; r0[7] -= r2[7] * m0;

        m0 = r0[1];                 /* now back substitute row 0 */
        s  = 1.0/r0[0];
        r0[4] = s * (r0[4] - r1[4] * m0); r0[5] = s * (r0[5] - r1[5] * m0);
        r0[6] = s * (r0[6] - r1[6] * m0); r0[7] = s * (r0[7] - r1[7] * m0);

        M(0,0) = r0[4]; M(0,1) = r0[5];
        M(0,2) = r0[6]; M(0,3) = r0[7];
        M(1,0) = r1[4]; M(1,1) = r1[5];
        M(1,2) = r1[6]; M(1,3) = r1[7];
        M(2,0) = r2[4]; M(2,1) = r2[5];
        M(2,2) = r2[6]; M(2,3) = r2[7];
        M(3,0) = r3[4]; M(3,1) = r3[5];
        M(3,2) = r3[6]; M(3,3) = r3[7];

        return true;
#undef SWAP_ROWS
    }

    // taken from OpenFlipper
    T determinant()
    {
        return  entries[12] * entries[9] * entries[6] * entries[3] - entries[8] * entries[13] * entries[6] * entries[3] -
                entries[12] * entries[5] * entries[10] * entries[3] + entries[4] * entries[13] * entries[10] * entries[3] +
                entries[8] * entries[5] * entries[14] * entries[3] - entries[4] * entries[9] * entries[14] * entries[3] -
                entries[12] * entries[9] * entries[2] * entries[7] + entries[8] * entries[13] * entries[2] * entries[7] +
                entries[12] * entries[1] * entries[10] * entries[7] - entries[0] * entries[13] * entries[10] * entries[7] -
                entries[8] * entries[1] * entries[14] * entries[7] + entries[0] * entries[9] * entries[14] * entries[7] +
                entries[12] * entries[5] * entries[2] * entries[11] - entries[4] * entries[13] * entries[2] * entries[11] -
                entries[12] * entries[1] * entries[6] * entries[11] + entries[0] * entries[13] * entries[6] * entries[11] +
                entries[4] * entries[1] * entries[14] * entries[11] - entries[0] * entries[5] * entries[14] * entries[11] -
                entries[8] * entries[5] * entries[2] * entries[15] + entries[4] * entries[9] * entries[2] * entries[15] +
                entries[8] * entries[1] * entries[6] * entries[15] - entries[0] * entries[9] * entries[6] * entries[15] -
                entries[4] * entries[1] * entries[10] * entries[15] + entries[0] * entries[5] * entries[10] * entries[15];
    }

private:
    std::array<T, 16> entries;
};

template <typename T>
bool operator==(const Matrix4x4T<T>& _mat, const Matrix4x4T<T>& _other)
{
  return _mat.entries == _other.entries;
}


template <typename T>
std::ostream& operator<<(std::ostream& _os, const Matrix4x4T<T>& _mat)
{
    for (unsigned int i = 0; i < 4; ++i)
    {
        _os << _mat(i,0);
        for (unsigned int j = 1; j < 4; ++j)
            _os << " " << _mat(i,j);
        _os << std::endl;
    }
    return _os;
}



} // namespace HexEx
