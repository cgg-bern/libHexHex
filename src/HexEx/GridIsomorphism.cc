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


#include <vector>

#include "GridIsomorphism.hh"
#include "Typedefs.hh"

namespace HexEx
{

std::array<char,24*24> RestrictedRotation::mulitplicationMap;
std::array<Matrix4x4d,24> RestrictedRotation::matrices;
std::array<char, 24> RestrictedRotation::inverseMap;
RestrictedRotation::MultiplicationMapInitializer RestrictedRotation::initializer;

RestrictedRotation::MultiplicationMapInitializer::MultiplicationMapInitializer()
{

    for (auto i = 0; i < 24; ++i)
        matrices[i] = convertToMatrix(i);


    for (auto i = 0; i < 24; ++i)
        for (auto j = 0; j < 24; ++j)
        {
            auto lhs = RestrictedRotation(i);
            auto rhs = RestrictedRotation(j);
            auto res = RestrictedRotation(lhs.toMatrix()*rhs.toMatrix());
            mulitplicationMap[24*i+j] = res.i;
        }

    inverseMap[ 0] =  0;
    inverseMap[ 1] =  1;
    inverseMap[ 2] =  3;
    inverseMap[ 3] =  2;
    inverseMap[ 4] =  4;
    inverseMap[ 5] =  5;
    inverseMap[ 6] =  6;
    inverseMap[ 7] =  7;
    inverseMap[ 8] = 16;
    inverseMap[ 9] = 20;
    inverseMap[10] = 10;
    inverseMap[11] = 14;
    inverseMap[12] = 21;
    inverseMap[13] = 17;
    inverseMap[14] = 11;
    inverseMap[15] = 15;
    inverseMap[16] =  8;
    inverseMap[17] = 13;
    inverseMap[18] = 22;
    inverseMap[19] = 19;
    inverseMap[20] =  9;
    inverseMap[21] = 12;
    inverseMap[22] = 18;
    inverseMap[23] = 23;

}

RestrictedRotation::RestrictedRotation(char i)
    :
      i(i)
{

}

RestrictedRotation::RestrictedRotation(const Matrix4x4d m)
{
    auto p = 0;
    if (m(1,0) != 0)
        p = 1;
    else if (m(2,0) != 0)
        p = 2;

    auto nb = m(p,0) < 0 ? 1 : 0;

    auto p2 = (p+1)%3;
    auto p2b = 0;

    if (m(p2,1) == 0)
    {
        p2 = ((p-1)+3)%3;
        p2b = 1;
    }

    auto n2b = 0;

    if (m(p2,1) < 0)
        n2b = 1;

    i = (p << 3) + (nb << 2) + (p2b << 1) + n2b;
}

RestrictedRotation RestrictedRotation::operator*(const RestrictedRotation& rhs) const
{
    return RestrictedRotation(mulitplicationMap[24*i+rhs.i]);
}


Vec3d RestrictedRotation::transform(const Vec3d& v) const
{
    return transformT(i, v);
}

Vec3i RestrictedRotation::transform(const Vec3i& v) const
{
    return transformT(i, v);
}

RestrictedRotation RestrictedRotation::inverted() const
{
    auto res = RestrictedRotation(i);
    res.invert();
    return res;
}

void RestrictedRotation::invert()
{
    i = inverseMap[i];

//    switch(i)
//    {
//    case 0:  i =  0; break;
//    case 1:  i =  1; break;
//    case 2:  i =  3; break;
//    case 3:  i =  2; break;
//    case 4:  i =  4; break;
//    case 5:  i =  5; break;
//    case 6:  i =  6; break;
//    case 7:  i =  7; break;
//    case 8:  i = 16; break;
//    case 9:  i = 20; break;
//    case 10: i = 10; break;
//    case 11: i = 14; break;
//    case 12: i = 21; break;
//    case 13: i = 17; break;
//    case 14: i = 11; break;
//    case 15: i = 15; break;
//    case 16: i =  8; break;
//    case 17: i = 13; break;
//    case 18: i = 22; break;
//    case 19: i = 19; break;
//    case 20: i =  9; break;
//    case 21: i = 12; break;
//    case 22: i = 18; break;
//    case 23: i = 23; break;
//    default: assert(false);
//    }
}

Matrix4x4d RestrictedRotation::toMatrix() const
{
    return matrices[i];
}

Matrix4x4d RestrictedRotation::convertToMatrix(char i)
{
    auto m = Matrix4x4d();

    auto vecs = std::vector<Vec3d>(3, Vec3d(0,0,0));

    auto p = i >> 3;
    vecs[0][p] = 1;

    if (i & 1<<2)
        vecs[0] *= -1.0;

    auto p2 = 0;
    if (i & 1<<1)
        p2 = p-1;
    else
        p2 = p+1;
    p2 = (p2+3)%3;

    vecs[1][p2] = 1;
    if (i & 1<<0)
        vecs[1] *= -1.0;

    vecs[2] = vecs[0] % vecs[1];

    for (auto i = 0u; i < 3; ++i)
    {
        for (auto j = 0u; j < 3; ++j)
            m(i,j) = vecs[j][i];
        m(3,i) = 0;
        m(i,3) = 0;
    }
    m(3,3) = 1;

    return m;
}

GridIsomorphism::GridIsomorphism(int i, const Translation& t)
    :
      r(i),
      t(t)
{}

GridIsomorphism::GridIsomorphism(RestrictedRotation r, const Translation& t)
    :
      r(r),
      t(t)
{}

GridIsomorphism::GridIsomorphism(const Matrix4x4d& m)
    :
      r(m),
      t(Vec3i(0,0,0))
{
    for (auto i = 0u; i < 3; ++i)
        t[i] = (int)round(m(i,3));
}

GridIsomorphism GridIsomorphism::operator*(const GridIsomorphism& rhs)
{    
    auto res = GridIsomorphism(0);

    res.r = r * rhs.r;
    res.t = r.transform(rhs.t) + t;

    return res;
}

Vec3d GridIsomorphism::transform_point(const Vec3d& v) const
{
    return r.transform(v) + Vec3d(t[0], t[1], t[2]);
}

Vec3d GridIsomorphism::transform_vector(const Vec3d& v) const
{
    return r.transform(v);
}

GridIsomorphism GridIsomorphism::inverted() const
{
    auto res = *this;
    res.invert();
    return res;
}

void GridIsomorphism::invert()
{
    r.invert();
    t = r.transform(-t);
}

Matrix4x4d GridIsomorphism::toMatrix() const
{
    auto res = r.toMatrix();
    for (auto i = 0; i < 3; ++i)
        res(i,3) = t[i];
    return res;
}

std::ostream& operator<<(std::ostream& os, const GridIsomorphism& gi)
{
    os << gi.toMatrix() << std::endl;
    return os;
}


}
