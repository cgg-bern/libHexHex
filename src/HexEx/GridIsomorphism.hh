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
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include "Matrix4x4T.hh"

namespace HexEx
{

using Vec3d = OpenVolumeMesh::Geometry::Vec3d;
using Vec3i = OpenVolumeMesh::Geometry::Vec3i;
using Matrix4x4d = Matrix4x4T<double>;

class RestrictedRotation
{
private:
    class MultiplicationMapInitializer
    {
    public:
        MultiplicationMapInitializer();
    };

public:
    explicit RestrictedRotation(char i);
    explicit RestrictedRotation(const Matrix4x4d m);

    bool operator==(const RestrictedRotation& other) const { return i == other.i; }
    bool operator!=(const RestrictedRotation& other) const { return i != other.i; }

    RestrictedRotation operator*(const RestrictedRotation& rhs) const;


    Vec3d transform(const Vec3d& v) const;
    Vec3i transform(const Vec3i& v) const;

    RestrictedRotation inverted() const;
    void invert();

    Matrix4x4d toMatrix() const;

private:


    template <typename Vec>
    Vec transformT(char i, const Vec& v) const
    {
//        return RestrictedRotation::matrices[i].transform_point(v);

        switch(i)
        {
        case 0:  return Vec( v[0],  v[1],  v[2]);
        case 1:  return Vec( v[0], -v[1], -v[2]);
        case 2:  return Vec( v[0], -v[2],  v[1]);
        case 3:  return Vec( v[0],  v[2], -v[1]);
        case 4:  return Vec(-v[0],  v[1], -v[2]);
        case 5:  return Vec(-v[0], -v[1],  v[2]);
        case 6:  return Vec(-v[0],  v[2],  v[1]);
        case 7:  return Vec(-v[0], -v[2], -v[1]);
        case 8:  return Vec( v[2],  v[0],  v[1]);
        case 9:  return Vec(-v[2],  v[0], -v[1]);
        case 10: return Vec( v[1],  v[0], -v[2]);
        case 11: return Vec(-v[1],  v[0],  v[2]);
        case 12: return Vec(-v[2], -v[0],  v[1]);
        case 13: return Vec( v[2], -v[0], -v[1]);
        case 14: return Vec( v[1], -v[0],  v[2]);
        case 15: return Vec(-v[1], -v[0], -v[2]);
        case 16: return Vec( v[1],  v[2],  v[0]);
        case 17: return Vec(-v[1], -v[2],  v[0]);
        case 18: return Vec(-v[2],  v[1],  v[0]);
        case 19: return Vec( v[2], -v[1],  v[0]);
        case 20: return Vec( v[1], -v[2], -v[0]);
        case 21: return Vec(-v[1],  v[2], -v[0]);
        case 22: return Vec( v[2],  v[1], -v[0]);
        case 23: return Vec(-v[2], -v[1], -v[0]);
        default: assert(false); return v;
        }
    }

    char i;

    static Matrix4x4d convertToMatrix(char i);
    static std::array<char, 24*24> mulitplicationMap;
    static std::array<Matrix4x4d, 24> matrices;
    static std::array<char, 24> inverseMap;
    static MultiplicationMapInitializer initializer;

};

class GridIsomorphism
{
public:
    using Translation = Vec3d;

    explicit GridIsomorphism(int i = 0, const Translation& t = Translation(0,0,0));
    explicit GridIsomorphism(RestrictedRotation r, const Translation& t = Translation(0,0,0));
    explicit GridIsomorphism(const Matrix4x4d& m);

    bool operator==(const GridIsomorphism& other) const { return r == other.r && t == other.t; }
    bool operator!=(const GridIsomorphism& other) const { return !operator==(other); }

    GridIsomorphism operator* (const GridIsomorphism& rhs);

    Vec3d transform_point(const Vec3d& v) const;
    Vec3d transform_vector(const Vec3d& v) const;

    void setTranslation(const Translation& translation) { t = translation; }

    GridIsomorphism inverted() const;
    void invert();

    Matrix4x4d toMatrix() const;

private:

    RestrictedRotation r;
    Translation t;
};


std::ostream& operator<<(std::ostream& os, const GridIsomorphism& gi);

}
