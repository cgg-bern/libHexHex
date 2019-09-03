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


#include "HexExtractor.hh"

#include "Utils.hh"


namespace HexEx
{

Position HexExtractor::getMergePosition(unsigned int equivalenceClass, bool optimizedMergePosition)
{
    assert(!equivalenceClasses[equivalenceClass].empty());

    auto& vertices = equivalenceClasses[equivalenceClass];

    if (vertices.size() == 1)
        return intermediateHexMesh.vertex(VertexHandle(vertices.front()));

    bool containsBoundaryVertex = false;

    for (auto vh : vertices)
        if (isBoundaryHexVertex(VertexHandle(vh)))
        {
            containsBoundaryVertex = true;
            break;
        }


    if (!containsBoundaryVertex || !optimizedMergePosition)
    {
        // all vertices inside
        // use center of gravity as merge position
        auto pos = Position(0,0,0);
        for (auto vh : vertices)
            pos += intermediateHexMesh.vertex(VertexHandle(vh));
        pos /= vertices.size();
        return pos;
    }
    else
    {
        return getComplicatedMergePosition(equivalenceClass);

    }

}

Position HexExtractor::getComplicatedMergePosition(unsigned int equivalenceClass)
{
    // only use boundary vertices

    auto vertices = equivalenceClasses[equivalenceClass];

    auto isBoundaryHV = [&](int vh){ return isBoundaryHexVertex(VertexHandle(vh)); };
    vertices.erase(std::partition(vertices.begin(), vertices.end(), isBoundaryHV), vertices.end());

    auto cog = Position(0,0,0);
    for (auto vh : vertices)
        cog += intermediateHexMesh.vertex(VertexHandle(vh));
    cog /= vertices.size();

    auto sixteenzeros = std::vector<double>(16, 0.0);
    auto quadric = Matrix4x4d(sixteenzeros.data());
//    for (unsigned int i = 0; i < 4; i++)
//        quadric(i,i) = 1;
    for (auto hexVh : vertices)
    {

        auto boundaryHalffaces = getBoundaryHalffacesOfHexVertex(VertexHandle(hexVh));

        for (auto hfh : boundaryHalffaces)
        {
            Vec3d n = getNormal(hfh);
//            auto a = getArea(hfh);
            auto bla = (intermediateHexMesh.vertex(VertexHandle(hexVh))-cog);
            auto d = -1.0* (n | bla);

            auto tmpquadric = getQuadric(n, d);

            quadric += tmpquadric;//*a;
        }
    }


    double minDist = 100000;
    auto bestVh = VertexHandle();

    for (auto vh : vertices)
    {
        auto pos = intermediateHexMesh.vertex(VertexHandle(vh)) - cog;
        //auto pos2 = Vec4d(pos[0], pos[1], pos[2], 1);
        //auto dist = (quadric * pos2) | pos2;
        auto pos3 = quadric.transform_point(pos);
        auto dist = pos3 | pos;
        if (dist < minDist)
        {
            minDist = dist;
            bestVh = VertexHandle(vh);
        }

    }

    return intermediateHexMesh.vertex(bestVh);
}

Matrix4x4d HexExtractor::getQuadric(Vec3d n, double d)
{
    auto sixteenzeros = std::vector<double>(16, 0.0);
    auto quadric = Matrix4x4d(sixteenzeros.data());
    for (auto i = 0u; i < 3; ++i)
        for (auto j = 0u; j < 3; ++j)
            quadric(i,j) = n[i]*n[j];
    for (auto i = 0u; i < 3; ++i)
    {
        quadric(3,i) = n[i] * d;
        quadric(i,3) = n[i] * d;
    }
    quadric(3,3) = d*d;

    return quadric;
}

std::vector<HalfFaceHandle> HexExtractor::getBoundaryHalffacesOfHexVertex(VertexHandle hexVh)
{

    auto halffaces = std::vector<HalfFaceHandle>();

    auto vertexType = vertexTypes[hexVh];
    if (vertexType == FHVertex)
    {
        auto fh = FaceHandle(incidentElementId[hexVh]);
        halffaces.push_back(getBoundaryHalfface(inputMesh, fh));
    }
    else if (vertexType == EHVertex)
    {
        auto eh = EdgeHandle(incidentElementId[hexVh]);
        auto heh = inputMesh.halfedge_handle(eh, 0);
        halffaces = getBoundaryHalffaces(inputMesh, heh);
    }
    else if (vertexType == VHVertex)
    {
        auto vh = VertexHandle(incidentElementId[hexVh]);
        halffaces = getBoundaryHalffaces(inputMesh, vh);
    }

    return halffaces;
}

Vec3d HexExtractor::getNormal(HalfFaceHandle hfh)
{

    auto vertices = inputMesh.get_halfface_vertices(hfh);

    auto p0 = inputPosition(vertices[0]);
    auto p1 = inputPosition(vertices[1]);
    auto p2 = inputPosition(vertices[2]);

    auto d1 = p1 - p0;
    auto d2 = p2 - p0;
    auto n = d1 % d2;

    assert(n != Vec3d(0,0,0));

    n.normalize();

    return n;
}

double HexExtractor::getArea(HalfFaceHandle hfh)
{
  auto vertices = inputMesh.get_halfface_vertices(hfh);

   auto p0 = inputPosition(vertices[0]);
   auto p1 = inputPosition(vertices[1]);
   auto p2 = inputPosition(vertices[2]);

   auto d1 = p1 - p0;
   auto d2 = p2 - p0;
   auto n = d1 % d2;

   assert(n != Vec3d(0,0,0));

   return n.length();
}




bool HexExtractor::isBoundaryHexVertex(VertexHandle hexVh)
{
    auto vertexType = vertexTypes[hexVh];
    if (vertexType == FHVertex)
    {
        auto fh = FaceHandle(incidentElementId[hexVh]);
        if (inputMesh.is_boundary(fh))
            return true;
    }
    else if (vertexType == EHVertex)
    {
        auto eh = EdgeHandle(incidentElementId[hexVh]);
        if (inputMesh.is_boundary(eh))
            return true;
    }
    else if (vertexType == VHVertex)
    {
        auto vh = VertexHandle(incidentElementId[hexVh]);
        if (inputMesh.is_boundary(vh))
            return true;
    }

    return false;
}

}
