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

#include "Typedefs.hh"

namespace HexEx {

template <typename Vec>
Vec3d toVec3d(const Vec& vec)
{
    return Vec3d(vec[0], vec[1], vec[2]);
}

template <typename Vec>
Vec toVec(const Vec3d& vec)
{
    return Vec(vec[0], vec[1], vec[2]);
}

bool isTetCell(const PolyhedralMesh& mesh, OpenVolumeMesh::CellHandle ch);
bool isTetMesh(const PolyhedralMesh& mesh);

bool isFlipped(const TetrahedralMesh& mesh, OpenVolumeMesh::CellHandle ch);
bool containsFlippedTets(const TetrahedralMesh& mesh);

bool isIrregularEdge(bool boundary, int valence);
Vec3d getValenceColor(bool boundary, int valence);

unsigned int getMaxCoord(Vec3d n);

std::vector<OpenVolumeMesh::EdgeHandle> getCellEdges(const TetrahedralMesh& mesh, OpenVolumeMesh::CellHandle ch);

template <typename PolyMeshT>
std::vector<OpenVolumeMesh::VertexHandle> getVertices(PolyMeshT& mesh, OpenVolumeMesh::HalfFaceHandle hfh)
{
    std::vector<VertexHandle> res;
    Face f = mesh.halfface(hfh);
    for (auto heh : f.halfedges())
        res.push_back(mesh.halfedge(heh).from_vertex());
    return res;
}

std::vector<OpenVolumeMesh::HalfFaceHandle> getIncidentHalfFaces(TetrahedralMesh& mesh, OpenVolumeMesh::CellHandle ch, OpenVolumeMesh::EdgeHandle eh);
std::vector<OpenVolumeMesh::HalfFaceHandle> getIncidentHalfFaces(TetrahedralMesh& mesh, OpenVolumeMesh::CellHandle ch, OpenVolumeMesh::VertexHandle vh);
OpenVolumeMesh::HalfFaceHandle getIncidentHalfFaceInCell(TetrahedralMesh& mesh, OpenVolumeMesh::CellHandle ch, OpenVolumeMesh::HalfEdgeHandle heh);
OpenVolumeMesh::HalfFaceHandle getIncidentHalfFaceInCell(TetrahedralMesh& mesh, OpenVolumeMesh::CellHandle ch, OpenVolumeMesh::FaceHandle fh);
std::vector<OpenVolumeMesh::HalfEdgeHandle> getIncidentHalfEdges(TetrahedralMesh& mesh, OpenVolumeMesh::CellHandle ch, OpenVolumeMesh::VertexHandle vh);

std::vector<OpenVolumeMesh::VertexHandle> getNeighborVertices(TetrahedralMesh& mesh, OpenVolumeMesh::VertexHandle vh);

OpenVolumeMesh::HalfEdgeHandle getHalfEdgeInHalfFace(TetrahedralMesh& mesh, OpenVolumeMesh::EdgeHandle eh, OpenVolumeMesh::HalfFaceHandle hfh);

template <typename MeshT>
OpenVolumeMesh::HalfFaceHandle getBoundaryHalfface(const MeshT& mesh, OpenVolumeMesh::FaceHandle fh)
{
    auto hfh = mesh.halfface_handle(fh, 0);
    if (mesh.is_boundary(hfh))
        return hfh;
    else
        return mesh.opposite_halfface_handle(hfh);
}

template <typename MeshT>
std::vector<OpenVolumeMesh::HalfFaceHandle> getBoundaryHalffaces(const MeshT& mesh, OpenVolumeMesh::HalfEdgeHandle heh)
{
    auto halffaces = std::vector<OpenVolumeMesh::HalfFaceHandle>();
    for (auto hehf_it = mesh.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
    {
        auto fh = mesh.face_handle(*hehf_it);
        if (mesh.is_boundary(fh))
            halffaces.push_back(getBoundaryHalfface(mesh, fh));
    }
    return halffaces;
}

template <typename MeshT>
std::vector<OpenVolumeMesh::HalfFaceHandle> getBoundaryHalffaces(const MeshT& mesh, OpenVolumeMesh::VertexHandle vh)
{
    auto halffaces = std::vector<OpenVolumeMesh::HalfFaceHandle>();
    for (auto voh_it = mesh.voh_iter(vh); voh_it.valid(); ++voh_it)
        for (auto hehf_it = mesh.hehf_iter(*voh_it); hehf_it.valid(); ++hehf_it)
        {
            auto fh = mesh.face_handle(*hehf_it);
            if (mesh.is_boundary(fh))
                halffaces.push_back(getBoundaryHalfface(mesh, fh));
        }

    std::sort(halffaces.begin(), halffaces.end());
    halffaces.erase(std::unique(halffaces.begin(), halffaces.end()), halffaces.end());

    return halffaces;
}




void roundMatrix(Matrix4x4d& m);

template <typename VectorT>
VectorT roundVector(VectorT vec)
{
    return VectorT(round(vec[0]), round(vec[1]), round(vec[2]));
}

bool containsVertex(const std::vector<OpenVolumeMesh::VertexHandle> vertices, OpenVolumeMesh::VertexHandle vh);


template <typename MeshT>
void checkHalfFaceValidity(MeshT& mesh)
{
    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    {
        for (auto vohe_it = mesh.voh_iter(*v_it); vohe_it.valid(); ++vohe_it)
        {
            for (auto hehf_it = mesh.hehf_iter(*vohe_it); hehf_it.valid(); ++hehf_it)
            {
                auto nextHeh = mesh.next_halfedge_in_halfface(*vohe_it, *hehf_it);
                assert(nextHeh.is_valid());
            }
        }
    }
}

template <typename MeshT>
void checkFaceValidity(MeshT& mesh)
{
    for (auto f : mesh.faces())
    {
        /*assert(mesh.incident_cell(mesh.halfface_handle(f,0)).is_valid() ||
                   mesh.incident_cell(mesh.halfface_handle(f,1)).is_valid());*/
        if (!(mesh.incident_cell(mesh.halfface_handle(f,0)).is_valid() ||
              mesh.incident_cell(mesh.halfface_handle(f,1)).is_valid()))
            mesh.delete_face(f);
    }
}

template <typename MeshT>
void checkHalfEdgeValidity(MeshT& mesh)
{

    for (auto he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end(); ++he_it)
    {
        for (auto hehf_it = mesh.hehf_iter(*he_it); hehf_it.valid(); ++hehf_it)
        {
            auto hf = mesh.halfface(*hehf_it);
            assert(std::find(hf.halfedges().begin(), hf.halfedges().end(), *he_it) != hf.halfedges().end());
            assert(mesh.next_halfedge_in_halfface(*he_it, *hehf_it).is_valid());
        }
    }
}


template <typename MeshT>
OpenVolumeMesh::HalfFaceHandle getHalfFace(MeshT& mesh, OpenVolumeMesh::VertexHandle vh0, OpenVolumeMesh::VertexHandle vh1, OpenVolumeMesh::VertexHandle vh2)
{

    auto heh = OpenVolumeMesh::HalfEdgeHandle();
    for (auto vohe_it = mesh.voh_iter(vh0); vohe_it.valid(); ++vohe_it)
    {
        auto he = mesh.halfedge(*vohe_it);
        if (he.to_vertex() == vh1)
        {
            heh = *vohe_it;
            break;
        }
    }

    if (!heh.is_valid())
        // no edge found
        return OpenVolumeMesh::HalfFaceHandle();

    auto hfh = OpenVolumeMesh::HalfFaceHandle();
    for (auto hehf_it = mesh.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
    {
        auto nextHeh = mesh.next_halfedge_in_halfface(heh, *hehf_it);
        if (!nextHeh.is_valid())
        {
#ifdef DEBUG
            checkHalfFaceValidity(mesh);
#endif
            getHalfFace(mesh, vh0,vh1,vh2);
        }
        auto he = mesh.halfedge(nextHeh);
        if (he.to_vertex() == vh2)
        {
            hfh = *hehf_it;
            break;
        }
    }

    return hfh;
}

template <typename MeshT>
OpenVolumeMesh::HalfFaceHandle getHalfFace(MeshT& mesh, OpenVolumeMesh::HalfEdgeHandle heh0, OpenVolumeMesh::HalfEdgeHandle heh1)
{

    auto heh = heh0;

    auto hfh = OpenVolumeMesh::HalfFaceHandle();
    for (auto hehf_it = mesh.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
    {
        auto nextHeh = mesh.next_halfedge_in_halfface(heh, *hehf_it);
        if (nextHeh == heh1)
        {
            hfh = *hehf_it;
            break;
        }
    }

    return hfh;
}

template <typename MeshT>
OpenVolumeMesh::HalfFaceHandle getHalfFace(MeshT& mesh, const std::vector<OpenVolumeMesh::HalfEdgeHandle> hehs)
{
    if (hehs.empty())
        return OpenVolumeMesh::HalfFaceHandle();

    auto heh = hehs.front();

    for (auto hehf_it = mesh.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
    {
        auto hf = mesh.halfface(*hehf_it);
        auto& hfedges = hf.halfedges();
        if (hfedges.size() != hehs.size())
            continue;
        auto allFound = true;
        for (auto tmpheh : hehs)
            if (std::find(hfedges.begin(), hfedges.end(), tmpheh) == hfedges.end())
                allFound = false;
        if (allFound)
            return *hehf_it;
    }

    return OpenVolumeMesh::HalfFaceHandle();;
}

template <typename MeshT>
std::vector<OpenVolumeMesh::HalfEdgeHandle> getHalfedges(MeshT& mesh, OpenVolumeMesh::VertexHandle fvh, OpenVolumeMesh::VertexHandle tvh)
{
    auto res = std::vector<OpenVolumeMesh::HalfEdgeHandle>();
    for (auto voh_it = mesh.voh_iter(fvh); voh_it.valid(); ++voh_it)
        if (mesh.halfedge(*voh_it).to_vertex() == tvh)
            res.push_back(*voh_it);
    return res;
}

OpenVolumeMesh::HalfFaceHandle getOppositeHalfFace(TetrahedralMesh& mesh, OpenVolumeMesh::CellHandle ch, OpenVolumeMesh::VertexHandle vh);


template <typename Vec>
Vec deCasteljau(std::vector<Vec> vecs, double alpha)
{
    for (auto i = 0u; i < vecs.size()-1; ++i)
        for (auto j = 0u; j < vecs.size()-1-i; ++j)
            vecs[j] = (1.0-alpha) * vecs[j] + alpha * vecs[j+1];
    return vecs[0];
}

Vec3d getRandomVector(double size);

}
