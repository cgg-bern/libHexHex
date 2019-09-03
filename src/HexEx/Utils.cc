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


#include "Utils.hh"

#include <random>

using namespace OpenVolumeMesh;

typedef OpenVolumeMeshEdge Edge;
typedef OpenVolumeMeshFace Face;
typedef OpenVolumeMeshCell Cell;
typedef HexEx::PolyhedralMesh::PointT Vec3;

bool HexEx::isTetCell(const PolyhedralMesh& mesh, CellHandle ch)
{
    if (mesh.n_vertices_in_cell(ch) == 4)
        return true;
    else
        return false;
}


bool HexEx::isTetMesh(const PolyhedralMesh& mesh)
{
    for (CellIter c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it)
        if (!isTetCell(mesh, *c_it))
            return false;
    return true;
}


bool HexEx::isFlipped(const TetrahedralMesh& mesh, CellHandle ch)
{
    std::vector<VertexHandle> vertices = mesh.get_cell_vertices(ch);
    std::vector<PolyhedralMesh::PointT> positions;
    for (unsigned int i = 0; i < vertices.size(); ++i)
        positions.push_back(mesh.vertex(vertices[i]));
    Vec3 d1 = positions[1] - positions[0];
    Vec3 d2 = positions[2] - positions[0];
    Vec3 d3 = positions[3] - positions[0];
    return (((d1%d2)|d3) < 0.0);
}

bool HexEx::containsFlippedTets(const TetrahedralMesh& mesh)
{
    for (CellIter c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it)
        if (isFlipped(mesh, *c_it))
            return true;
    return false;
}

bool HexEx::isIrregularEdge(bool boundary, int valence)
{
    if (boundary)
        return valence != 2;
    else
        return valence != 4;
}

HexEx::Vec3d HexEx::getValenceColor(bool boundary, int valence)
{
    if(!boundary && valence == 3) {
        return HexEx::Vec3d(0.0f, 1.0f, 1.0f);
    } else if(!boundary && valence == 5) {
        return HexEx::Vec3d(1.0f, 0.0f, 1.0f);
    } else if(boundary && valence > 2) {
        return HexEx::Vec3d(0.0f, 1.0f, 0.0f);
    } else if(boundary && valence == 1) {
        return HexEx::Vec3d(1.0f, 1.0f, 0.0f);
    } else if(!boundary && valence > 5) {
        return HexEx::Vec3d(0.5f, 1.0f, 0.5f);
    }
    return HexEx::Vec3d(0.0f, 0.0f, 0.0f);

    if (valence == 0)
        return HexEx::Vec3d(0, 0, 0);
    else if (valence == 1)
        return HexEx::Vec3d(1, 0, 0);
    else if (valence == 2)
        return HexEx::Vec3d(0, 1, 0);
    else if (valence == 3)
        return HexEx::Vec3d(0, 0, 1);
    else if (valence == 4)
        return HexEx::Vec3d(1, 1, 0);
    else if (valence == 5)
        return HexEx::Vec3d(1, 0, 1);
    else if (valence == 6)
        return HexEx::Vec3d(0, 1, 1);
    else
        return HexEx::Vec3d(1.0, 1.0, 1.0);
}


std::vector<EdgeHandle> HexEx::getCellEdges(const TetrahedralMesh& mesh, CellHandle ch)
{
    auto res = std::vector<EdgeHandle>();

    auto c = mesh.cell(ch);
    for (auto hfh : c.halffaces())
    {
        auto hf = mesh.halfface(hfh);
        for (auto heh : hf.halfedges())
            res.push_back(mesh.edge_handle(heh));
    }


    std::sort( res.begin(), res.end() );
    res.erase( std::unique( res.begin(), res.end() ), res.end() );

    return res;
}


std::vector<HalfFaceHandle> HexEx::getIncidentHalfFaces(TetrahedralMesh& mesh, CellHandle ch, EdgeHandle eh)
{
    Cell c = mesh.cell(ch);
    Edge e = mesh.edge(eh);
    std::vector<HalfFaceHandle> res;
    for (auto hfh : c.halffaces())
        if (containsVertex(mesh.get_halfface_vertices(hfh), e.from_vertex()) &&
            containsVertex(mesh.get_halfface_vertices(hfh), e.to_vertex()))
            res.push_back(hfh);
    return res;
}

std::vector<HalfFaceHandle> HexEx::getIncidentHalfFaces(TetrahedralMesh& mesh, CellHandle ch, VertexHandle vh)
{
    Cell c = mesh.cell(ch);
    std::vector<HalfFaceHandle> res;
    for (auto hfh : c.halffaces())
        if (containsVertex(mesh.get_halfface_vertices(hfh), vh))
            res.push_back(hfh);
    return res;
}

HalfFaceHandle HexEx::getIncidentHalfFaceInCell(TetrahedralMesh& mesh, CellHandle ch, HalfEdgeHandle heh)
{
    Cell c = mesh.cell(ch);
    for (auto hfh : c.halffaces())
    {
        Face f = mesh.halfface(hfh);
        if (std::find(f.halfedges().begin(), f.halfedges().end(), heh) != f.halfedges().end())
            return hfh;
    }

    return HalfFaceHandle();
}

std::vector<HalfEdgeHandle> HexEx::getIncidentHalfEdges(TetrahedralMesh& mesh, CellHandle ch, VertexHandle vh)
{
    assert(containsVertex(mesh.get_cell_vertices(ch), vh));

    std::vector<HalfEdgeHandle> res;
    for (auto hfh : getIncidentHalfFaces(mesh, ch, vh))
    {
        Face hf = mesh.halfface(hfh);
        for (auto heh : hf.halfedges())
        {
            Edge he = mesh.halfedge(heh);
            if (he.from_vertex() == vh)
                res.push_back(heh);
        }
    }
    return res;
}


HalfEdgeHandle HexEx::getHalfEdgeInHalfFace(TetrahedralMesh& mesh, EdgeHandle eh, HalfFaceHandle hfh)
{
    auto hf = mesh.halfface(hfh);
    auto heh = mesh.halfedge_handle(eh, 0);
    if (std::find(hf.halfedges().begin(), hf.halfedges().end(), heh) == hf.halfedges().end())
        heh = mesh.opposite_halfedge_handle(heh);

    return heh;
}

void HexEx::roundMatrix(Matrix4x4d& m)
{
    for (int row = 0; row < 4; ++row)
        for (int col = 0; col < 4; ++col)
        {
            if (fabs(m(row,col) - round(m(row,col))) > 0.01)
            {
#ifdef DEBUG
                std::cout << "WAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAh!!!!!!!!!!!!!!!!! " << fabs(m(row,col) - round(m(row,col)));
                std::cout << std::endl;
#endif

            }
            m(row, col) = round(m(row,col));
        }
}


unsigned int HexEx::getMaxCoord(Vec3d n)
{
    if ((fabs(n[0]) > fabs(n[1])) && (fabs(n[0]) > fabs(n[2])))
        return 0u;
    else if (fabs(n[1]) > fabs(n[2]))
        return 1u;
    else
        return 2u;
}


bool HexEx::containsVertex(const std::vector<VertexHandle> vertices, VertexHandle vh)
{
    return std::find(vertices.cbegin(), vertices.cend(), vh) != vertices.cend();
}


HalfFaceHandle HexEx::getIncidentHalfFaceInCell(TetrahedralMesh& mesh, CellHandle ch, FaceHandle fh)
{
    auto hfh = mesh.halfface_handle(fh, 0);
    if (mesh.incident_cell(hfh) != ch)
        hfh = mesh.halfface_handle(fh, 1);
    if (mesh.incident_cell(hfh) != ch)
        hfh = HalfFaceHandle();
    return hfh;
}



std::vector<VertexHandle> HexEx::getNeighborVertices(TetrahedralMesh& mesh, VertexHandle vh)
{
    auto res = std::vector<VertexHandle>();
    for (auto voh_it = mesh.voh_iter(vh); voh_it.valid(); ++voh_it)
    {
        auto he = mesh.halfedge(*voh_it);
        res.push_back(he.to_vertex());
    }
    std::sort(res.begin(), res.end());
    return res;
}


HalfFaceHandle HexEx::getOppositeHalfFace(TetrahedralMesh& mesh, CellHandle ch, VertexHandle vh)
{
    auto vertices = mesh.get_cell_vertices(ch, vh);

    assert(vertices.front() == vh);

    return getHalfFace(mesh, vertices[1], vertices[3], vertices[2]);
}


HexEx::Vec3d HexEx::getRandomVector(double size)
{
    static std::default_random_engine e1(3);
    static std::uniform_int_distribution<int> uniform_dist1(-100, 100);

    auto res = HexEx::Vec3d(0, 0, 0);
    for (unsigned int i = 0; i < 3; ++i)
        res[i] = uniform_dist1(e1)*size/100;

    return res;
}

