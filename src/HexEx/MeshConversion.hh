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


template <typename TetMeshT>
void convertToHexExTetrahedralMesh(const TetMeshT& _mesh, TetrahedralMesh& _tetMesh)
{
    using inPoint = typename TetMeshT::PointT;

    auto toVec3d = [&](const inPoint& point)
    {
        return Vec3d(point[0], point[1], point[2]);
    };

    // add vertices
    _tetMesh.clear(false);

    for (auto v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); ++v_it)
        _tetMesh.add_vertex(toVec3d(_mesh.vertex(*v_it)));

    // add tets
    for (auto c_it = _mesh.cells_begin(); c_it != _mesh.cells_end(); ++c_it)
    {
        auto vertices = _mesh.get_cell_vertices(*c_it);
        _tetMesh.add_cell(vertices);
    }

    // copy feature tags if available
    copy_feature_tags(_mesh, _tetMesh);
}


template <typename TetMeshT>
void copy_feature_tags(const TetMeshT& _mesh, TetrahedralMesh& _tetMesh)
{
  copy_vertex_feature_tags(_mesh, _tetMesh);
  copy_edge_feature_tags(_mesh, _tetMesh);
  copy_face_feature_tags(_mesh, _tetMesh);
}

template <typename TetMeshT>
void copy_vertex_feature_tags(const TetMeshT& _mesh, TetrahedralMesh& _tetMesh)
{
  // vertex ordering is identical --> simply copy tags if available
  if(_mesh.template vertex_property_exists<int>("AlgoHex::FeatureVertices"))
  {
    VertexProperty<int> input_vfeature;
    _mesh.request_vertex_property("AlgoHex::FeatureVertices", input_vfeature);

    VertexProperty<int> output_vfeature;
    _tetMesh.request_vertex_property("AlgoHex::FeatureVertices", output_vfeature);
    _tetMesh.set_persistent(output_vfeature,true);

    for (auto vh : _mesh.vertices())
    {
      output_vfeature[vh] = input_vfeature[vh];
    }
  }
}

template <typename TetMeshT>
void copy_edge_feature_tags(const TetMeshT& _mesh, TetrahedralMesh& _tetMesh)
{
}

template <typename TetMeshT>
void copy_face_feature_tags(const TetMeshT& _mesh, TetrahedralMesh& _tetMesh)
{
}

} // namespace HexEx
