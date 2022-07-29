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

// forward declarations
template <typename TetMeshT> void copy_feature_tags( TetMeshT& _mesh, TetrahedralMesh& _tetMesh);
template <typename TetMeshT> void copy_vertex_feature_tags( TetMeshT& _mesh, TetrahedralMesh& _tetMesh);
template <typename TetMeshT> void copy_edge_feature_tags( TetMeshT& _mesh, TetrahedralMesh& _tetMesh);
template <typename TetMeshT> void copy_face_feature_tags( TetMeshT& _mesh, TetrahedralMesh& _tetMesh);

template <typename TetMeshT>
void convertToHexExTetrahedralMesh( TetMeshT& _mesh, TetrahedralMesh& _tetMesh)
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
void copy_feature_tags( TetMeshT& _mesh, TetrahedralMesh& _tetMesh)
{
  copy_vertex_feature_tags(_mesh, _tetMesh);
  copy_edge_feature_tags(_mesh, _tetMesh);
  copy_face_feature_tags(_mesh, _tetMesh);
}

template <typename TetMeshT>
void copy_vertex_feature_tags( TetMeshT& _mesh, TetrahedralMesh& _tetMesh)
{
  // vertex ordering is identical --> simply copy tags if available
  if(_mesh.template vertex_property_exists<int>("AlgoHex::FeatureVertices"))
  {
    auto input_vfeature = _mesh.template request_vertex_property<int>("AlgoHex::FeatureVertices");
    auto output_vfeature = _tetMesh.request_vertex_property<int>("AlgoHex::FeatureVertices");;
    _tetMesh.set_persistent(output_vfeature,true);

    for (auto vh : _mesh.vertices())
    {
      output_vfeature[VertexHandle(vh.idx())] = input_vfeature[vh];
    }
  }
}

template <typename TetMeshT>
void copy_edge_feature_tags( TetMeshT& _mesh, TetrahedralMesh& _tetMesh)
{
  // vertex ordering is identical --> identify edge correspondences via vertices (multiple edges between vertex pairs not allowed!!)
  if(_mesh.template edge_property_exists<int>("AlgoHex::FeatureEdges"))
  {
    auto input_efeature = _mesh.template request_edge_property<int>("AlgoHex::FeatureEdges");
    auto output_efeature = _tetMesh.request_edge_property<int>("AlgoHex::FeatureEdges");;
    _tetMesh.set_persistent(output_efeature,true);

    for (auto eh : _mesh.edges())
    {
      auto vh0 = _mesh.halfedge(_mesh.halfedge_handle(eh, 0)).to_vertex();
      auto vh1 = _mesh.halfedge(_mesh.halfedge_handle(eh, 1)).to_vertex();

      auto output_heh = _tetMesh.find_halfedge(vh0,vh1);

      if(output_heh.is_valid())
      {
        auto output_eh = _tetMesh.edge_handle(output_heh);
        output_efeature[output_eh] = input_efeature[eh];
      }
      else
        std::cerr << "ERROR: copy_edge_feature_tags failed to find corresponding edge_handle" << std::endl;
    }
  }
}

template <typename TetMeshT>
void copy_face_feature_tags( TetMeshT& _mesh, TetrahedralMesh& _tetMesh)
{
  // vertex ordering is identical --> identify face correspondences via vertices)
  if(_mesh.template face_property_exists<int>("AlgoHex::FeatureFaces"))
  {
    auto input_ffeature = _mesh.template request_face_property<int>("AlgoHex::FeatureFaces");
    auto output_ffeature = _tetMesh.request_face_property<int>("AlgoHex::FeatureFaces");;
    _tetMesh.set_persistent(output_ffeature,true);

    for (auto fh : _mesh.faces())
    {
      auto hfh0 = _mesh.halfface_handle(fh,0);
      auto f0 = _mesh.halfface(hfh0);
      std::vector<VertexHandle> vhs0;
      vhs0.push_back(_mesh.halfedge(f0.halfedges()[0]).to_vertex());
      vhs0.push_back(_mesh.halfedge(f0.halfedges()[1]).to_vertex());
      vhs0.push_back(_mesh.halfedge(f0.halfedges()[2]).to_vertex());

      // get corresponding halfface in original mesh
      auto output_hfh0 = _tetMesh.find_halfface(vhs0);
      if(output_hfh0.is_valid())
      {
        auto output_fh = _tetMesh.face_handle(output_hfh0);
        output_ffeature[output_fh] = input_ffeature[fh];
      }
      else
        std::cerr << "ERROR: copy_face_feature_tags failed to find corresponding face_handle" << std::endl;
    }
  }
}

} // namespace HexEx
