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


#include "FileAccessor.hh"

#include <fstream>
#include <iomanip>
#include <iostream>

bool HexEx::writeToFile(std::string fileName, TetrahedralMesh& mesh, PerCellVertexProperty<Parameter>& parameters)
{
    std::ofstream filestream(fileName, std::ofstream::out);
    return writeToStream(filestream, mesh, parameters);
//    return writeToStreamBinary(filestream, mesh, parameters);
}


bool HexEx::readFromFile(std::string fileName, TetrahedralMesh& mesh, PerCellVertexProperty<Parameter>& parameters)
{
    std::ifstream filestream(fileName, std::ifstream::in);
    if (!filestream.is_open())
    {
        std::cout << "could not open " << fileName << std::endl;
        return false;
    }
    return readFromStream(filestream, mesh, parameters);
//    return readFromStreamBinary(filestream, mesh, parameters);
}


bool HexEx::writeToStream(std::ostream& os, TetrahedralMesh& mesh, PerCellVertexProperty<Parameter>& parameters)
{
    os << std::setprecision(100);

    os << mesh.n_vertices() << std::endl;

    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
        os << mesh.vertex(*v_it) << std::endl;

    os << mesh.n_cells() << std::endl;

    for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it)
    {
        auto vh = *mesh.cv_iter(*c_it);
        auto vertices = mesh.get_cell_vertices(*c_it, vh);
        for (auto v : vertices)
            os << v << " ";
        for (auto v : vertices)
            os << parameters[*c_it][v] << " ";
        os << std::endl;
    }

    // store feature information if available
    std::vector<OpenVolumeMesh::VertexHandle> ft_vertices;
    std::vector<OpenVolumeMesh::EdgeHandle>   ft_edges;
    std::vector<OpenVolumeMesh::FaceHandle>   ft_faces;

    // collect
    if(mesh.template vertex_property_exists<bool>("AlgoHex::FeatureVertices")) {
      auto feature_vprop = mesh.template request_vertex_property<bool>("AlgoHex::FeatureVertices", false);
      for (const auto &vh: mesh.vertices())
        if (feature_vprop[vh])
          ft_vertices.push_back(vh);
    }

    if(mesh.template edge_property_exists<bool>("AlgoHex::FeatureEdges")) {
    auto feature_eprop = mesh.template request_edge_property<bool>("AlgoHex::FeatureEdges", false);
    for (const auto &eh: mesh.edges())
      if (feature_eprop[eh])
        ft_edges.push_back(eh);
    }

    if(mesh.template face_property_exists<bool>("AlgoHex::FeatureFaces")) {
      auto feature_fprop = mesh.template request_face_property<bool>("AlgoHex::FeatureFaces", false);
      for (const auto &fh: mesh.faces())
        if (feature_fprop[fh])
          ft_faces.push_back(fh);
    }

    std::cerr << "store #feature_vertices = " << ft_vertices.size()
              << "store #feature_edges = " << ft_edges.size()
              << "store #feature_vertices = " << ft_faces.size() << std::endl;

    os << ft_vertices.size() << " " << ft_edges.size() << ft_faces.size() << std::endl;

    for(const auto& vh : ft_vertices)
      os << vh.idx() << std::endl;

    for(const auto& eh : ft_edges)
      os << mesh.edge(eh).from_vertex().idx() << " " << mesh.edge(eh).to_vertex().idx() << std::endl;

    for(const auto& fh : ft_faces)
    {
      for(int i=0; i<3; ++i)
        os << mesh.edge(mesh.face(fh).halfedges()[i]).to_vertex().idx() << " ";
      os << std::endl;
    }

  return true;
}

bool HexEx::readFromStream(std::istream& is, TetrahedralMesh& mesh, PerCellVertexProperty<Parameter>& parameters)
{
    auto n_vertices = 0u;
    is >> n_vertices;

    for (unsigned int i = 0; i < n_vertices; ++i)
    {
        auto pos = Position();
        is >> pos;
        mesh.add_vertex(pos);
    }

    auto n_cells = 0u;
    is >> n_cells;

    for (auto i = 0u; i < n_cells; ++i)
    {
        VertexHandle vh0, vh1, vh2, vh3;
        is >> vh0;
        is >> vh1;
        is >> vh2;
        is >> vh3;
        auto ch = mesh.add_cell(vh0, vh1, vh2, vh3, true);

        auto param = Parameter();
        is >> param;
        parameters[ch][vh0] = param;
        is >> param;
        parameters[ch][vh1] = param;
        is >> param;
        parameters[ch][vh2] = param;
        is >> param;
        parameters[ch][vh3] = param;
    }

    // number of feature vertices/edges/faces stored in file
    int n_ftv(0), n_fte(0), n_ftf(0);

    is >> n_ftv >> n_fte >> n_ftf;

    std::cerr << "read #feature_vertices = " << n_ftv
              << "read #feature_edges = " << n_fte
              << "read #feature_vertices = " << n_ftf << std::endl;

    // request/add feature properties
    auto feature_vprop = mesh.template request_vertex_property<bool>("AlgoHex::FeatureVertices", false);
    auto feature_eprop  = mesh.template request_edge_property<bool>("AlgoHex::FeatureEdges", false);
    auto feature_fprop  = mesh.template request_face_property<bool>("AlgoHex::FeatureFaces", false);

    mesh.set_persistent(feature_vprop, true);
    mesh.set_persistent(feature_eprop, true);
    mesh.set_persistent(feature_fprop, true);

    for(int i=0; i<n_ftv; ++i)
    {
      int vidx;
      is >> vidx;
      feature_vprop[OpenVolumeMesh::VertexHandle(vidx)] = true;
    }

    for(int i=0; i<n_fte; ++i)
    {
      int v0idx, v1idx;
      is >> v0idx >> v1idx;

      OpenVolumeMesh::HalfEdgeHandle heh = mesh.halfedge(OpenVolumeMesh::VertexHandle(v0idx),OpenVolumeMesh::VertexHandle(v1idx));
      if(!heh.is_valid()) std::cerr << "ERROR: could not obtain feature edge stored in .hexex file " << v0idx << " " << v1idx << std::endl;
      else
      {
        auto eh = mesh.edge_handle(heh);
        feature_eprop[eh] = true;
      }
    }

  for(int i=0; i<n_ftf; ++i)
  {
    int v0idx, v1idx, v2idx;
    is >> v0idx >> v1idx >> v2idx;

    // map vertex indices
    std::vector<OpenVolumeMesh::VertexHandle> vhs;
    vhs.push_back(OpenVolumeMesh::VertexHandle(v0idx));
    vhs.push_back(OpenVolumeMesh::VertexHandle(v1idx));
    vhs.push_back(OpenVolumeMesh::VertexHandle(v2idx));

    // get corresponding halfface in original mesh
    OpenVolumeMesh::HalfFaceHandle hfh = mesh.halfface(vhs);
    if(!hfh.is_valid()) std::cerr << "ERROR: could not obtain feature face stored in .hexex file" << std::endl;
    else
    {
      OpenVolumeMesh::FaceHandle fh = mesh.face_handle(hfh);
      feature_fprop[fh] = true;
    }

  return true;
}


bool HexEx::writeToStreamBinary(std::ostream& os, TetrahedralMesh& mesh, PerCellVertexProperty<HexEx::Parameter>& parameters)
{
    uint32_t nVertices = mesh.n_vertices();
    os.write(reinterpret_cast<char*>(&nVertices), sizeof(decltype(nVertices)));

    auto buffer = std::vector<Position>();
    buffer.reserve(nVertices);
    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
        buffer.push_back(mesh.vertex(*v_it));

    os.write(reinterpret_cast<char*>(buffer.data()), nVertices*sizeof(Position));

    auto nCells = mesh.n_cells();
    os.write(reinterpret_cast<char*>(&nCells), sizeof(decltype(nCells)));

    for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it)
    {
        auto vh = *mesh.cv_iter(*c_it);
        auto vertices = mesh.get_cell_vertices(*c_it, vh);

        const auto verticesPerCell = 4u;

        auto handles = std::vector<VertexHandle>();
        handles.reserve(verticesPerCell);
        for (auto v : vertices)
            handles.push_back(v);
        os.write(reinterpret_cast<char*>(handles.data()), handles.size()*sizeof(decltype(handles)::value_type));

        auto params = std::vector<Parameter>();
        params.reserve(verticesPerCell);
        for (auto v : vertices)
            params.push_back(parameters[*c_it][v]);
        os.write(reinterpret_cast<char*>(params.data()), params.size()*sizeof(decltype(params)::value_type));

    }

    return true;

}


bool HexEx::readFromStreamBinary(std::istream& is, TetrahedralMesh& mesh, PerCellVertexProperty<HexEx::Parameter>& parameters)
{

    auto n_vertices = uint32_t();
    is.read(reinterpret_cast<char*>(&n_vertices), sizeof(decltype(n_vertices)));

    for (unsigned int i = 0; i < n_vertices; ++i)
    {
        auto pos = Position();
        is.read(reinterpret_cast<char*>(&pos[0]), sizeof(decltype(pos)));
        mesh.add_vertex(pos);
    }

    auto n_cells = mesh.n_cells();
    is.read(reinterpret_cast<char*>(&n_cells), sizeof(decltype(n_cells)));

    for (auto i = 0u; i < n_cells; ++i)
    {
        const auto verticesPerCell = 4u;

        auto vertices = std::vector<VertexHandle>();
        for (auto j = 0u; j < verticesPerCell; ++j)
        {
            auto vh = VertexHandle();
            is.read(reinterpret_cast<char*>(&vh), sizeof(decltype(vh)));
            vertices.push_back(vh);
        }

        auto ch = mesh.add_cell(vertices, true);

        auto params = std::vector<Parameter>();
        for (auto j = 0u; j < verticesPerCell; ++j)
        {
            auto param = Parameter();
            is.read(reinterpret_cast<char*>(&param), sizeof(decltype(param)));
            params.push_back(param);
        }

        for (auto j = 0u; j < verticesPerCell; ++j)
        {
            parameters[ch][vertices[j]] = 1.0*params[j];
        }
    }

    return true;
}
