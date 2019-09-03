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
