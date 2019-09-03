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

#include <string>

#include "Typedefs.hh"
#include "HexExtractor.hh"

namespace HexEx
{

void extractHexMesh(const std::string& inFileName, const std::string& outFileName);

template <typename TetMeshT, typename HexMeshT, typename ParameterT>
void extractHexMesh(const TetMeshT& tetMesh, OpenVolumeMesh::CellPropertyT<std::map<OpenVolumeMesh::VertexHandle, ParameterT> >& parameters, HexMeshT& hexMesh)
{
    HexEx::HexExtractor he(tetMesh, parameters);
    he.extract();
    he.getHexMesh(hexMesh);
}

} // namespace HexEx
