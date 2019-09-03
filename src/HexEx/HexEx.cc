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


#include "HexEx.hh"

#include <OpenVolumeMesh/FileManager/FileManager.hh>

namespace HexEx {

void extractHexMesh(const std::string& inFileName, const std::string& outFileName)
{
    HexEx::HexExtractor he(inFileName);
    he.extract();

    HexEx::HexahedralMesh hexMesh;
    he.getHexMesh(hexMesh);

    OpenVolumeMesh::IO::FileManager fm;

    fm.writeFile(outFileName, hexMesh);
}

} // namespace HexEx
