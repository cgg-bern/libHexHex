#pragma once

/*
 * Copyright 2024 Computer Graphics Group, University of Bern - Tobias Kohler <tobias.kohler@unibe.ch>
 * Copyright 2016 Computer Graphics Group, RWTH Aachen University - Max Lyon <lyon@cs.rwth-aachen.de>
 *
 * This file is part of HexHex.
 */


#include <HexHex/Config/Export.hh>
#include <HexHex/Utils/Typedefs.hh>
#include <HexHex/Utils/ParametrizedMesh.hh>
#include <filesystem>

namespace HexHex
{

std::optional<ParametrizedMesh> OVM_EXPORT loadInputFromFile(const std::filesystem::path& filename);

bool OVM_EXPORT saveInputToHEXEX(const std::filesystem::path& filename, const TetrahedralMesh& tetmesh, const OVM::HalfFacePropertyT<Vec3d>& igm);

bool OVM_EXPORT saveOutputToFile(const std::filesystem::path& filename, const HexahedralMesh& mesh);

}

