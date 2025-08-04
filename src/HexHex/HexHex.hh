/*
 * Copyright 2024 Computer Graphics Group, University of Bern - Tobias Kohler <tobias.kohler@unibe.ch>
 * Copyright 2016 Computer Graphics Group, RWTH Aachen University - Max Lyon <lyon@cs.rwth-aachen.de>
 *
 * This file is part of HexHex.
 */


#pragma once

#include <HexHex/Config/Export.hh>
#include <HexHex/Config.hh>
#include <HexHex/Utils/Typedefs.hh>
#include <nlohmann/json.hpp>

namespace HexHex
{

struct HexHexOutput
{
    std::unique_ptr<HexahedralMesh> hex_mesh;
    std::unique_ptr<PolyhedralMesh> piecewise_linear_mesh;
    nlohmann::json report;
    bool success;
};

HexHexOutput HEXHEX_EXPORT extractHexMesh(const TetrahedralMesh& tetmesh,
                                          const OVM::HalfFacePropertyT<Vec3d>& igm,
                                          const Config& config = Config());

}
