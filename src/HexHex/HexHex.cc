/*
 * Copyright 2024 Computer Graphics Group, University of Bern - Tobias Kohler <tobias.kohler@unibe.ch>
 * Copyright 2016 Computer Graphics Group, RWTH Aachen University - Max Lyon <lyon@cs.rwth-aachen.de>
 *
 * This file is part of HexHex.
 */


#include <HexHex/HexHex.hh>
#include <HexHex/HexExtractor.hh>

#include <nlohmann/json.hpp>

namespace HexHex {

HexHexOutput HEXHEX_EXPORT extractHexMesh(const TetrahedralMesh& tetmesh, const OVM::HalfFacePropertyT<Vec3d>& igm, const Config& config)
{
    ScopedStopWatch _{sw::root};
    //auto prop = tetmesh.create_private_halfface_property<double>("test");
    HexHex::HexExtractor he(tetmesh, igm);
    HexHexOutput res;
    res.success = he.extract(config);
    res.report = he.getHexHexReport();
    res.hex_mesh = he.takeOutputMesh();
    if (config.extract_piecewise_linear_edges || config.extract_piecewise_linear_faces) {
        res.piecewise_linear_mesh = std::make_unique<PolyhedralMesh>(std::move(he.getPiecewiseLinearMesh()));
    }
    return res;
}

} // namespace HexEx
