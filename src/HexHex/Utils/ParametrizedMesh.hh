#pragma once

#include <HexHex/Utils/Typedefs.hh>
namespace HexHex
{

struct ParametrizedMesh {
    ParametrizedMesh() {
        mesh.set_persistent(igm);
    };
    ParametrizedMesh(TetrahedralMesh _mesh, HFParam _igm)
        : mesh(std::move(_mesh))
        , igm(std::move(_igm))
    {
        mesh.set_persistent(igm);
    };
    TetrahedralMesh mesh;
    HFParam igm = mesh.request_halfface_property<Vec3d>("HexHex::Parametrization");
};

} // namespace HexHex
