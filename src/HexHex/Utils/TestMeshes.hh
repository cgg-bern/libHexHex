#pragma once

#include <HexHex/Config/Export.hh>
#include <HexHex/Utils/Utils.hh>

namespace HexHex
{

struct TestMesh {
    TestMesh() {
        mesh.set_persistent(igm);
    };
    TetrahedralMesh mesh;
    HFParam igm = mesh.request_halfface_property<Vec3d>("HexHex::Parametrization");
};

//TestMesh createHexCube();
TestMesh HEXHEX_EXPORT createCube(unsigned int tetLength=1, unsigned int hexLength=1);
TestMesh HEXHEX_EXPORT createCylinder(unsigned int valence=5, unsigned int hexScale=1);
void HEXHEX_EXPORT randomizeTransitions(TestMesh&);

}

