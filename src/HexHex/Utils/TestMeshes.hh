#pragma once

#include <HexHex/Config/Export.hh>
#include <HexHex/Utils/Utils.hh>
#include <HexHex/Utils/ParametrizedMesh.hh>

namespace HexHex
{

//TestMesh createHexCube();
ParametrizedMesh HEXHEX_EXPORT createCube(unsigned int tetLength=1, unsigned int hexLength=1);
ParametrizedMesh HEXHEX_EXPORT createCylinder(unsigned int valence=5, unsigned int hexScale=1);
void HEXHEX_EXPORT randomizeTransitions(ParametrizedMesh&);

}

