#pragma once
#include <HexHex/Utils/Typedefs.hh>
#include <HexHex/Commons/MeshElement.hh>

namespace HexHex {

// Stores info about a picewise linear point
struct PLEPoint
{
    Position pos;
    double t;
    MeshElement generator;
    CellHandle ch;
};

} // namespace HexHex
