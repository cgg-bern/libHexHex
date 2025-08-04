#pragma once
#include <HexHex/Utils/Typedefs.hh>
#include <HexHex/Commons/MeshElement.hh>

namespace HexHex {

struct PLFPoint
{
    Vec3d pos;
    MeshElement generator;
    VertexHandle ple_vertex; // Vertex preexists from Edge Extraction
    Vec2d uv;
};

inline std::ostream& operator<<(std::ostream& os, const PLFPoint& point)
{
    return os << "PLFPoint("
              << point.pos << ", "
              << point.generator << ", "
              << point.ple_vertex << ", "
              << point.uv << ")";
}


} // namespace HexHex
