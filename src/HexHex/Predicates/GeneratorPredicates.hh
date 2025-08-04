#pragma once

#include <HexHex/Utils/TetMeshCache.hh>
#include <HexHex/Utils/Typedefs.hh>

namespace HexHex
{

// Forward declaration
class MeshElement;
class Direction;

// Returns the Element of the cell on the interior of which param lies.
MeshElement computeVertexGeneratorElement(const TetMeshCache& mesh,
                                          const CellHandle ch,
                                          const std::array<VertexHandle,4>& vhs,
                                          const std::array<Parameter,4>& params, const Parameter& param);


MeshElement computePropellerHolderElement(
    const TetMeshCache& mesh,
    const CellHandle ch, const std::array<VertexHandle,4>& vhs, const std::array<Parameter,4>& params,
    const MeshElement& generator,
    const Parameter& from, const Direction dir);

MeshElement computePropellerBladeCasingElement(
    const TetMeshCache& mesh,
    const CellHandle ch, const std::array<VertexHandle,4>& vhs, const std::array<Parameter,4>& params,
    const MeshElement& holder, const Parameter &from,
    const Direction dir1, const Direction dir2);

}

