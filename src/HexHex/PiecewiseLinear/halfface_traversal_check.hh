#pragma once

#include <HexHex/Utils/Typedefs.hh>

namespace HexHex
{

// Flood Fill for Piecewise Linear Face Extraction
bool shouldVisitOppositeCellOverHalfface(HexExtractor& he, const HalfFaceHandle& hfh,
                            const Transition& trans2Ref, const std::array<Parameter, 4>& refCorners3d);

}

