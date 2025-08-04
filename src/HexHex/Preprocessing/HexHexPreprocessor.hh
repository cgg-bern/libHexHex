#pragma once

#include <HexHex/Utils/Typedefs.hh>

namespace HexHex
{

class HexExtractor;
class HexHexPreprocessor
{
public:
    HexHexPreprocessor(HexExtractor& hexEx);

    void preprocess();

private:
    void extractTransitionFunctions();

    void calculateSingularities();

    void sanitize();

    void scaleParameters();

    void checkIGMValidity();

    HexExtractor& hexEx;
    OVM::VertexPropertyT<bool> vertexHasIncidentNonIdentityTransitions;
    OVM::VertexPropertyT<int> nIncidentSingularEdges;
};
}

