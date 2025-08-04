#pragma once

#include <HexHex/Commons/MeshElement.hh>
#include <HexHex/Commons/Permutation.hh>
#include <HexHex/Utils/Typedefs.hh>

namespace HexHex {

class HexExtractor;

class VertexExtractor
{
public:
    VertexExtractor(HexExtractor& hexEx);

    bool extract();

private:
    HexExtractor& hexEx;
    size_t num_generators;
    size_t num_hex_vertices;

    void addIntegerGridPoint(const CellHandle& ch, const MeshElement& elem, const Parameter &param);

    OVM::VertexPropertyT<std::vector<Parameter>> v_igps;
    OVM::EdgePropertyT<std::vector<Parameter>> e_igps;
    OVM::FacePropertyT<std::vector<Parameter>> f_igps;
    OVM::CellPropertyT<std::vector<Parameter>> c_igps;

    OVM::CellPropertyT<Permutation> piInv;

    OVM::CellPropertyT<std::vector<Parameter>> tet_debug_props;
};

}

