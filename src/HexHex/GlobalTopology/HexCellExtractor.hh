#pragma once

#include <HexHex/GlobalTopology/cell_connectivity.hh>
#include <HexHex/LocalTopology/Propellers.hh>

namespace HexHex
{

class HexExtractor;

class HexCellExtractor
{
public:
    explicit HexCellExtractor(HexExtractor& hexEx) : hexEx(hexEx) {}

    //void extractHexFaces();

    void extractHexCells();

private:
    HexExtractor& hexEx;

    bool extractHexCell(const GlobalCornerHandle& gch, std::vector<HexCell>& hex_cells);

    std::vector<HalfFaceHandle> hfh_per_halfface; // identified by unique halfface index (4x bigger than n_halffaces, but faster than map)
    std::vector<bool> hex_cell_added; // identified by unique cell index (8x bigger than n_cells, but faster than hashmap)

    bool addHexCell(const HexCell& hex_cell);
};

}

