#pragma once

#include <HexHex/LocalTopology/Propellers.hh>
#include <HexHex/PiecewiseLinear/PLEPoint.hh>
#include <HexHex/Utils/Typedefs.hh>

namespace HexHex
{

class HexEdgeExtractor
{
public:
    explicit HexEdgeExtractor(HexExtractor& hexEx) : hexEx(hexEx) {}

    void tracePropellers();

private:
    HexExtractor& hexEx;

    std::vector<GlobalPropellerOpposite> global_propellers_opposites;

    bool tracePropeller(const GlobalPropellerHandle& gph, size_t thread_id = 0);

    // For each thread store the result
    typedef struct
    {
        CellHandle startCh;

        GlobalPropellerHandle gph1;
        GlobalPropellerHandle gph2;
        int offset;

        std::vector<PLEPoint> ple_points; // points do not include start and end! ch[i] is cell of segment between points i and i+1

    } HexHalfEdge;
    size_t longest_ple_segment;
    std::vector<std::vector<HexHalfEdge>> thread_hex_edges;

    std::vector<std::atomic<bool>> global_propellers_visited_in_parallel_edge_extraction;

    EdgeHandle addHexEdge( HexHalfEdge& hex_edge);
};

}

