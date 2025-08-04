#pragma once

#include <HexHex/LocalTopology/Propellers.hh>
#include <absl/container/flat_hash_map.h>
#include <vector>

namespace HexHex
{

static const std::vector<LocalCorner> CELL_GENERATOR_LOCAL_CORNERS = {{
    {0, 2, 4},
    {0, 4, 3},
    {0, 3, 5},
    {0, 5, 2},
    {1, 2, 5},
    {1, 5, 3},
    {1, 3, 4},
    {1, 4, 2}
}};

static const absl::flat_hash_map<std::pair<LPH,LPH>,LCH> CELL_GENERATORS_LOCAL_CORNERS_FROM_PROPELLER_PAIR = {{
    {{0,2}, 0},
    {{2,4}, 0},
    {{4,0}, 0},

    {{0,4}, 1},
    {{4,3}, 1},
    {{3,0}, 1},

    {{0,3}, 2},
    {{3,5}, 2},
    {{5,0}, 2},

    {{0,5}, 3},
    {{5,2}, 3},
    {{2,0}, 3},

    {{1,2}, 4},
    {{2,5}, 4},
    {{5,1}, 4},

    {{1,5}, 5},
    {{5,3}, 5},
    {{3,1}, 5},

    {{1,3}, 6},
    {{3,4}, 6},
    {{4,1}, 6},

    {{1,4}, 7},
    {{4,2}, 7},
    {{2,1}, 7},
}};

static const std::vector<LocalCorner> BOUNDARY_FACE_GENERATOR_LOCAL_CORNERS = {{
    {0, 1, 2},
    {0, 2, 3},
    {0, 3, 4},
    {0, 4, 1}
}};

// For a bounadry face the local propeller handle with index 0 points into the cell and the propellers with indices
// 1,2,3,4 are orthogonal along the (integer aligned) face.
static const std::array<std::vector<LPH>, 5> BOUNDARY_FACE_GENERATOR_BLADE_LPHS_PER_LPH_CCW = {{
    {1, 2, 3, 4},
    {2, 0, 4},
    {3, 0, 1},
    {4, 0, 2},
    {1, 0, 3}
}};

static const std::array<std::array<int, 5>, 5> BOUNDARY_FACE_GENERATOR_BLADE_INDICES_PER_LPH_CCW = {{
    {-1, 0, 1, 2, 3},
    {1, -1, 0, -1, 2},
    {1, 2, -1, 0, -1},
    {1, -1, 2, -1, 0},
    {1, 0, -2, 2, -1},
}};


}

