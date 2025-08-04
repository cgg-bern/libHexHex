#pragma once

#include <HexHex/HexExtractor.hh>
#include <HexHex/LocalTopology/Propellers.hh>
#include <array>

namespace HexHex {

// A quarter of a halfface (essentially a dart, lol)
typedef struct
{
    VertexHandle vh;
    LocalCornerAxis lca;
} HexHalfFaceQuarter;

static const HexHalfFaceQuarter INVALID_HALFFACE_QUARTER = HexHalfFaceQuarter{
    .vh = VertexHandle(-1), .lca = LocalCornerAxis{.lch = LCH(-1), .forward = LPH(-1)}
};

inline GPH getGPH(const HexHalfFaceQuarter& hfq)
{
    return GPH{.vh = hfq.vh, .lph = hfq.lca.forward};
}

// Compare two HexHalfFaceQuarters
inline bool operator==(const HexHalfFaceQuarter& hfq1, const HexHalfFaceQuarter& hfq2)
{
    return ((hfq1.vh == hfq2.vh) && (hfq1.lca.lch == hfq2.lca.lch) && (hfq1.lca.forward == hfq2.lca.forward));
};

// a halfface
typedef struct
{
    std::array<HexHalfFaceQuarter, 4> quarters;
} HexHalfFace;

static const HexHalfFace INVALID_HALFFACE = HexHalfFace{.quarters = {
    INVALID_HALFFACE_QUARTER,
    INVALID_HALFFACE_QUARTER,
    INVALID_HALFFACE_QUARTER,
    INVALID_HALFFACE_QUARTER
}};

// a face
typedef struct
{
    HexHalfFace halfface;
    HexHalfFace opposite_halfface;
} HexHalfFaceWithOpposite;

// Get the four ordered global propeller handles of a hex halfface
inline std::array<GPH,4> getGPHs(const HexHalfFace& hf)
{
    return {getGPH(hf.quarters[0]),getGPH(hf.quarters[1]),getGPH(hf.quarters[2]),getGPH(hf.quarters[3])};
}

// Compare two HexHalfFaces, since they are in canonical form, the order of quarters must match
inline bool operator==(const HexHalfFace& f1, const HexHalfFace& f2)
{
    return f1.quarters == f2.quarters;
};

inline HexHalfFaceQuarter create_halfface_quarter(const VertexHandle vh, const LCH lch, const LPH forward)
{
    return HexHalfFaceQuarter{.vh = vh, .lca = LocalCornerAxis{.lch = lch, .forward = forward}};
}

// Create a canonical (uniquely identifiable) hex-half-face by starting with the smallest vertex index
HexHalfFace create_canonical_halfface(
    const VertexHandle vh1, const LCH lch1, const LPH forward1,
    const VertexHandle vh2, const LCH lch2, const LPH forward2,
    const VertexHandle vh3, const LCH lch3, const LPH forward3,
    const VertexHandle vh4, const LCH lch4, const LPH forward4
    );

HexHalfFace create_opposite_halfface(const HexExtractor& hexEx, const HexHalfFace& hf);

inline HexHalfFaceWithOpposite create_canonical_halfface_with_opposite(const HexExtractor& hexEx,
    const VertexHandle vh1, const LCH lch1, const LPH forward1,
    const VertexHandle vh2, const LCH lch2, const LPH forward2,
    const VertexHandle vh3, const LCH lch3, const LPH forward3,
    const VertexHandle vh4, const LCH lch4, const LPH forward4
    )
{
    HexHalfFace hf = create_canonical_halfface(
        vh1,lch1,forward1,
        vh2,lch2,forward2,
        vh3,lch3,forward3,
        vh4,lch4,forward4
    );
    return HexHalfFaceWithOpposite{.halfface = hf, .opposite_halfface = create_opposite_halfface(hexEx, hf)};
}

// cell
typedef struct
{
    std::array<HexHalfFaceWithOpposite,6> faces; // 1st halfface of a face is the incident one
    std::vector<uint32> global_corner_indices; // 8 for a complete cell
    uint32 min_corner_idx;
} HexCell;
}

