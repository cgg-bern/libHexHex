#include <HexHex/GlobalTopology/cell_connectivity.hh>
#include <HexHex/LocalTopology/HexVertexGenerator.hh>

namespace HexHex {

HexHalfFace create_opposite_halfface(const HexExtractor& hexEx, const HexHalfFace& hf)
{
    // Get the four generators
    const auto& gen1 = hexEx.getHexVertexGenerator(hf.quarters[0].vh);
    const auto& gen2 = hexEx.getHexVertexGenerator(hf.quarters[1].vh);
    const auto& gen3 = hexEx.getHexVertexGenerator(hf.quarters[2].vh);
    const auto& gen4 = hexEx.getHexVertexGenerator(hf.quarters[3].vh);

    // Find the opposite corners
    const auto o1 = gen1.findLocalCornerOpposite(hf.quarters[0].lca);
    if (!o1.lch.is_valid()) return INVALID_HALFFACE;
    const auto o2 = gen2.findLocalCornerOpposite(hf.quarters[1].lca);
    if (!o2.lch.is_valid()) return INVALID_HALFFACE;
    const auto o3 = gen3.findLocalCornerOpposite(hf.quarters[2].lca);
    if (!o3.lch.is_valid()) return INVALID_HALFFACE;
    const auto o4 = gen4.findLocalCornerOpposite(hf.quarters[3].lca);
    if (!o4.lch.is_valid()) return INVALID_HALFFACE;

    // Create the opposite face (ensure reversed order)
    return create_canonical_halfface(
        hf.quarters[3].vh,o4.lch,o4.forward,
        hf.quarters[2].vh,o3.lch,o3.forward,
        hf.quarters[1].vh,o2.lch,o2.forward,
        hf.quarters[0].vh,o1.lch,o1.forward
    );
};

HexHalfFace create_canonical_halfface(
    const VertexHandle vh1, const LCH lch1, const LPH forward1,
    const VertexHandle vh2, const LCH lch2, const LPH forward2,
    const VertexHandle vh3, const LCH lch3, const LPH forward3,
    const VertexHandle vh4, const LCH lch4, const LPH forward4
    )
{
    // A face must have four different vertices
    assert(vh1 != vh2 && vh1 != vh3 && vh1 != vh4 && vh2 != vh3 && vh2 != vh4 && vh3 != vh4);

    // Create the canonical order
    if (vh1 < vh2 && vh1 < vh3 && vh1 < vh4)
    {
        return HexHalfFace{.quarters = {
        create_halfface_quarter(vh1, lch1, forward1),
        create_halfface_quarter(vh2, lch2, forward2),
        create_halfface_quarter(vh3, lch3, forward3),
        create_halfface_quarter(vh4, lch4, forward4)
        }};
    }
    if (vh2 < vh3 && vh2 < vh4)
    {
        return HexHalfFace{.quarters = {
        create_halfface_quarter(vh2, lch2, forward2),
        create_halfface_quarter(vh3, lch3, forward3),
        create_halfface_quarter(vh4, lch4, forward4),
        create_halfface_quarter(vh1, lch1, forward1)
       }};
    }
    if (vh3 < vh4)
    {
        return HexHalfFace{.quarters = {
           create_halfface_quarter(vh3, lch3, forward3),
           create_halfface_quarter(vh4, lch4, forward4),
           create_halfface_quarter(vh1, lch1, forward1),
           create_halfface_quarter(vh2, lch2, forward2)
       }};
    }
    return HexHalfFace{.quarters = {
       create_halfface_quarter(vh4, lch4, forward4),
       create_halfface_quarter(vh1, lch1, forward1),
       create_halfface_quarter(vh2, lch2, forward2),
    create_halfface_quarter(vh3, lch3, forward3)
   }};
}

} // namespace HexHex
