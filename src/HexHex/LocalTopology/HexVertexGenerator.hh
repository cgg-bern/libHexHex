#pragma once

#include <HexHex/Commons/Direction.hh>
#include <HexHex/Commons/MeshElement.hh>
#include <HexHex/Utils/Typedefs.hh>
#include <HexHex/LocalTopology/Propellers.hh>

namespace HexHex
{

class HexExtractor;

typedef struct
{
    uint32 hex_vertices_offset;
    uint32 global_propellers_offset;
    uint32 global_corners_offset;
} GeneratorOffsets; // How many Elements in total in the previous generators

class HexVertexGenerator
{
public:
    friend class HexEdgeExtractor;
    friend class HexCellExtractor;
    friend class PLHexFaceExtractor;

    HexExtractor& hexEx;
    MeshElement generator;
    std::pair<VertexHandle, VertexHandle> vh_range;
    bool integer_aligned;
    std::vector<LocalPropeller> localPropellers;
    std::vector<LocalCorner> localCorners;
    GeneratorOffsets global_offsets;

    HexVertexGenerator(HexExtractor& hexEx, VertexHandle generator, VertexHandle hexVertex);

    HexVertexGenerator(HexExtractor& hexEx, EdgeHandle generator, std::pair<VertexHandle, VertexHandle> vh_range);

    HexVertexGenerator(HexExtractor& hexEx, HalfFaceHandle generator, std::pair<VertexHandle, VertexHandle> vh_range);

    HexVertexGenerator(HexExtractor& hexEx, CellHandle generator, std::pair<VertexHandle, VertexHandle> vh_range);

    std::array<size_t,3> extractLocalTopology();

    void enumerateLocalPropellers();

    void enumerateLocalPropellerBlades();

    void enumerateLocalCorners();

    EdgeHandle getHolderEdge(const HalfEdgeHandle& hexHeh) const;
    FaceHandle getCasingFace(const HalfEdgeHandle& heh1, const HalfEdgeHandle& heh2) const;
    MeshElement getHolder(const LPH& lph) const;
    size_t getValence(const LPH& lph) const;

    inline uint nVertices() const {return vh_range.second.idx()-vh_range.first.idx()+1u;}
    inline uint nLocalPropellers() const {return generator.is_cell()? 6u : (generator.is_face()&&isBoundary())? 5u : localPropellers.size();}
    inline uint nGlobalPropellers() const {return nVertices() * nLocalPropellers();}
    inline uint nLocalCorners() const {return localCorners.size();}
    inline uint nGlobalCorners() const {return nVertices() * nLocalCorners();}
    inline bool hasVertex(VertexHandle vh) const {return vh.idx() >= vh_range.first.idx() && vh.idx() <= vh_range.second.idx();}

    bool isBoundary() const;

    LocalPropellerHandle addLocalPropeller(MeshElement holder, CellHandle ch, Direction dir);
    LocalCornerHandle addLocalCorner(LocalPropellerHandle lph1, LocalPropellerHandle lph2, LocalPropellerHandle lph3);

    // Returns the generator specific index of the global propeller
    inline uint32 toLocalIndex(const GPH& gph) const
    {
        return gph.lph.idx()*nVertices()+(gph.vh.idx()-vh_range.first.idx());
    }

    inline uint32 toLocalIndex(const VertexHandle& vh, const LocalCornerHandle& lch) const
    {
        return lch.idx()*nVertices()+(vh.idx()-vh_range.first.idx());
    }

    // Type soecific

    void getDirections(const LPH& lph1, const LPH& lph2, CellHandle& ch, Direction& d1, Direction& d2);
    int nBlades(const LPH& lph) const;
    int findBladeIndex(const LPH& lph1, const LPH& lph2) const;
    int getBladeIndex(const LPH& lph1, const CellHandle& cch, const Direction& d1, const Direction& d2) const;
    LocalPropellerHandle getBlade(const LPH& lph, const int& j) const;

    std::tuple<LocalPropellerHandle, CellHandle, Transition> findBlade(const LocalPropeller& lp, CellHandle ch, Direction d1, Direction d2);
    HalfFaceHandle findBladeNextHalfFace(const LocalPropeller& lp, CellHandle ch, HalfFaceHandle hfh, Direction d1, Direction d2);

    LocalPropellerHandle findLocalPropeller(const CellHandle& cch, const Direction& dir) const;

    LocalCornerHandle findLocalCorner(const LocalPropellerHandle& lph1, const LocalPropellerHandle& lph2) const;

    // Finds the local corner opposite to the implied halfface quarter and returns it with the equivalent forward
    LocalCornerAxis findLocalCornerOpposite(const LocalCornerAxis& lca) const;

    inline const LocalCorner& getLocalCorner(const GCH& gch) const
    {
        return localCorners[gch.lch.idx()];
    }

    inline const LocalCorner& getLocalCorner(const LCH& lch) const
    {
        return localCorners[lch.idx()];
    }

private:
    enum GeneratorType : unsigned char
    {
        UNDEFINED = 0,
        CELL = 1,
        BOUNDARY_FACE = 2,
        GENERAL = 3
    };
    GeneratorType type;

    uint n_local_blades;
    uint n_local_connected_blades;

    Direction boundary_face_direction_into_cell = Direction(0);

    // Tracing specific
    void pickTracingStart(const GPH& gph, CellHandle& ch, Direction& d1, Direction& d2, int& bladeIndex1);
};

}

