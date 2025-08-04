#pragma once

#include <HexHex/Commons/CellIGM.hh>
#include <HexHex/Commons/FlatMap.hh>
#include <HexHex/Commons/Transition.hh>
#include <HexHex/Config.hh>
#include <HexHex/Utils/TetMeshCache.hh>
#include <HexHex/Utils/Typedefs.hh>
#include <HexHex/Utils/Utils.hh>
#include <absl/container/flat_hash_map.h>

// Abseil: Parameter Hash
namespace OpenVolumeMesh::Geometry
{
template <typename H>
H AbslHashValue(H h, const Vec3d& v) {
    return H::combine(std::move(h), v[0], v[1], v[2]);
}
}  // namespace OpenVolumeMesh::Geometry

namespace HexHex
{

class Parametrization
{

friend class HexHexPreprocessor;
//friend class NavigationTest;

public:
    Parametrization(TetMeshCache& cache, HexahedralMesh& hexmesh);

    template<typename VecT>
    void transferFromHalfFaceProp(const TetrahedralMesh& tetmesh, const OVM::HalfFacePropertyT<VecT>& igm)
    {
#pragma omp parallel for //default(none) shared(igm,tetmesh,tetVertexCellParameters)
        for (uint i = 0; i < tetmesh.n_cells(); ++i) {
            CellHandle ch(i);
            for (HalfFaceHandle hfh : tetmesh.cell(ch).halffaces()) {
                VertexHandle vh = tetmesh.halfface_opposite_vertex(hfh);
                tetVertexCellParameters[ch][vh] = toVec3d(igm[hfh]);
            }
        }
    }

    template <typename VecT, typename MapT>
    void transferFromCellVertexMapProp(const TetrahedralMeshT<VecT>& tetmesh, const OVM::CellPropertyT<MapT>& igm)
    {
#pragma omp parallel for // default(none) shared(igm,tetmesh,tetVertexCellParameters)
        for (uint i = 0; i < tetmesh.n_cells(); ++i) {
            CellHandle ch(i);
            const auto& vhs = tetmesh.get_cell_vertices(ch);
            for (VertexHandle vh : vhs) {
                tetVertexCellParameters[ch][vh] = toVec3d(igm[ch][vh]);
            }
        }
    }

    inline void set(const CellHandle ch, const VertexHandle vh, const Parameter& param)
    {
        tetVertexCellParameters[ch][vh] = param;
    }

    inline const Parameter& get(const CellHandle ch, const VertexHandle vh) const
    {
        return tetVertexCellParameters[ch].at(vh);
    }

    inline const std::array<Parameter, 3> get(const CellHandle& ch, const VertexHandle& vh1, const VertexHandle& vh2, const VertexHandle& vh3)
    {
        return {get(ch,vh1),get(ch,vh2),get(ch,vh3)};
    }

    inline const std::array<Parameter, 4> get(const CellHandle& ch) const
    {
        const auto& vhs = cache.cellVertices[ch];
        return {get(ch, vhs[0]), get(ch, vhs[1]), get(ch, vhs[2]), get(ch, vhs[3])};
    }

    inline const std::array<Parameter, 3> get(const HalfFaceHandle& hfh) const
    {
        const auto& ch = cache.tetmesh.incident_cell(hfh);
        const auto& vhs = cache.tetmesh.get_halfface_vertices(hfh);
        return {get(ch, vhs[0]), get(ch, vhs[1]), get(ch, vhs[2])};
    }

    inline Transition& transition(const HalfFaceHandle& hfh)
    {
        return transitions[hfh];
    }

    Transition getTransition(const CellHandle& fromCell, const CellHandle& toCell, const EdgeHandle& eh);

    Transition getTransition(const CellHandle& fromCell, const CellHandle& toCell);

    inline Parameter getHexVertexParameter(const CellHandle& ch, const VertexHandle& vh) const
    {
        return parameterInCellPerHexVertex[vh].at(ch);
    }

    void setHexVertexParameters(const CellHandle& ch, const std::pair<VertexHandle, VertexHandle>& vh_range, const std::vector<Parameter>& us);

    VertexHandle getHexVertexInCellWithParameter(const CellHandle& ch, const Parameter& u, const Config& config) const;

    // Parameter to World
    inline const Matrix4x4d& inverseMapping(const CellHandle& ch)
    {
        return igmsParam2World[ch].second? igmsParam2World[ch].first : (igmsParam2World[ch] = {calculateInverseMap(ch), true}).first;
    }

    // World to Parameter
    inline const Matrix4x4d& mapping(const CellHandle& ch)
    {
        return igmsWorld2Param[ch].second? igmsWorld2Param[ch].first : (igmsWorld2Param[ch] = {calculateMap(ch), true}).first;
    }

    inline bool edgeIsSingular(const EdgeHandle& eh)
    {
        return edgeSingularity[eh];
    }

private:
    TetMeshCache& cache;

public:
    OVM::CellPropertyT<CellIGM> tetVertexCellParameters;
    OVM::HalfFacePropertyT<Transition> transitions;


private:
    Matrix4x4d calculateMap(const CellHandle& ch);
    Matrix4x4d calculateInverseMap(const CellHandle& ch);

    //=========================
    // Transitions
    //=========================

    HalfFaceHandle rotateAroundHalfedge(const CellHandle& startCell, HalfEdgeHandle currentEdge, bool ccw = true);

    void doTransition(HalfFaceHandle hfh, CellHandle& ch);
    void doTransition(HalfFaceHandle hfh, Parameter& parameter);
    void doTransition(HalfFaceHandle hfh, Direction& dir);
    void doTransition(HalfFaceHandle hfh, Transition& tranFun);

    template <typename T, typename... Rest>
    void doTransition(HalfFaceHandle hfh, T& target, Rest&... rest)
    {
        doTransition(hfh, target);
        doTransition(hfh, rest...);
    }

private:


    OVM::EdgePropertyT<bool> edgeSingularity;
    OVM::EdgePropertyT<bool> edgeHasIncidentNonIdentityTransitions;

    OVM::CellPropertyT<std::pair<Matrix4x4d,bool>> igmsParam2World; // inverse and whether it has already been computed
    OVM::CellPropertyT<std::pair<Matrix4x4d,bool>> igmsWorld2Param;

    OVM::CellPropertyT<absl::flat_hash_map<Parameter, VertexHandle>> hexVertexPerParameterInCell;
    OVM::VertexPropertyT<FlatMap<CellHandle, Parameter>> parameterInCellPerHexVertex;
};

}

