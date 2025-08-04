#include <HexHex/Preprocessing/Parametrization.hh>
#include <HexHex/Utils/TetMeshCache.hh>
#include <HexHex/Commons/Direction.hh>
#include <HexHex/HexExtractor.hh>

namespace HexHex
{

Parametrization::Parametrization(TetMeshCache &cache, HexahedralMesh &hexmesh) :
    cache(cache),

    tetVertexCellParameters(cache.tetmesh.create_private_cell_property<CellIGM>()),

    transitions(cache.tetmesh.create_private_halfface_property<Transition>()),

    edgeSingularity(cache.tetmesh.create_private_edge_property<bool>()),

    edgeHasIncidentNonIdentityTransitions(cache.tetmesh.create_private_edge_property<bool>()),

    igmsParam2World(cache.tetmesh.create_private_cell_property<std::pair<Matrix4x4d,bool>>()),
    igmsWorld2Param(cache.tetmesh.create_private_cell_property<std::pair<Matrix4x4d,bool>>()),

    hexVertexPerParameterInCell(cache.tetmesh.create_private_cell_property<absl::flat_hash_map<Parameter, VertexHandle>>()),
    parameterInCellPerHexVertex(hexmesh.request_vertex_property<FlatMap<CellHandle, Parameter>>())
{}

Matrix4x4d Parametrization::calculateMap(const CellHandle& ch)
{
    const auto& vhs = cache.cellVertices[ch];

    const auto& p = cache.tetmesh.vertex(vhs[0]);
    const auto& q = cache.tetmesh.vertex(vhs[1]);
    const auto& r = cache.tetmesh.vertex(vhs[2]);
    const auto& s = cache.tetmesh.vertex(vhs[3]);

    const auto& u = get(ch, vhs[0]);
    const auto& v = get(ch, vhs[1]);
    const auto& w = get(ch, vhs[2]);
    const auto& t = get(ch, vhs[3]);

    Matrix4x4d mat1;
    mat1(0,0) = q[0]-p[0]; mat1(0,1) = r[0]-p[0]; mat1(0,2) = s[0]-p[0]; mat1(0,3) = p[0];
    mat1(1,0) = q[1]-p[1]; mat1(1,1) = r[1]-p[1]; mat1(1,2) = s[1]-p[1]; mat1(1,3) = p[1];
    mat1(2,0) = q[2]-p[2]; mat1(2,1) = r[2]-p[2]; mat1(2,2) = s[2]-p[2]; mat1(2,3) = p[2];
    mat1(3,0) = 0;         mat1(3,1) = 0;         mat1(3,2) = 0;         mat1(3,3) = 1;

    Matrix4x4d mat2;
    mat2(0,0) = v[0]-u[0]; mat2(0,1) = w[0]-u[0]; mat2(0,2) = t[0]-u[0]; mat2(0,3) = u[0];
    mat2(1,0) = v[1]-u[1]; mat2(1,1) = w[1]-u[1]; mat2(1,2) = t[1]-u[1]; mat2(1,3) = u[1];
    mat2(2,0) = v[2]-u[2]; mat2(2,1) = w[2]-u[2]; mat2(2,2) = t[2]-u[2]; mat2(2,3) = u[2];
    mat2(3,0) = 0;         mat2(3,1) = 0;         mat2(3,2) = 0;         mat2(3,3) = 1;

    mat1.invert();

    return mat2*mat1;
}

Matrix4x4d Parametrization::calculateInverseMap(const CellHandle& ch)
{
    const auto& vhs = cache.cellVertices[ch];// getCellVertices(ch);

    const auto& p = cache.tetmesh.vertex(vhs[0]);
    const auto& q = cache.tetmesh.vertex(vhs[1]);
    const auto& r = cache.tetmesh.vertex(vhs[2]);
    const auto& s = cache.tetmesh.vertex(vhs[3]);

    const auto& u = get(ch, vhs[0]);
    const auto& v = get(ch, vhs[1]);
    const auto& w = get(ch, vhs[2]);
    const auto& t = get(ch, vhs[3]);

    Matrix4x4d mat1;
    mat1(0,0) = p[0]; mat1(0,1) = q[0]; mat1(0,2) = r[0]; mat1(0,3) = s[0];
    mat1(1,0) = p[1]; mat1(1,1) = q[1]; mat1(1,2) = r[1]; mat1(1,3) = s[1];
    mat1(2,0) = p[2]; mat1(2,1) = q[2]; mat1(2,2) = r[2]; mat1(2,3) = s[2];
    mat1(3,0) = 1;    mat1(3,1) = 1;    mat1(3,2) = 1;    mat1(3,3) = 1;

    Matrix4x4d mat2;
    mat2(0,0) = u[0]; mat2(0,1) = v[0]; mat2(0,2) = w[0]; mat2(0,3) = t[0];
    mat2(1,0) = u[1]; mat2(1,1) = v[1]; mat2(1,2) = w[1]; mat2(1,3) = t[1];
    mat2(2,0) = u[2]; mat2(2,1) = v[2]; mat2(2,2) = w[2]; mat2(2,3) = t[2];
    mat2(3,0) = 1;    mat2(3,1) = 1;    mat2(3,2) = 1;    mat2(3,3) = 1;

    mat2.invert();

    return mat1*mat2;
}

void Parametrization::setHexVertexParameters(const CellHandle& ch, const std::pair<VertexHandle, VertexHandle>& vh_range, const std::vector<Parameter> &us)
{
    size_t n = us.size();
    assert(static_cast<size_t>(vh_range.second.idx() - vh_range.first.idx() + 1) == n);
    for (auto i = 0u; i < n; ++i) {
        auto vh = VertexHandle(vh_range.first.idx()+i);

        parameterInCellPerHexVertex[vh].set(ch, us[i]);
        hexVertexPerParameterInCell[ch][us[i]] = vh;
    }
}

VertexHandle Parametrization::getHexVertexInCellWithParameter(const CellHandle& ch, const Parameter &u, const Config &config) const
{
#if 0
    // Search for Parameter with linear search (only for experimental purposes)
    for (const auto& m : hexVertexPerParameterInCell[ch])
    {
        if (m.first == u)
        {
            return m.second;
        }
    }
    return VertexHandle(-1);
#else
    // Search for Parameter via map
    if (hexVertexPerParameterInCell[ch].contains(u)) {
        return hexVertexPerParameterInCell[ch].at(u);
    }
    return VertexHandle(-1);
#endif
}

HalfFaceHandle Parametrization::rotateAroundHalfedge(const CellHandle& startCell, HalfEdgeHandle currentEdge, bool ccw)
{
    if (ccw) currentEdge = cache.tetmesh.opposite_halfedge_handle(currentEdge);
    for (auto hehf_it = cache.tetmesh.hehf_iter(currentEdge); hehf_it.valid(); ++hehf_it)
        if (startCell == cache.tetmesh.incident_cell(*hehf_it))
            return *hehf_it;
    return HalfFaceHandle();
}

Transition Parametrization::getTransition(const CellHandle& fromCell, const CellHandle& toCell, const EdgeHandle& eh)
{
    assert(MeshElement(fromCell).is_incident(cache, MeshElement(eh)) && MeshElement(toCell).is_incident(cache, MeshElement(eh)));
    if  (!edgeHasIncidentNonIdentityTransitions[eh] || fromCell == toCell) return Transition::IDENTITY;

    auto heh = cache.tetmesh.halfedge_handle(eh, 0);
    if (cache.tetmesh.is_boundary(eh))
    {
        for (auto ccw : {true, false})
        {
            auto tranFun = Transition::IDENTITY;
            auto currentCell = fromCell;
            while ((currentCell != toCell) && currentCell.is_valid())
            {
                auto transitionFace =  rotateAroundHalfedge(currentCell, heh, ccw);
                doTransition(transitionFace, currentCell, tranFun);
            }
            if (currentCell.is_valid())
                return tranFun;
        }
        assert(false); // multiple connected components of incident tets
        return Transition::IDENTITY;
    }
    else
    {
        auto tranFun = Transition::IDENTITY;
        auto currentCell = fromCell;
        while ((currentCell != toCell))
        {
            auto transitionFace =  rotateAroundHalfedge(currentCell, heh, true);
            doTransition(transitionFace, currentCell, tranFun);
        }
        return tranFun;
    }
}

Transition Parametrization::getTransition(const CellHandle& fromCell, const CellHandle& toCell)
{
    for (const auto& hfh : cache.tetmesh.cell(fromCell).halffaces())
        if (cache.tetmesh.incident_cell(hfh.opposite_handle()) == toCell)
            return transitions[hfh];
    throw std::runtime_error("Can't get transition between non-adjacent cells!");
    return Transition::IDENTITY;
}

void Parametrization::doTransition(HalfFaceHandle hfh, CellHandle& ch)
{
    ch = cache.tetmesh.incident_cell(cache.tetmesh.opposite_halfface_handle(hfh));
}

void Parametrization::doTransition(HalfFaceHandle hfh, Parameter& parameter)
{
    parameter = transition(hfh).transform_point(parameter);
}

void Parametrization::doTransition(HalfFaceHandle hfh, Direction& dir)
{
    dir = transition(hfh).transform(dir);
}

void Parametrization::doTransition(HalfFaceHandle hfh, Transition& tranFun)
{
    tranFun = transition(hfh) * tranFun;
}

}
