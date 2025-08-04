#include <HexHex/PiecewiseLinear/PLMesh.hh>
#include <HexHex/Utils/TetMeshCache.hh>
#include <HexHex/GlobalTopology/HexEdgeExtractor.hh>
#include <HexHex/PiecewiseLinear/PLHexFaceExtractor.hh>

namespace HexHex {

PLMesh::PLMesh(TetMeshCache& cache, HexahedralMesh& hexmesh) :
    cache(cache),
    hexmesh(hexmesh),
    mesh(),

    linear_hex_edge_indices(mesh.create_persistent_edge_property<int>("e_linear_hex_edge", -1).value()),
    linear_hex_face_indices(mesh.create_persistent_face_property<int>("f_linear_hex_face", -1).value()),

    vertex_tet_generator_indices(mesh.create_persistent_vertex_property<int>("v_tet_generator_index", -1).value()),
    vertex_tet_generator_types(mesh.create_persistent_vertex_property<int>("v_tet_generator_type", -1).value()),

    vertex_types(mesh.create_persistent_vertex_property<int>("v_type").value()),

    uv_prop(mesh.create_persistent_vertex_property<Vec2d>("v_igm_uv").value()),

    edge_tet_generator_indices(mesh.create_persistent_edge_property<int>("e_tet_generator_index", -1).value()),
    edge_tet_generator_types(mesh.create_persistent_edge_property<int>("e_tet_generator_type", -1).value()),

    face_tet_generator_indices(mesh.create_persistent_face_property<int>("f_intersected_tet_index").value()),

    hex_edge_ple_segments(hexmesh.request_edge_property<absl::flat_hash_map<CellHandle, MeshElement>>(""))
{
    // Some speedup
    mesh.enable_bottom_up_incidences(false);
    mesh.enable_vertex_bottom_up_incidences(true);
}


VertexHandle PLMesh::addPLEVertex(const EdgeHandle hex_eh, const PLEPoint& point) {
    VertexHandle vh = addVertex(point.pos, point.generator, Vec2d(point.t, 0), PLVertexType::ARC);

    // Remember Vertex in incident tets
    for (CellHandle incident_tet : point.generator.get_incident_cells(cache))
    {
        if (!hex_edge_ple_segments[hex_eh][incident_tet].is_valid())
        {
            hex_edge_ple_segments[hex_eh][incident_tet] = vh;
        }
    }

    return vh;
}

VertexHandle PLMesh::addPLFVertex(const PLFPoint& point)
{
    if (point.ple_vertex.is_valid()) return point.ple_vertex;
    return addVertex(point.pos, point.generator, point.uv, PLVertexType::PATCH);
}
VertexHandle PLMesh::addVertex(const Vec3d& pos, const MeshElement& generator, const Vec2d &uv, PLVertexType type)
{
    assert(generator.is_valid());
    VertexHandle vh = mesh.add_vertex(pos);
    vertex_tet_generator_indices[vh] = generator.idx();
    vertex_tet_generator_types[vh] = generator.type();
    vertex_types[vh] = type;
    uv_prop[vh] = uv;
    return vh;
}
MeshElement PLMesh::getVertexTetGenerator(VertexHandle vh) const
{
    if (vertex_tet_generator_types[vh]==MeshElementType::ELEMENT_VERTEX_TYPE)
        return VertexHandle(vertex_tet_generator_indices[vh]);
    if (vertex_tet_generator_types[vh]==MeshElementType::ELEMENT_EDGE_TYPE)
        return HalfEdgeHandle(vertex_tet_generator_indices[vh]);
    if (vertex_tet_generator_types[vh]==MeshElementType::ELEMENT_FACE_TYPE)
        return HalfFaceHandle(vertex_tet_generator_indices[vh]);
    if (vertex_tet_generator_types[vh]==MeshElementType::ELEMENT_CELL_TYPE)
        return CellHandle(vertex_tet_generator_indices[vh]);
    return MeshElement();
}
EdgeHandle PLMesh::addPLEdgeSegment(const EdgeHandle hex_eh,
                                    const VertexHandle vh1, const VertexHandle vh2,
                                    const CellHandle ch)
{
    EdgeHandle eh = mesh.add_edge(vh1, vh2);
    linear_hex_edge_indices[eh] = hex_eh.idx();

    // Get the Generator of the Edge Segment
    MeshElement eh_gen = getElementBetweenInCell(cache, ch,
        getVertexTetGenerator(vh1), getVertexTetGenerator(vh2));

    edge_tet_generator_types[eh] = eh_gen.type();
    edge_tet_generator_indices[eh] = eh_gen.idx();

    // Remember Edge in incident tets
    auto incident_chs = eh_gen.get_incident_cells(cache);
    for (auto incident_ch : incident_chs)
    {
        hex_edge_ple_segments[hex_eh][incident_ch] = eh;
    }

    return eh;
}

FaceHandle PLMesh::addPLFacePatch(const FaceHandle hex_fh, const std::vector<VertexHandle>& vhs, int gen_tet_idx)
{
    FaceHandle fh = mesh.add_face(vhs);

    linear_hex_face_indices[fh] = hex_fh.idx();
    face_tet_generator_indices[fh] = gen_tet_idx;

    return fh;
}

std::vector<VertexHandle> PLMesh::findPVerticesInCell(const FaceHandle hex_fh, const CellHandle tet_ch)
{
    assert(cache.tetmesh.is_valid(tet_ch));

    std::vector<VertexHandle> res;
    res.reserve(8);

    // For each of the 4 hex edges
    for (HalfEdgeHandle hex_heh : hexmesh.halfface(hex_fh.halfface_handle(0)).halfedges())
    {
        if (hex_edge_ple_segments[hex_heh.edge_handle()].contains(tet_ch))
        {
            MeshElement segEnt = hex_edge_ple_segments[hex_heh.edge_handle()][tet_ch];
            if (segEnt.is_edge())
            {
                EdgeHandle eh = segEnt.eh();
                VertexHandle vh1 = mesh.from_vertex_handle(eh.halfedge_handle(0));
                VertexHandle vh2 = mesh.from_vertex_handle(eh.halfedge_handle(1));
                res.push_back(vh1);
                res.push_back(vh2);
            } else if (segEnt.is_vertex())
            {
                res.push_back(segEnt.vh());
            }
        }
    }
    return res;
}


} // namespace HexHex
