#pragma once

#include <HexHex/Commons/MeshElement.hh>
#include <HexHex/Utils/Typedefs.hh>
#include <HexHex/Utils/Utils.hh>
#include <HexHex/PiecewiseLinear/PLEPoint.hh>
#include <HexHex/PiecewiseLinear/PLFPoint.hh>
#include <absl/container/flat_hash_map.h>


namespace HexHex
{

enum PLVertexType : int
{
    VERTEX, ARC, PATCH
};

struct PLMesh
{
    PLMesh(TetMeshCache& cache, HexahedralMesh& hexmesh);
    TetMeshCache& cache;
    HexahedralMesh& hexmesh;
    PolyhedralMesh mesh;

    // Piecewise linear mesh props
    OVM::EdgePropertyT<int> linear_hex_edge_indices; // -1 if inner part of face patch
    OVM::FacePropertyT<int> linear_hex_face_indices;

    OVM::VertexPropertyT<int> vertex_tet_generator_indices;
    OVM::VertexPropertyT<int> vertex_tet_generator_types; //0=vertex, 1=edge, 2=face, 3=cell

    OVM::VertexPropertyT<int> vertex_types;

    /*
     * Linear Hex-Vertex -> uv = (0,0)
     * PLE Vertex -> v = 0, u = w.r.t. first hex-halfedge
     * PLF Vertex -> uv = w.r.t. first hex-halfface
     */
    OVM::VertexPropertyT<Vec2d> uv_prop;

    OVM::EdgePropertyT<int> edge_tet_generator_indices;
    OVM::EdgePropertyT<int> edge_tet_generator_types; // -1=none, 1=edge, 2=face, 3=cell

    OVM::FacePropertyT<int> face_tet_generator_indices; // Tet index of the tet, the patch intersects

    // Hexmesh props

    // Per Tet, remember the pl edge segments or vertices extracted there
    OVM::EdgePropertyT<absl::flat_hash_map<CellHandle, MeshElement>> hex_edge_ple_segments;

    // Add Vertex on Edge Arc
    VertexHandle addPLEVertex(const EdgeHandle hex_eh, const PLEPoint& point);

    // Add Vertex within Face Patch
    VertexHandle addPLFVertex(const PLFPoint& point);

    VertexHandle addVertex(const Vec3d& pos, const MeshElement& generator, const Vec2d &uv = Vec2d(0,0), PLVertexType type = VERTEX);
    MeshElement getVertexTetGenerator(VertexHandle vh) const;

    inline const Vec2d& getVertexParameterUV(VertexHandle vh) const
    {
        return uv_prop[vh];
    }

    EdgeHandle addPLEdgeSegment(const EdgeHandle hex_eh,
                                       const VertexHandle vh1, const VertexHandle vh2,
                                       const CellHandle ch);

    FaceHandle addPLFacePatch(const FaceHandle hex_fh, const std::vector<VertexHandle>& vhs, int gen_tet_idx);

    std::vector<VertexHandle> findPVerticesInCell(const FaceHandle hex_fh, const CellHandle tet_ch);

    template <typename PolyMeshT>
    void copyMesh(PolyMeshT& polyMesh)
    {
        auto tovec = toVec<typename PolyMeshT::PointT>;

        polyMesh.clear(false);

        // Add vertices
        for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
            polyMesh.add_vertex(tovec(mesh.vertex(*v_it)));

        // Add Edges
        for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
            auto vh1 = mesh.edge(*e_it).from_vertex();
            auto vh2 = mesh.edge(*e_it).to_vertex();
            polyMesh.add_edge(vh1, vh2);
        }

        // Add Faces
        for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
            polyMesh.add_face(mesh.face(*f_it).halfedges());;
        }

        // Copy properties
        copy_edge_property<int>(mesh, polyMesh, "hex_eh");
        copy_face_property<int>(mesh, polyMesh, "hex_fh");
    }
};

}

