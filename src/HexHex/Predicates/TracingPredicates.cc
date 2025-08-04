#include <HexHex/Commons/Direction.hh>
#include <HexHex/Predicates/ExactPredicates.hh>
#include <HexHex/Predicates/GeneratorPredicates.hh>
#include <HexHex/Predicates/TracingPredicates.hh>
#include <HexHex/Utils/Utils.hh>

HexHex::HalfFaceHandle HexHex::findOppositeNextHalfFace (
    const TetMeshCache& mesh,
    const CellHandle ch, const std::array<VertexHandle,4>& cell_vhs, const std::array<Parameter,4>& cell_params,
    const HalfFaceHandle hfh, const Parameter& u,
    const Direction d1, const Direction d2, MeshElement& traceThrough, bool computeTraceThrough)
{

    assert(!hfh.is_valid() || mesh.tetmesh.incident_cell(hfh) == ch);
    assert((d1 | d2) == 0);

    assert(!computeVertexGeneratorElement(mesh, ch, cell_vhs, cell_params, u+d1).is_valid());

    auto VERTEX_PARAM = [&](VertexHandle vh) -> const Parameter&
    {
        for (int i = 0; i < 4; ++i)
            if (cell_vhs[i] == vh)
                return cell_params[i];
        assert(false);
        throw std::runtime_error("IMPOSSIBLE VERTEX");
    };

    auto VERTEX_PARAMS = [&](VertexHandle vh1, VertexHandle vh2, VertexHandle vh3) -> std::array<Parameter,3>
    {
        return {VERTEX_PARAM(vh1),VERTEX_PARAM(vh2),VERTEX_PARAM(vh3)};
    };

    for (const auto& checkHfh : mesh.tetmesh.cell(ch).halffaces()) {
        if (checkHfh == hfh) continue;

        const auto& vertices = mesh.get_halfface_vertices(checkHfh);
        auto params = VERTEX_PARAMS(vertices[0], vertices[1], vertices[2]);

        // Primary direction must cut through the plane given by the face
        if (ori3d(params[0], params[1], params[2], u) != ORI_ABOVE || ori3d(params[0], params[1], params[2], u + d1) != ORI_BELOW)
            continue;

        // Get the orientations of the 3 tets formed around u, u+d1
        std::array<ORIENTATION, 3> oris;
        oris[0] = ori3d(params[0], params[1], u, u + d1);
        oris[1] = ori3d(params[1], params[2], u, u + d1);
        oris[2] = ori3d(params[2], params[0], u, u + d1);

        // Check if primary direction cuts through interior/center of triangle
        if (oris[0] == oris[1] && oris[0] == oris[2]) {
            assert(oris[0] == ORI_BELOW);
            if (computeTraceThrough) {traceThrough = checkHfh.face_handle();}
            return checkHfh; // case 1
        }

        // Check if primary direction does not cut through triangle at all
        //if (!std::any_of(oris.begin(), oris.end(), [](ORIENTATION ori){return ori == ORI_ZERO;}) {
        if (std::any_of(oris.begin(), oris.end(), [](ORIENTATION ori){return ori == ORI_BELOW;})
            && std::any_of(oris.begin(), oris.end(), [](ORIENTATION ori){return ori == ORI_ABOVE;})) {

            continue;
        }

        // Otherwise the primary direction cuts through the triangles boundary (edge or vertex). First get all edge intersection (amount must be 1 or 2)
        std::vector<int> edgeIntersections;
        for (unsigned char i = 0; i < oris.size(); ++i) if (oris[i] == ORI_ZERO) edgeIntersections.push_back(i);

        if (edgeIntersections.size() == 1) {  // Pointing through edge - case 2

            // Get vertices in order s.t. first and second correspond to the intersected edge
            auto index = edgeIntersections[0];
            auto from_vertex = vertices[index];
            auto to_vertex = vertices[(index+1)%3];
            auto other_vertex = vertices[(index+2)%3];
            params = VERTEX_PARAMS(from_vertex, to_vertex, other_vertex);
            if (computeTraceThrough)
            {traceThrough = mesh.tetmesh.find_halfedge_in_cell(from_vertex, to_vertex, ch).edge_handle();}

            // Get on which side of the plane (from_vertex, to_vertex, u + d1) u + d2 is.
            // Either it's on the plane in which case we have two halffaces that we can pick (doesn't matter which but might as well pick the current one, ey?)
            // Or it is left or right in which case it has to coinside with the direction other_vertex is on.
            auto ori = ori3d(params[0], params[1], u + d1, u + d2);

            if (ori == ORI_ZERO) {
                return checkHfh; // case 2b
            }

            if (ori == ori3d(params[0], params[1], u + d1, params[2]))
                return checkHfh; // case 2a
            else
                return mesh.get_other_incident_halfface_in_cell(
                    mesh.tetmesh.find_halfedge_in_cell(from_vertex, to_vertex, ch).edge_handle(),
                    checkHfh); // case 2a (pick the one halfface which IS on the correct side)

        } else { // Pointing through vertex - case 3
            assert(edgeIntersections.size() == 2);

            // Get vertices s.t. first vertex is the intersected one
            auto index = 0;
            if (edgeIntersections[0] == 1) index = 2; else if (edgeIntersections[1] == 1) index = 1;

            auto intersectedVertex = vertices[index];
            auto oneOtherVertex = vertices[(index+1)%3];
            auto otherOtherVertex = vertices[(index+2)%3];
            params = VERTEX_PARAMS(intersectedVertex, oneOtherVertex, otherOtherVertex);
            if (computeTraceThrough) {traceThrough = intersectedVertex;}

            // the secondary direction must be on the correct side
            if (ori3d(params[1], params[0], u, u + d2) == ORI_BELOW || ori3d(params[0], params[2], u, u + d2) == ORI_BELOW)
                continue;

            return checkHfh;
        }
    }

    //std::cerr << std::hexfloat;
    std::cerr << "Failed to trace through face:" << std::endl;
    std::cerr << "From " << u << " to " << (u+d1) << std::endl;
    std::cerr << "Cell params: " << ch << std::endl;
    std::cerr << cell_params[0] << std::endl;
    std::cerr << cell_params[1] << std::endl;
    std::cerr << cell_params[2] << std::endl;
    std::cerr << cell_params[3] << std::endl;
    assert(false);
    return HalfFaceHandle(-1); // fail :(

}

HexHex::HalfFaceHandle HexHex::findBladeNextHalfFace (
    const TetMeshCache& mesh,
    const CellHandle ch, const std::array<VertexHandle,4>& cell_vhs, const std::array<Parameter,4>& cell_params,
    const MeshElement& generator, const MeshElement& holder,
    const HalfFaceHandle hfh, const Parameter& u,
    const Direction d1, const Direction d2)
{
    assert(!hfh.is_valid() || mesh.tetmesh.incident_cell(hfh) == ch);
    assert((d1 | d2) == 0);

    auto VERTEX_PARAM = [&](VertexHandle vh) -> const Parameter&
    {
        for (int i = 0; i < 4; ++i)
            if (cell_vhs[i] == vh)
                return cell_params[i];
        assert(false);
        throw std::runtime_error("IMPOSSIBLE VERTEX");
    };

    auto VERTEX_PARAMS = [&](VertexHandle vh1, VertexHandle vh2, VertexHandle vh3) -> std::array<Parameter,3>
    {
        return {VERTEX_PARAM(vh1),VERTEX_PARAM(vh2),VERTEX_PARAM(vh3)};
    };

    auto HALFFACE_PARAMS = [&](HalfFaceHandle hfh) -> std::array<Parameter,3>
    {
        auto vhs = mesh.get_halfface_vertices(hfh);
        return {VERTEX_PARAM(vhs[0]),VERTEX_PARAM(vhs[1]),VERTEX_PARAM(vhs[2])};
    };

    // For a vertex extracted on a face there is only one other cell to rotate into
    if (generator.is_face()) return mesh.get_incident_halfface_in_cell(ch, generator.fh());

    // There is also only one other option when not in the first cell and rotating around an edge
    // -> The only other incident face to the edge (that is not equal to the entering face)
    if (hfh.is_valid() && generator.is_edge()) {return mesh.get_other_incident_halfface_in_cell(generator.eh(), hfh);}

    for (auto checkHfh : mesh.tetmesh.cell(ch).halffaces()) {
        if (checkHfh == hfh) continue;

        if (generator.is_edge()) {
            assert(!hfh.is_valid());
            if (holder.is_face()) {
                // Rotating from edge from holder halfface to other incident halfface
                return mesh.get_other_incident_halfface_in_cell(
                        generator.eh(),
                        mesh.get_incident_halfface_in_cell(ch, holder.fh()));
            } else {
                // Rotating from edge from cell
                assert(holder.is_cell());

                // Halfface must be incident to edge holder
                if (!mesh.is_incident(generator.eh(), checkHfh.face_handle())) continue;

                auto e = mesh.tetmesh.edge(generator.eh());
                std::array<VertexHandle, 3> vhs = {
                    e.from_vertex(),
                    e.to_vertex(),
                    mesh.get_nonincident_vertex_in_face(checkHfh.face_handle(), generator.eh())
                };
                std::array<Parameter, 3> params = VERTEX_PARAMS(vhs[0],vhs[1],vhs[2]);

                auto dir2Side = ori3d(params[0], params[1], u + d1, u + d2);
                if (dir2Side == ORI_ZERO) return checkHfh; // u + d1 -> u + d2 goes through generator edge

                auto faceSide = ori3d(params[0], params[1], u + d1, params[2]);

                if (faceSide == dir2Side) {return checkHfh;}
                else {return mesh.get_other_incident_halfface_in_cell(generator.eh(), checkHfh);}
            }
        } else {
            // Rotating around vertex
            assert(generator.is_vertex());

            auto vh = generator.vh();
            assert(u == VERTEX_PARAM(vh));

            auto vertices = mesh.get_halfface_vertices(checkHfh);
            auto params = HALFFACE_PARAMS(checkHfh);

            // u + d1 must be strictly above the face and u + d2 must be strictly below
            if (ori3d(params, u + d1) != ORI_ABOVE) continue;
            if (ori3d(params, u + d2) != ORI_BELOW) continue;

            // Get v, w s.t. (u,v,w) keep the face orientation and u is the parameter of this vertex
            Parameter v, w;
            bool faceIsIncidentToVertex = false;
            for (auto i = 0u; i < 3; ++i) {
                if (vertices[i] == vh) {
                    faceIsIncidentToVertex = true;
                    v = params[(i+1)%3];
                    w = params[(i+2)%3];
                    break;
                }
            }

            // Don't pick the halfface if it is the one not incident to the generator vertex
            if (!faceIsIncidentToVertex) continue;

            auto ori1 = ori3d(u, v, u + d1, u + d2);
            auto ori2 = ori3d(w, u, u + d1, u + d2);

            if (ori1 == ORI_BELOW && ori2 == ORI_BELOW) return checkHfh; // Rotation goes through inner part of triangle
            if (ori1 == ORI_ZERO && ori2 == ORI_BELOW) return checkHfh; // Rotation goes through edge uv
            if (ori1 == ORI_BELOW && ori2 == ORI_ZERO) return checkHfh; // Rotation goes through edge wu
        }
    }

    // Fail
    return HalfFaceHandle(-1);
}
