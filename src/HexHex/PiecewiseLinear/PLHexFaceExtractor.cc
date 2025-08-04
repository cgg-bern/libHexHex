
#include <HexHex/PiecewiseLinear/PLMesh.hh>
#include <HexHex/PiecewiseLinear/halfface_traversal_check.hh>
#include <HexHex/Predicates/ExactPredicates.hh>
#include <HexHex/PiecewiseLinear/PLHexFaceExtractor.hh>
#include <HexHex/LocalTopology/HexVertexGenerator.hh>
#include <HexHex/HexExtractor.hh>

namespace HexHex
{

void PLHexFaceExtractor::computeReference(FaceHandle fh)
{
    const auto& hexMesh = he.getOutputMesh();
    //auto& tetMesh = he.inputMesh;

    // Get hex-halfedges of hex-face
    auto hfh = fh.halfface_handle(0);
    const auto& fh_hehs = hexMesh.halfface(hfh).halfedges();
    assert(fh_hehs.size() == 4);

    // Reference vertex is first hex-vertex, p1, p2 are propellers incident to that vertex
    // get reference tet and reference directions and integer grid face parameters
    VertexHandle refVh = hexMesh.from_vertex_handle(fh_hehs[0]);
    auto& gen = he.getHexVertexGenerator(refVh);
    auto ph1 = he.hexHalfEdgeLocalPropellers[fh_hehs[0]];
    auto ph2 = he.hexHalfEdgeLocalPropellers[fh_hehs[3].opposite_handle()];
    Direction refDir1(0), refDir2(0);
    gen.getDirections(ph1, ph2, reference.tet, refDir1, refDir2);

    // Parameter of corner vertex
    Parameter refCorner = he.getParametrization().getHexVertexParameter(reference.tet, refVh);

    // Axis is the one along which the face does not extend
    reference.axis = argMin(absVector(refDir1.vector() - refDir2.vector()));
    reference.axisValue = refCorner[reference.axis];

    // Sort corner points so first is minimum
    auto u = reference.parameter3dTo2d(refCorner);
    auto v = reference.parameter3dTo2d(refCorner + refDir1 + refDir2);
    auto uMin = Vec2d(std::min(u[0],v[0]), std::min(u[1],v[1]));
    reference.quad2d = {uMin, uMin+Vec2d(1,0), uMin+Vec2d(1,1), uMin+Vec2d(0,1)};
    for (unsigned char i = 0u; i < 4u; ++i) reference.quad3d[i] = reference.parameter2dTo3d(reference.quad2d[i]);

    assert(refCorner[reference.axis] == (refCorner+refDir1+refDir2)[reference.axis]);
}

void PLHexFaceExtractor::extract(FaceHandle fh)
{
    assert(unprocessed.empty() && cellsProps.empty() && mesh.points.empty());

    // Reset
    mesh.points.clear();
    mesh.polygons.clear();
    hex_fh = fh;

    // Compute the Reference System
    computeReference(fh);

    // Set Reference Properties
    cellsProps[reference.tet].transitionToRef = Transition::IDENTITY;
    cellsProps[reference.tet].distanceToRef = 0;

    // Start with reference tet and flood fill to get other tets that intersect the integer grid face
    unprocessed.insert(reference.tet);

    const auto &inputMesh = he.getInputMesh();
    //PiecewiseLinearHexFaceExtractor::nVisitedTets = 0;
    while (!unprocessed.empty()) {
        CellHandle ch = unprocessed.extract(unprocessed.begin()).value(); // pop

        assert(cellsProps.contains(ch));
        assert(!cellsProps[ch].visited);

        // Extract part of triangle mesh - intersection of tet ch with integer grid face
        uint nPts = extract(ch);
        if (nPts==0) {continue;}

        // Check adjacent tets
        for (const auto& hfh : inputMesh.cell(ch).halffaces()) {

            auto ch2 = inputMesh.incident_cell(hfh.opposite_handle());

            if (!ch2.is_valid() || cellsProps[ch2].visited || unprocessed.contains(ch2)) {continue;}

            // Visit adjacent tet if projection of triangle face in between intersects interior of integer grid face
            if (shouldVisitOppositeCellOverHalfface(he, hfh, cellsProps[ch].transitionToRef, reference.quad3d)) {
                // Set Transition and add adjacent tet to list to be visited
                cellsProps[ch2].transitionToRef = cellsProps[ch].transitionToRef
                                                * he.getParametrization().transition(hfh.opposite_handle());
                cellsProps[ch2].distanceToRef = cellsProps[ch].distanceToRef + 1;
                unprocessed.insert(ch2);
            }
            else
            {
            }
        }

    }

}

uint PLHexFaceExtractor::extract(const CellHandle& ch)
{
    assert(!cellsProps[ch].visited);
    cellsProps[ch].visited = true;

    const auto& param2World = he.getParametrization().inverseMapping(ch);
    const auto& ref2Tet = cellsProps.at(ch).transitionToRef.inverted();

    auto points = tetIntegerGridFaceIntersection(ch); // in reference system
    const uint n = points.size();

    // Only add at least triangles
    if (n <= 2) return n;

    // Add Poly Vertices
    auto worldPoints = std::vector<PLFPoint>();
    for (auto& pt : points)
    {
        // Transform from param reference to world
        pt.pos = ref2Tet.transform_point(pt.pos);
        pt.pos = param2World.transform_point(pt.pos); // world

        worldPoints.push_back(pt);
    }

    // Add Poly Faces
    mesh.addPolygon(ch, worldPoints);

    return n;
}

std::vector<PLFPoint> PLHexFaceExtractor::tetIntegerGridPlaneIntersection(const CellHandle& ch)
{
    std::vector<PLFPoint> points;
    points.reserve(4);

    auto trans2Ref = cellsProps.at(ch).transitionToRef;

    // Check which vertices lie on the plane
    for (auto vh : he.getCache().cellVertices[ch]) {
        auto param = trans2Ref.transform_point(he.getParametrization().get(ch, vh));

        if (param[reference.axis] == reference.axisValue)
        {
            points.push_back(PLFPoint{
                .pos = param,
                .generator = vh,
                .ple_vertex = VertexHandle(-1),
                .uv = (reference.parameter3dTo2d(param) - reference.quad2d[0])
            });
        }
    }

    const auto &inputMesh = he.getInputMesh();
    // Check which edges intersect the plane
    for (auto ce_it = inputMesh.ce_iter(ch); ce_it.is_valid(); ++ce_it) {
        auto eh = *ce_it;
        auto heh = eh.halfedge_handle(0);
        auto vh1 = inputMesh.from_vertex_handle(heh);
        auto vh2 = inputMesh.to_vertex_handle(heh);
        auto param1 = trans2Ref.transform_point(he.getParametrization().get(ch, vh1));
        auto param2 = trans2Ref.transform_point(he.getParametrization().get(ch, vh2));

        double p = reference.axisValue;
        double a = param1[reference.axis];
        double b = param2[reference.axis];

        if ((a<p&&p<b) || (b<p&&p<a)) {

            auto t = (p - a) / (b - a);
            auto param = param1 + t*(param2 - param1);
            Vec2d uv = (reference.parameter3dTo2d(param) - reference.quad2d[0]);
            points.push_back(PLFPoint{
                .pos = param,
                .generator = eh,
                .ple_vertex = VertexHandle(-1),
                .uv = uv
            });
        }

    }

    return points;
};

void PLHexFaceExtractor::sortPointsCCW(std::vector<PLFPoint>& res)
{
    // Sort the points ccw around centroid


    // Compute Centroid
    Vec2d centroid = Vec2d(0,0);
    for (const auto& p : res) centroid += reference.parameter3dTo2d(p.pos);
    centroid /= res.size();

    auto computeAngle = [&centroid](const Vec2d& p)
    {
        return std::atan2(p[1] - centroid[1], p[0] - centroid[0]);
    };

    std::sort(res.begin(), res.end(), [&](const PLFPoint& a, const PLFPoint& b)
  {
      return computeAngle(reference.parameter3dTo2d(a.pos)) < computeAngle(reference.parameter3dTo2d(b.pos));
  });
}

std::vector<PLFPoint> PLHexFaceExtractor::tetIntegerGridFaceIntersection(const CellHandle& ch)
{
    auto poly = tetIntegerGridPlaneIntersection(ch);

    auto &piecewiseLinearMesh = he.getPLMesh();

    // Find the (up to 8) piecewise linear edge vertices overlapping with the cell.
    std::vector<VertexHandle> segment_vhs = piecewiseLinearMesh.findPVerticesInCell(hex_fh, ch);
    assert(segment_vhs.size() <= 8);

    // No segments?
    if (segment_vhs.size() == 0)
    {

        // Either the entire polygon must be inside the quad or outside the quad.
        // Take the approx. majority
        uint n_outside = 0;
        for (const auto& poly_pt : poly)
        {
            if (!reference.isUVApproxInQuad(poly_pt.uv)) n_outside++;
        }

        if (poly.size() >= 3)
        {
            // Entire poly is outside
            if (n_outside >= poly.size()/2) return {};

            // Entire poly is inside, dont forget to order points
            sortPointsCCW(poly);
        }

        // Return Poly (might also be a single point ot edge)
        return poly;
    }

    // Get the map to map piecewise linear edge segments to parameter space
    const auto& world2Param = he.getParametrization().mapping(ch);

    // Get the clipline points in the reference system
    std::vector<PLFPoint> res;
    res.reserve(8);

    // Add the vertices that were already extracted in the edge or vertex extraction while keeping track of duplicates
    absl::flat_hash_set<VertexHandle> has_vertices;
    for (VertexHandle vh : segment_vhs)
    {
        assert(vh.is_valid());
        if (has_vertices.contains(vh)) continue;

        // Get World Position
        Vec3d pos = piecewiseLinearMesh.mesh.vertex(vh);

        // Transform World to tet Parameter
        pos = world2Param.transform_point(pos);

        // Transform to Reference Tet Parameter
        pos = cellsProps.at(ch).transitionToRef.transform_point(pos);

        // Add the point
        Vec2d uv = piecewiseLinearMesh.getVertexParameterUV(vh);
        res.push_back({
            .pos = pos,
            .generator = piecewiseLinearMesh.getVertexTetGenerator(vh),
            .ple_vertex = vh,
            .uv = uv
        });
        has_vertices.insert(vh);
    }

    // Add the polygon points inside the quad
    for (const auto& polyPoint : poly)
    {
        // Check if inside Quad: TODO: Evaluate with Exact Predicates
        if (reference.isUVApproxInQuad(polyPoint.uv))
        {
            // Check if the point is already on an edge segment (meaning, same generator)
            bool found_on_hex_face_boundary = false;
            for (const auto segment_vh : has_vertices)
            {
                if (piecewiseLinearMesh.getVertexTetGenerator(segment_vh) == polyPoint.generator)
                {
                    found_on_hex_face_boundary = true;
                    break;
                }
            }
            if (!found_on_hex_face_boundary)
            {
                res.push_back(polyPoint);
            }
        }
    }

    // If we have 2 or less points, there is no need to sort them
    const uint8_t n = res.size();
    if (n <= 2) return res;

    sortPointsCCW(res);

    return res;
}

} // namespace HexHex
