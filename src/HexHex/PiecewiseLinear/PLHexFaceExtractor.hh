#pragma once

#include <HexHex/Commons/MeshElement.hh>
#include <HexHex/Commons/Transition.hh>
#include <HexHex/Utils/Typedefs.hh>
#include <absl/container/flat_hash_set.h>
#include <absl/container/flat_hash_map.h>
#include <HexHex/PiecewiseLinear/PLFPoint.hh>

namespace HexHex
{

// Forward Declaration
class HexExtractor;


class PLHexFaceExtractor
{
public:

    struct PolyMesh
    {
        std::vector<PLFPoint> points;
        std::vector<std::vector<uint32>> polygons;
        std::vector<int> tet_per_polygon;
        absl::flat_hash_map<MeshElement, uint32> generator_poly_vertices;

        uint32 addPoint(const PLFPoint& point)
        {
            // Check if we have poly point already (poly point with same generator)
            const bool is_poly_point = !point.ple_vertex.is_valid();
            if (is_poly_point && generator_poly_vertices.contains(point.generator))
            {
                uint32 idx = generator_poly_vertices[point.generator];
                return idx;
            }

            // Add new vertex
            points.push_back(point);
            uint32 idx = points.size()-1;
            if (is_poly_point) {generator_poly_vertices[point.generator] = idx;}
            return idx;
        }

        void addPolygon(CellHandle tet, const std::vector<PLFPoint>& polyPts)
        {
            uint poly_idx = polygons.size();
            polygons.push_back({});
            tet_per_polygon.push_back(tet.idx());
            polygons[poly_idx].reserve(polyPts.size());
            for (const auto& polyPt : polyPts)
            {
                uint32 pt_idx = addPoint(polyPt);
                polygons[poly_idx].push_back(pt_idx);
            }
        }
    };

public:
    explicit PLHexFaceExtractor(HexExtractor& he)
        :
        he(he)
    {
    }

    void extract(FaceHandle fh);

    inline PolyMesh& getHexFaceMesh() {return mesh;}

private:
    HexExtractor& he;
    PolyMesh mesh;
    FaceHandle hex_fh;

    struct CellProps
    {
        Transition transitionToRef; // transition from tet to reference system
        int distanceToRef; // DEBUG: how many transitions from the reference tet to this tet
        bool visited; // computed the integer grid face intersection
    };

    // Our Calculations are done with respect to some reference coordinate system.
    struct ReferenceSystem
    {
        CellHandle tet; // Reference Tet
        std::array<Vec2d, 4> quad2d; // Parameter Corners in 2d
        std::array<Vec3d, 4> quad3d; // Parameter Corners in 3d
        uint8_t axis; // Axis which with which the Quad is aligned (0, 1 or 2)
        int32 axisValue; // The height of the Quad (i.e. = quad3d[i][axis])

        // Returns the given 2d point on the integer grid plane (w.r.t. the reference system) as a 3d parameter.
        inline Parameter parameter2dTo3d(const Vec2d& u) const
        {
            Parameter u3d;
            u3d[axis] = axisValue;
            u3d[(axis+1)%3] = u[0];
            u3d[(axis+2)%3] = u[1];
            return u3d;
        }

        // Projects a 3d parameter w.r.t. the reference sytem onto the integer grid plane and returns the 2d projection.
        inline Vec2d parameter3dTo2d(const Parameter& u) const
        {
            return Vec2d(u[(axis+1)%3], u[(axis+2)%3]);
        }

        inline bool isUVApproxInQuad(const Vec2d& uv) const
        {
            constexpr double EPS = 1e-9;
            return (uv[0] >= -EPS && uv[1] >= -EPS && uv[0] <= 1+EPS && uv[1] <= 1+EPS);
        }

    } reference;
    void computeReference(FaceHandle fh);


    absl::flat_hash_set<CellHandle> unprocessed; // tet cells that are not yet visited but should be checked
    absl::flat_hash_map<CellHandle, CellProps> cellsProps;

    uint extract(const CellHandle& ch);

    void sortPointsCCW(std::vector<PLFPoint>& res);

    // Computes the intersection points of the infinite integer grid plane with the given tet and
    // returns the intersection points w.r.t. the reference system ordered ccw.
    std::vector<PLFPoint> tetIntegerGridPlaneIntersection(const CellHandle& ch);

    // Computes the intersection points of the integer grid face with the given tet and
    // returns the intersection points w.r.t. the reference system ordered ccw.
    std::vector<PLFPoint> tetIntegerGridFaceIntersection(const CellHandle& ch);
};
}

