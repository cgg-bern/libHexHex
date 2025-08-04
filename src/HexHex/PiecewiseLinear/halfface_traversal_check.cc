#include <HexHex/PiecewiseLinear/halfface_traversal_check.hh>
#include <HexHex/HexExtractor.hh>
#include <HexHex/Predicates/ExactPredicates.hh>

namespace HexHex
{

bool shouldVisitOppositeCellOverHalfface(HexExtractor& he, const HalfFaceHandle& hfh,
                            const Transition& trans2Ref, const std::array<Parameter, 4>& refCorners3d)
{
    using BoundingBox2d = std::pair<Vec2d, Vec2d>;
    const auto& ch = he.getInputMesh().incident_cell(hfh);
    const auto& vhs = he.getInputMesh().get_halfface_vertices(hfh);

    // Get halfface vertex parameters in reference system
    auto T = std::array<Parameter,3>();
    for (unsigned char i = 0u; i < 3u; ++i) {
        auto param = he.getParametrization().get(ch, vhs[i]);
        param = trans2Ref.transform_point(param);
        T[i] = param;
    }

    // Does Triangle intersect quad interior?
    // (The case where the triangle and quad are parallel is excluded
    // as there are other cell faces which then satisfy the predicate)
    {
        const auto& Q = refCorners3d;
        const char c3d = (Q[0][0]==Q[1][0]&&Q[0][0]==Q[2][0])? 0 : (Q[0][1]==Q[1][1]&&Q[0][1]==Q[2][1])? 1 : 2;

        auto computeTriangleBoundingBox = [](const std::array<Vec2d,3>& vs) -> BoundingBox2d
        {
            Vec2d min, max;
            for (unsigned char i = 0u; i < 2u; ++i)
            {
                min[i] = min3(vs[0][i], vs[1][i], vs[2][i]);
                max[i] = max3(vs[0][i], vs[1][i], vs[2][i]);
            }
            return {min, max};
        };

        // auto computeQuadBoundingBox = [](const std::array<Vec2d,4>& vs) -> BoundingBox2d
        // {
        //     Vec2d min, max;
        //     for (unsigned char i = 0u; i < 2u; ++i)
        //     {
        //         min[i] = min4(vs[0][i], vs[1][i], vs[2][i], vs[3][i]);
        //         max[i] = max4(vs[0][i], vs[1][i], vs[2][i], vs[3][i]);
        //     }
        //     return {min, max};
        // };

        // auto isStrictlyInBoundingBox = [](const Vec2d& P, const BoundingBox2d& bb) -> bool
        // {
        //     return (P[0]>bb.first[0]) && (P[1]>bb.first[1]) && (P[0]<bb.second[0]) && (P[1]<bb.second[1]);
        // };

        auto projectPoint = [](const Vec3d& P, const char& p) -> Vec2d {
            return Vec2d(P[(p+1)%3], P[(p+2)%3]);
        };
        auto projectTriangleCCW = [&projectPoint, &computeTriangleBoundingBox](const std::array<Vec3d,3>& T, const char& p) -> std::vector<Vec2d>
        {
            //std::vector<Vec2d> res; res.reserve(3);
            const Vec2d A = projectPoint(T[0], p);
            const Vec2d B = projectPoint(T[1], p);
            const Vec2d C = projectPoint(T[2], p);
            auto ori = ori2d(A, B, C);
            if (ori == ORI_CCW) return {A, B, C};
            if (ori == ORI_CW) return {A, C, B};
            const auto bb = computeTriangleBoundingBox({A,B,C});
            const auto c = argAbsMax(bb.second-bb.first);
            if (isInBetween(A[c], B[c], C[c], true)) return {B, C};
            if (isInBetween(B[c], A[c], C[c], true)) return {A, C};
            assert(isInBetween(C[c], A[c], B[c], true));
            return {A, B};
        };

        auto isInTriangle = [](const Vec2d& P, const Vec2d& T1, const Vec2d& T2, const Vec2d& T3) -> bool
        {
            assert(ori2d(T1, T2, T3) == ORI_CCW);
            return (ori2d(T1, T2, P) != ORI_CW && ori2d(T2, T3, P) != ORI_CW && ori2d(T3, T1, P) != ORI_CW);
        };

        // auto isStrictlyInTriangle = [](const Vec2d& P, const Vec2d& T1, const Vec2d& T2, const Vec2d& T3) -> bool
        // {
        //     assert(ori2d(T1, T2, T3) == ORI_CCW);
        //     return (ori2d(T1, T2, P) == ORI_CCW && ori2d(T2, T3, P) == ORI_CCW && ori2d(T3, T1, P) == ORI_CCW);
        // };

        // auto edgesStrictlyIntersect = [](const Vec2d& A1, const Vec2d& A2, const Vec2d& B1, const Vec2d& B2) -> bool
        // {
        //     return (ori2d(A1, A2, B1) * ori2d(A1, A2, B2) < 0) && (ori2d(B1, B2, A1) * ori2d(B1, B2, A2) < 0);
        // };

        // 1. Check if triangle projected to 2d planes intersects projected quad
        for (const char p : {(c3d+1)%3, (c3d+2)%3})
        {
            const std::vector<Vec2d> T2d = projectTriangleCCW(T, p);
            const std::array<Vec2d, 2> Q2d = {projectPoint(Q[0], p), projectPoint(Q[2], p)};
            const char c2d = (Q2d[0][0]==Q2d[1][0])? 0 : 1;

            if (T2d.size() == 3) // Projected Triangle is Triangle
            {
                // a) Triangle must intersect the infinite quad line
                if (min3(T2d[0][c2d], T2d[1][c2d], T2d[2][c2d]) > Q2d[0][c2d] || max3(T2d[0][c2d], T2d[1][c2d], T2d[2][c2d]) < Q2d[0][c2d]) return false;

                // b1) Quad Line Segment must pierce one of the infinite triangle edge lines
                for (const char i : {0,1,2})
                    if (ori2d(T2d[(i+0)%3], T2d[(i+1)%3], Q2d[0]) * ori2d(T2d[(i+0)%3], T2d[(i+1)%3], Q2d[1]) < 0)
                        goto L1;

                // b2) Or any point (<-> all points) strictly on the quad line segment must be inside the triangle or on its boundary
                if (!isInTriangle((Q2d[0]+Q2d[1])/2, T2d[0], T2d[1], T2d[2])) return false;

            L1: continue;
            }
            else // Projected Triangle is Line Segment
            {
                assert(T2d.size() == 2);

                // a) Line Segment must intersect the infinite quad line
                if (std::min(T2d[0][c2d], T2d[1][c2d]) > Q2d[0][c2d] || std::max(T2d[0][c2d], T2d[1][c2d]) < Q2d[0][c2d]) return false;

                // b1) Quad Line Segment must pierce the infinite triangle edge line
                if (ori2d(T2d[0], T2d[1], Q2d[0]) * ori2d(T2d[0], T2d[1], Q2d[1]) < 0) continue;

                // b2) Or one of the 2 projected triangle points is strictly on the quad line segment
                const auto c = (c2d+1)%2;
                for (const char i : {0,1})
                    if (T2d[i][c2d] == Q2d[0][c2d])
                        if (isInBetween(T2d[i][c], Q2d[0][c], Q2d[1][c], true))
                            goto L2;
                return false;

            L2: continue;
            }
        }

        // Check if any quad edge pierces through the triangle plane
        for (char i = 0; i < 4; ++i)
            if (ori3d(T[0], T[1], T[2], Q[(i+0)%4]) * ori3d(T[0], T[1], T[2], Q[(i+1)%4]) < 0)
                return true;

        // Or if any quad diagonal pierces through the triangle plane
        for (char i = 0; i < 2; ++i)
            if (ori3d(T[0], T[1], T[2], Q[(i+0)%4]) * ori3d(T[0], T[1], T[2], Q[(i+2)%4]) < 0)
                return true;

        return false;
    }
}
}
