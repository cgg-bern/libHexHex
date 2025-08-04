#include <HexHex/VertexExtraction/HexVertexExtractor.hh>
#include <HexHex/VertexExtraction/epsilon_geometry.hh>
#include <HexHex/LocalTopology/HexVertexGenerator.hh>
#include <HexHex/Utils/Stopwatches.hh>
#include <HexHex/HexExtractor.hh>
#include <HexHex/Commons/Permutation.hh>

//#define TEST_VERTEX_EXTRACTION_CORRECTNESS

namespace HexHex
{

    struct IntegerBoundingBox3d
    {
        Vec3i min;
        Vec3i max;
        inline bool empty()
        {
            return (max[0]<min[0])||(max[1]<min[1])||(max[2]<min[2]);
        }
        inline Vec3i size()
        {
            return Vec3i(max[0]-min[0],max[1]-min[1],max[2]-min[2]);
        }
        inline bool isInside(const Vec3d& p)
        {
            return (p[0]>=min[0]&&p[0]<=max[0])&&(p[1]>=min[1]&&p[1]<=max[1])&&(p[2]>=min[2]&&p[2]<=max[2]);
        }
        inline double clamp_x(const double x) {return (x<min[0])? min[0] : (x>max[0])? max[0] : x;}
        inline double clamp_y(const double y) {return (y<min[1])? min[1] : (y>max[1])? max[1] : y;}
    };

    IntegerBoundingBox3d computeTetIntegerBoundingBox(const std::array<Vec3d,4>& vs)
    {
        Vec3i min, max;
        for (unsigned char i = 0u; i < 3u; ++i)
        {
            min[i] = std::ceil(min4(vs[0][i], vs[1][i], vs[2][i], vs[3][i]));
            max[i] = std::floor(max4(vs[0][i], vs[1][i], vs[2][i], vs[3][i]));
        }
        return IntegerBoundingBox3d{.min = min, .max = max};
    }

    VertexExtractor::VertexExtractor(HexExtractor& hexEx) :
        hexEx(hexEx),
        v_igps(hexEx.inputMesh.create_private_vertex_property<std::vector<Parameter>>()),
        e_igps(hexEx.inputMesh.create_private_edge_property<std::vector<Parameter>>()),
        f_igps(hexEx.inputMesh.create_private_face_property<std::vector<Parameter>>()),
        c_igps(hexEx.inputMesh.create_private_cell_property<std::vector<Parameter>>()),
        piInv(hexEx.inputMesh.create_private_cell_property<Permutation>()),
        tet_debug_props(hexEx.inputMesh.create_private_cell_property<std::vector<Parameter>>())
    {}

    void VertexExtractor::addIntegerGridPoint(const CellHandle& ch, const MeshElement& elem, const Parameter& param)
    {
#ifdef TEST_VERTEX_EXTRACTION_CORRECTNESS
        {
            MeshElement assElem;
            assert(hexEx.computeVertexGeneratorElement(ch, param, assElem));
            assert(elem.type() == assElem.type());
            tet_debug_props[ch].push_back(param);
        }
#endif

        switch (elem.type()) {
        case ELEMENT_VERTEX_TYPE:
            if (ch.idx() != hexEx.getCache().v_owners[elem.vh()]) return;
            assert(v_igps[elem.vh()].empty());
            #pragma omp atomic
            num_hex_vertices++;
            #pragma omp atomic
            num_generators++;
            v_igps[elem.vh()].push_back(param);
            break;
        case ELEMENT_EDGE_TYPE:
            if (ch.idx() != hexEx.getCache().e_owners[elem.eh()]) return;
            e_igps[elem.eh()].push_back(param);
            #pragma omp atomic
            num_hex_vertices++;
            if (e_igps[elem.eh()].size()==1) {
                #pragma omp atomic
                num_generators++;
            }
            break;
        case ELEMENT_FACE_TYPE:
            if (ch.idx() != hexEx.getCache().f_owners[elem.fh()]) return;
            f_igps[elem.fh()].push_back(param);
            #pragma omp atomic
            num_hex_vertices++;
            if (f_igps[elem.fh()].size()==1) {
                #pragma omp atomic
                num_generators++;
            }
            break;
        case ELEMENT_CELL_TYPE:
            c_igps[elem.ch()].push_back(param);
            #pragma omp atomic
            num_hex_vertices++;
            if (c_igps[elem.ch()].size()==1) {
                #pragma omp atomic
                num_generators++;
            }
            break;
        default: break;
        }
    }

    bool VertexExtractor::extract()
    {
        ScopedStopWatch _{sw::extractHexVertices};
        auto& inputMesh = hexEx.inputMesh;
        num_hex_vertices = 0;
        num_generators = 0;

        // Prepare variables
        const double EPSILON = hexEx.config.rasterization_epsilon;

        tet_debug_props.fill(tet_debug_props.def());

        auto bb_per_cell_permuted = inputMesh.create_private_cell_property<IntegerBoundingBox3d>();

        // =========================
        // Rasterization Functions
        //==========================

        const Vec3d ZERO_3D(0,0,0);
        const TetEdge INVALID_TET_EDGE(ZERO_3D, ZERO_3D, 1);

        /*
         * Add integer-grid points on the (permuted) line segment (x,y,ceil(xld)), (x,y,floor(xrd))
         */
        auto rasterizeSegment1dExp = [&](const CellHandle& ch, const int32 z, const int32 y, const double xld, const double xrd) -> void
        {

            if (xld > xrd) return;
            assert(xld <= xrd);

            // Safety measure, and clamp to bounding box
            int32 xl = ceil(bb_per_cell_permuted[ch].clamp_x(xld));
            int32 xr = floor(bb_per_cell_permuted[ch].clamp_x(xrd));
            if (xl > xr) return;

            // check from left and from right
            MeshElement eleml, elemr;
            Parameter paraml, paramr;
            for (; xl <= xr; ++xl) {paraml = piInv[ch].permuted(Parameter(xl, y, z)); if (hexEx.computeVertexGeneratorElement(ch, paraml, eleml)) break;}
            for (; xr > xl; --xr) {paramr = piInv[ch].permuted(Parameter(xr, y, z)); if (hexEx.computeVertexGeneratorElement(ch, paramr, elemr)) break;}
            if (!eleml.is_valid()) return;
            if (!elemr.is_valid()) elemr = eleml;
            int n = xr-xl+1;

            // add left and right points
            addIntegerGridPoint(ch, eleml, paraml);
            if (n > 1) addIntegerGridPoint(ch, elemr, paramr);

            if (n > 2)
            {
                // see, if we can trivially get the in between generator
                MeshElement elemm;
                if (eleml.is_cell() || eleml == elemr) elemm = eleml;
                else if (elemr.is_cell()) elemm = elemr;

                // add in between points
                for (int x = xl+1; x <= xr-1; ++x)
                {
                    auto param = piInv[ch].permuted(Parameter(x, y, z)); assert(isInteger(param));
                    if (!elemm.is_valid()) hexEx.computeVertexGeneratorElement(ch, param, elemm);
                    assert(elemm.is_valid());
                    addIntegerGridPoint(ch, elemm, param);
                }
            }
        };



        auto rasterizeTriangleExp = [&](const CellHandle& ch, const int z, BB2d A, BB2d B, BB2d C) -> void
        {

            // Get the y range
            const double y_min_d = std::ceil(min3(C.b(),B.b(),A.b()));
            const double y_max_d = std::floor(max3(A.t(),B.t(),C.t()));
            if (y_max_d < y_min_d) return;
            const int y_min_int = bb_per_cell_permuted[ch].clamp_y(y_min_d);
            const int y_max_int = bb_per_cell_permuted[ch].clamp_y(y_max_d);
            assert(y_max_int >= y_min_int);

            // Sort by  s.t. A.b >= B.b >= C.b
            if (C.b() > B.b()) std::swap(C, B);
            if (B.b() > A.b()) std::swap(B, A);
            if (C.b() > B.b()) std::swap(C, B);
            assert(A.b() >= B.b() && B.b() >= C.b());

            // Cache the Poly Edges
            PolyEdge CA(C, A, EPSILON);
            PolyEdge CB(C, B, EPSILON);
            PolyEdge BA(B, A, EPSILON);

            // Sweep along the y-axis
            int y = y_min_int;
            std::pair<double, double> x_range;

            // 1. In C epsilon range
            for (; y <= C.t(); ++y)
            {
                if (y > y_max_int) break;

                // Consider all edges
                x_range = poly_edges_x_range(y, {CA,CB,BA});
                rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
            }

            // 2. Between C and B
            for (; y < B.b(); ++y)
            {
                if (y > y_max_int) break;

                // Consider edges [CA][CB]
                x_range = poly_edges_x_range(y, {CA,CB});
                rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
            }

            // 3. In B epsilon range
            for (; y <= B.t(); ++y)
            {
                if (y > y_max_int) break;

                // Consider all edges
                x_range = poly_edges_x_range(y, {CA,CB,BA});
                rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
            }

            // 4. Between B and A
            for (; y < A.b(); ++y)
            {
                if (y > y_max_int) break;

                // Consider Edges [CA][BA]
                x_range = poly_edges_x_range(y, {CA,BA});
                rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
            }

            // 5. In A epsilon range
            for (; y <= y_max_int; ++y)
            {
                // Consider all edges
                x_range = poly_edges_x_range(y, {CA,CB,BA});
                rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
            }
        };

        auto rasterizeQuadrilateralExp = [&](const CellHandle& ch, const int z, BB2d A, BB2d B, BB2d C, BB2d D) -> void
        {
            // Get the y range
            const double min_b = min4(D.b(),C.b(),B.b(),A.b());
            const double max_t = max4(A.t(),B.t(),C.t(),D.t());
            const double y_min_d = std::ceil(min_b);
            const double y_max_d = std::floor(max_t);
            if (y_max_d < y_min_d) return;
            const int y_min_int = bb_per_cell_permuted[ch].clamp_y(y_min_d);
            const int y_max_int = bb_per_cell_permuted[ch].clamp_y(y_max_d);
            assert(y_max_int >= y_min_int);

            // Rotate the points such that D bottom is the bottomest
            while (D.b() > min_b) {
                std::swap(A, B);
                std::swap(C, D);
                std::swap(B, D);
            }
            assert(D.b() == min_b);
            std::pair<double, double> x_range;
            std::vector<PolyEdge> edges; edges.reserve(4);

            // We have 3 cases, A is the top or B or C
            if (B.t() == max_t)
            {
                // B is top, we have the edges [DA][DC][AB][CB]
                PolyEdge DA(D, A, EPSILON);
                PolyEdge DC(D, C, EPSILON);
                PolyEdge AB(A, B, EPSILON);
                PolyEdge CB(C, B, EPSILON);

                // Sweep along y-axis
                int y = y_min_int;

                for (; y <= y_max_int; ++y)
                {
                    x_range = poly_edges_x_range(y, {DA, DC, AB, CB});
                    rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
                }
            }
            else if (A.t() == max_t)
            {
                // A is top, we have the edges [DC][CB][BA][DA]
                PolyEdge DC(D, C, EPSILON);
                PolyEdge CB(C, B, EPSILON);
                PolyEdge BA(B, A, EPSILON);
                PolyEdge DA(D, A, EPSILON);

                // Sweep along y-axis
                int y = y_min_int;

                // 1. D epsilon range
                for (; y <= D.t(); ++y)
                {
                    if (y > y_max_int) break;

                    // Consider all edges
                    x_range = poly_edges_x_range(y, {DC, CB, BA, DA});
                    rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
                }

                // 2. Between D and C
                for (; y < C.b(); ++y)
                {
                    if (y > y_max_int) break;

                    // Consider edges [DA][DC]
                    x_range = poly_edges_x_range(y, {DA,DC});
                    rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
                }

                // 3. C epsilon range
                for (; y <= C.t(); ++y)
                {
                    if (y > y_max_int) break;

                    // Consider all edges
                    x_range = poly_edges_x_range(y, {DC, CB, BA, DA});
                    rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
                }

                // 4. Between C and B
                for (; y < B.b(); ++y)
                {
                    if (y > y_max_int) break;

                    // Consider edges [DA][CB]
                    x_range = poly_edges_x_range(y, {DA,CB});
                    rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
                }

                // 5. B epsilon range
                for (; y < B.t(); ++y)
                {
                    if (y > y_max_int) break;

                    // Consider all edges
                    x_range = poly_edges_x_range(y, {DC, CB, BA, DA});
                    rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
                }

                // 6. Between B and A
                for (; y < A.b(); ++y)
                {
                    if (y > y_max_int) break;

                    // Consider edges [DA][BA]
                    x_range = poly_edges_x_range(y, {DA,BA});
                    rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
                }

                // 7. A epsilon range
                for (; y <= y_max_int; ++y)
                {
                    // Consider all edges
                    x_range = poly_edges_x_range(y, {DC, CB, BA, DA});
                    rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
                }
            }
            else if (C.t() == max_t)
            {

                // C is top, we have the edges [DC][DA][AB][BC]
                PolyEdge DC(D, C, EPSILON);
                PolyEdge DA(D, A, EPSILON);
                PolyEdge AB(A, B, EPSILON);
                PolyEdge BC(B, C, EPSILON);

                // Sweep along y-axis
                int y = y_min_int;


                // 1. D epsilon range
                for (; y <= D.t(); ++y)
                {
                    if (y > y_max_int) break;

                    // Consider all edges
                    x_range = poly_edges_x_range(y, {DC, DA, AB, BC});
                    rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
                }

                // 2. Between D and A
                for (; y < A.b(); ++y)
                {
                    if (y > y_max_int) break;

                    // Consider edges [DA][DC]
                    x_range = poly_edges_x_range(y, {DA,DC});
                    rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
                }

                // 3. A epsilon range
                for (; y <= A.t(); ++y)
                {
                    if (y > y_max_int) break;

                    // Consider all edges
                    x_range = poly_edges_x_range(y, {DC, DA, AB, BC});
                    rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
                }

                // 4. Between A and B
                for (; y < B.b(); ++y)
                {
                    if (y > y_max_int) break;

                    // Consider edges [DC][AB]
                    x_range = poly_edges_x_range(y, {DC,AB});
                    rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
                }

                // 5. B epsilon range
                for (; y < B.t(); ++y)
                {
                    if (y > y_max_int) break;

                    // Consider all edges
                    x_range = poly_edges_x_range(y, {DC, DA, AB, BC});
                    rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
                }

                // 6. Between B and C
                for (; y < C.b(); ++y)
                {
                    if (y > y_max_int) break;

                    // Consider edges [DC][BC]
                    x_range = poly_edges_x_range(y, {DC,BC});
                    rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
                }

                // 7. C epsilon range
                for (; y <= y_max_int; ++y)
                {
                    // Consider all edges
                    x_range = poly_edges_x_range(y, {DC, DA, AB, BC});
                    rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
                }
            }
            else
            {
                // D is both bottom and top, just rasterize rectangle
                assert(D.t() == max_t);
                x_range.first = min4(D.l(), C.l(), B.l(), A.l());
                x_range.second = max4(D.r(), C.r(), B.r(), A.r());

                // Sweep along y-axis
                int y = y_min_int;
                for (; y <= y_max_int; ++y)
                {
                    rasterizeSegment1dExp(ch, z, y, x_range.first, x_range.second);
                }
            }
        };

        /*
         * Checks the z-aligned (permuted) parameter of a tet-vertex
         */
        auto checkTetVertexForIGP = [&](const CellHandle& ch, const int z, const Vec2d& u) -> void
        {
            if (u[0] != (int)u[0] || u[1] != (int)u[1]) return; // If not integer, return
            Parameter U = piInv[ch].permuted(Parameter(u[0], u[1], z));
            for (const auto& vh : hexEx.getCache().cellVertices[ch]) // Find which vertex it is
            {
                if (hexEx.getParametrization().get(ch, vh) == U)
                {
                    addIntegerGridPoint(ch, vh, U);
                    return;
                }
            }
        };

        /*
         * Tests a (permuted) 2d line segment. This will only be called if a tet is touching the z-plane with one or two of its vertices.
         * In particular, the endpoints are exact vertices!
         */
        auto checkTetEdgeForIGPs = [&](const CellHandle& ch, const int z, Vec2d A, Vec2d B) -> void
        {

            // Check if the segment is aligned in the y dimension -> Call the standard 1d line algorithm
            if (A[1] == B[1]) {
                const int y = (int)A[1];
                if (A[1] == y) {
                    rasterizeSegment1dExp(ch, z, y, std::min(A[0],B[0]), std::max(A[0],B[0]));
                }
                return;
            }

            // Check if the segment is aligned in the x dimension
            if (false && A[0] == B[0])
            {
                const int x = (int)A[0];
                if (A[0] == x)
                {
                    // Segment goes from (x, yb, z) to (x, yt, z)
                    int yb = std::ceil(std::min(A[1], B[1]));
                    int yt = std::floor(std::max(A[1], B[1]));
                    if (yb > yt) return;

                    // Get top and bottom generators
                    MeshElement elemb, elemt;
                    Parameter paramb, paramt;
                    for (; yb <= yt; ++yb) {paramb = piInv[ch].permuted(Parameter(x, yb, z)); if (hexEx.computeVertexGeneratorElement(ch, paramb, elemb)) break;}
                    for (; yt > yb; --yt) {paramt = piInv[ch].permuted(Parameter(x, yt, z)); if (hexEx.computeVertexGeneratorElement(ch, paramt, elemt)) break;}
                    assert(elemb.is_vertex() || elemb.is_edge());
                    if (!elemt.is_valid()) elemt = elemb;
                    const int n = yt-yb+1;
                    assert(n >= 1);

                    // add bottom and top points
                    addIntegerGridPoint(ch, elemb, paramb);
                    if (n > 1) addIntegerGridPoint(ch, elemt, paramt);

                    if (n > 2)
                    {
                        // The inbetween generator must be the edge!
                        MeshElement elemm;
                        if (elemb.is_edge()) {elemm = elemb;} else if (elemt.is_edge()) {elemm = elemt;}

                        // Add between points
                        for (int y = yb+1; y <= yt-1; ++y)
                        {
                            auto param = piInv[ch].permuted(Parameter(x, y, z)); assert(isInteger(param));
                            if (!elemm.is_valid()) hexEx.computeVertexGeneratorElement(ch, param, elemm);
                            assert(elemm.is_edge());
                            addIntegerGridPoint(ch, elemm, param);
                        }
                    }

                }
                return;
            }

            // Check 2d line segment by sweeping along y axis and rounding interpolated values
            // The generator element will either be a tet-vertex or a tet-edge
            if (A[1] > B[1]) std::swap(A, B);
            assert(A[1] < B[1]);
            const int y_min = ceil(A[1]);
            const int y_max = floor(B[1]);
            for (int y = y_min; y <= y_max; ++y)
            {
                const double t = (y - A[1]) / (B[1] - A[1]);
                const double x = std::round(A[0] + t * (B[0] - A[0]));

                auto param = piInv[ch].permuted(Parameter(x,y,z));
                MeshElement elem;
                if (hexEx.computeVertexGeneratorElement(ch, param, elem))
                {
                    assert(elem.is_vertex() || elem.is_edge());
                    addIntegerGridPoint(ch, elem, param);
                }
            }
        };

        auto rasterizeTetrahedronExp = [&](const CellHandle& ch, Parameter A, Parameter B, Parameter C, Parameter D)
        {
            // Sort by z st Dz <= Cz <= Bz <= Az
            if(D[2] > C[2]) std::swap(D, C);
            if(B[2] > A[2]) std::swap(B, A);
            if(D[2] > B[2]) std::swap(D, B);
            if(C[2] > A[2]) std::swap(C, A);
            if(C[2] > B[2]) std::swap(C, B);
            assert(A[2] >= B[2] && B[2] >= C[2] && C[2] >= D[2]);

            // Cache Tet Edges if needed (an edge is needed if it strictly intersects any integer z-plane)
            TetEdge DC = (std::floor(D[2])+1 < std::ceil(C[2]))? TetEdge(D,C,EPSILON) : INVALID_TET_EDGE;
            TetEdge DB = (std::floor(D[2])+1 < std::ceil(B[2]))? TetEdge(D,B,EPSILON) : INVALID_TET_EDGE;
            TetEdge DA = (std::floor(D[2])+1 < std::ceil(A[2]))? TetEdge(D,A,EPSILON) : INVALID_TET_EDGE;
            TetEdge CB = (std::floor(C[2])+1 < std::ceil(B[2]))? TetEdge(C,B,EPSILON) : INVALID_TET_EDGE;
            TetEdge CA = (std::floor(C[2])+1 < std::ceil(A[2]))? TetEdge(C,A,EPSILON) : INVALID_TET_EDGE;
            TetEdge BA = (std::floor(B[2])+1 < std::ceil(A[2]))? TetEdge(B,A,EPSILON) : INVALID_TET_EDGE;

            // Sweep along z-axis
            int z = ceil(D[2]);

            // 1. z = D_z
            if (D[2] == z)
            {

                if (B[2] == z) // Face DCB
                {
                    rasterizeTriangleExp(ch, z,
                        tet_vertex_intersection_z(z, D),
                        tet_vertex_intersection_z(z, C),
                        tet_vertex_intersection_z(z, B)
                    );
                }
                else if (C[2] == z) // Edge DC
                {
                    checkTetEdgeForIGPs(ch, z, {D[0], D[1]}, {C[0], C[1]});
                }
                else // Vertex D
                {
                    checkTetVertexForIGP(ch, z, {D[0], D[1]});
                }
                ++z; // increment z
            }

            // 2. D_z < z < C_z
            for (; z < C[2]; ++z)
            {
                assert(z > D[2] && z < C[2]);

                // Triangle [DC][DB][DA]
                rasterizeTriangleExp(ch, z,
                    tet_edge_intersection_z(z, DC),
                    tet_edge_intersection_z(z, DB),
                    tet_edge_intersection_z(z, DA)
                );
            }

            // 3. z = C_z
            if (C[2] == z)
            {
                assert(C[2] > D[2]);

                if (A[2] == z) // Triangle CBA
                {
                    rasterizeTriangleExp(ch, z,
                        tet_vertex_intersection_z(z, C),
                        tet_vertex_intersection_z(z, B),
                        tet_vertex_intersection_z(z, A)
                    );
                }
                else if (B[2] == z) // Triangle CB[DA]
                {
                    rasterizeTriangleExp(ch, z,
                        tet_vertex_intersection_z(z, C),
                        tet_vertex_intersection_z(z, B),
                        tet_edge_intersection_z(z, DA)
                    );
                }
                else // Triangle C[DB][DA]
                {
                    rasterizeTriangleExp(ch, z,
                        tet_vertex_intersection_z(z, C),
                       tet_edge_intersection_z(z, DB),
                       tet_edge_intersection_z(z, DA)
                    );
                }
                ++z; // increment z
            }

            // 4. C_z < z < B_z
            for (; z < B[2]; ++z)
            {
                assert(z > C[2] && z < B[2]);

                // Quad [DA][DB][CB][CA] (consider order to prevent self intersecting quad)
                rasterizeQuadrilateralExp(ch, z,
                    tet_edge_intersection_z(z, DA),
                    tet_edge_intersection_z(z, DB),
                    tet_edge_intersection_z(z, CB),
                    tet_edge_intersection_z(z, CA)
                );
            }

            // 5. z = B_z
            if (B[2] == z)
            {
                assert(B[2] > C[2]);

                if (A[2] == z) // Edge BA
                {
                    checkTetEdgeForIGPs(ch, z, {B[0], B[1]}, {A[0], A[1]});
                }
                else // Triangle B[DA][CA]
                {
                    rasterizeTriangleExp(ch, z,
                        tet_vertex_intersection_z(z, B),
                       tet_edge_intersection_z(z, DA),
                       tet_edge_intersection_z(z, CA)
                   );
                }
                ++z; // increment z
            }

            // 6. B_z < z < A_z
            for (; z < A[2]; ++z)
            {
                assert(z > B[2] && z < A[2]);

                // Triangle [DA][CA][BA]
                rasterizeTriangleExp(ch, z,
                   tet_edge_intersection_z(z, DA),
                   tet_edge_intersection_z(z, CA),
                   tet_edge_intersection_z(z, BA)
               );
            }

            // 7. z = A_z
            if (A[2] == z)
            {
                assert(A[2] > B[2]);

                // Vertex A
                checkTetVertexForIGP(ch, z, {A[0], A[1]});
            }
        };

        // =================================
        // Main Hex Vertex Extraction Loop
        //==================================
#ifdef TEST_VERTEX_EXTRACTION_CORRECTNESS
        bool is_correct = true;
#endif
#pragma omp parallel for schedule(guided)
        for (int ch_idx = 0; ch_idx < (int)(inputMesh.n_cells()); ++ch_idx)
        {
            auto ch = CellHandle(ch_idx);

            // Permute coordinate system s.t. dz <= dy <= dx
            auto us = hexEx.getParametrization().get(ch);
            auto bb = computeTetIntegerBoundingBox(us);
            if (bb.empty()) continue;

            auto bbs = bb.size();
            auto pi = Permutation::descending(bbs);
            piInv[ch] = pi.inverse();
            for (auto i = 0u; i < 4; ++i) us[i] = pi.permuted(us[i]);
            bb.min = pi.permuted(bb.min);
            bb.max = pi.permuted(bb.max);
            bb_per_cell_permuted[ch] = bb;

#if 1

            rasterizeTetrahedronExp(ch, us[0], us[1], us[2], us[3]);

#else
        // Vertex Extraction via naive Bounding Box check (only for experimental purposes)
        MeshElement elem;
        for (int x = bb.min[0]; x <= bb.max[0]; ++x)
            for (int y = bb.min[1]; y <= bb.max[1]; ++y)
                for (int z = bb.min[2]; z <= bb.max[2]; ++z)
                {
                    Parameter param = piInv[ch].permuted(Parameter(x,y,z));
                    if (hexEx.computeVertexGeneratorElement(ch, param, elem))
                        addIntegerGridPoint(ch, elem, param);
                }
#endif

#ifdef TEST_VERTEX_EXTRACTION_CORRECTNESS
            // Assert Correctness

            std::vector<std::pair<MeshElement, Parameter>> igps;
            auto params = hexEx.getParametrization().get(ch);
            bb = computeTetIntegerBoundingBox(params);
            MeshElement elem;
            for (int x = bb.min[0]; x <= bb.max[0]; ++x)
                for (int y = bb.min[1]; y <= bb.max[1]; ++y)
                    for (int z = bb.min[2]; z <= bb.max[2]; ++z)
                        if (hexEx.computeVertexGeneratorElement(ch, Parameter(x,y,z), elem))
                            igps.push_back({elem, Parameter(x,y,z)});

            if (tet_debug_props[ch].size() != igps.size())
            {
                std::cerr << "We should have " << igps.size() << " vertices in cell " << ch << " but we found only " << tet_debug_props[ch].size() << std::endl;

                std::cerr << "Cell Parameters:" << std::endl;
                //std::cerr << std::hexfloat;
                std::cerr << us[0] << std::endl;
                std::cerr << us[1] << std::endl;
                std::cerr << us[2] << std::endl;
                std::cerr << us[3] << std::endl;

                std::cerr << "Missing Vertices:" << std::endl;

                for (const auto& igp : igps)
                {
                    bool found = false;
                    for (const auto& found_igp : tet_debug_props[ch])
                    {
                        if (igp.second == toVec3d(found_igp))
                        {
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                    {
                        std::cerr << pi.permuted(igp.second) << " (" << igp.first << ")" << std::endl;
                    }
                }

                hexEx.addMessage(HHM_MissingHexVertex());
                is_correct = false;
            }
#endif
        }

#ifdef TEST_VERTEX_EXTRACTION_CORRECTNESS
        hexEx.printLog("HexHex: Asserted Vertex Extraction Correctness by Checking Bounding Boxes: ", is_correct);
        if (!is_correct) {return false;}
#endif

        // Reserve space for generators and hex vertices
        hexEx.generators.reserve(num_generators);
        hexEx.getOutputMesh().reserve_vertices(num_hex_vertices);
        hexEx.printLog("HexHex: Detected ", num_hex_vertices, " Vertices on ",
                       num_generators, " Generators.");

        // Add Vertex Hex-Vertices
        for (size_t v_i = 0; v_i < inputMesh.n_vertices(); ++v_i)
        {
            VertexHandle vh(v_i);
            hexEx.addHexVertices(v_igps[vh], vh, CellHandle(hexEx.getCache().v_owners[vh]));
        }

        // Add Edge Hex-Vertices
        for (size_t e_i = 0; e_i < inputMesh.n_edges(); ++e_i)
        {
            EdgeHandle eh(e_i);
            hexEx.addHexVertices(e_igps[eh], eh, CellHandle(hexEx.getCache().e_owners[eh]));
        }

        // Add Face Hex-Vertices
        for (size_t f_i = 0; f_i < inputMesh.n_faces(); ++f_i)
        {
            FaceHandle fh(f_i);
            CellHandle ch(hexEx.getCache().f_owners[fh]);
            HalfFaceHandle hfh = hexEx.getCache().get_incident_halfface_in_cell(ch, fh);
            hexEx.addHexVertices(f_igps[fh], hfh, ch);
        }

        // Add Cell Hex-Vertices
        for (size_t c_i = 0; c_i < inputMesh.n_cells(); ++c_i)
        {
            CellHandle ch(c_i);
            hexEx.addHexVertices(c_igps[ch], ch, ch);
        }

        assert(num_generators == hexEx.generators.size());

        hexEx.printLog("HexHex: Extracted ", hexEx.getOutputMesh().n_vertices(), " Vertices on ",
                       hexEx.generators.size(), " Generators. ",
        " (User Epsilon = ", EPSILON, ")");

        return true;
    }
}

