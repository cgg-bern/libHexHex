#pragma once

#include <HexHex/Utils/Typedefs.hh>

namespace HexHex
{

/*
 * Helper Struct to represent a 2d point with some x and y uncertainty.
 */
struct BB2d
{
    Vec2d p; // center
    Vec2d hs; // half size
    inline bool is_valid() const {return !(hs[0] < 0.0);}
    inline Vec2d tl() const {return p + Vec2d{-hs[0], +hs[1]};} // top left
    inline Vec2d tr() const {return p + Vec2d{+hs[0], +hs[1]};} // top right
    inline Vec2d bl() const {return p + Vec2d{-hs[0], -hs[1]};} // bottom left
    inline Vec2d br() const {return p + Vec2d{+hs[0], -hs[1]};} // bottom right
    inline double t() const {return p[1] + hs[1];} // top
    inline double b() const {return p[1] - hs[1];} // bottom
    inline double l() const {return p[0] - hs[0];} // left
    inline double r() const {return p[0] + hs[0];} // right
};

inline double compute_epsilon(const double GLOBAL_EPSILON, const double dx, const double one_over_dy)
{
    return GLOBAL_EPSILON + GLOBAL_EPSILON * abs(dx * one_over_dy);
}

/*
 * Helper Struct to represent a value with some uncertainty.
 */
struct BB1d
{
    double p; // center
    double hs; // half size
    inline bool is_valid() const {return !(hs < 0.0);}
    inline double l() const {return p - hs;} // left
    inline double r() const {return p + hs;} // right
};

/*
 * Helper Struct to represent an edge (left and right) between two uncertain points.
 * We must have A.y <= B.y (within epsilon range)
 */
struct PolyEdge
{
    /*
     * Helper Struct to represent a (nonparallel) left edge between two uncertain points.
     */
    struct Left
    {
        Left() : dx(0), dy(0), dyInv(0), epsilon(0) {}

        Left(const BB2d& A_, const BB2d& B_, const double EPSILON) :
            A((B_.l() < A_.l())? A_.bl() : A_.tl()),
            B((B_.l() < A_.l())? B_.bl() : B_.tl()),
            dx(B[0] - A[0]),
            dy(B[1] - A[1]),
            dyInv(1.0 / dy),
            epsilon(compute_epsilon(EPSILON, dx, dyInv))
            //epsilon(EPSILON * (std::abs(dx*dyInv) + 1))
        {
            assert(A_.b() <= B_.b());
            assert(epsilon >= 0);
        }

        const Vec2d A;
        const Vec2d B;
        const double dx;
        const double dy;
        const double dyInv;
        const double epsilon;
    };

    /*
     * Helper Struct to represent a (nonparallel) right edge between two uncertain points.
     */
    struct Right
    {
        Right() : dx(0), dy(0), dyInv(0), epsilon(0) {}

        Right(const BB2d& A_, const BB2d& B_, const double EPSILON) :
            A((B_.r() < A_.r())? A_.tr() : A_.br()),
            B((B_.r() < A_.r())? B_.tr() : B_.br()),
            dx(B[0] - A[0]),
            dy(B[1] - A[1]),
            dyInv(1.0 / dy),
            epsilon(EPSILON * (std::abs(dx*dyInv) + 1))
        {
            assert(A_.b() <= B_.b());
            assert(epsilon >= 0);
        }

        const Vec2d A;
        const Vec2d B;
        const double dx;
        const double dy;
        const double dyInv;
        const double epsilon;
    };

    PolyEdge(const BB2d& A_, const BB2d& B_, const double EPSILON) :
        A(A_), B(B_),
        L((A.t() < B.b())? Left{A_, B_, EPSILON} : Left{}), // ignore left & right edge if parallel
        R((A.t() < B.b())? Right{A_, B_, EPSILON} : Right{})
    {
        assert(A_.b()-B_.b()<=10*(A_.hs[1]+B_.hs[1]));
    }

    inline bool is_y_parallel() const {return !(A.t() < B.b());}

    const BB2d& A;
    const BB2d& B;

    const Left L;
    const Right R;
};

inline std::ostream& operator<<(std::ostream &os, const PolyEdge &E)
{
    os << "PolyEdge: A=(" << E.A.p << ", " << E.A.hs << "), B = " << E.B.p << ", " << E.B.hs << ")";
    return os;
}

/*
 * Computes the intersection x (range) of the given Poly Edge with a y-line
 */
std::pair<double, double> poly_edge_x_range(const int y, const PolyEdge& E)
{
    // Return empty range if Edge does not cut the y-line
    if (y < E.A.b() || y > std::max(E.A.t(), E.B.t())) return {INFINITY, -INFINITY};

    // Return max range if epsilon-parallel
    if (E.is_y_parallel()) return {std::min(E.A.l(),E.B.l()),std::max(E.A.r(),E.B.r())};

    // Left
    const double t_l = (y<=E.L.A[1])? 0 : (y>=E.L.B[1])? 1 : (y - E.L.A[1]) * E.L.dyInv;
    const double x_l = E.L.A[0] + t_l*E.L.dx;
    const double left = x_l - E.L.epsilon;

    // Right
    const double t_r = (y<=E.R.A[1])? 0 : (y>=E.R.B[1])? 1 : (y - E.R.A[1]) * E.R.dyInv;
    const double x_r = E.R.A[0] + t_r*E.R.dx;
    const double right = x_r + E.R.epsilon;

    return {left, right};
};

std::pair<double, double> poly_edges_x_range(const int32 y, const std::vector<PolyEdge>& Es)
{
    double x_min = INFINITY;
    double x_max = -INFINITY;

    for (const auto& E : Es)
    {
        const auto x_r = poly_edge_x_range(y, E);

        if (x_r.first < x_min) x_min = x_r.first;
        if (x_r.second > x_max) x_max = x_r.second;
    }

    return {x_min, x_max};
};

/*
 * Helper Struct to represent a (nonparallel) edge of a rasterized tet
 * We must have A.z < B.z!
 */
struct TetEdge
{
    TetEdge(const Vec3d& A, const Vec3d& B, const double EPSILON) :
        A(A), B(B),
        dx(B[0] - A[0]), dy(B[1] - A[1]), dz(B[2] - A[2]), dzInv(1.0 / dz),
        epsilon(
            compute_epsilon(EPSILON, dx, dzInv),
            compute_epsilon(EPSILON, dy, dzInv)
        )
        // epsilon(
        //     EPSILON * (std::abs(dx)*dzInv + 1),
        //     EPSILON * (std::abs(dy)*dzInv + 1)
        //     )
    {
        assert((B==Vec3d(0,0,0)&&A==Vec3d(0,0,0)) || B[2] > A[2]);
    }

    const Vec3d& A;
    const Vec3d& B;
    const double dx;
    const double dy;
    const double dz;
    const double dzInv;
    const Vec2d epsilon;
};



/*
 * "Computes" the intersection of a point A which lies on the z-plane as an uncertainty point in xy
 */
BB2d tet_vertex_intersection_z(const int z, const Vec3d& A)
{
    assert(A[2] == z);
    return BB2d{.p = {A[0],A[1]}, .hs = {0.0,0.0}};
};

/*
 * Computes the intersection of a line segment AB which pierces the z-plane as an uncertainty point in xy
 */
BB2d tet_edge_intersection_z(const int z, const TetEdge& E)
{
    // Assert strict piercing
    assert((E.A[2] < z && E.B[2] > z));

    const double t = (z<=E.A[2])? 0 : (z>=E.B[2])? 1 : (z - E.A[2]) * E.dzInv;
    const double x = E.A[0] + t * E.dx;
    const double y = E.A[1] + t * E.dy;

    return BB2d{.p = Vec2d{x, y}, .hs = E.epsilon};
};

}

