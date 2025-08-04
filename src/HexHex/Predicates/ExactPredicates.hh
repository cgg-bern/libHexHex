#pragma once

/** @file
 * @brief Defines the the exact predicate ::orient2d and ::orient3d, and
 * enums, types and functions it needs.
 */


#include <cassert>
#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Initializes the exact predicates.
 *
 * This function has to be called before calling ::orient2d.
 */
void exactinit();

/**
 * @brief Computes the orientation of the supplied points.
 *
 * This function only returns correct values if exactinit() was called before.
 *
 * @param pa, pb, pc Arrays containing the two coordinates of a point, each.
 * @return A value the sign of which reflects the orientation of the three points.
 */
double orient2d(const double* pa, const double* pb, const double* pc);

/**
 * @brief Computes the orientation of the supplied points.
 *
 * This function only returns correct values if exactinit() was called before.
 *
 * @param pa, pb, pc, pd Arrays containing the three coordinates of a point, each.
 * @return A value the sign of which reflects the orientation of the four points.
 */
double orient3d(const double* pa, const double* pb, const double* pc, const double*  pd);

/**
 * @brief Specifies the orientation of a tuple of three 2D points.
 *
 * There is a bijection between the orientation of points (a, b, c)
 * and the sign of the determinant of the matrix [b-a c-a] which
 * is reflected by the aliases for ::ORI_CW, ::ORI_CCW and ::ORI_COLLINEAR.
 */
typedef enum {
    ORI_CW = -1,      //!< The points are ordered in clockwise order.
    ORI_CCW = 1,      //!< The points are ordered in counterclockwise order.
    ORI_COLLINEAR = 0,//!< The points are collinear.
    ORI_NEGATIVE = -1,//!< The determinant/area of the triangle is negative, equivalent to ::ORI_CW.
    ORI_ZERO = 0,     //!< The determinant/area of the triangle is zero, equivalent to ::ORI_COLLINEAR.
    ORI_POSITIVE = 1, //!< The determinant/area of the triangle is positive, equivalent to ::ORI_CCW.
    ORI_RIGHT = -1,   //!< Connecting the points one has to make a right, equivalent to ::ORI_CW.
    ORI_LEFT = 1,     //!< Connecting the points one has to make a left, equivalent to ::ORI_CCW.
    ORI_BELOW = 1,   //!< The forth point lies below the plane defined by the first three ::ORI_CW.
    ORI_ABOVE = -1     //!< The forth point lies above the plane defined by the first three ::ORI_CCW.
} ORIENTATION;

/// Maps an ::ORIENTATION to a string. Useful for debugging output.
inline const char *Orientation2Str(ORIENTATION value) {
    static const char *strs[] = {
        "ORI_NEGATIVE", "ORI_ZERO", "ORI_POSITIVE"
    };

    assert(value >= -1 && value <= 1);

    return strs[value + 1];
}

/// Wrapper around ::orient2d. Returns the result as an ::ORIENTATION.
static inline ORIENTATION sign_orient2d(const double* pa, const double* pb, const double* pc) {
    const double result = orient2d(pa, pb, pc);
    // A little convoluted but branchless.
    return (ORIENTATION) ((result > 0.0) - (result < 0.0));
}

/// Wrapper around ::orient3d. Returns the result as an ::ORIENTATION.
static inline ORIENTATION sign_orient3d(const double* pa, const double* pb, const double* pc, const double* pd) {
    const double result = orient3d(pa, pb, pc, pd);
    // A little convoluted but branchless.
    return (ORIENTATION) ((result > 0.0) - (result < 0.0));
}

static inline double abs_orient2d(const double* pa, const double* pb, const double* pc) {
    const double result = orient2d(pa, pb, pc);
    return (result >= 0.0)? result : -result;
}

static inline double abs_orient3d(const double* pa, const double* pb, const double* pc, const double* pd) {
    const double result = orient3d(pa, pb, pc, pd);
    return (result >= 0.0)? result : -result;
}

static inline char sign(const double& x) {return ((x > 0.0) - (x < 0.0));}


#ifdef __cplusplus
} // extern "C"
#endif

#include <HexHex/Utils/Typedefs.hh>
namespace HexHex
{

bool isInsideTet(Parameter u, Parameter v, Parameter w, Parameter x,  Parameter y);
bool isInsideTriangle(Parameter u, Parameter v, Parameter w, Parameter x);
bool isInsideEdge(Parameter u, Parameter v, Parameter w);
bool isInsideEdgeOrOnBoundary(Parameter u, Parameter v, Parameter w);


inline ORIENTATION ori2d(const Vec2d& u, const Vec2d& v, const Vec2d& w) {return sign_orient2d(u.data(), v.data(), w.data());}
inline ORIENTATION ori3d(const Vec3d& u, const Vec3d& v, const Vec3d& w, const Vec3d& x) {return sign_orient3d(u.data(), v.data(), w.data(), x.data());}
inline ORIENTATION ori3d(const Vec3i& u, const Vec3d& v, const Vec3d& w, const Vec3d& x) {return ori3d(Vec3d(u), v, w, x);}
inline ORIENTATION ori3d(const Vec3d& u, const Vec3i& v, const Vec3d& w, const Vec3d& x) {return ori3d(u, Vec3d(v), w, x);}
inline ORIENTATION ori3d(const Vec3d& u, const Vec3d& v, const Vec3i& w, const Vec3d& x) {return ori3d(u, v, Vec3d(w), x);}
inline ORIENTATION ori3d(const Vec3d& u, const Vec3d& v, const Vec3d& w, const Vec3i& x) {return ori3d(u, v, w, Vec3d(x));}
inline ORIENTATION ori3d(const Vec3d& u, const Vec3d& v, const Vec3i& w, const Vec3i& x) {return ori3d(u, v, Vec3d(w), Vec3d(x));}
inline ORIENTATION ori3d(const Vec3i& u, const Vec3d& v, const Vec3i& w, const Vec3i& x) {return ori3d(Vec3d(u), v, Vec3d(w), Vec3d(x));}
inline ORIENTATION ori3d(const std::array<Vec3d,3>& us, const Vec3d& x) {return sign_orient3d(us[0].data(), us[1].data(), us[2].data(), x.data());}
inline ORIENTATION ori3d(const std::array<Vec3d,3>& us, const Vec3i& x) {return ori3d(us[0], us[1], us[2], x);}
inline ORIENTATION ori3d(const std::array<Vec3d,4>& us) {return sign_orient3d(us[0].data(), us[1].data(), us[2].data(), us[3].data());}

inline bool areCCW(const std::vector<Vec2d>& pts)
{
    auto n = pts.size();
    if (n <= 2) return true;

    // each point must be to the left of each non-incident edge
    for (auto i1 = 0u; i1 < n; ++i1) {
        auto i2 = (i1+1)%n;
        for (auto j = 2u; j < n; ++j) {
            auto i3 = (i1+j)%n;
            if (ori2d(pts[i1], pts[i2], pts[i3]) == ORI_CW) {
                return false;
            }
        }
    }
    return true;
}

class PredicatesInitalizer
{
    PredicatesInitalizer();
    static PredicatesInitalizer instance;
};
}

