#pragma once
/*
 * Copyright 2023 Computer Graphics Group, University of Bern - Tobias Kohler <tobias.kohler@unibe.ch>
 * Copyright 2016 Computer Graphics Group, RWTH Aachen University - Max Lyon <lyon@cs.rwth-aachen.de>
 *
 * This file is part of HexHex.
 */

#include <HexHex/Utils/Typedefs.hh>
#include <libTimekeeper/Duration.hh>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace HexHex {

// Returns the index of the largest absolute element
inline unsigned char argAbsMax(const Vec3d& n)
{
    if ((fabs(n[0]) > fabs(n[1])) && (fabs(n[0]) > fabs(n[2]))) return 0u;
    else if (fabs(n[1]) > fabs(n[2])) return 1u;
    else return 2u;
}

inline unsigned char argAbsMax(const Vec2d& v)
{
    return (fabs(v[0]) > fabs(v[1]))? 0u : 1u;
}

// Returns the index of the largest element
inline unsigned char argMax(const Vec3d& v)
{
    if (v[0] > v[1] && v[0] > v[2]) return 0u;
    else if (v[1] > v[2]) return 1u;
    else return 2u;
}

// Returns the index of the smallest element
inline unsigned char argMin(const Vec3d& v)
{
    if (v[0] < v[1] && v[0] < v[2]) return 0u;
    else if (v[1] < v[2]) return 1u;
    else return 2u;
}

inline std::string toLower(const std::string& s)
{
    std::string r = s;
    std::transform(r.begin(), r.end(), r.begin(),[](unsigned char c){return std::tolower(c);});
    return r;
}

inline Vec3d roundVector(const Vec3d& v) {return Vec3d(round(v[0]), round(v[1]), round(v[2]));}

inline bool isInteger(const Vec3d& v) {return (v == roundVector(v));}

inline Vec3d absVector(const Vec3d& vec) {return Vec3d((vec[0]>=0)? vec[0] : -vec[0], (vec[1]>=0)? vec[1] : -vec[1], (vec[2]>=0)? vec[2] : -vec[2]);}

inline double max3(const double& a, const double& b, const double& c) {return std::max(std::max(a,b),c);}

inline double min3(const double& a, const double& b, const double& c) {return std::min(std::min(a,b),c);}

inline double max4(const double& a, const double& b, const double& c, const double& d) {return std::max(std::max(a,b),std::max(c,d));}

inline double min4(const double& a, const double& b, const double& c, const double& d) {return std::min(std::min(a,b),std::min(c,d));}

inline uint32 min7u32(
    const uint32 a, const uint32 b, const uint32 c, const uint32 d,
    const uint32 e, const uint32 f, const uint32 g
)
{
    return std::min(std::min(std::min(a,b),std::min(c,d)),std::min(std::min(e,f),g));
}

inline char argmin4(const double& a, const double& b, const double& c, const double& d)
{
    if (a <= b && a <= c && a <= d) return 0;
    if (b <= c && b <= d) return 1;
    if (c <= d) return 2;
    return 3;
}

inline char argmin4i(const int a, const int b, const int c, const int d)
{
    if (a <= b && a <= c && a <= d) return 0;
    if (b <= c && b <= d) return 1;
    if (c <= d) return 2;
    return 3;
}

inline char argmin3(const double& a, const double& b, const double& c)
{
    if (a <= b && a <= c) return 0;
    if (b <= c) return 1;
    return 2;
}

inline bool isnan(const Vec3d& v) {return std::isnan(v[0]) || std::isnan(v[1]) || std::isnan(v[2]);}

inline bool isnan(const Vec2d& v) {return std::isnan(v[0]) || std::isnan(v[1]);}

inline double sceil(const double& x) {return (floor(x) + 1);}

inline double sfloor(const double& x) {return (ceil(x) - 1);}

inline bool isInBetween(const double& x, const double& a, const double& b, bool allowEquals) {return allowEquals? ((x>=a&&x<=b)||(x>=b&&x<=a)) : ((x>a&&x<b)||(x>b&&x<a));}

inline bool isApproxEqual(const Vec2d& v1, const Vec2d& v2, const double EPSILON) {return (v1-v2).sqrnorm() < (EPSILON*EPSILON);}

template <typename Vec>
Vec3d toVec3d(const Vec& vec)
{
    return Vec3d(vec[0], vec[1], vec[2]);
}

template <typename Vec>
Vec toVec(const Vec3d& vec)
{
    return Vec(vec[0], vec[1], vec[2]);
}

inline double toSeconds(const Timekeeper::Duration& dur)
{
    return (double)(dur.wallclock.count())*1e-6;
}

} // namespace hexhex

