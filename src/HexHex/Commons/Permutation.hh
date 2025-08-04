#pragma once


#include <HexHex/Utils/Typedefs.hh>
namespace HexHex
{
class Permutation
{
public:
    Permutation(char idx) : idx(idx) {assert(idx >= 0 && idx <= 5);};
    Permutation(const Vec3i& p) : idx(vecToChar(p)) {}
    Permutation(int p0, int p1, int p2) : Permutation(Vec3i(p0, p1, p2)) {}
    Permutation() : idx(0) {}

    template<typename Vec3T>
    static Permutation descending(const Vec3T& v)
    {
        auto p = Vec3i(0, 1, 2);
        if (v[p[2]] > v[p[1]]) {std::swap(p[2], p[1]);}
        if (v[p[1]] > v[p[0]]) {std::swap(p[1], p[0]);}
        if (v[p[2]] > v[p[1]]) {std::swap(p[2], p[1]);}
        return Permutation(p);
    }

    static Permutation ascending(const Vec3d& v)
    {
        auto p = Vec3i(0, 1, 2);
        if (v[p[2]] < v[p[1]]) {std::swap(p[2], p[1]);}
        if (v[p[1]] < v[p[0]]) {std::swap(p[1], p[0]);}
        if (v[p[2]] < v[p[1]]) {std::swap(p[2], p[1]);}
        return Permutation(p);
    }

    inline const Vec3d permuted(const Vec3d& v) const {const auto& p = vector(); return Vec3d(v[p[0]], v[p[1]], v[p[2]]);};
    inline const Vec3i permuted(const Vec3i& v) const {const auto& p = vector(); return Vec3i(v[p[0]], v[p[1]], v[p[2]]);};
    inline const Permutation permuted(const Permutation& v) const {const auto& p = vector(); return Permutation(Vec3i(v[p[0]], v[p[1]], v[p[2]]));};

    inline const Permutation inverse() const {return Permutation(inverses[idx]);};

    inline const Vec3i& vector() const {assert(idx >= 0 && idx <= 5); return permutations[idx]; }

    inline int operator[] (int i) const {return vector()[i];}

    friend const Permutation operator*(const Permutation& lhs, const Permutation& rhs) {return Permutation(rhs.permuted(lhs));}

    friend bool operator==(const Permutation& p1, const Permutation& p2) {return p1.vector() == p2.vector();}
    friend bool operator!=(const Permutation& p1, const Permutation& p2) {return p1.vector() != p2.vector();}

    friend std::ostream& operator<<(std::ostream& os, const Permutation& p) { os << p.vector(); return os;}

private:
    char idx;

    char vecToChar(const Vec3i& v)
    {
        if (v[0] == 0) return (v[1] == 1)? 0 : 3;
        if (v[0] == 1) return (v[1] == 2)? 2 : 5;
        return (v[1] == 0)? 1 : 4;
    }

    constexpr const static std::array<Vec3i, 6> permutations = {Vec3i(0,1,2),
                                                                Vec3i(2,0,1),
                                                                Vec3i(1,2,0),
                                                                Vec3i(0,2,1),
                                                                Vec3i(2,1,0),
                                                                Vec3i(1,0,2)};

    constexpr const static std::array<char, 6> inverses = {0, 2, 1, 3, 4, 5};
}; // Permutation
}

