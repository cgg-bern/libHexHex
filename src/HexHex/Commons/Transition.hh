#pragma once

#include <HexHex/Commons/Direction.hh>
#include <HexHex/Utils/Typedefs.hh>
#include <HexHex/Commons/Matrix4x4T.hh>
#include <random>

namespace HexHex
{

/*
Rotation Index is different from AlgoHex' Quaternions to HexHex' Restricted Rotations
0 -> 0
1 -> 1
2 -> 4
3 -> 5
4 -> 3
5 -> 2
6 -> 18
7 -> 22
8 -> 14
9 -> 11
10 -> 10
11 -> 6
12 -> 19
13 -> 15
14 -> 7
15 -> 23
16 -> 16
17 -> 9
18 -> 13
19 -> 21
20 -> 12
21 -> 17
22 -> 20
23 -> 8
*/
const std::array<char, 24> ALGOHEX_TO_HEXEX_TRANSITION_ROTATION_INDICES = {0,1,4,5,3,2,18,22,14,11,10,6,19,15,7,23,16,9,13,21,12,17,20,8};

class RestrictedRotation
{
private:
    class MultiplicationMapInitializer
    {
    public:
        MultiplicationMapInitializer()
        {

            for (auto i = 0; i < 24; ++i)
                matrices[i] = convertToMatrix(i);


            for (auto i = 0; i < 24; ++i)
                for (auto j = 0; j < 24; ++j)
                {
                    auto lhs = RestrictedRotation(i);
                    auto rhs = RestrictedRotation(j);
                    auto res = RestrictedRotation(lhs.toMatrix()*rhs.toMatrix());
                    multiplicationTable[24*i+j] = res.i;
                }

            inverseMap[ 0] =  0;
            inverseMap[ 1] =  1;
            inverseMap[ 2] =  3;
            inverseMap[ 3] =  2;
            inverseMap[ 4] =  4;
            inverseMap[ 5] =  5;
            inverseMap[ 6] =  6;
            inverseMap[ 7] =  7;
            inverseMap[ 8] = 16;
            inverseMap[ 9] = 20;
            inverseMap[10] = 10;
            inverseMap[11] = 14;
            inverseMap[12] = 21;
            inverseMap[13] = 17;
            inverseMap[14] = 11;
            inverseMap[15] = 15;
            inverseMap[16] =  8;
            inverseMap[17] = 13;
            inverseMap[18] = 22;
            inverseMap[19] = 19;
            inverseMap[20] =  9;
            inverseMap[21] = 12;
            inverseMap[22] = 18;
            inverseMap[23] = 23;

        }
    };

public:
    static const RestrictedRotation IDENTITY;

    explicit RestrictedRotation(char i) : i(i) {}

    RestrictedRotation() : i(0) {}

    explicit RestrictedRotation(const Matrix4x4d m)
    {
        auto p = 0;
        if (m(1,0) != 0)
            p = 1;
        else if (m(2,0) != 0)
            p = 2;

        auto nb = m(p,0) < 0 ? 1 : 0;

        auto p2 = (p+1)%3;
        auto p2b = 0;

        if (m(p2,1) == 0)
        {
            p2 = ((p-1)+3)%3;
            p2b = 1;
        }

        auto n2b = 0;

        if (m(p2,1) < 0)
            n2b = 1;

        i = (p << 3) + (nb << 2) + (p2b << 1) + n2b;
    }

    bool operator==(const RestrictedRotation& other) const { return i == other.i; }
    bool operator!=(const RestrictedRotation& other) const { return i != other.i; }

    RestrictedRotation operator*(const RestrictedRotation& rhs) const {return RestrictedRotation(multiplicationTable[24*i+rhs.i]);}

    inline int index() const {return i;}

    inline Vec3d transform(const Vec3d& v) const
    {
        return transformT(i, v);
    }

    inline Vec3i transform(const Vec3i& v) const
    {
        return transformT(i, v);
    }

    template<typename Vec3T>
    inline Vec3T transformVec3T(const Vec3T& v) const
    {
        return transformT(i, v);
    }

    template<typename Mat3x3T>
    inline Mat3x3T transformMat3x3T(const Mat3x3T& m) const
    {
        Mat3x3T res;
        for (unsigned char i = 0; i < 3; ++i)
        {
            Vec3d w = transform(Vec3d(m(0,i), m(1,i), m(2,i)));
            res(0,i) = w[0];
            res(1,i) = w[1];
            res(2,i) = w[2];
        }
        return res;
    }

    inline Direction transform(const Direction& d) const
    {
        return Direction(transformT(i, d.vector()));
    }

    inline std::vector<Vec3d> transform(const std::vector<Vec3d>& vs) const
    {
        std::vector<Vec3d> res;
        res.reserve(vs.size());
        for (const auto& v : vs) res.push_back(transform(v));
        return res;
    }

    inline RestrictedRotation inverted() const
    {
        return RestrictedRotation(inverseMap[i]);
    }

    inline void invert()
    {
        i = inverseMap[i];
    }

    inline Matrix4x4d toMatrix() const
    {
        return matrices[i];
    }

private:
    template <typename Vec>
    static Vec transformT(char i, const Vec& v)
    {
        switch(i)
        {
        case 0:  return Vec( v[0],  v[1],  v[2]);
        case 1:  return Vec( v[0], -v[1], -v[2]);
        case 2:  return Vec( v[0], -v[2],  v[1]);
        case 3:  return Vec( v[0],  v[2], -v[1]);
        case 4:  return Vec(-v[0],  v[1], -v[2]);
        case 5:  return Vec(-v[0], -v[1],  v[2]);
        case 6:  return Vec(-v[0],  v[2],  v[1]);
        case 7:  return Vec(-v[0], -v[2], -v[1]);
        case 8:  return Vec( v[2],  v[0],  v[1]);
        case 9:  return Vec(-v[2],  v[0], -v[1]);
        case 10: return Vec( v[1],  v[0], -v[2]);
        case 11: return Vec(-v[1],  v[0],  v[2]);
        case 12: return Vec(-v[2], -v[0],  v[1]);
        case 13: return Vec( v[2], -v[0], -v[1]);
        case 14: return Vec( v[1], -v[0],  v[2]);
        case 15: return Vec(-v[1], -v[0], -v[2]);
        case 16: return Vec( v[1],  v[2],  v[0]);
        case 17: return Vec(-v[1], -v[2],  v[0]);
        case 18: return Vec(-v[2],  v[1],  v[0]);
        case 19: return Vec( v[2], -v[1],  v[0]);
        case 20: return Vec( v[1], -v[2], -v[0]);
        case 21: return Vec(-v[1],  v[2], -v[0]);
        case 22: return Vec( v[2],  v[1], -v[0]);
        case 23: return Vec(-v[2], -v[1], -v[0]);
        default: assert(false); return v;
        }
    }

    char i;

    static Matrix4x4d convertToMatrix(char i)
    {
        auto m = Matrix4x4d();

        std::array<Vec3d, 3> vecs {Vec3d(0.), Vec3d(0.), Vec3d(0.)};

        auto p = i >> 3;
        vecs[0][p] = 1;

        if (i & 1<<2)
            vecs[0] *= -1.0;

        auto p2 = 0;
        if (i & 1<<1)
            p2 = p-1;
        else
            p2 = p+1;
        p2 = (p2+3)%3;

        vecs[1][p2] = 1;
        if (i & 1<<0)
            vecs[1] *= -1.0;

        vecs[2] = vecs[0] % vecs[1];

        for (auto i = 0u; i < 3; ++i)
        {
            for (auto j = 0u; j < 3; ++j)
                m(i,j) = vecs[j][i];
            m(3,i) = 0;
            m(i,3) = 0;
        }
        m(3,3) = 1;

        return m;
    }

    static std::array<char,24*24> multiplicationTable;
    static std::array<Matrix4x4d,24> matrices;
    static std::array<char, 24> inverseMap;
    static MultiplicationMapInitializer initializer;

}; // class RestrictedRotation

class Transition
{
public:
    using Translation = Vec3d;
    static const Transition IDENTITY;

    explicit Transition(int i = 0, const Translation& t = Translation(0,0,0)) : r(i), t(t) {}
    explicit Transition(RestrictedRotation r, const Translation& t = Translation(0,0,0)) : r(r), t(t) {}
    // explicit Transition(const Matrix4x4d& m) : r(m), t(Vec3i(0,0,0))
    // {
    //     for (auto i = 0u; i < 3; ++i)
    //         t[i] = (int)round(m(i,3));
    // }


    bool operator==(const Transition& other) const { return r == other.r && t == other.t; }
    bool operator!=(const Transition& other) const { return !operator==(other); }

    inline Transition operator*(const Transition& rhs)
    {
        return Transition(r * rhs.r, r.transform(rhs.t) + t);
    }

    inline Vec3d transform_point(const Vec3d& v) const {return r.transform(v) + Vec3d(t[0], t[1], t[2]);}
    inline Vec3d transform_vector(const Vec3d& v) const {return r.transform(v);}
    inline Direction transform(const Direction& d) const {return r.transform(d);}
    inline std::vector<Vec3d> transform_points(const std::vector<Vec3d>& vs) const
    {
        std::vector<Vec3d> res;
        res.reserve(vs.size());
        for (const auto& v : vs) res.push_back(transform_point(v));
        return res;
    }

    Transition inverted() const
    {
        auto res = *this;
        res.invert();
        return res;
    }

    void invert()
    {
        r.invert();
        t = r.transform(-t);
    }

    inline void setTranslation(const Translation& translation) { t = translation; }

    inline RestrictedRotation rotation() const {return r;}

    Matrix4x4d toMatrix() const
    {
        auto res = r.toMatrix();
        for (auto i = 0; i < 3; ++i)
            res(i,3) = t[i];
        return res;
    }

    /// Pass a RandomNumberEngine
    static Transition random(auto &rng)
    {
        int r = std::uniform_int_distribution(0, 24)(rng);
        auto pos_dist = std::uniform_int_distribution(-10, 10);
        int tx = pos_dist(rng);
        int ty = pos_dist(rng);
        int tz = pos_dist(rng);
        return Transition(r, Vec3d(tx,ty,tz));
    }

    friend std::ostream& operator<<(std::ostream& os, const Transition& gi)
    {
        os << gi.toMatrix() << std::endl;
        return os;
    }

private:

    RestrictedRotation r;
    Translation t;
};

} // namespace hexhex

