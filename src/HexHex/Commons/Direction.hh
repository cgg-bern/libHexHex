#pragma once

#include <HexHex/Utils/Typedefs.hh>

namespace HexHex
{
class Direction
{
public:
    Direction(Vec3d direction) : dir(vecToChar(direction)) {}
    Direction(char dir) : dir(dir) {assert(dir >= 0 && dir <= 5);}
    Direction(int x, int y, int z) : Direction(Vec3d(x,y,z)) {}

    const Vec3d vector() const { assert(dir >= 0 && dir <= 5); return directions[dir]; }
    inline char idx() const {return dir;}
    inline char axis() {return dir>>1;}
    Direction flipped() const {return Direction(dir^1);}

    static inline std::array<std::array<char,4>,6> ORTHOGONAL_DIRECTIONS_PER_DIRECTION_CCW = {{
        {2, 4, 3, 5},
        {2, 5, 3, 4},
        {0, 5, 1, 4},
        {0, 4, 1, 5},
        {0, 2, 1, 3},
        {0, 3, 1, 2}
    }};
    static inline std::array<std::array<char,6>,6> ORTHOGONAL_INDICES_PER_DIRECTION_CCW = {{
        {-1, -1, 0, 2, 1, 3},
        {-1, -1, 0, 2, 3, 1},
        {0, 2, -1, -1, 3, 1},
        {0, 2, -1, -1, 1, 3},
        {0, 2, 1, 3, -1, -1},
        {0, 2, 3, 1, -1, -1}
    }};

    friend const Vec3d operator+(const Vec3d& vec, Direction dir) {return vec + dir.vector();}
    friend const Direction operator-(Direction dir) {return Direction(dir.idx()%2==0? dir.idx()+1 : dir.idx()-1);}
    friend const Vec3d operator+(Direction dir, const Vec3d& vec) { return vec + dir; }
    friend const Direction operator%(Direction dir1, Direction dir2) { return Direction(dir1.vector() % dir2.vector()); }
    friend const Vec3d operator-(const Vec3d& vec, Direction dir) {return vec + dir.flipped(); }
    friend const Vec3d operator*(Direction dir, double d) { return dir.vector() * d; }
    friend const Vec3d operator*(double d, Direction dir) { return dir * d; }

    friend double operator|(const Direction& dir1, Direction dir2) { return dir1.vector() | dir2.vector(); }
    friend double operator|(const Vec3d& vec, Direction dir) { return vec | dir.vector(); }
    friend double operator|(Direction dir, const Vec3d& vec) { return vec | dir; }
    friend bool operator==(const Direction& dir1, const Direction& dir2) { return dir1.dir == dir2.dir; }
    friend bool operator!=(const Direction& dir1, const Direction& dir2) { return dir1.dir != dir2.dir; }
    friend std::ostream& operator<<(std::ostream& os, const Direction& dir) { os << dir.vector(); return os; }

private:

    template<typename Vec3T>
    char vecToChar(const Vec3T& direction)
    {
        if (direction[0] != 0) return (direction[0] > 0)? 0 : 1;
        if (direction[1] != 0) return (direction[1] > 0)? 2 : 3;
        if (direction[2] != 0) return (direction[2] > 0)? 4 : 5;
        assert(false);
        return -1;
    }

    char dir;

    constexpr const static std::array<Vec3d, 6> directions = {
        Vec3d( 1, 0, 0),
        Vec3d(-1, 0, 0),
        Vec3d( 0, 1, 0),
        Vec3d( 0,-1, 0),
        Vec3d( 0, 0, 1),
        Vec3d( 0, 0,-1)
    };
}; // Direction

const inline std::vector<Direction> getAll6Directions() {return {0, 1, 2, 3, 4, 5};}
const inline std::vector<Direction> getAll3PositiveDirections() {return {0, 2, 4};}
const inline std::vector<Direction> getAll3NegativeDirections() {return {1, 3, 5};}

inline std::vector<Direction> getAllOrthogonalDirectionsCCW(Direction dir)
{
    const auto& dirs = Direction::ORTHOGONAL_DIRECTIONS_PER_DIRECTION_CCW[dir.idx()];
    return {dirs[0], dirs[1], dirs[2], dirs[3]};
}

// inline std::vector<Direction> getAllOrthogonalDirectionsCW(Direction dir)
// {
//     const auto& dirs = Direction::ORTHOGONAL_DIRECTIONS_PER_DIRECTION_CCW[dir.index()];
//     return {dirs[3], dirs[2], dirs[1], dirs[0]};
// }
}

