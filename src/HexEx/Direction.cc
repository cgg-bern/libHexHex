/*
 * Copyright 2016 Computer Graphics Group, RWTH Aachen University
 * Author: Max Lyon <lyon@cs.rwth-aachen.de>
 *
 * This file is part of HexEx.
 *
 * HexEx is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * HexEx is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HexEx.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "Direction.hh"
#include "Utils.hh"

namespace HexEx
{

    const std::array<Vec3d, 6> Direction::directions = { Vec3d( 1, 0, 0),
                                                         Vec3d(-1, 0, 0),
                                                         Vec3d( 0, 1, 0),
                                                         Vec3d( 0,-1, 0),
                                                         Vec3d( 0, 0, 1),
                                                         Vec3d( 0, 0,-1) };

std::vector<Direction> getAll6Directions()
{
    return { Direction( 1, 0, 0),
             Direction(-1, 0, 0),
             Direction( 0, 1, 0),
             Direction( 0,-1, 0),
             Direction( 0, 0, 1),
             Direction( 0, 0,-1) };
}

std::vector<Direction> getAllOrthogonalDirections(Direction dir, bool reverse)
{
    auto res = std::vector<Direction>();
    if ((dir | Direction(0,0,1)) != 0)
        if (dir.vector()[2] > 0)
             res = { Direction(1,0,0), Direction(0,1,0), Direction(-1,0,0), Direction(0,-1,0) };
        else
             res = { Direction(1,0,0), Direction(0,-1,0), Direction(-1,0,0), Direction(0,1,0) };
    else if ((dir | Direction(0,1,0)) != 0)
        if (dir.vector()[1] > 0)
             res = { Direction(1,0,0), Direction(0,0,-1), Direction(-1,0,0), Direction(0,0,1) };
        else
             res = { Direction(1,0,0), Direction(0,0,1), Direction(-1,0,0), Direction(0,0,-1) };
    else if ((dir | Direction(1,0,0)) != 0)
        if (dir.vector()[0] > 0)
             res = { Direction(0,1,0), Direction(0,0,1), Direction(0,-1,0), Direction(0,0,-1) };
        else
             res = { Direction(0,1,0), Direction(0,0,-1), Direction(0,-1,0), Direction(0,0,1) };
    else
        assert(false); // should not happen

    if (reverse)
        std::reverse(res.begin(), res.end());

    return res;
}

Direction::Direction(Vec3d direction)
    :
      dir(vecToChar(direction))
{
    assert(isValidDirection(direction.normalized()));
    assert(dir >= 0 && dir <= 5);
}

Direction::Direction(char dir)
    :
      dir(dir)
{
    assert(dir >= 0 && dir <= 5);
}

Direction::Direction(int x, int y, int z)
    :
      Direction(Vec3d(x,y,z))
{
    assert(dir >= 0 && dir <= 5);
}

void Direction::transform(Transition tranFun)
{
    dir = vecToChar(tranFun.transform_vector(vector()).normalized());
}

const Direction Direction::transformed(Transition tranFun) const
{
    auto newDir = *this;
    newDir.transform(tranFun);
    return newDir;
}

char Direction::vecToChar(Vec3d direction)
{
    assert(isValidDirection(direction));

    if (direction[0] != 0)
    {
        if (direction[0] > 0)
            return 0;
        else
            return 1;
    }
    else if (direction[1] != 0)
    {
        if (direction[1] > 0)
            return 2;
        else
            return 3;
    }
    else if (direction[2] != 0)
    {
        if (direction[2] > 0)
            return 4;
        else
            return 5;
    }
    else
    {
        assert(false);
        return -1;
    }
}


bool isValidDirection(Vec3d dir)
{
    bool b0 = (fabs(dir[0]) == 0) || (fabs(dir[0]) == 1);
    bool b1 = (fabs(dir[1]) == 0) || (fabs(dir[1]) == 1);
    bool b2 = (fabs(dir[2]) == 0) || (fabs(dir[2]) == 1);
    bool sum = fabs(dir[0]) + fabs(dir[1]) + fabs(dir[2]) == 1;
    return b0 && b1 && b2 && sum;
}

bool isInteger(Vec3d dir)
{
    return dir == roundVector(dir);
}

const Vec3d operator+(const Vec3d& vec, Direction dir)
{
    assert(isInteger(vec));
    return vec + dir.vector();
}

const Vec3d operator+(Direction dir, const Vec3d& vec) { return vec + dir; }
const Direction operator%(Direction dir1, Direction dir2) { return Direction(dir1.vector() % dir2.vector()); }
const Vec3d operator-(const Vec3d& vec, Direction dir) { dir.flip(); return vec + dir; }
const Vec3d operator*(Direction dir, double d) { return dir.vector() * d; }
const Vec3d operator*(double d, Direction dir) { return dir * d; }

double operator|(const Direction& dir1, Direction dir2) { return dir1.vector() | dir2.vector(); }
double operator|(const Vec3d& vec, Direction dir) { return vec | dir.vector(); }
double operator|(Direction dir, const Vec3d& vec) { return vec | dir; }
bool operator==(const Direction& dir1, const Direction& dir2) { return dir1.vector() == dir2.vector(); }
bool operator!=(const Direction& dir1, const Direction& dir2) { return dir1.vector() != dir2.vector(); }
std::ostream& operator<<(std::ostream& os, const Direction& dir) { os << dir.vector(); return os; }



}
