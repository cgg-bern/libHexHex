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


#pragma once


#include <ostream>
#include <array>

#include "Typedefs.hh"

namespace HexEx
{

class Direction
{
public:
    explicit Direction(Vec3d direction);
    explicit Direction(char direction);
    Direction(int x, int y, int z);
    const Vec3d vector() const { assert(dir >= 0 && dir <= 5); return directions[dir]; }
    void transform(Transition tranFun);
    const Direction transformed(Transition tranFun) const;
    void flip() { dir += ((dir % 2) == 0 ? 1 : -1); }
    Direction flipped() const { return Direction(dir + ((dir % 2) == 0 ? 1 : -1)); }

private:

    char vecToChar(Vec3d direction);

    char dir;

    const static std::array<Vec3d, 6> directions;
};

bool isValidDirection(Vec3d dir);
bool isInteger(Vec3d dir);

const Vec3d operator+(const Vec3d& vec, Direction dir);
const Vec3d operator+(Direction dir, const Vec3d& vec);
const Vec3d operator-(const Vec3d& vec, Direction dir);
const Direction operator%(Direction dir1, Direction dir2);
const Vec3d operator*(Direction dir, double d);
const Vec3d operator*(double d, Direction dir);

double operator|(const Direction& dir1, Direction dir2);
double operator|(const Vec3d& vec, Direction dir);
double operator|(Direction dir, const Vec3d& vec);
bool operator==(const Direction& dir1, const Direction& dir2);
bool operator!=(const Direction& dir1, const Direction& dir2);
std::ostream& operator<<(std::ostream& os, const Direction& dir);

std::vector<Direction> getAll6Directions();
std::vector<Direction> getAllOrthogonalDirections(Direction dir, bool reverse);

}
