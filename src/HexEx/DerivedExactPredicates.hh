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

#include "Typedefs.hh"

namespace HexEx {

class Direction;

class PredicatesInitalizer
{
    PredicatesInitalizer();
    static PredicatesInitalizer instance;
};

bool isInside(Parameter u, Parameter v, Parameter w, Parameter x);
bool isInsideOrOnBoundary(Parameter u, Parameter v, Parameter w, Parameter x);
bool isInside(Parameter u, Parameter v, Parameter w, Parameter x,  Parameter y);
bool isInsideOrOnBoundary(Parameter u, Parameter v, Parameter w, Parameter x,  Parameter y);
bool isInPlane(Parameter u, Parameter v, Parameter w, Parameter x);
bool isOnLine(Parameter u, Parameter v, Parameter w);
bool isOnLineBetween(Parameter u, Parameter v, Parameter w);
bool isFlipped(Parameter u, Parameter v, Parameter w, Parameter x);

/// only works if dir is one of the six main directions
bool areCoLinear(Parameter u, Parameter v, Direction dir);

/// Checks if vector y points from u into tet uvwx. False if y lies in a face.
bool pointsIntoTet(Parameter u, Parameter v, Parameter w, Parameter x,  Direction y);
/// Checks if vector y points from u into tet uvwx. True if y lies in a face or an edge.
bool pointsIntoTetRelaxed(Parameter u, Parameter v, Parameter w, Parameter x,  Direction y);
bool pointsIntoTriangle(Parameter u, Parameter v, Parameter w, Direction x);

bool pointsIntoProjectedTriangleFromEdge(Parameter u, Parameter v, Parameter w, Parameter x,  Direction dir);

bool isDegenerate(Parameter u, Parameter v, Parameter w, Parameter x);
bool isRegular(Parameter u, Parameter v, Parameter w, Parameter x);
bool intersectsTriangle(Parameter u, Parameter v, Parameter w, Parameter x, Parameter y);
bool intersectsTriangleRelaxed(Parameter u, Parameter v, Parameter w, Parameter x, Parameter y);

bool isRegular(Vec2d u, Vec2d v, Vec2d w);

}
