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


#include "DerivedExactPredicates.hh"
#include "Direction.hh"
#include "ExactPredicates.hh"
#include "Typedefs.hh"


namespace  HexEx {

bool isInside(Parameter u, Parameter v, Parameter w, Parameter x)
{
    // for the 2d checks we pretend we are in 2d
    // if u, v and w are collinear in the first two dimension
    // we use the second two instead

    //Todo: geiler machen
    int offset = 0;
    if (sign_orient2d(u.data(),v.data(),w.data()) == ORI_COLLINEAR)
        offset = 1;
    if (sign_orient2d(u.data()+offset,v.data()+offset,w.data()+offset) == ORI_COLLINEAR)
    {
        offset = 0;
        std::swap(u[1], u[2]);
        std::swap(v[1], v[2]);
        std::swap(w[1], w[2]);
        std::swap(x[1], x[2]);
    }

    ORIENTATION side = ORI_LEFT;
    if (sign_orient2d(u.data()+offset,v.data()+offset,w.data()+offset) == ORI_CW)
        side = ORI_RIGHT;

    return (sign_orient3d(u.data(), v.data(), w.data(), x.data()) == ORI_ZERO) &&
           (sign_orient2d(u.data()+offset, v.data()+offset, x.data()+offset)) == side &&
           (sign_orient2d(v.data()+offset, w.data()+offset, x.data()+offset)) == side &&
           (sign_orient2d(w.data()+offset, u.data()+offset, x.data()+offset)) == side;
}

bool isInsideOrOnBoundary(Parameter u, Parameter v, Parameter w, Parameter x)
{
    if (isInside(u,v,w,x))
        return true;
    if (isOnLine(u,v,x))
        return true;
    if (isOnLine(v,w,x))
        return true;
    if (isOnLine(w,u,x))
        return true;
    return false;
}

bool isInside(Parameter u, Parameter v, Parameter w, Parameter x, Parameter y)
{
    return (sign_orient3d(u.data(), v.data(), w.data(), y.data()) == ORI_ABOVE) &&
           (sign_orient3d(u.data(), w.data(), x.data(), y.data()) == ORI_ABOVE) &&
           (sign_orient3d(u.data(), x.data(), v.data(), y.data()) == ORI_ABOVE) &&
           (sign_orient3d(v.data(), x.data(), w.data(), y.data()) == ORI_ABOVE);

}

bool isInsideOrOnBoundary(Parameter u, Parameter v, Parameter w, Parameter x, Parameter y)
{

    return (sign_orient3d(u.data(), v.data(), w.data(), y.data()) != ORI_BELOW) &&
           (sign_orient3d(u.data(), w.data(), x.data(), y.data()) != ORI_BELOW) &&
           (sign_orient3d(u.data(), x.data(), v.data(), y.data()) != ORI_BELOW) &&
           (sign_orient3d(v.data(), x.data(), w.data(), y.data()) != ORI_BELOW);
}

bool isInPlane(Parameter u, Parameter v, Parameter w, Parameter x)
{
    return (sign_orient3d(u.data(), v.data(), w.data(), x.data()) == ORI_ZERO);
}

bool isOnLine(Parameter u, Parameter v, Parameter w)
{
    if (sign_orient2d(u.data(), v.data(), w.data()) != ORI_COLLINEAR)
        return false;
    if (sign_orient2d(u.data()+1, v.data()+1, w.data()+1) != ORI_COLLINEAR)
        return false;
    std::swap(u[1],u[2]);
    std::swap(v[1],v[2]);
    std::swap(w[1],w[2]);
    if (sign_orient2d(u.data(), v.data(), w.data()) != ORI_COLLINEAR)
        return false;
    return true;
}

// check if w is on line between u and v
bool isOnLineBetween(Parameter u, Parameter v, Parameter w)
{
//    return isOnLine(u,v,w); // todo: richtig machen

    if (isOnLine(u,v,w))
    {
        bool between0 = false;
        bool between1 = false;
        bool between2 = false;
        if (u[0] < w[0] && w[0] < v[0])
            between0 = true;
        if (u[1] < w[1] && w[1] < v[1])
            between1 = true;
        if (u[2] < w[2] && w[2] < v[2])
            between2 = true;
        if (u[0] > w[0] && w[0] > v[0])
            between0 = true;
        if (u[1] > w[1] && w[1] > v[1])
            between1 = true;
        if (u[2] > w[2] && w[2] > v[2])
            between2 = true;


//        if (!(between0 || between1 || between2))
//            std::cout << "Warning: isOnLineBetween different than before" << std::endl;
        return between0 || between1 || between2;
    }

    return false;
}

bool isFlipped(Parameter u, Parameter v, Parameter w, Parameter x)
{
    return sign_orient3d(u.data(), v.data(), w.data(), x.data()) == ORI_BELOW;
}

bool pointsIntoTet(Parameter u, Parameter v, Parameter w, Parameter x, Direction y)
{
    auto z = u + y;
    return isRegular(u,v,w,z) &&
           isRegular(u,x,v,z) &&
           isRegular(u,w,x,z);
}

bool pointsIntoTetRelaxed(Parameter u, Parameter v, Parameter w, Parameter x, Direction y)
{
    auto z = u + y;
    return !isFlipped(u,v,w,z) &&
           !isFlipped(u,x,v,z) &&
           !isFlipped(u,w,x,z);
}

bool pointsIntoTriangle(Parameter u, Parameter v, Parameter w, Direction x)
{
    auto y = u+x;

    // for the 2d checks we pretend we are in 2d
    // if u, v and w are collinear in the first two dimension
    // we use the second two instead

    //Todo: geiler machen
    int offset = 0;
    if (sign_orient2d(u.data(),v.data(),w.data()) == ORI_COLLINEAR)
        offset = 1;
    if (sign_orient2d(u.data()+offset,v.data()+offset,w.data()+offset) == ORI_COLLINEAR)
    {
        offset = 0;
        std::swap(u[1], u[2]);
        std::swap(v[1], v[2]);
        std::swap(w[1], w[2]);
        std::swap(y[1], y[2]);
    }

    ORIENTATION side = ORI_LEFT;
    if (sign_orient2d(u.data()+offset,v.data()+offset,w.data()+offset) == ORI_CW)
        side = ORI_RIGHT;

    return (sign_orient3d(u.data(), v.data(), w.data(), y.data()) == ORI_ZERO) &&
           (sign_orient2d(u.data()+offset, v.data()+offset, y.data()+offset)) == side &&
           (sign_orient2d(w.data()+offset, u.data()+offset, y.data()+offset)) == side;
}

/// x is the point on the edge uv
bool pointsIntoProjectedTriangleFromEdge(Parameter u, Parameter v, Parameter w, Parameter x, Direction dir)
{

    auto y = x+dir;

    // for the 2d checks we pretend we are in 2d
    // if u, v and w are collinear in the first two dimension
    // we use the second two instead

    //Todo: geiler machen
    int offset = 0;
    if (sign_orient2d(u.data(),v.data(),w.data()) == ORI_COLLINEAR)
        offset = 1;
    if (sign_orient2d(u.data()+offset,v.data()+offset,w.data()+offset) == ORI_COLLINEAR)
    {
        offset = 0;
        std::swap(u[1], u[2]);
        std::swap(v[1], v[2]);
        std::swap(w[1], w[2]);
        std::swap(x[1], x[2]);
        std::swap(y[1], y[2]);
    }
    if (sign_orient2d(u.data()+offset,v.data()+offset,w.data()+offset) == ORI_COLLINEAR)
    {
        // triangle is a line
        return false;
    }


    ORIENTATION ori1 = sign_orient2d(x.data()+offset, y.data()+offset, u.data()+offset);
    if (ori1 == ORI_ZERO)
        return false;
    ORIENTATION ori2 = sign_orient2d(x.data()+offset, v.data()+offset, y.data()+offset);
    if (ori2 == ORI_ZERO)
        return false;

    auto ori3 = sign_orient2d(x.data()+offset, w.data()+offset, u.data()+offset);
    auto ori4 = sign_orient2d(x.data()+offset, v.data()+offset, w.data()+offset);

    return (ori3 == ori1) &&
           (ori4 == ori2) &&
           (sign_orient3d(u.data(), v.data(), w.data(), y.data()) == ORI_ZERO);

}

bool isDegenerate(Parameter u, Parameter v, Parameter w, Parameter x)
{
    return sign_orient3d(u.data(), v.data(), w.data(), x.data()) == ORI_ZERO;
}

bool isRegular(Parameter u, Parameter v, Parameter w, Parameter x)
{
    return sign_orient3d(u.data(), v.data(), w.data(), x.data()) == ORI_ABOVE;
}

// false if intersection lies on edge, false if line in plane
bool intersectsTriangle(Parameter u, Parameter v, Parameter w, Parameter x, Parameter y)
{
    // intersection test depends on the direction. x should be above and y below the triangle uvw

    return  isRegular(u,v,w,x) &&  // check if x is on the one side
            isRegular(u,w,v,y) &&  // and y on the other
            isRegular(u,v,y,x) &&
            isRegular(v,w,y,x) &&
            isRegular(w,u,y,x);
}

// true if intersection lies on edge
bool intersectsTriangleRelaxed(Parameter u, Parameter v, Parameter w, Parameter x, Parameter y)
{
    // intersection test depends on the direction. x should be above and y below the triangle uvw

    return  (isRegular(u,v,w,x) || isDegenerate(u,v,w,x)) &&  // check if x is on the one side or on the face itself
            isRegular(u,w,v,y) &&  // and y on the other
            !isFlipped(u,v,y,x) &&
            !isFlipped(v,w,y,x) &&
            !isFlipped(w,u,y,x);
}

PredicatesInitalizer PredicatesInitalizer::instance;

PredicatesInitalizer::PredicatesInitalizer()
{
    exactinit();
}



bool areCoLinear(Parameter u, Parameter v, Direction dir)
{
    auto d = dir.vector();
    if (d[0] != 0)
    {
        if ((u[1] != v[1]) || (u[2] != v[2]))
            return false;
        if ((d[0] > 0) != (u[0] < v[0]))
            return false;
        return true;
    }
    else if (d[1] != 0)
    {
        if ((u[0] != v[0]) || (u[2] != v[2]))
            return false;
        if ((d[1] > 0) != (u[1] < v[1]))
            return false;
        return true;
    }
    if (d[2] != 0)
    {
        if ((u[0] != v[0]) || (u[1] != v[1]))
            return false;
        if ((d[2] > 0) != (u[2] < v[2]))
            return false;
        return true;
    }
    else
    {
        assert(false);
        return false;
    }
}

bool isRegular(Vec2d u, Vec2d v, Vec2d w)
{
    return sign_orient2d(u.data(), v.data(), w.data()) == ORI_CCW;
}


}
