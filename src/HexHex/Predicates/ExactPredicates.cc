#include <HexHex/Commons/Direction.hh>
#include <HexHex/Commons/MeshElement.hh>
#include <HexHex/Predicates/ExactPredicates.hh>

namespace HexHex
{
PredicatesInitalizer PredicatesInitalizer::instance;
PredicatesInitalizer::PredicatesInitalizer() {exactinit();} // This must be called before using the predicates, otherwise - BOOOM!!!

bool isInsideTet(Parameter u, Parameter v, Parameter w, Parameter x,  Parameter y)
{
    return (sign_orient3d(u.data(), v.data(), w.data(), y.data()) == ORI_ABOVE) &&
           (sign_orient3d(u.data(), w.data(), x.data(), y.data()) == ORI_ABOVE) &&
           (sign_orient3d(u.data(), x.data(), v.data(), y.data()) == ORI_ABOVE) &&
           (sign_orient3d(v.data(), x.data(), w.data(), y.data()) == ORI_ABOVE);

}

bool isInsideTriangle(Parameter u, Parameter v, Parameter w, Parameter x)
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

bool isInsideEdgeOrOnBoundary(Parameter u, Parameter v, Parameter w)
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


bool isInsideEdge(Parameter u, Parameter v, Parameter w)
{
    //    return isOnLine(u,v,w); // todo: richtig machen

    if (isInsideEdgeOrOnBoundary(u,v,w))
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

}
