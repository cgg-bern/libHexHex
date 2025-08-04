#include <HexHex/Commons/Transition.hh>

namespace HexHex
{
const Transition Transition::IDENTITY = Transition(0);
const RestrictedRotation RestrictedRotation::IDENTITY = RestrictedRotation(0);
std::array<char,24*24> RestrictedRotation::multiplicationTable;
std::array<Matrix4x4d,24> RestrictedRotation::matrices;
std::array<char, 24> RestrictedRotation::inverseMap;
RestrictedRotation::MultiplicationMapInitializer RestrictedRotation::initializer;
}
