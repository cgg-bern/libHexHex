#pragma once

#include <HexHex/Utils/Typedefs.hh>
#include <iostream>
#include <ostream>
#include <string>

namespace HexHex
{

class HexHexMessage
{
public:
    enum class Severity : char
    {
        Info,
        Error,
        Warning
    };

    HexHexMessage(const Severity& severity, const std::string& msg) : severity(severity), msg(msg) {}

    HexHexMessage() : HexHexMessage(Severity::Info, "") {}

    ~HexHexMessage() {}

    inline bool isOk() const {return severity== Severity::Info;}
    inline bool isError() const {return severity== Severity::Error;}

    inline const std::string& message() const {return msg;}

    inline void print() const
    {
        std::cerr << "HexHex: " << ((severity== Severity::Info)? "Info: " :
            (severity== Severity::Error)? "Error: " : "Warning: ") << msg << std::endl;
    }

protected:
    Severity severity;
    std::string msg;
};
using HHM = HexHexMessage;

struct HHM_FlippedOrDegenerateTet : public HHM
{
    HHM_FlippedOrDegenerateTet() : HHM(Severity::Error, "Invalid Input - Flipped or Degenerate Tet") {};
};

struct HHM_LargeSanitizationFix : public HHM
{
    HHM_LargeSanitizationFix() : HHM(Severity::Warning, "Potentially Invalid Input - Large Sanitization Fix!") {};
};

struct HHM_FoundNoOppositePropeller : public HHM {
    HHM_FoundNoOppositePropeller() : HHM(Severity::Warning, "Failure - FoundNoOppositePropeller!") {};
};

struct HHM_PropellerOppositesDifferentBlades : HHM {
    HHM_PropellerOppositesDifferentBlades(const size_t l1, const size_t l2) : HHM(Severity::Warning, "Failure - Two opposing Propellers have a different number of blades!") {};
};

struct HHM_FoundNoOppositeBlade : public HHM {
    HHM_FoundNoOppositeBlade() : HHM(Severity::Warning, "Failure - FoundNoOppositeBlade!") {};
};

struct HHM_IsolatedTetVertex : public HHM {
    HHM_IsolatedTetVertex() : HHM(Severity::Error, "Invalid Input - IsolatedTetVertex!") {};
};

struct HHM_MissingHexVertex : public HHM {
    HHM_MissingHexVertex() : HHM(Severity::Error, "Failure - Missed Point in Vertex Extraction!") {};
};

struct HHM_DuplicateHexVertex : public HHM
{
    HHM_DuplicateHexVertex(const CellHandle& ch)
    {
        this->severity = Severity::Error;
        std::stringstream ss;
        ss << "Failure - Cell " << ch << ": Duplicate Point in Vertex Extraction!";
        this->msg = ss.str();
    };
};

}

