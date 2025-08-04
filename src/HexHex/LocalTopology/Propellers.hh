#pragma once

#include <HexHex/Commons/FlatMap.hh>
#include <HexHex/Commons/MeshElement.hh>
#include <HexHex/Commons/Transition.hh>
#include <HexHex/Utils/Typedefs.hh>

namespace HexHex
{

class HexVertexGenerator;

struct LocalPropellerHandle
{
private:
    int idx_;
public:
    LocalPropellerHandle(int idx=-1) : idx_(idx) {}
    inline bool is_valid() const { return idx_ != -1; }
    constexpr const int& idx() const { return idx_; }
    int& idx_mutable() {return idx_;}
    LocalPropellerHandle& operator=(int idx) {idx_ = idx; return *this;}
    constexpr bool operator==(const LocalPropellerHandle& _h) const { return _h.idx_ == this->idx_; }
    constexpr bool operator!=(const LocalPropellerHandle& _h) const { return _h.idx_ != this->idx_; }

    template <typename H>
    friend H AbslHashValue(H h, const LocalPropellerHandle& lph) {
        return H::combine(std::move(h), lph.idx_);
    }
};

struct GeneratorHandle
{
private:
    int idx_;
public:
    GeneratorHandle(int idx=-1) : idx_(idx) {}
    inline bool is_valid() const { return idx_ != -1; }
    constexpr const int& idx() const { return idx_; }
    GeneratorHandle& operator=(int idx) {idx_ = idx; return *this;}
};

struct LocalCornerHandle
{
private:
    int idx_;
public:
    LocalCornerHandle(int idx=-1) : idx_(idx) {}
    inline bool is_valid() const { return idx_ != -1; }
    constexpr const int& idx() const { return idx_; }
    LocalCornerHandle& operator=(int idx) {idx_ = idx; return *this;}
    constexpr bool operator==(const LocalCornerHandle& _h) const { return _h.idx_ == this->idx_; }
    constexpr bool operator!=(const LocalCornerHandle& _h) const { return _h.idx_ != this->idx_; }
};

using PropellerImages = FlatMap<CellHandle, Direction>;

class LocalPropeller
{
public:
    MeshElement holder;

    //bool has_incident_non_identity_transitions;

    LocalPropeller(HexVertexGenerator& generator, MeshElement holder);

    void setDirection(CellHandle cch, Direction dir) {images.set(cch, dir);}

    Direction getDirection(CellHandle cch);

    bool hasDirection(CellHandle cch, Direction dir) const;

    void initBlade(HexVertexGenerator& generator, MeshElement casing, CellHandle ch, Direction dir);

    void setBlade(int32 i, LocalPropellerHandle lph, CellHandle to, Transition transition);

    int32 findBladeIndex(CellHandle ch, Direction dir) const;
    int32 findBladeIndex(LocalPropellerHandle lph2) const;

    inline const FlatMap<CellHandle, Direction>& getBladeDirections(int32 j) const {return blades_from[j];}
    inline int nBlades() const {return blades_casings.size();}
    inline LocalPropellerHandle getBlade(int32 j) const {return blades[j];}
    inline bool bladeIsConnected(int32 j) const {return j < (int32)(blades.size()) && blades[j].is_valid();}
    inline LocalPropellerHandle getNextBlade(LocalPropellerHandle lph) const {int j = findBladeIndex(lph); return (j==-1)? LocalPropellerHandle(-1) : getBlade((j+1)%nBlades());}

    inline MeshElement getCasing(const LocalPropellerHandle& lph2) const {return blades_casings[findBladeIndex(lph2)];}

private:
    PropellerImages images;

    std::vector<MeshElement> blades_casings; // face or cell
    std::vector<FlatMap<CellHandle, Direction>> blades_from;
    std::vector<std::pair<CellHandle, Transition>> blades_to; // transition goes from casing (1st hf if is face) to cell stored here
    std::vector<LocalPropellerHandle> blades;
};

using LocalCorner = std::array<LocalPropellerHandle,3>;

struct GlobalPropellerHandle
{
    VertexHandle vh;
    LocalPropellerHandle lph;
};

inline bool operator==(const GlobalPropellerHandle& gph1, const GlobalPropellerHandle& gph2)
{
    return gph1.vh == gph2.vh && gph1.lph == gph2.lph;
}

typedef struct
{
    VertexHandle vh;
    LocalCornerHandle lch;
} GlobalCornerHandle;

typedef struct
{
    GlobalPropellerHandle gph;
    int32 offset;
} GlobalPropellerOpposite;

typedef struct
{
    LocalCornerHandle lch;
    LocalPropellerHandle forward;
} LocalCornerAxis;

using LPH = LocalPropellerHandle;
using LCH = LocalCornerHandle;
using GPH = GlobalPropellerHandle;
using GCH = GlobalCornerHandle;
using GH = GeneratorHandle;


inline GPH getGPH(const GCH& gch, const LPH lph)
{
    return {gch.vh, lph};
}

}

