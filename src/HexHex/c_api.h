#pragma once

#include <HexHex/Config/Export.hh>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

struct HexMesh {
    size_t n_vertices;
    double *vpos;
    size_t n_hexes;
    uint32_t *hex_vertices; // TODO: specify order per hex (VTK?)
};
struct HexMesh HEXHEX_EXPORT hexhex_extract_simple(
        size_t n_vertices,
        const double *vpos,            // length: 3 * n_vertices
        size_t n_tets,
        const uint32_t *tet_vertices, // length: 4 * n_cells
        const double *tet_igm         // length: 3 * 4 * n_cells
        );

// Frees the buffers allocated by hexhex_extract_simple. Must be called on every
// HexMesh returned by that function to avoid memory leaks.
void HEXHEX_EXPORT hexhex_free_hex_mesh(struct HexMesh *hm);

#ifdef __cplusplus
}
#endif

