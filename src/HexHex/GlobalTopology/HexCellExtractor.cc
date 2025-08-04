#include <HexHex/GlobalTopology/HexCellExtractor.hh>
#include <HexHex/Utils/Stopwatches.hh>
#include <HexHex/Utils/Threading.hh>
#include <HexHex/HexExtractor.hh>
#include <HexHex/LocalTopology/HexVertexGenerator.hh>

namespace HexHex
{

void HexCellExtractor::extractHexCells()
{
    ScopedStopWatch _{sw::extractHexFacesAndCells};

    std::chrono::steady_clock::time_point t1, t2;
    t1 = std::chrono::steady_clock::now();

    // Threaded Cell Collection

    // Get Number of Threads
    size_t num_threads = Threading::getNumThreads();

    // Each thread keeps its hex_cells. When its a single thread, we store only one hex-cell at a time and add it to the mesh right away
    std::vector<std::vector<HexCell>> thread_hex_cells(num_threads);
    if (num_threads==1) {thread_hex_cells[0].resize(1);} // If using a single thread, we only store one hex-cell at a time
    std::vector<std::atomic<bool>> global_hex_corners_visited_in_parallel_cell_extraction(hexEx.numGlobalCorners());
    hex_cell_added.resize(hexEx.numGlobalCorners());
    hfh_per_halfface.resize(3 * hexEx.numGlobalCorners());

    size_t num_collected_cells = 0;
    t2 = std::chrono::steady_clock::now();
    hexEx.printLog("HexHex: Prepared for Cell Extraction in ", std::chrono::duration_cast<std::chrono::milliseconds> (t2 - t1).count(), " ms");

    t1 = std::chrono::steady_clock::now();

// Each Thread gets its vertices to process
#pragma omp parallel reduction(+:num_collected_cells)
    {
        // Get current thread id
        size_t t_id = Threading::getThreadID();

#pragma omp for nowait
        for (int i = 0; i < static_cast<int>(hexEx.getOutputMesh().n_vertices()); ++i)
        {
            auto n_local_corners = static_cast<int>(hexEx.getHexVertexGenerator(VertexHandle(i)).nLocalCorners());
            for (int j = 0; j < n_local_corners; ++j)
            {
                GCH gch{VertexHandle(i), LCH(j)};

                // Compute Unique corner index
                const uint32 corner_idx = hexEx.toGlobalIndex(gch);

                // Skip corner if already visited
                if (global_hex_corners_visited_in_parallel_cell_extraction[corner_idx].load(std::memory_order_acquire))
                {
                    continue;
                }

                // Collect Cell
                bool success = extractHexCell(gch, thread_hex_cells[t_id]);

                if (success)
                {
                    // Get Hex Cell
                    auto& hex_cell = thread_hex_cells[t_id].back();

                    // Is Hex Cell already processed? -> If not, claim ownership
                    bool f = false;
                    if (global_hex_corners_visited_in_parallel_cell_extraction[hex_cell.min_corner_idx]
                        .compare_exchange_strong(f, true, std::memory_order_acq_rel))
                    {
                        // Mark all corners as visited
                        for (uint32 idx : hex_cell.global_corner_indices)
                        {
                            if (idx != hex_cell.min_corner_idx)
                            {
                                global_hex_corners_visited_in_parallel_cell_extraction[idx]
                                    .store(true, std::memory_order_release);
                            }
                        }

                        // Add Hex Cell to the Mesh if using single thread (no danger)
                        if (num_threads == 1)
                        {
                            addHexCell(hex_cell);
                        }
                    }

                    num_collected_cells++;
                }
            }
        }
    }

    t2 = std::chrono::steady_clock::now();
    hexEx.printLog("HexHex: Successfully extracted ", num_collected_cells, " Cells in "
                   ,std::chrono::duration_cast<std::chrono::milliseconds> (t2 - t1).count()," ms ("
                   ,(num_collected_cells-hexEx.numGlobalCorners()/8), " were found multiple times)");


    // Add the collected cells and faces to the mesh if working multithreaded
    if (num_threads > 1)
    {
        t1 = std::chrono::steady_clock::now();
        uint32 num_cells_added = 0;
        uint32 num_cells_skipped = 0;
        for (size_t t_id = 0; t_id < num_threads; ++t_id)
        {
            for (size_t i = 0; i < thread_hex_cells[t_id].size(); ++i)
            {
                if (addHexCell(thread_hex_cells[t_id][i]))
                {
                    // we are happy c(^v^c)
                    num_cells_added++;
                }
                else
                {
                    // we aren't sad :I
                    num_cells_skipped++;
                }
            }
        }
        t2 = std::chrono::steady_clock::now();
        hexEx.printLog("HexHex: Combined per Thread Collection Results in ", std::chrono::duration_cast<std::chrono::milliseconds> (t2 - t1).count(), " ms");
        hexEx.printLog("HexHex: Added ", num_cells_added, " Cells, Skipped ", num_cells_skipped);
    }

    // If extracting a complete mesh and the IGM is valid, the number of corners should be 8 times the number of cells
    if ((hexEx.numGlobalCorners() != 8*hexEx.getOutputMesh().n_cells()))
    {
        hexEx.printErr("HexHex: The Number of Hex Corners implies ", hexEx.numGlobalCorners()/8, " Hex Cells but we have ", hexEx.getOutputMesh().n_cells());
    }

}

bool HexCellExtractor::extractHexCell(const GCH& gch, std::vector<HexCell>& hex_cells)
{
    HexCell hex_cell;
    hex_cell.min_corner_idx = UINT32_MAX;

    const VertexHandle vh0 = gch.vh;
    const LCH hch0 = gch.lch;
    auto& gen0 = hexEx.getHexVertexGenerator(vh0);

    const auto& hc0 = gen0.localCorners[hch0.idx()];

    // Propellers on vertex 0
    const LPH& ph0 = hc0[0];
    const LPH& ph1 = hc0[1];
    const LPH& ph2 = hc0[2];

    // Distance 1 connections
    const auto& conn01 = hexEx.getGlobalPropellerOpposite({vh0, ph0});
    const auto& conn03 = hexEx.getGlobalPropellerOpposite({vh0, ph1});
    const auto& conn04 = hexEx.getGlobalPropellerOpposite({vh0, ph2});
    const VertexHandle vh1 = conn01.gph.vh;
    const VertexHandle vh3 = conn03.gph.vh;
    const VertexHandle vh4 = conn04.gph.vh;
    if (!vh1.is_valid() || !vh3.is_valid() || !vh4.is_valid()) {
        {hexEx.printErr("HexHex: Skipping cell corner on invalid vertex.");}
        return false;
    }
    auto& gen1 = hexEx.getHexVertexGenerator(vh1);
    auto& gen3 = hexEx.getHexVertexGenerator(vh3);
    auto& gen4 = hexEx.getHexVertexGenerator(vh4);
    const LPH ph10 = conn01.gph.lph;
    const LPH ph12 = hexEx.getOppositeBlade(vh0, ph0, ph2).lph;
    const LPH ph31 = conn03.gph.lph;
    const LPH ph30 = hexEx.getOppositeBlade(vh0, ph1, ph0).lph;
    const LPH ph42 = conn04.gph.lph;
    const LPH ph41 = hexEx.getOppositeBlade(vh0, ph2, ph1).lph;
    const LCH hch1 = gen1.findLocalCorner(ph10, ph12);
    const LCH hch3 = gen3.findLocalCorner(ph31, ph30);
    const LCH hch4 = gen4.findLocalCorner(ph42, ph41);

    // Distance 2 connections
    const auto& conn17 = hexEx.getGlobalPropellerOpposite({vh1, ph12});
    const auto& conn32 = hexEx.getGlobalPropellerOpposite({vh3, ph30});
    const auto& conn45 = hexEx.getGlobalPropellerOpposite({vh4, ph41});
    const VertexHandle vh7 = conn17.gph.vh;
    const VertexHandle vh2 = conn32.gph.vh;
    const VertexHandle vh5 = conn45.gph.vh;
    if (!vh7.is_valid() || !vh2.is_valid() || !vh5.is_valid()) {
        {hexEx.printErr("HexHex: Skipping cell corner on invalid vertex.");}
        return false;
    }
    auto& gen7 = hexEx.getHexVertexGenerator(vh7);
    auto& gen2 = hexEx.getHexVertexGenerator(vh2);
    auto& gen5 = hexEx.getHexVertexGenerator(vh5);
    const LPH ph230 = conn32.gph.lph;
    const LPH ph231 = hexEx.getOppositeBlade(vh3, ph30, ph31).lph;
    const LPH ph541 = conn45.gph.lph;
    const LPH ph542 = hexEx.getOppositeBlade(vh4, ph41, ph42).lph;
    const LPH ph712 = conn17.gph.lph;
    const LPH ph710 = hexEx.getOppositeBlade(vh1, ph12, ph10).lph;
    const LCH hch2 = gen2.findLocalCorner(ph230, ph231);
    const LCH hch5 = gen5.findLocalCorner(ph541, ph542);
    const LCH hch7 = gen7.findLocalCorner(ph712, ph710);

    const LPH ph32 = hexEx.getOppositeBlade(vh0, ph1, ph2).lph;
    const LPH ph40 = hexEx.getGlobalPropellerOpposite({vh7, ph710}).gph.lph;
    const LPH ph11 = hexEx.getGlobalPropellerOpposite({vh2, ph231}).gph.lph;

    // Front HalfFace
    hex_cell.faces[4] = create_canonical_halfface_with_opposite(hexEx,
        vh0,hch0,ph2,
        vh4,hch4,ph40,
        vh7,hch7,ph712,
        vh1,hch1,ph10
    );

    // Left HalfFace
    hex_cell.faces[1] = create_canonical_halfface_with_opposite(hexEx,
        vh0,hch0,ph1,
        vh3,hch3,ph32,
        vh5,hch5,ph541,
        vh4,hch4,ph42
    );

    // Bottom HalfFace
    hex_cell.faces[2] = create_canonical_halfface_with_opposite(hexEx,
        vh0,hch0,ph0,
        vh1,hch1,ph11,
        vh2,hch2,ph230,
        vh3,hch3,ph31
    );

    // Store the corner indices of the cell
    hex_cell.global_corner_indices = {
        hexEx.toGlobalIndex({vh0,hch0}),
        hexEx.toGlobalIndex({vh1,hch1}),
        hexEx.toGlobalIndex({vh2,hch2}),
        hexEx.toGlobalIndex({vh3,hch3}),
        hexEx.toGlobalIndex({vh4,hch4}),
        hexEx.toGlobalIndex({vh5,hch5}),
        hexEx.toGlobalIndex({vh7,hch7})
    };

    // Compute the min. index as a unique identifier
    hex_cell.min_corner_idx = *std::min_element(hex_cell.global_corner_indices.begin(), hex_cell.global_corner_indices.end());

    // Distance 3 connection
    const LPH ph232 = hexEx.getOppositeBlade(vh3, ph30, ph32).lph;
    const auto& conn26 = hexEx.getGlobalPropellerOpposite({vh2, ph232});
    const VertexHandle vh6 = conn26.gph.vh;
    if (!vh6.is_valid()) {
        {hexEx.printErr("HexHex: Skipping cell corner on invalid vertex.");}
        return false;
    }
    auto& gen6 = hexEx.getHexVertexGenerator(vh6);
    const LPH ph6232 = conn26.gph.lph;
    const LPH ph6231 = hexEx.getOppositeBlade(vh2, ph232, ph231).lph;
    const LCH hch6 = gen6.findLocalCorner(ph6232, ph6231);

    // Get other propellers of hex cell
    const LPH ph711 = hexEx.getGlobalPropellerOpposite({vh6, ph6231}).gph.lph;
    const LPH ph6230 = hexEx.getOppositeBlade(vh2, ph232, ph230).lph;
    const LPH ph532 = hexEx.getGlobalPropellerOpposite({vh6, ph6230}).gph.lph;

    assert(ph40 == hexEx.getOppositeBlade(vh0, ph2, ph0).lph);
    assert(hexEx.getGlobalPropellerOpposite({vh5,ph532}).gph.lph == ph6230);

    // Top HalfFace
    hex_cell.faces[3] = create_canonical_halfface_with_opposite(hexEx,
        vh4,hch4,ph41,
        vh5,hch5,ph532,
        vh6,hch6,ph6231,
        vh7,hch7,ph710
    );

    // Right HalfFace
    hex_cell.faces[0] = create_canonical_halfface_with_opposite(hexEx,
        vh1,hch1,ph12,
        vh7,hch7,ph711,
        vh6,hch6,ph6232,
        vh2,hch2,ph231
    );

    // Back HalfFace
    hex_cell.faces[5] = create_canonical_halfface_with_opposite(hexEx,
        vh3,hch3,ph30,
        vh2,hch2,ph232,
        vh6,hch6,ph6230,
        vh5,hch5,ph542
    );

    // Insert the last corner index
    const uint32 idx6 = hexEx.toGlobalIndex({vh6,hch6});
    hex_cell.global_corner_indices.insert(hex_cell.global_corner_indices.begin()+6, idx6);
    assert(hex_cell.global_corner_indices[6] == idx6);

    // Compute the min. corner index
    hex_cell.min_corner_idx = std::min(hex_cell.min_corner_idx, idx6);

    // Add the hex-cell to the per thread list if using multiple threads
    if (Threading::getNumThreads() > 1)
    {
        hex_cells.emplace_back(hex_cell);
    }
    else
    {
    // If  using a single thread, we always override the first and only element
        hex_cells.clear();
        hex_cells.emplace_back(hex_cell);
    }

    return true;
}

bool HexCellExtractor::addHexCell(const HexCell& hex_cell)
{
    if (hex_cell_added[hex_cell.min_corner_idx]) return false; // cell already exists

    auto compute_unique_halfface_index = [&](const HexHalfFace& hf) -> uint32
    {
        const auto& q = hf.quarters[0]; // Since ordering is canonical, there is a 1:1 between hf quarter and hf
        const auto& gen = hexEx.getHexVertexGenerator(q.vh);
        const auto& lc = gen.getLocalCorner(q.lca.lch);
        const uint8_t axis = (lc[0]== q.lca.forward)? 0 : (lc[1]== q.lca.forward)? 1 : 2;
        assert(lc[axis] == q.lca.forward);
        return 3 * hexEx.toGlobalIndex(GCH{.vh = q.vh, .lch = q.lca.lch}) + axis;
    };

    auto get_hfh = [&](const HexHalfFace& hf) -> HalfFaceHandle
    {
        if (!hf.quarters[0].vh.is_valid()) return HalfFaceHandle(-1);
        return hfh_per_halfface[compute_unique_halfface_index(hf)];
    };

    auto set_hfh = [&](const HexHalfFace& hf, const HalfFaceHandle hfh) -> void
    {
        if (!hf.quarters[0].vh.is_valid()) return;
        hfh_per_halfface[compute_unique_halfface_index(hf)] = hfh;
    };

    // Collect the halfface handles
    std::vector<HalfFaceHandle> hfhs;
    hfhs.reserve(6);
    for (const auto& hex_face : hex_cell.faces)
    {
        // If Face is invalid, skip it (only happens for an incomplete cell which should never happen if the igm is valid)
        if (!hex_face.halfface.quarters[0].vh.is_valid())
        {
            continue;
        }

        // Assert canonical order
        assert(hex_face.halfface.quarters[0].vh < hex_face.halfface.quarters[1].vh);
        assert(hex_face.halfface.quarters[0].vh < hex_face.halfface.quarters[2].vh);
        assert(hex_face.halfface.quarters[0].vh < hex_face.halfface.quarters[3].vh);

        HalfFaceHandle hfh = get_hfh(hex_face.halfface);
        if (hfh.is_valid())
        {
            // Face already exists
            hfhs.push_back(hfh);
        }
        else
        {
            // Face does not yet exist, add it to the mesh

            // Collect the four halfedges
            const auto gphs = getGPHs(hex_face.halfface);
            std::vector<HalfEdgeHandle> hehs = {
                hexEx.getHexHalfEdge(gphs[0]),
                hexEx.getHexHalfEdge(gphs[1]),
                hexEx.getHexHalfEdge(gphs[2]),
                hexEx.getHexHalfEdge(gphs[3])
            };

            // Add the face
            const FaceHandle fh = hexEx.getOutputMesh().add_face(std::move(hehs));
            hfh = fh.halfface_handle(0);
            hfhs.push_back(hfh);

            // Store it in the map as well as the opposite
            set_hfh(hex_face.halfface, hfh);
            set_hfh(hex_face.opposite_halfface, hfh.opposite_handle());
        }
    }

    // Mark the cell as visited
    hex_cell_added[hex_cell.min_corner_idx] = true;

    // Add the cell if we have 6 halffaces
    if (hfhs.size() == 6u)
    {
        hexEx.getOutputMesh().add_cell(hfhs);

        return true;
    }
    return false; // incomplete cell
}

}
