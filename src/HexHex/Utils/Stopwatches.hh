#pragma once

#include <libTimekeeper/StopWatch.hh>
#include <libTimekeeper/StopWatchPrinting.hh>

namespace HexHex {
    using Timekeeper::ScopedStopWatch;
}

namespace HexHex::sw {
    using HSW = Timekeeper::HierarchicalStopWatch;

    inline HSW root("hexhex");
        inline HSW readInputFiles("read_input_files", root);
        inline HSW writeOutputFiles("write_output_files", root);

        inline HSW extract("extract", root);
        inline HSW preprocessSanitization("preprocess_parametrization", extract);
            inline HSW cachePreSanitizeProperties("cache_pre_sanitize_properties", preprocessSanitization);
            inline HSW extractTransitionFunctions("extract_transitions", preprocessSanitization);
            inline HSW calculateSingularities("calculate_singularities", preprocessSanitization);
            inline HSW sanitizeParametrization("sanitize", preprocessSanitization);
            inline HSW snapBoundaryFaces("snap_boundary_faces", preprocessSanitization);
            inline HSW checkIGMValidity("check_igm_validity", preprocessSanitization);

            inline HSW extractHexVertices("extract_hex_vertices", extract);

            inline HSW extractLocalTopology("extract_local_topology", extract);

            inline HSW traceHexEdges("trace_hex_edges", extract);

            inline HSW extractHexFacesAndCells("extract_hex_faces_and_cells", extract);

            inline HSW transferFeatureTags("transfer_feature_tags", extract);

            inline HSW extractPiecewiseLinearFaces("extract_piecewise_linear_faces", extract);
}
