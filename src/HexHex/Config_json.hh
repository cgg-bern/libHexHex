#pragma once
#include <HexHex/Config.hh>
#include <nlohmann/json.hpp>

namespace HexHex {

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Config,
   extract_piecewise_linear_edges,
   extract_piecewise_linear_faces,
   verbose,
   hf_transition_prop_name,
   e_valence_prop_name,
   sanitize_igm,
   num_threads,
   igm_scaling_factor,
   sanitization_epsilon,
   assert_igm_validity,
   rasterization_epsilon);
} // namespace HexHex

