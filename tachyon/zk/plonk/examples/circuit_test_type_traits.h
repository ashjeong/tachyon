#ifndef TACHYON_ZK_PLONK_EXAMPLES_CIRCUIT_TEST_TYPE_TRAITS_H_
#define TACHYON_ZK_PLONK_EXAMPLES_CIRCUIT_TEST_TYPE_TRAITS_H_

#include <type_traits>

#include "tachyon/math/elliptic_curves/bn/bn254/bn254.h"
#include "tachyon/zk/base/commitments/gwc_extension.h"
#include "tachyon/zk/base/commitments/shplonk_extension.h"
#include "tachyon/zk/lookup/halo2/scheme.h"
#include "tachyon/zk/plonk/layout/floor_planner/simple_floor_planner.h"
#include "tachyon/zk/plonk/layout/floor_planner/v1/v1_floor_planner.h"

namespace tachyon::zk::plonk::halo2 {

namespace {

using SHPlonk =
    SHPlonkExtension<math::bn254::BN254Curve, kMaxDegree, kMaxExtendedDegree,
                     math::bn254::G1AffinePoint>;

using GWC = GWCExtension<math::bn254::BN254Curve, kMaxDegree,
                         kMaxExtendedDegree, math::bn254::G1AffinePoint>;
using Halo2LS =
    lookup::halo2::Scheme<typename SHPlonk::Poly, typename SHPlonk::Evals,
                          typename SHPlonk::Commitment>;

}  // namespace

template <typename Circuit>
constexpr bool IsSimpleFloorPlanner =
    std::is_same_v<typename Circuit::FloorPlanner, SimpleFloorPlanner<Circuit>>;

template <typename Circuit>
constexpr bool IsV1FloorPlanner =
    std::is_same_v<typename Circuit::FloorPlanner, V1FloorPlanner<Circuit>>;

template <typename PCS>
constexpr bool IsSHPlonk = std::is_same_v<PCS, SHPlonk>;

template <typename PCS>
constexpr bool IsGWC = std::is_same_v<PCS, GWC>;

template <typename LS>
constexpr bool IsHalo2LS = std::is_same_v<LS, Halo2LS>;

}  // namespace tachyon::zk::plonk::halo2

#endif  // TACHYON_ZK_PLONK_EXAMPLES_CIRCUIT_TEST_TYPE_TRAITS_H_
