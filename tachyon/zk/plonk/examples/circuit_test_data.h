#ifndef TACHYON_ZK_PLONK_EXAMPLES_CIRCUIT_TEST_DATA_H_
#define TACHYON_ZK_PLONK_EXAMPLES_CIRCUIT_TEST_DATA_H_

#include <stdint.h>

#include <string_view>
#include <utility>
#include <vector>

#include "tachyon/base/range.h"
#include "tachyon/zk/plonk/examples/circuit_test.h"

namespace tachyon::zk::plonk::halo2 {

template <typename Circuit, typename PCS, typename LS>
class CircuitTestData {
 public:
  constexpr static size_t kN = 1;
  constexpr static std::string_view kPinnedConstraintSystem = "";
  constexpr static std::string_view kFixedColumnsOther[0][0] = {};
  constexpr static AnyColumnKey kColumns[0] = {};
  constexpr static Label kMapping[0][0] = {};
  constexpr static Label kAux[0][0] = {};
  constexpr static size_t kSizes[0][0] = {};
  constexpr static bool kSelectors[0][0] = {};
  constexpr static base::Range<RowIndex> kUsableRows =
      base::Range<RowIndex>::Until(10);
  constexpr static std::string_view kPinnedVerifyingKey = "";
  constexpr static std::string_view kTranscriptRepr = "";
  constexpr static std::string_view kLFirst[0] = {};
  constexpr static std::string_view kLLast[0] = {};
  constexpr static std::string_view kLActiveRow[0] = {};
  constexpr static std::string_view kFixedColumns[0][0] = {};
  constexpr static std::string_view kFixedPolys[0][0] = {};
  constexpr static std::string_view kPermutationsColumns[0][0] = {};
  constexpr static std::string_view kPermutationsPolys[0][0] = {};
  constexpr static uint8_t kProof[0] = {};
  constexpr static
      typename CircuitTest<PCS, LS>::Point kAdviceCommitments[0][0] = {};
  constexpr static std::string_view kChallenges[0] = {};
  constexpr static std::string_view kTheta = "";
  constexpr static typename CircuitTest<PCS, LS>::Point
      kPermutationProductCommitmentsInputPoints[0][0] = {};
  constexpr static typename CircuitTest<PCS, LS>::Point
      kPermutationProductCommitmentsTablePoints[0][0] = {};
  constexpr static std::string_view kBeta = "";
  constexpr static std::string_view kGamma = "";
  constexpr static typename CircuitTest<PCS, LS>::Point
      kPermutationProductCommitments[0][0] = {};
  constexpr static
      typename CircuitTest<PCS, LS>::Point kLookupProductCommitments[0][0] = {};
  constexpr static std::string_view kY = "";
  constexpr static
      typename CircuitTest<PCS, LS>::Point kVanishingHPolyCommitments[0] = {};
  constexpr static std::string_view kX = "";
  constexpr static std::string_view kAdviceEvals[0][0] = {};
  constexpr static std::string_view kFixedEvals[0] = {};
  constexpr static std::string_view kCommonPermutationEvals[0] = {};
  constexpr static std::string_view kPermutationProductEvals[0][0] = {};
  constexpr static std::string_view kPermutationProductNextEvals[0][0] = {};
  constexpr static std::string_view kPermutationProductLastEvals[0][0] = {};
  constexpr static std::string_view kLookupProductEvals[0][0] = {};
  constexpr static std::string_view kLookupProductNextEvals[0][0] = {};
  constexpr static std::string_view kLookupPermutedInputEvals[0][0] = {};
  constexpr static std::string_view kLookupPermutedInputPrevEvals[0][0] = {};
  constexpr static std::string_view kLookupPermutedTableEvals[0][0] = {};
  constexpr static std::string_view kHEval = "";

  constexpr static
      typename CircuitTest<PCS, LS>::Point kVanishingRandomPolyCommitment{
          "0x0000000000000000000000000000000000000000000000000000000000000001",
          "0x0000000000000000000000000000000000000000000000000000000000000002"};

  constexpr static std::string_view kVanishingRandomEval =
      "0x0000000000000000000000000000000000000000000000000000000000000001";

  static Circuit GetCircuit() { return Circuit(); }

  static std::vector<Circuit> GetCircuits() {
    Circuit circuit = GetCircuit();
    return {circuit, std::move(circuit)};
  }

  static std::vector<typename PCS::RationalEvals> GetFixedColumns() {
    return {};
  }

  static std::vector<typename PCS::Evals> GetInstanceColumns() { return {}; }
};

}  // namespace tachyon::zk::plonk::halo2

#endif  // TACHYON_ZK_PLONK_EXAMPLES_CIRCUIT_TEST_DATA_H_
