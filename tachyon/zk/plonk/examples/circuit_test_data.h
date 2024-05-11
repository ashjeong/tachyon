#ifndef TACHYON_ZK_PLONK_EXAMPLES_CIRCUIT_TEST_DATA_H_
#define TACHYON_ZK_PLONK_EXAMPLES_CIRCUIT_TEST_DATA_H_

#include <string_view>
#include <utility>
#include <vector>

#include "tachyon/base/range.h"
#include "tachyon/zk/plonk/examples/circuit_test.h"

namespace tachyon::zk::plonk::halo2 {

template <typename Circuit, typename PCS, typename LS>
class CircuitTestData {
 public:
  // constexpr static size_t kN; // REQUIRED
  // constexpr static bool kPinnedConstraintSystemFlag = false; // REQUIRED
  constexpr static bool kFixedColumnsOtherFlag = false;
  constexpr static bool kColumnsFlag = false;
  constexpr static bool kMappingFlag = false;
  constexpr static bool kAuxFlag = false;
  constexpr static bool kSizesFlag = false;
  // constexpr static bool kSelectorsFlag = false; // REQUIRED
  // constexpr static bool kUsableRowsFlag = false;
  // constexpr static bool kPinnedVerifyingKeyFlag = false; // REQUIRED
  // constexpr static bool kTranscriptReprFlag = false; // REQUIRED
  constexpr static bool kLFirstFlag = false;
  constexpr static bool kLLastFlag = false;
  constexpr static bool kLActiveRowFlag = false;
  constexpr static bool kFixedColumnsFlag = false;
  constexpr static bool kFixedPolysFlag = false;
  constexpr static bool kPermutationsColumnsFlag = false;
  constexpr static bool kPermutationsPolysFlag = false;
  // constexpr static bool kProofFlag = false; // REQUIRED
  constexpr static bool kAdviceCommitmentsFlag = false;
  constexpr static bool kChallengesFlag = false;
  // constexpr static bool kThetaFlag = false; // REQUIRED
  constexpr static bool kPermutationProductCommitmentsPointsFlag = false;
  // constexpr static bool kBetaFlag = false; // REQUIRED
  // constexpr static bool kGammaFlag = false; // REQUIRED
  constexpr static bool kPermutationProductCommitmentsFlag = false;
  constexpr static bool kLookupProductCommitmentsFlag = false;
  // constexpr static bool kYFlag = false; // REQUIRED
  constexpr static bool kVanishingHPolyCommitmentsFlag = false;
  // constexpr static bool kXFlag = false; // REQUIRED
  constexpr static bool kAdviceEvalsFlag = false;
  constexpr static bool kFixedEvalsFlag = false;
  constexpr static bool kCommonPermutationEvalsFlag = false;
  constexpr static bool kPermutationProductEvalsFlag = false;
  constexpr static bool kPermutationProductNextEvalsFlag = false;
  constexpr static bool kPermutationProductLastEvalsFlag = false;
  constexpr static bool kLookupProductEvalsFlag = false;
  constexpr static bool kLookupProductNextEvalsFlag = false;
  constexpr static bool kLookupPermutedInputEvalsFlag = false;
  constexpr static bool kLookupPermutedInputPrevEvalsFlag = false;
  constexpr static bool kLookupPermutedTableEvalsFlag = false;
  // constexpr static bool kHEvalFlag = false; // REQUIRED

  constexpr static base::Range<RowIndex> kUsableRows =
      base::Range<RowIndex>::Until(10);
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
