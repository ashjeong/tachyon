#include <utility>

#include "gtest/gtest.h"

#include "tachyon/base/array_to_vector.h"
#include "tachyon/math/elliptic_curves/bn/bn254/bn254.h"
#include "tachyon/zk/base/commitments/gwc_extension.h"
#include "tachyon/zk/base/commitments/shplonk_extension.h"
#include "tachyon/zk/lookup/halo2/scheme.h"
#include "tachyon/zk/plonk/examples/circuit_test.h"
#include "tachyon/zk/plonk/examples/circuit_test_type_traits.h"
#include "tachyon/zk/plonk/examples/shuffle_circuit.h"
#include "tachyon/zk/plonk/examples/shuffle_circuit_test_data.h"
#include "tachyon/zk/plonk/examples/simple_circuit.h"
#include "tachyon/zk/plonk/examples/simple_circuit_test_data.h"
#include "tachyon/zk/plonk/examples/simple_lookup_circuit.h"
#include "tachyon/zk/plonk/examples/simple_lookup_circuit_test_data.h"
#include "tachyon/zk/plonk/halo2/pinned_constraint_system.h"
#include "tachyon/zk/plonk/halo2/pinned_verifying_key.h"
#include "tachyon/zk/plonk/keys/proving_key.h"
#include "tachyon/zk/plonk/layout/floor_planner/simple_floor_planner.h"
#include "tachyon/zk/plonk/layout/floor_planner/v1/v1_floor_planner.h"

namespace tachyon::zk::plonk::halo2 {

namespace {

using SHPlonk =
    SHPlonkExtension<math::bn254::BN254Curve, kMaxDegree, kMaxExtendedDegree,
                     math::bn254::G1AffinePoint>;

using GWC = GWCExtension<math::bn254::BN254Curve, kMaxDegree,
                         kMaxExtendedDegree, math::bn254::G1AffinePoint>;

// Note(ashjeong): |Poly|, |Evals|, |Commitment|, and |Field| are currently the
// same across all types
using Halo2LS =
    lookup::halo2::Scheme<typename SHPlonk::Poly, typename SHPlonk::Evals,
                          typename SHPlonk::Commitment>;

const size_t kBits = 3;
const size_t kW = 2;
const size_t kH = 8;

template <typename _Circuit, typename _PCS, typename _LS>
struct TestingArguments {
  using Circuit = _Circuit;
  using PCS = _PCS;
  using LS = _LS;
};

template <typename Circuit, typename PCS, typename LS>
struct ExampleTestDataSelector {
  using F = typename PCS::Field;

  // clang-format off
  using Type =
      std::conditional_t<
        std::is_same_v<Circuit, SimpleCircuit<F, SimpleFloorPlanner>>,
          SimpleTestData<SimpleCircuit<F, SimpleFloorPlanner>, PCS, LS>,
      std::conditional_t<
        std::is_same_v<Circuit, SimpleCircuit<F, V1FloorPlanner>>,
          SimpleTestData<SimpleCircuit<F, V1FloorPlanner>, PCS, LS>,
      std::conditional_t<
        std::is_same_v<Circuit, SimpleLookupCircuit<F, kBits, SimpleFloorPlanner>>,
          SimpleLookupTestData<SimpleLookupCircuit<F, kBits, SimpleFloorPlanner>, PCS, LS>,
      std::conditional_t<
        std::is_same_v<Circuit, SimpleLookupCircuit<F, kBits, V1FloorPlanner>>,
          SimpleLookupTestData<SimpleLookupCircuit<F, kBits, V1FloorPlanner>, PCS, LS>,
      std::conditional_t<
        std::is_same_v<Circuit, ShuffleCircuit<F, kW, kH, SimpleFloorPlanner>>,
          ShuffleTestData<ShuffleCircuit<F, kW, kH, SimpleFloorPlanner>, PCS, LS>,
          ShuffleTestData<ShuffleCircuit<F, kW, kH, V1FloorPlanner>, PCS, LS>
      >>>>>;
  // clang-format on
};
template <typename TestingType>
class ExampleCircuitTest
    : public CircuitTest<typename TestingType::PCS, typename TestingType::LS> {
 public:
  static void SetUpTestSuite() { math::bn254::BN254Curve::Init(); }
};

using TestingTypes = testing::Types<
    TestingArguments<SimpleCircuit<SHPlonk::Field, SimpleFloorPlanner>, SHPlonk,
                     Halo2LS>,
    TestingArguments<SimpleCircuit<SHPlonk::Field, V1FloorPlanner>, SHPlonk,
                     Halo2LS>,
    TestingArguments<
        SimpleLookupCircuit<SHPlonk::Field, kBits, SimpleFloorPlanner>, SHPlonk,
        Halo2LS>,
    TestingArguments<SimpleLookupCircuit<SHPlonk::Field, kBits, V1FloorPlanner>,
                     SHPlonk, Halo2LS>,
    TestingArguments<ShuffleCircuit<SHPlonk::Field, kW, kH, SimpleFloorPlanner>,
                     SHPlonk, Halo2LS>,
    TestingArguments<ShuffleCircuit<SHPlonk::Field, kW, kH, V1FloorPlanner>,
                     SHPlonk, Halo2LS>,
    TestingArguments<ShuffleCircuit<GWC::Field, kW, kH, SimpleFloorPlanner>,
                     GWC, Halo2LS>,
    TestingArguments<ShuffleCircuit<GWC::Field, kW, kH, V1FloorPlanner>, GWC,
                     Halo2LS>>;
}  // namespace

TYPED_TEST_SUITE(ExampleCircuitTest, TestingTypes);

TYPED_TEST(ExampleCircuitTest, Configure) {
  using TestingType = TypeParam;
  using Circuit = typename TestingType::Circuit;
  using PCS = typename TestingType::PCS;
  using LS = typename TestingType::LS;
  using F = typename PCS::Field;
  using ExampleTestData =
      typename ExampleTestDataSelector<Circuit, PCS, LS>::Type;

  ConstraintSystem<F> constraint_system;
  auto config = Circuit::Configure(constraint_system);

  ExampleTestData::TestConfig(config);

  halo2::PinnedConstraintSystem<F> pinned_constraint_system(constraint_system);
  EXPECT_EQ(ExampleTestData::kPinnedConstraintSystem,
            base::ToRustDebugString(pinned_constraint_system));

  EXPECT_TRUE(constraint_system.selector_map().empty());
  EXPECT_TRUE(constraint_system.general_column_annotations().empty());
}

TYPED_TEST(ExampleCircuitTest, Synthesize) {
  using TestingType = TypeParam;
  using Circuit = typename TestingType::Circuit;
  using PCS = typename TestingType::PCS;
  using LS = typename TestingType::LS;
  using F = typename PCS::Field;
  using Domain = typename PCS::Domain;
  using Commitment = typename PCS::Commitment;
  using RationalEvals = typename PCS::RationalEvals;
  using ExampleTestData =
      typename ExampleTestDataSelector<Circuit, PCS, LS>::Type;

  CHECK(this->prover_->pcs().UnsafeSetup(ExampleTestData::kN, F(2)));
  this->prover_->set_domain(Domain::Create(ExampleTestData::kN));
  const Domain* domain = this->prover_->domain();

  ConstraintSystem<F> constraint_system;
  auto config = Circuit::Configure(constraint_system);
  Assembly<RationalEvals> assembly =
      VerifyingKey<F, Commitment>::template CreateAssembly<RationalEvals>(
          domain, constraint_system);

  Circuit circuit = ExampleTestData::GetCircuit();
  typename Circuit::FloorPlanner floor_planner;
  floor_planner.Synthesize(&assembly, circuit, std::move(config),
                           constraint_system.constants());

  {
    std::vector<RationalEvals> expected_fixed_columns =
        this->CreateRationalColumns(
            base::Array2DToVector2D(ExampleTestData::kFixedColumnsOther));
    EXPECT_EQ(assembly.fixed_columns(), expected_fixed_columns);
  }
  {
    std::vector<AnyColumnKey> expected_columns =
        base::ArrayToVector(ExampleTestData::kColumns);
    EXPECT_EQ(assembly.permutation().columns(), expected_columns);
  }

  const CycleStore& cycle_store = assembly.permutation().cycle_store();

  {
    CycleStore::Table<Label> expected_mapping(
        base::Array2DToVector2D(ExampleTestData::kMapping));
    EXPECT_EQ(cycle_store.mapping(), expected_mapping);
  }
  {
    CycleStore::Table<Label> expected_aux(
        base::Array2DToVector2D(ExampleTestData::kAux));
    EXPECT_EQ(cycle_store.aux(), expected_aux);
  }
  {
    CycleStore::Table<size_t> expected_sizes(
        base::Array2DToVector2D(ExampleTestData::kSizes));
    EXPECT_EQ(cycle_store.sizes(), expected_sizes);
  }
  {
    std::vector<std::vector<bool>> expected_selectors =
        base::Array2DToVector2D(ExampleTestData::kSelectors);
    EXPECT_EQ(assembly.selectors(), expected_selectors);
  }
  EXPECT_EQ(assembly.usable_rows(), ExampleTestData::kUsableRows);
}

TYPED_TEST(ExampleCircuitTest, LoadVerifyingKey) {
  using TestingType = TypeParam;
  using Circuit = typename TestingType::Circuit;
  using PCS = typename TestingType::PCS;
  using LS = typename TestingType::LS;
  using F = typename PCS::Field;
  using Domain = typename PCS::Domain;
  using Commitment = typename PCS::Commitment;
  using ExampleTestData =
      typename ExampleTestDataSelector<Circuit, PCS, LS>::Type;

  // TODO(ashjeong): Implement test for |ShuffleCircuit| GWC version
  if constexpr ((std::is_same_v<
                     Circuit, ShuffleCircuit<F, kW, kH, SimpleFloorPlanner>> ||
                 std::is_same_v<
                     Circuit,
                     ShuffleCircuit<F, kW, kH, V1FloorPlanner>>)&&IsGWC<PCS>) {
    GTEST_SKIP() << "LoadVerifyingKey test skipped for ShuffleCircuit with GWC";
    return;
  }

  CHECK(this->prover_->pcs().UnsafeSetup(ExampleTestData::kN, F(2)));
  this->prover_->set_domain(Domain::Create(ExampleTestData::kN));

  Circuit circuit = ExampleTestData::GetCircuit();

  VerifyingKey<F, Commitment> vkey;
  ASSERT_TRUE(vkey.Load(this->prover_.get(), circuit));

  {
    halo2::PinnedVerifyingKey pinned_vkey(this->prover_.get(), vkey);
    EXPECT_EQ(base::ToRustDebugString(pinned_vkey),
              ExampleTestData::kPinnedVerifyingKey);
  }
  {
    F expected_transcript_repr =
        *F::FromHexString(ExampleTestData::kTranscriptRepr);
    EXPECT_EQ(vkey.transcript_repr(), expected_transcript_repr);
  }
}

TYPED_TEST(ExampleCircuitTest, LoadProvingKey) {
  using TestingType = TypeParam;
  using Circuit = typename TestingType::Circuit;
  using PCS = typename TestingType::PCS;
  using LS = typename TestingType::LS;
  using F = typename PCS::Field;
  using Poly = typename PCS::Poly;
  using Evals = typename PCS::Evals;
  using Domain = typename PCS::Domain;
  using Commitment = typename PCS::Commitment;
  using ExampleTestData =
      typename ExampleTestDataSelector<Circuit, PCS, LS>::Type;

  // TODO(ashjeong): Implement test for |ShuffleCircuit| GWC version
  if constexpr ((std::is_same_v<
                     Circuit, ShuffleCircuit<F, kW, kH, SimpleFloorPlanner>> ||
                 std::is_same_v<
                     Circuit,
                     ShuffleCircuit<F, kW, kH, V1FloorPlanner>>)&&IsGWC<PCS>) {
    GTEST_SKIP() << "LoadProvingKey test skipped for ShuffleCircuit with GWC";
    return;
  }

  CHECK(this->prover_->pcs().UnsafeSetup(ExampleTestData::kN, F(2)));
  this->prover_->set_domain(Domain::Create(ExampleTestData::kN));

  Circuit circuit = ExampleTestData::GetCircuit();

  for (size_t i = 0; i < 2; ++i) {
    ProvingKey<LS> pkey;
    bool load_verifying_key = i == 0;
    SCOPED_TRACE(
        absl::Substitute("load_verifying_key: $0", load_verifying_key));
    if (load_verifying_key) {
      VerifyingKey<F, Commitment> vkey;
      ASSERT_TRUE(vkey.Load(this->prover_.get(), circuit));
      ASSERT_TRUE(pkey.LoadWithVerifyingKey(this->prover_.get(), circuit,
                                            std::move(vkey)));
    } else {
      ASSERT_TRUE(pkey.Load(this->prover_.get(), circuit));
    }

    {
      Poly expected_l_first =
          this->CreatePoly(base::ArrayToVector(ExampleTestData::kLFirst));
      EXPECT_EQ(pkey.l_first(), expected_l_first);
    }
    {
      Poly expected_l_last =
          this->CreatePoly(base::ArrayToVector(ExampleTestData::kLLast));
      EXPECT_EQ(pkey.l_last(), expected_l_last);
    }
    {
      Poly expected_l_active_row =
          this->CreatePoly(base::ArrayToVector(ExampleTestData::kLActiveRow));
      EXPECT_EQ(pkey.l_active_row(), expected_l_active_row);
    }
    {
      std::vector<Evals> expected_fixed_columns = this->CreateColumns(
          base::Array2DToVector2D(ExampleTestData::kFixedColumns));
      EXPECT_EQ(pkey.fixed_columns(), expected_fixed_columns);
    }
    {
      std::vector<Poly> expected_fixed_polys = this->CreatePolys(
          base::Array2DToVector2D(ExampleTestData::kFixedPolys));
      EXPECT_EQ(pkey.fixed_polys(), expected_fixed_polys);
    }
    {
      std::vector<Evals> expected_permutations_columns = this->CreateColumns(
          base::Array2DToVector2D(ExampleTestData::kPermutationsColumns));
      EXPECT_EQ(pkey.permutation_proving_key().permutations(),
                expected_permutations_columns);
    }
    {
      std::vector<Poly> expected_fixed_polys = this->CreatePolys(
          base::Array2DToVector2D(ExampleTestData::kPermutationsPolys));
      EXPECT_EQ(pkey.permutation_proving_key().polys(), expected_fixed_polys);
    }
  }
}

TYPED_TEST(ExampleCircuitTest, CreateProof) {
  using TestingType = TypeParam;
  using Circuit = typename TestingType::Circuit;
  using PCS = typename TestingType::PCS;
  using LS = typename TestingType::LS;
  using F = typename PCS::Field;
  using Evals = typename PCS::Evals;
  using Domain = typename PCS::Domain;
  using ExampleTestData =
      typename ExampleTestDataSelector<Circuit, PCS, LS>::Type;

  CHECK(this->prover_->pcs().UnsafeSetup(ExampleTestData::kN, F(2)));
  this->prover_->set_domain(Domain::Create(ExampleTestData::kN));

  std::vector<Circuit> circuits = ExampleTestData::GetCircuits();

  std::vector<Evals> instance_columns = ExampleTestData::GetInstanceColumns();
  std::vector<std::vector<Evals>> instance_columns_vec = {
      instance_columns, std::move(instance_columns)};

  ProvingKey<LS> pkey;
  ASSERT_TRUE(pkey.Load(this->prover_.get(), circuits[0]));
  this->prover_->CreateProof(pkey, std::move(instance_columns_vec), circuits);

  {
    std::vector<uint8_t> proof =
        this->prover_->GetWriter()->buffer().owned_buffer();
    EXPECT_EQ(proof, base::ArrayToVector(ExampleTestData::kProof));
  }
}

TYPED_TEST(ExampleCircuitTest, VerifyProof) {
  using TestingType = TypeParam;
  using Circuit = typename TestingType::Circuit;
  using PCS = typename TestingType::PCS;
  using LS = typename TestingType::LS;
  using F = typename PCS::Field;
  using Evals = typename PCS::Evals;
  using Domain = typename PCS::Domain;
  using Commitment = typename PCS::Commitment;
  using ExampleTestData =
      typename ExampleTestDataSelector<Circuit, PCS, LS>::Type;

  CHECK(this->prover_->pcs().UnsafeSetup(ExampleTestData::kN, F(2)));
  this->prover_->set_domain(Domain::Create(ExampleTestData::kN));

  Circuit circuit = ExampleTestData::GetCircuit();

  VerifyingKey<F, Commitment> vkey;
  ASSERT_TRUE(vkey.Load(this->prover_.get(), circuit));

  std::vector<uint8_t> owned_proof =
      base::ArrayToVector(ExampleTestData::kProof);

  Verifier<PCS, LS> verifier = this->CreateVerifier(
      this->CreateBufferWithProof(absl::MakeSpan(owned_proof)));

  std::vector<Evals> instance_columns = ExampleTestData::GetInstanceColumns();
  std::vector<std::vector<Evals>> instance_columns_vec = {
      instance_columns, std::move(instance_columns)};

  Proof<F, Commitment> proof;
  F h_eval;
  ASSERT_TRUE(verifier.VerifyProofForTesting(vkey, instance_columns_vec, &proof,
                                             &h_eval));

  // TODO(ashjeong): Implement rest of test for |ShuffleCircuit| GWC version
  if constexpr ((std::is_same_v<
                     Circuit, ShuffleCircuit<F, kW, kH, SimpleFloorPlanner>> ||
                 std::is_same_v<
                     Circuit,
                     ShuffleCircuit<F, kW, kH, V1FloorPlanner>>)&&IsGWC<PCS>) {
    GTEST_SKIP()
        << "Remaining VerifyProof testing skipped for ShuffleCircuit with GWC";
    return;
  }

  {
    std::vector<std::vector<Commitment>> expected_advice_commitments_vec{
        this->CreateCommitments(
            base::ArrayToVector(ExampleTestData::kAdviceCommitments[0])),
        this->CreateCommitments(
            base::ArrayToVector(ExampleTestData::kAdviceCommitments[1])),
    };
    EXPECT_EQ(proof.advices_commitments_vec, expected_advice_commitments_vec);
  }
  {
    std::vector<F> expected_challenges =
        this->CreateEvals(base::ArrayToVector(ExampleTestData::kChallenges));
    EXPECT_EQ(proof.challenges, expected_challenges);
  }
  {
    F expected_theta = *F::FromHexString(ExampleTestData::kTheta);
    EXPECT_EQ(proof.theta, expected_theta);
  }
  {
    std::vector<std::vector<lookup::Pair<Commitment>>>
        expected_lookup_permuted_commitments_vec{
            this->CreateLookupPermutedCommitments(
                base::ArrayToVector(
                    ExampleTestData::kPermutationProductCommitmentsInputPoints
                        [0]),
                base::ArrayToVector(
                    ExampleTestData::kPermutationProductCommitmentsTablePoints
                        [0])),
            this->CreateLookupPermutedCommitments(
                base::ArrayToVector(
                    ExampleTestData::kPermutationProductCommitmentsInputPoints
                        [1]),
                base::ArrayToVector(
                    ExampleTestData::kPermutationProductCommitmentsTablePoints
                        [1])),
        };
    EXPECT_EQ(proof.lookup_permuted_commitments_vec,
              expected_lookup_permuted_commitments_vec);
  }
  {
    F expected_beta = *F::FromHexString(ExampleTestData::kBeta);
    EXPECT_EQ(proof.beta, expected_beta);
  }
  {
    F expected_gamma = *F::FromHexString(ExampleTestData::kGamma);
    EXPECT_EQ(proof.gamma, expected_gamma);
  }
  {
    std::vector<std::vector<Commitment>>
        expected_permutation_product_commitments_vec{
            this->CreateCommitments(base::ArrayToVector(
                ExampleTestData::kPermutationProductCommitments[0])),
            this->CreateCommitments(base::ArrayToVector(
                ExampleTestData::kPermutationProductCommitments[1])),
        };
    EXPECT_EQ(proof.permutation_product_commitments_vec,
              expected_permutation_product_commitments_vec);
  }
  {
    std::vector<std::vector<Commitment>>
        expected_lookup_product_commitments_vec{
            this->CreateCommitments(base::ArrayToVector(
                ExampleTestData::kLookupProductCommitments[0])),
            this->CreateCommitments(base::ArrayToVector(
                ExampleTestData::kLookupProductCommitments[1])),
        };
    EXPECT_EQ(proof.lookup_product_commitments_vec,
              expected_lookup_product_commitments_vec);
  }
  {
    Commitment expected_vanishing_random_poly_commitment =
        this->CreateCommitment(ExampleTestData::kVanishingRandomPolyCommitment);
    EXPECT_EQ(proof.vanishing_random_poly_commitment,
              expected_vanishing_random_poly_commitment);
  }
  {
    F expected_y = *F::FromHexString(ExampleTestData::kY);
    EXPECT_EQ(proof.y, expected_y);
  }
  {
    std::vector<Commitment> expected_vanishing_h_poly_commitments =
        this->CreateCommitments(
            base::ArrayToVector(ExampleTestData::kVanishingHPolyCommitments));
    EXPECT_EQ(proof.vanishing_h_poly_commitments,
              expected_vanishing_h_poly_commitments);
  }
  {
    F expected_x = *F::FromHexString(ExampleTestData::kX);
    EXPECT_EQ(proof.x, expected_x);
  }
  {
    std::vector<std::vector<F>> expected_advice_evals_vec{
        this->CreateEvals(
            base::ArrayToVector(ExampleTestData::kAdviceEvals[0])),
        this->CreateEvals(
            base::ArrayToVector(ExampleTestData::kAdviceEvals[1])),
    };
    EXPECT_EQ(proof.advice_evals_vec, expected_advice_evals_vec);
  }
  {
    std::vector<F> expected_fixed_evals =
        this->CreateEvals(base::ArrayToVector(ExampleTestData::kFixedEvals));
    EXPECT_EQ(proof.fixed_evals, expected_fixed_evals);
  }
  {
    F expected_vanishing_random_eval =
        *F::FromHexString(ExampleTestData::kVanishingRandomEval);
    EXPECT_EQ(proof.vanishing_random_eval, expected_vanishing_random_eval);
  }
  {
    std::vector<F> expected_common_permutation_evals = this->CreateEvals(
        base::ArrayToVector(ExampleTestData::kCommonPermutationEvals));
    EXPECT_EQ(proof.common_permutation_evals,
              expected_common_permutation_evals);
  }
  {
    std::vector<std::vector<F>> expected_permutation_product_evals_vec{
        this->CreateEvals(
            base::ArrayToVector(ExampleTestData::kPermutationProductEvals[0])),
        this->CreateEvals(
            base::ArrayToVector(ExampleTestData::kPermutationProductEvals[1])),
    };
    EXPECT_EQ(proof.permutation_product_evals_vec,
              expected_permutation_product_evals_vec);
  }
  {
    std::vector<std::vector<F>> expected_permutation_product_next_evals_vec{
        this->CreateEvals(base::ArrayToVector(
            ExampleTestData::kPermutationProductNextEvals[0])),
        this->CreateEvals(base::ArrayToVector(
            ExampleTestData::kPermutationProductNextEvals[1])),
    };
    EXPECT_EQ(proof.permutation_product_next_evals_vec,
              expected_permutation_product_next_evals_vec);
  }
  {
    std::vector<std::vector<std::optional<F>>>
        expected_permutation_product_last_evals_vec{
            this->CreateOptionalEvals(base::ArrayToVector(
                ExampleTestData::kPermutationProductLastEvals[0])),
            this->CreateOptionalEvals(base::ArrayToVector(
                ExampleTestData::kPermutationProductLastEvals[1])),
        };
    EXPECT_EQ(proof.permutation_product_last_evals_vec,
              expected_permutation_product_last_evals_vec);
  }
  {
    std::vector<std::vector<F>> expected_lookup_product_evals_vec{
        this->CreateEvals(
            base::ArrayToVector(ExampleTestData::kLookupProductEvals[0])),
        this->CreateEvals(
            base::ArrayToVector(ExampleTestData::kLookupProductEvals[1])),
    };
    EXPECT_EQ(proof.lookup_product_evals_vec,
              expected_lookup_product_evals_vec);
  }
  {
    std::vector<std::vector<F>> expected_lookup_product_next_evals_vec{
        this->CreateEvals(
            base::ArrayToVector(ExampleTestData::kLookupProductNextEvals[0])),
        this->CreateEvals(
            base::ArrayToVector(ExampleTestData::kLookupProductNextEvals[1])),
    };
    EXPECT_EQ(proof.lookup_product_next_evals_vec,
              expected_lookup_product_next_evals_vec);
  }
  {
    std::vector<std::vector<F>> expected_lookup_permuted_input_evals_vec{
        this->CreateEvals(
            base::ArrayToVector(ExampleTestData::kLookupPermutedInputEvals[0])),
        this->CreateEvals(
            base::ArrayToVector(ExampleTestData::kLookupPermutedInputEvals[1])),
    };
    EXPECT_EQ(proof.lookup_permuted_input_evals_vec,
              expected_lookup_permuted_input_evals_vec);
  }
  {
    std::vector<std::vector<F>> expected_lookup_permuted_input_prev_evals_vec{
        this->CreateEvals(base::ArrayToVector(
            ExampleTestData::kLookupPermutedInputPrevEvals[0])),
        this->CreateEvals(base::ArrayToVector(
            ExampleTestData::kLookupPermutedInputPrevEvals[1])),
    };
    EXPECT_EQ(proof.lookup_permuted_input_prev_evals_vec,
              expected_lookup_permuted_input_prev_evals_vec);
  }
  {
    std::vector<std::vector<F>> expected_lookup_permuted_table_evals_vec{
        this->CreateEvals(
            base::ArrayToVector(ExampleTestData::kLookupPermutedTableEvals[0])),
        this->CreateEvals(
            base::ArrayToVector(ExampleTestData::kLookupPermutedTableEvals[1])),
    };
    EXPECT_EQ(proof.lookup_permuted_table_evals_vec,
              expected_lookup_permuted_table_evals_vec);
  }
  {
    F expected_h_eval = *F::FromHexString(ExampleTestData::kHEval);
    EXPECT_EQ(h_eval, expected_h_eval);
  }
}

}  // namespace tachyon::zk::plonk::halo2
