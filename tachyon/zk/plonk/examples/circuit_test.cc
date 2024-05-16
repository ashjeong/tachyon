#include "tachyon/zk/plonk/examples/circuit_test.h"

#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "tachyon/base/array_to_vector.h"
#include "tachyon/zk/plonk/examples/simple_circuit.h"
#include "tachyon/zk/plonk/examples/simple_circuit_test_data.h"
#include "tachyon/zk/plonk/examples/simple_lookup_circuit.h"
#include "tachyon/zk/plonk/examples/simple_lookup_circuit_test_data.h"
#include "tachyon/zk/plonk/halo2/pinned_constraint_system.h"
#include "tachyon/zk/plonk/halo2/pinned_verifying_key.h"
#include "tachyon/zk/plonk/keys/proving_key.h"
#include "tachyon/zk/plonk/layout/floor_planner/simple_floor_planner.h"
#include "tachyon/zk/plonk/layout/floor_planner/v1/v1_floor_planner.h"

namespace tachyon::zk::plonk {

template <typename TestArguments, typename TestData>
void CircuitTest<TestArguments, TestData>::ConfigureTest() {
  ConstraintSystem<F> constraint_system;
  auto config = Circuit::Configure(constraint_system);

  TestData::TestConfig(config);

  halo2::PinnedConstraintSystem<F> pinned_constraint_system(constraint_system);
  EXPECT_EQ(TestData::kPinnedConstraintSystem,
            base::ToRustDebugString(pinned_constraint_system));

  EXPECT_TRUE(constraint_system.selector_map().empty());
  EXPECT_TRUE(constraint_system.general_column_annotations().empty());
}

template <typename TestArguments, typename TestData>
void CircuitTest<TestArguments, TestData>::SynthesizeTest() {
  CHECK(this->prover_->pcs().UnsafeSetup(TestData::kN, F(2)));
  this->prover_->set_domain(Domain::Create(TestData::kN));
  const Domain* domain = this->prover_->domain();

  ConstraintSystem<F> constraint_system;
  auto config = Circuit::Configure(constraint_system);
  Assembly<RationalEvals> assembly =
      VerifyingKey<F, Commitment>::template CreateAssembly<RationalEvals>(
          domain, constraint_system);

  Circuit circuit = TestData::GetCircuit();
  typename Circuit::FloorPlanner floor_planner;
  floor_planner.Synthesize(&assembly, circuit, std::move(config),
                           constraint_system.constants());

  if constexpr (TestData::kAssemblyFixedColumnsFlag) {
    std::vector<RationalEvals> expected_fixed_columns =
        this->CreateRationalColumns(
            base::Array2DToVector2D(TestData::kAssemblyFixedColumns));
    EXPECT_EQ(assembly.fixed_columns(), expected_fixed_columns);
  } else {
    EXPECT_TRUE(assembly.fixed_columns().empty());
  }

  if constexpr (TestData::kAssemblyPermutationColumnsFlag) {
    std::vector<AnyColumnKey> expected_columns =
        base::ArrayToVector(TestData::kAssemblyPermutationColumns);
    EXPECT_EQ(assembly.permutation().columns(), expected_columns);
  } else {
    EXPECT_TRUE(assembly.permutation().columns().empty());
  }

  const CycleStore& cycle_store = assembly.permutation().cycle_store();

  if constexpr (TestData::kCycleStoreMappingFlag) {
    CycleStore::Table<Label> expected_mapping(
        base::Array2DToVector2D(TestData::kCycleStoreMapping));
    EXPECT_EQ(cycle_store.mapping(), expected_mapping);
  } else {
    EXPECT_TRUE(cycle_store.mapping().IsEmpty());
  }

  if constexpr (TestData::kCycleStoreAuxFlag) {
    CycleStore::Table<Label> expected_aux(
        base::Array2DToVector2D(TestData::kCycleStoreAux));
    EXPECT_EQ(cycle_store.aux(), expected_aux);
  } else {
    EXPECT_TRUE(cycle_store.aux().IsEmpty());
  }

  if constexpr (TestData::kCycleStoreSizesFlag) {
    CycleStore::Table<size_t> expected_sizes(
        base::Array2DToVector2D(TestData::kCycleStoreSizes));
    EXPECT_EQ(cycle_store.sizes(), expected_sizes);
  } else {
    EXPECT_TRUE(cycle_store.sizes().IsEmpty());
  }

  std::vector<std::vector<bool>> expected_selectors =
      base::Array2DToVector2D(TestData::kCycleStoreSelectors);
  EXPECT_EQ(assembly.selectors(), expected_selectors);

  EXPECT_EQ(assembly.usable_rows(), TestData::kUsableRows);
}

template <typename TestArguments, typename TestData>
void CircuitTest<TestArguments, TestData>::LoadVerifyingKeyTest() {
  CHECK(this->prover_->pcs().UnsafeSetup(TestData::kN, F(2)));
  this->prover_->set_domain(Domain::Create(TestData::kN));

  Circuit circuit = TestData::GetCircuit();

  VerifyingKey<F, Commitment> vkey;
  ASSERT_TRUE(vkey.Load(this->prover_.get(), circuit));

  halo2::PinnedVerifyingKey pinned_vkey(this->prover_.get(), vkey);
  EXPECT_EQ(base::ToRustDebugString(pinned_vkey),
            TestData::kPinnedVerifyingKey);

  F expected_transcript_repr = *F::FromHexString(TestData::kTranscriptRepr);
  EXPECT_EQ(vkey.transcript_repr(), expected_transcript_repr);
}

template <typename TestArguments, typename TestData>
void CircuitTest<TestArguments, TestData>::LoadProvingKeyTest() {
  CHECK(this->prover_->pcs().UnsafeSetup(TestData::kN, F(2)));
  this->prover_->set_domain(Domain::Create(TestData::kN));

  Circuit circuit = TestData::GetCircuit();

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

    if constexpr (TestData::kLFirstFlag) {
      Poly expected_l_first =
          this->CreatePoly(base::ArrayToVector(TestData::kLFirst));
      EXPECT_EQ(pkey.l_first(), expected_l_first);
    } else {
      EXPECT_TRUE(pkey.l_first().empty());
    }

    if constexpr (TestData::kLLastFlag) {
      Poly expected_l_last =
          this->CreatePoly(base::ArrayToVector(TestData::kLLast));
      EXPECT_EQ(pkey.l_last(), expected_l_last);
    } else {
      EXPECT_TRUE(pkey.l_last().empty());
    }

    if constexpr (TestData::kLActiveRowFlag) {
      Poly expected_l_active_row =
          this->CreatePoly(base::ArrayToVector(TestData::kLActiveRow));
      EXPECT_EQ(pkey.l_active_row(), expected_l_active_row);
    } else {
      EXPECT_TRUE(pkey.l_active_row().empty());
    }

    if constexpr (TestData::kFixedColumnsFlag) {
      std::vector<Evals> expected_fixed_columns =
          this->CreateColumns(base::Array2DToVector2D(TestData::kFixedColumns));
      EXPECT_EQ(pkey.fixed_columns(), expected_fixed_columns);
    } else {
      EXPECT_TRUE(pkey.fixed_columns().empty());
    }

    if constexpr (TestData::kFixedPolysFlag) {
      std::vector<Poly> expected_fixed_polys =
          this->CreatePolys(base::Array2DToVector2D(TestData::kFixedPolys));
      EXPECT_EQ(pkey.fixed_polys(), expected_fixed_polys);
    } else {
      EXPECT_TRUE(pkey.fixed_polys().empty());
    }

    if constexpr (TestData::kPermutationsColumnsFlag) {
      std::vector<Evals> expected_permutations_columns = this->CreateColumns(
          base::Array2DToVector2D(TestData::kPermutationsColumns));
      EXPECT_EQ(pkey.permutation_proving_key().permutations(),
                expected_permutations_columns);
    } else {
      EXPECT_TRUE(pkey.permutation_proving_key().permutations().empty());
    }

    if constexpr (TestData::kPermutationsPolysFlag) {
      std::vector<Poly> expected_fixed_polys = this->CreatePolys(
          base::Array2DToVector2D(TestData::kPermutationsPolys));
      EXPECT_EQ(pkey.permutation_proving_key().polys(), expected_fixed_polys);
    } else {
      EXPECT_TRUE(pkey.permutation_proving_key().polys().empty());
    }
  }
}

template <typename TestArguments, typename TestData>
void CircuitTest<TestArguments, TestData>::CreateProofTest() {
  CHECK(this->prover_->pcs().UnsafeSetup(TestData::kN, F(2)));
  this->prover_->set_domain(Domain::Create(TestData::kN));

  std::vector<Circuit> circuits = TestData::Get2Circuits();

  std::vector<Evals> instance_columns = TestData::GetInstanceColumns();
  std::vector<std::vector<Evals>> instance_columns_vec = {
      instance_columns, std::move(instance_columns)};

  ProvingKey<LS> pkey;
  ASSERT_TRUE(pkey.Load(this->prover_.get(), circuits[0]));
  this->prover_->CreateProof(pkey, std::move(instance_columns_vec), circuits);

  std::vector<uint8_t> proof =
      this->prover_->GetWriter()->buffer().owned_buffer();
  EXPECT_EQ(proof, base::ArrayToVector(TestData::kProof));
}

template <typename TestArguments, typename TestData>
void CircuitTest<TestArguments, TestData>::VerifyProofTest() {
  CHECK(this->prover_->pcs().UnsafeSetup(TestData::kN, F(2)));
  this->prover_->set_domain(Domain::Create(TestData::kN));

  Circuit circuit = TestData::GetCircuit();

  VerifyingKey<F, Commitment> vkey;
  ASSERT_TRUE(vkey.Load(this->prover_.get(), circuit));

  std::vector<uint8_t> owned_proof = base::ArrayToVector(TestData::kProof);

  halo2::Verifier<PCS, LS> verifier = this->CreateVerifier(
      this->CreateBufferWithProof(absl::MakeSpan(owned_proof)));

  std::vector<Evals> instance_columns = TestData::GetInstanceColumns();
  std::vector<std::vector<Evals>> instance_columns_vec = {
      instance_columns, std::move(instance_columns)};

  halo2::Proof<F, Commitment> proof;
  F h_eval;
  ASSERT_TRUE(verifier.VerifyProofForTesting(vkey, instance_columns_vec, &proof,
                                             &h_eval));

  if constexpr (TestData::kAdviceCommitmentsFlag) {
    std::vector<std::vector<Commitment>> expected_advice_commitments_vec{
        this->CreateCommitments(
            base::ArrayToVector(TestData::kAdviceCommitments[0])),
        this->CreateCommitments(
            base::ArrayToVector(TestData::kAdviceCommitments[1])),
    };
    EXPECT_EQ(proof.advices_commitments_vec, expected_advice_commitments_vec);
  } else {
    EXPECT_TRUE(proof.advices_commitments_vec[0].empty());
  }

  if constexpr (TestData::kChallengesFlag) {
    std::vector<F> expected_challenges =
        this->CreateEvals(base::ArrayToVector(TestData::kChallenges));
    EXPECT_EQ(proof.challenges, expected_challenges);
  } else {
    EXPECT_TRUE(proof.challenges.empty());
  }

  F expected_theta = *F::FromHexString(TestData::kTheta);
  EXPECT_EQ(proof.theta, expected_theta);

  if constexpr (TestData::kPermutationProductCommitmentsPointsFlag) {
    std::vector<std::vector<lookup::Pair<Commitment>>>
        expected_lookup_permuted_commitments_vec{
            this->CreateLookupPermutedCommitments(
                base::ArrayToVector(
                    TestData::kPermutationProductCommitmentsInputPoints[0]),
                base::ArrayToVector(
                    TestData::kPermutationProductCommitmentsTablePoints[0])),
            this->CreateLookupPermutedCommitments(
                base::ArrayToVector(
                    TestData::kPermutationProductCommitmentsInputPoints[1]),
                base::ArrayToVector(
                    TestData::kPermutationProductCommitmentsTablePoints[1])),
        };
    EXPECT_EQ(proof.lookup_permuted_commitments_vec,
              expected_lookup_permuted_commitments_vec);
  } else {
    EXPECT_TRUE(proof.lookup_permuted_commitments_vec[0].empty());
  }

  F expected_beta = *F::FromHexString(TestData::kBeta);
  EXPECT_EQ(proof.beta, expected_beta);

  F expected_gamma = *F::FromHexString(TestData::kGamma);
  EXPECT_EQ(proof.gamma, expected_gamma);

  if constexpr (TestData::kPermutationProductCommitmentsFlag) {
    std::vector<std::vector<Commitment>>
        expected_permutation_product_commitments_vec{
            this->CreateCommitments(base::ArrayToVector(
                TestData::kPermutationProductCommitments[0])),
            this->CreateCommitments(base::ArrayToVector(
                TestData::kPermutationProductCommitments[1])),
        };
    EXPECT_EQ(proof.permutation_product_commitments_vec,
              expected_permutation_product_commitments_vec);
  } else {
    EXPECT_TRUE(proof.permutation_product_commitments_vec[0].empty());
  }

  if constexpr (TestData::kLookupProductCommitmentsFlag) {
    std::vector<std::vector<Commitment>>
        expected_lookup_product_commitments_vec{
            this->CreateCommitments(
                base::ArrayToVector(TestData::kLookupProductCommitments[0])),
            this->CreateCommitments(
                base::ArrayToVector(TestData::kLookupProductCommitments[1])),
        };
    EXPECT_EQ(proof.lookup_product_commitments_vec,
              expected_lookup_product_commitments_vec);
  } else {
    EXPECT_TRUE(proof.lookup_product_commitments_vec[0].empty());
  }

  Commitment expected_vanishing_random_poly_commitment =
      this->CreateCommitment(TestData::kVanishingRandomPolyCommitment);
  EXPECT_EQ(proof.vanishing_random_poly_commitment,
            expected_vanishing_random_poly_commitment);

  F expected_y = *F::FromHexString(TestData::kY);
  EXPECT_EQ(proof.y, expected_y);

  if constexpr (TestData::kVanishingHPolyCommitmentsFlag) {
    std::vector<Commitment> expected_vanishing_h_poly_commitments =
        this->CreateCommitments(
            base::ArrayToVector(TestData::kVanishingHPolyCommitments));
    EXPECT_EQ(proof.vanishing_h_poly_commitments,
              expected_vanishing_h_poly_commitments);
  } else {
    EXPECT_TRUE(proof.vanishing_h_poly_commitments.empty());
  }

  F expected_x = *F::FromHexString(TestData::kX);
  EXPECT_EQ(proof.x, expected_x);

  if constexpr (TestData::kAdviceEvalsFlag) {
    std::vector<std::vector<F>> expected_advice_evals_vec{
        this->CreateEvals(base::ArrayToVector(TestData::kAdviceEvals[0])),
        this->CreateEvals(base::ArrayToVector(TestData::kAdviceEvals[1])),
    };
    EXPECT_EQ(proof.advice_evals_vec, expected_advice_evals_vec);
  } else {
    EXPECT_TRUE(proof.advice_evals_vec[0].empty());
  }

  if constexpr (TestData::kFixedEvalsFlag) {
    std::vector<F> expected_fixed_evals =
        this->CreateEvals(base::ArrayToVector(TestData::kFixedEvals));
    EXPECT_EQ(proof.fixed_evals, expected_fixed_evals);
  } else {
    EXPECT_TRUE(proof.fixed_evals.empty());
  }

  F expected_vanishing_random_eval =
      *F::FromHexString(TestData::kVanishingRandomEval);
  EXPECT_EQ(proof.vanishing_random_eval, expected_vanishing_random_eval);

  if constexpr (TestData::kCommonPermutationEvalsFlag) {
    std::vector<F> expected_common_permutation_evals = this->CreateEvals(
        base::ArrayToVector(TestData::kCommonPermutationEvals));
    EXPECT_EQ(proof.common_permutation_evals,
              expected_common_permutation_evals);
  } else {
    EXPECT_TRUE(proof.common_permutation_evals.empty());
  }

  if constexpr (TestData::kPermutationProductEvalsFlag) {
    std::vector<std::vector<F>> expected_permutation_product_evals_vec{
        this->CreateEvals(
            base::ArrayToVector(TestData::kPermutationProductEvals[0])),
        this->CreateEvals(
            base::ArrayToVector(TestData::kPermutationProductEvals[1])),
    };
    EXPECT_EQ(proof.permutation_product_evals_vec,
              expected_permutation_product_evals_vec);
  } else {
    EXPECT_TRUE(proof.permutation_product_evals_vec[0].empty());
  }

  if constexpr (TestData::kPermutationProductNextEvalsFlag) {
    std::vector<std::vector<F>> expected_permutation_product_next_evals_vec{
        this->CreateEvals(
            base::ArrayToVector(TestData::kPermutationProductNextEvals[0])),
        this->CreateEvals(
            base::ArrayToVector(TestData::kPermutationProductNextEvals[1])),
    };
    EXPECT_EQ(proof.permutation_product_next_evals_vec,
              expected_permutation_product_next_evals_vec);
  } else {
    EXPECT_TRUE(proof.permutation_product_next_evals_vec[0].empty());
  }

  if constexpr (TestData::kPermutationProductLastEvalsFlag) {
    std::vector<std::vector<std::optional<F>>>
        expected_permutation_product_last_evals_vec{
            this->CreateOptionalEvals(
                base::ArrayToVector(TestData::kPermutationProductLastEvals[0])),
            this->CreateOptionalEvals(
                base::ArrayToVector(TestData::kPermutationProductLastEvals[1])),
        };
    EXPECT_EQ(proof.permutation_product_last_evals_vec,
              expected_permutation_product_last_evals_vec);
  } else {
    EXPECT_TRUE(proof.permutation_product_last_evals_vec[0].empty());
  }

  if constexpr (TestData::kLookupProductEvalsFlag) {
    std::vector<std::vector<F>> expected_lookup_product_evals_vec{
        this->CreateEvals(
            base::ArrayToVector(TestData::kLookupProductEvals[0])),
        this->CreateEvals(
            base::ArrayToVector(TestData::kLookupProductEvals[1])),
    };
    EXPECT_EQ(proof.lookup_product_evals_vec,
              expected_lookup_product_evals_vec);
  } else {
    EXPECT_TRUE(proof.lookup_product_evals_vec[0].empty());
  }

  if constexpr (TestData::kLookupProductNextEvalsFlag) {
    std::vector<std::vector<F>> expected_lookup_product_next_evals_vec{
        this->CreateEvals(
            base::ArrayToVector(TestData::kLookupProductNextEvals[0])),
        this->CreateEvals(
            base::ArrayToVector(TestData::kLookupProductNextEvals[1])),
    };
    EXPECT_EQ(proof.lookup_product_next_evals_vec,
              expected_lookup_product_next_evals_vec);
  } else {
    EXPECT_TRUE(proof.lookup_product_next_evals_vec[0].empty());
  }

  if constexpr (TestData::kLookupPermutedInputEvalsFlag) {
    std::vector<std::vector<F>> expected_lookup_permuted_input_evals_vec{
        this->CreateEvals(
            base::ArrayToVector(TestData::kLookupPermutedInputEvals[0])),
        this->CreateEvals(
            base::ArrayToVector(TestData::kLookupPermutedInputEvals[1])),
    };
    EXPECT_EQ(proof.lookup_permuted_input_evals_vec,
              expected_lookup_permuted_input_evals_vec);
  } else {
    EXPECT_TRUE(proof.lookup_permuted_input_evals_vec[0].empty());
  }

  if constexpr (TestData::kLookupPermutedInputPrevEvalsFlag) {
    std::vector<std::vector<F>> expected_lookup_permuted_input_prev_evals_vec{
        this->CreateEvals(
            base::ArrayToVector(TestData::kLookupPermutedInputPrevEvals[0])),
        this->CreateEvals(
            base::ArrayToVector(TestData::kLookupPermutedInputPrevEvals[1])),
    };
    EXPECT_EQ(proof.lookup_permuted_input_prev_evals_vec,
              expected_lookup_permuted_input_prev_evals_vec);
  } else {
    EXPECT_TRUE(proof.lookup_permuted_input_prev_evals_vec[0].empty());
  }

  if constexpr (TestData::kLookupPermutedTableEvalsFlag) {
    std::vector<std::vector<F>> expected_lookup_permuted_table_evals_vec{
        this->CreateEvals(
            base::ArrayToVector(TestData::kLookupPermutedTableEvals[0])),
        this->CreateEvals(
            base::ArrayToVector(TestData::kLookupPermutedTableEvals[1])),
    };
    EXPECT_EQ(proof.lookup_permuted_table_evals_vec,
              expected_lookup_permuted_table_evals_vec);
  } else {
    EXPECT_TRUE(proof.lookup_permuted_table_evals_vec[0].empty());
  }

  F expected_h_eval = *F::FromHexString(TestData::kHEval);
  EXPECT_EQ(h_eval, expected_h_eval);
}

namespace {

const size_t kBits = 3;

}  // namespace

template class CircuitTest<
    TestArguments<SimpleCircuit<BN254SHPlonk::Field, SimpleFloorPlanner>,
                  BN254SHPlonk, BN254Halo2LS>,
    SimpleTestData<SimpleCircuit<BN254SHPlonk::Field, SimpleFloorPlanner>,
                   BN254SHPlonk, BN254Halo2LS>>;
template class CircuitTest<
    TestArguments<SimpleCircuit<BN254SHPlonk::Field, V1FloorPlanner>,
                  BN254SHPlonk, BN254Halo2LS>,
    SimpleTestData<SimpleCircuit<BN254SHPlonk::Field, V1FloorPlanner>,
                   BN254SHPlonk, BN254Halo2LS>>;

template class CircuitTest<
    TestArguments<
        SimpleLookupCircuit<BN254SHPlonk::Field, kBits, SimpleFloorPlanner>,
        BN254SHPlonk, BN254Halo2LS>,
    SimpleLookupTestData<
        SimpleLookupCircuit<BN254SHPlonk::Field, kBits, SimpleFloorPlanner>,
        BN254SHPlonk, BN254Halo2LS>>;
template class CircuitTest<
    TestArguments<
        SimpleLookupCircuit<BN254SHPlonk::Field, kBits, V1FloorPlanner>,
        BN254SHPlonk, BN254Halo2LS>,
    SimpleLookupTestData<
        SimpleLookupCircuit<BN254SHPlonk::Field, kBits, V1FloorPlanner>,
        BN254SHPlonk, BN254Halo2LS>>;

}  // namespace tachyon::zk::plonk
