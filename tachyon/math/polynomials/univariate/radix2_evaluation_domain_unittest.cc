#include "tachyon/math/polynomials/univariate/radix2_evaluation_domain.h"

#include <memory>

#include "gtest/gtest.h"

#include "tachyon/math/finite_fields/baby_bear/packed_baby_bear.h"
#include "tachyon/math/finite_fields/koala_bear/packed_koala_bear.h"
#include "tachyon/math/finite_fields/test/finite_field_test.h"
#include "tachyon/math/matrix/matrix_types.h"
#include "tachyon/math/polynomials/univariate/naive_batch_fft.h"

namespace tachyon::math {

namespace {

template <typename F>
class Radix2EvaluationDomainTest
    : public FiniteFieldTest<
          typename PackedPrimeFieldTraits<F>::PackedPrimeField> {};

}  // namespace

using PrimeFieldTypes = testing::Types<BabyBear, KoalaBear>;
TYPED_TEST_SUITE(Radix2EvaluationDomainTest, PrimeFieldTypes);

TYPED_TEST(Radix2EvaluationDomainTest, FFTBatch) {
  using F = TypeParam;
  for (size_t log_r = 0; log_r < 5; ++log_r) {
    RowMajorMatrix<F> expected = RowMajorMatrix<F>::Random(1 << log_r, 3);
    RowMajorMatrix<F> result = expected;
    NaiveBatchFFT<F> naive;
    naive.FFTBatch(expected);
    std::unique_ptr<Radix2EvaluationDomain<F>> domain =
        Radix2EvaluationDomain<F>::Create(1 << log_r);
    domain->FFTBatch(result);
    EXPECT_EQ(expected, result);
  }
}

TYPED_TEST(Radix2EvaluationDomainTest, CosetLDEBatch) {
  using F = TypeParam;
  for (size_t log_r = 0; log_r < 5; ++log_r) {
    RowMajorMatrix<F> expected = RowMajorMatrix<F>::Random(1 << log_r, 3);
    RowMajorMatrix<F> result = expected;
    NaiveBatchFFT<F> naive;
    naive.CosetLDEBatch(expected, 1);
    std::unique_ptr<Radix2EvaluationDomain<F>> domain =
        Radix2EvaluationDomain<F>::Create(1 << log_r);
    domain->CosetLDEBatch(result, 1);
    EXPECT_EQ(expected, result);
  }
}

}  // namespace tachyon::math
