#include <optional>

#include "gtest/gtest.h"

#include "tachyon/c/math/elliptic_curves/bn/bn254/fr_traits.h"
#include "tachyon/cc/math/elliptic_curves/bn/bn254/fr.h"
#include "tachyon/math/elliptic_curves/bn/bn254/fr.h"

namespace tachyon::cc::math {

namespace {

class PrimeFieldTest : public testing::Test {
 public:
  void SetUp() override {
    a_ = tachyon::math::bn254::Fr::Random();
    b_ = tachyon::math::bn254::Fr::Random();

    cc_a_ = bn254::Fr(c::base::c_cast(a_));
    cc_b_ = bn254::Fr(c::base::c_cast(b_));
  }

 protected:
  tachyon::math::bn254::Fr a_;
  tachyon::math::bn254::Fr b_;
  bn254::Fr cc_a_;
  bn254::Fr cc_b_;
};

}  // namespace

TEST_F(PrimeFieldTest, Zero) {
  bn254::Fr zero = bn254::Fr::Zero();
  EXPECT_TRUE(c::base::native_cast(zero.value()).IsZero());
}

TEST_F(PrimeFieldTest, One) {
  bn254::Fr one = bn254::Fr::One();
  EXPECT_TRUE(c::base::native_cast(one.value()).IsOne());
}

TEST_F(PrimeFieldTest, Random) {
  bn254::Fr random = bn254::Fr::Random();
  EXPECT_NE(c::base::native_cast(random.value()), a_);
}

TEST_F(PrimeFieldTest, Add) {
  EXPECT_EQ(c::base::native_cast((cc_a_ + cc_b_).value()), a_ + b_);
  EXPECT_EQ(c::base::native_cast((cc_a_ += cc_b_).value()), a_ += b_);
}

TEST_F(PrimeFieldTest, Sub) {
  EXPECT_EQ(c::base::native_cast((cc_a_ - cc_b_).value()), a_ - b_);
  EXPECT_EQ(c::base::native_cast((cc_a_ -= cc_b_).value()), a_ -= b_);
}

TEST_F(PrimeFieldTest, Mul) {
  EXPECT_EQ(c::base::native_cast((cc_a_ * cc_b_).value()), a_ * b_);
  EXPECT_EQ(c::base::native_cast((cc_a_ *= cc_b_).value()), a_ *= b_);
}

TEST_F(PrimeFieldTest, Div) {
  std::optional<bn254::Fr> a_div_b = cc_a_ / cc_b_;
  CHECK(a_div_b);
  EXPECT_EQ(c::base::native_cast((*a_div_b).value()), a_ / b_);
  CHECK(cc_a_ /= cc_b_);
  CHECK(a_ /= b_);
  EXPECT_EQ(c::base::native_cast(cc_a_.value()), a_);
}

TEST_F(PrimeFieldTest, Negate) {
  EXPECT_EQ(c::base::native_cast((-cc_a_).value()), -a_);
}

TEST_F(PrimeFieldTest, Double) {
  EXPECT_EQ(c::base::native_cast(cc_a_.Double().value()), a_.Double());
}

TEST_F(PrimeFieldTest, Square) {
  EXPECT_EQ(c::base::native_cast(cc_a_.Square().value()), a_.Square());
}

TEST_F(PrimeFieldTest, Inverse) {
  std::optional<bn254::Fr> cc_a_inv = cc_a_.Inverse();
  std::optional<tachyon::math::bn254::Fr> a_inv = a_.Inverse();
  EXPECT_TRUE(cc_a_inv && a_inv);
  EXPECT_EQ(c::base::native_cast((*cc_a_inv).value()), *a_inv);
}

TEST_F(PrimeFieldTest, Eq) { EXPECT_EQ(cc_a_ == cc_b_, a_ == b_); }

TEST_F(PrimeFieldTest, Ne) { EXPECT_EQ(cc_a_ != cc_b_, a_ != b_); }

TEST_F(PrimeFieldTest, Gt) { EXPECT_EQ(cc_a_ > cc_b_, a_ > b_); }

TEST_F(PrimeFieldTest, Ge) { EXPECT_EQ(cc_a_ >= cc_b_, a_ >= b_); }

TEST_F(PrimeFieldTest, Lt) { EXPECT_EQ(cc_a_ < cc_b_, a_ < b_); }

TEST_F(PrimeFieldTest, Le) { EXPECT_EQ(cc_a_ <= cc_b_, a_ <= b_); }

}  // namespace tachyon::cc::math
