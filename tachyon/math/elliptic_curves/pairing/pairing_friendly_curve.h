#ifndef TACHYON_MATH_ELLIPTIC_CURVES_PAIRING_PAIRING_FRIENDLY_CURVE_H_
#define TACHYON_MATH_ELLIPTIC_CURVES_PAIRING_PAIRING_FRIENDLY_CURVE_H_

#include <vector>

#include "tachyon/math/elliptic_curves/pairing/ell_coeff.h"
#include "tachyon/math/elliptic_curves/pairing/twist_type.h"

namespace tachyon::math {

template <typename _Config>
class PairingFriendlyCurve {
 public:
  using Config = _Config;
  using G1Curve = typename Config::G1Curve;
  using G2Curve = typename Config::G2Curve;
  using Fp2 = typename G2Curve::BaseField;
  using Fp12 = typename Config::Fp12;
  using G1AffinePoint = typename G1Curve::AffinePoint;

  static void Init() {
    Config::Init();
    Fp12::Init();
  }

 protected:
  class Pair {
   public:
    Pair() = default;
    Pair(const G1AffinePoint* g1, const std::vector<EllCoeff<Fp2>>* ell_coeffs)
        : g1_(g1), ell_coeffs_(ell_coeffs) {}

    const G1AffinePoint& g1() const { return *g1_; }
    const EllCoeff<Fp2>& NextEllCoeff() const { return (*ell_coeffs_)[idx_++]; }

   private:
    const G1AffinePoint* g1_ = nullptr;
    const std::vector<EllCoeff<Fp2>>* ell_coeffs_ = nullptr;
    mutable size_t idx_ = 0;
  };

  static Fp12 PowByX(const Fp12& f_in) {
    Fp12 f = f_in.CyclotomicPow(Config::kX);
    if constexpr (Config::kXIsNegative) {
      CHECK(f.CyclotomicInverseInPlace());
    }
    return f;
  }

  static Fp12 PowByNegX(const Fp12& f_in) {
    Fp12 f = f_in.CyclotomicPow(Config::kX);
    if constexpr (!Config::kXIsNegative) {
      CHECK(f.CyclotomicInverseInPlace());
    }
    return f;
  }

  // TODO(chokobole): Leave a comment to help understand readers.
  // Evaluates the line function at point |p|.
  static void Ell(Fp12& f, const EllCoeff<Fp2>& coeffs,
                  const G1AffinePoint& p) {
    if constexpr (Config::kTwistType == TwistType::kM) {
      f.MulInPlaceBy014(coeffs.c0(), coeffs.c1() * p.x(), coeffs.c2() * p.y());
    } else {
      f.MulInPlaceBy034(coeffs.c0() * p.y(), coeffs.c1() * p.x(), coeffs.c2());
    }
  }

  template <typename G1AffinePointContainer, typename G2PreparedContainer>
  static std::vector<Pair> CreatePairs(const G1AffinePointContainer& a,
                                       const G2PreparedContainer& b) {
    size_t size = std::size(a);
    CHECK_EQ(size, std::size(b));
    std::vector<Pair> pairs;
    pairs.reserve(size);
    for (size_t i = 0; i < size; ++i) {
      if (!a[i].IsZero() && !b[i].infinity()) {
        pairs.emplace_back(&a[i], &b[i].ell_coeffs());
      }
    }
    return pairs;
  }
};

}  // namespace tachyon::math

#endif  // TACHYON_MATH_ELLIPTIC_CURVES_PAIRING_PAIRING_FRIENDLY_CURVE_H_
