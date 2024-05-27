#ifndef TACHYON_MATH_FINITE_FIELDS_PACKED_PRIME_FIELD_BASE_H_
#define TACHYON_MATH_FINITE_FIELDS_PACKED_PRIME_FIELD_BASE_H_

#include <stddef.h>

#include <array>
#include <optional>
#include <string>

#include "tachyon/base/containers/container_util.h"
#include "tachyon/base/strings/string_util.h"
#include "tachyon/math/finite_fields/packed_prime_field_traits_forward.h"

namespace tachyon::math {

template <typename Derived>
class PackedPrimeFieldBase : public Field<Derived> {
 public:
  using PrimeField = typename PackedPrimeFieldTraits<Derived>::PrimeField;

  constexpr static size_t N = PackedPrimeFieldTraits<Derived>::N;

  static Derived Random() {
    Derived ret;
    for (size_t i = 0; i < N; ++i) {
      ret.values_[i] = PrimeField::Random();
    }
    return ret;
  }

  const std::array<PrimeField, N>& values() const { return values_; }
  std::array<PrimeField, N>& values() { return values_; }

  constexpr bool IsZero() const {
    for (size_t i = 0; i < N; ++i) {
      if (!values_[i].IsZero()) return false;
    }
    return true;
  }

  constexpr bool IsOne() const {
    for (size_t i = 0; i < N; ++i) {
      if (!values_[i].IsOne()) return false;
    }
    return true;
  }

  std::string ToString() const { return base::ContainerToString(values_); }
  std::string ToHexString(bool pad_zero = false) const {
    return base::ContainerToString(
        base::Map(values_, [pad_zero](const PrimeField& value) {
          return value.ToHexString(pad_zero);
        }));
  }

  constexpr PrimeField& operator[](size_t i) { return values_[i]; }
  constexpr const PrimeField& operator[](size_t i) const { return values_[i]; }

  constexpr bool operator==(const Derived& other) const {
    return values_ == other.values_;
  }
  constexpr bool operator!=(const Derived& other) const {
    return values_ != other.values_;
  }

  // AdditiveSemigroup methods
  Derived& AddInPlace(const Derived& other) {
    Derived& self = static_cast<Derived&>(*this);
    return self = self + other;
  }

  // AdditiveGroup methods
  Derived& SubInPlace(const Derived& other) {
    Derived& self = static_cast<Derived&>(*this);
    return self = self - other;
  }

  Derived& NegateInPlace() {
    Derived& self = static_cast<Derived&>(*this);
    return self = self.Negate();
  }

  // MultiplicativeSemigroup methods
  Derived& MulInPlace(const Derived& other) {
    Derived& self = static_cast<Derived&>(*this);
    return self = self * other;
  }

  // MultiplicativeGroup methods
  std::optional<Derived> Inverse() const {
    NOTIMPLEMENTED();
    return std::nullopt;
  }

  std::optional<Derived*> InverseInPlace() {
    NOTIMPLEMENTED();
    return std::nullopt;
  }

 protected:
  std::array<PrimeField, N> values_;
};

}  // namespace tachyon::math

#endif  // TACHYON_MATH_FINITE_FIELDS_PACKED_PRIME_FIELD_BASE_H_
