// Copyright (c) 2022 The Plonky3 Authors
// Use of this source code is governed by a MIT/Apache-2.0 style license that
// can be found in the LICENSE-MIT.plonky3 and the LICENCE-APACHE.plonky3
// file.

#include "tachyon/math/finite_fields/koala_bear/packed_koala_bear_avx2.h"

#include <immintrin.h>

#include "tachyon/math/finite_fields/packed_prime_field32_avx2.h"

namespace tachyon::math {

namespace {

__m256i kP;
__m256i kInv;
__m256i kZero;
__m256i kOne;

__m256i ToVector(const PackedKoalaBearAVX2& packed) {
  return _mm256_loadu_si256(
      reinterpret_cast<const __m256i_u*>(packed.values().data()));
}

PackedKoalaBearAVX2 FromVector(__m256i vector) {
  PackedKoalaBearAVX2 ret;
  _mm256_storeu_si256(reinterpret_cast<__m256i_u*>(ret.values().data()),
                      vector);
  return ret;
}

__m256i Add(__m256i lhs, __m256i rhs) { return AddMod32(lhs, rhs, kP); }

__m256i Sub(__m256i lhs, __m256i rhs) { return SubMod32(lhs, rhs, kP); }

__m256i Negate(__m256i val) { return NegateMod32(val, kP); }

__m256i Mul(__m256i lhs, __m256i rhs) {
  return MontMulMod32(lhs, rhs, kP, kInv);
}

}  // namespace

// static
void PackedKoalaBearAVX2::Init() {
  kP = _mm256_set1_epi32(KoalaBear::Config::kModulus);
  kInv = _mm256_set1_epi32(KoalaBear::Config::kInverse32);
  kZero = _mm256_set1_epi32(0);
  kOne = _mm256_set1_epi32(KoalaBear::Config::kOne);
}

// static
PackedKoalaBearAVX2 PackedKoalaBearAVX2::Zero() { return FromVector(kZero); }

// static
PackedKoalaBearAVX2 PackedKoalaBearAVX2::One() { return FromVector(kOne); }

// static
PackedKoalaBearAVX2 PackedKoalaBearAVX2::Broadcast(const PrimeField& value) {
  return FromVector(_mm256_set1_epi32(value.value()));
}

PackedKoalaBearAVX2 PackedKoalaBearAVX2::Add(
    const PackedKoalaBearAVX2& other) const {
  return FromVector(tachyon::math::Add(ToVector(*this), ToVector(other)));
}

PackedKoalaBearAVX2 PackedKoalaBearAVX2::Sub(
    const PackedKoalaBearAVX2& other) const {
  return FromVector(tachyon::math::Sub(ToVector(*this), ToVector(other)));
}

PackedKoalaBearAVX2 PackedKoalaBearAVX2::Negate() const {
  return FromVector(tachyon::math::Negate(ToVector(*this)));
}

PackedKoalaBearAVX2 PackedKoalaBearAVX2::Mul(
    const PackedKoalaBearAVX2& other) const {
  return FromVector(tachyon::math::Mul(ToVector(*this), ToVector(other)));
}

}  // namespace tachyon::math
