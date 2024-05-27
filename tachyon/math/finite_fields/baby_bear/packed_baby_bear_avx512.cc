// Copyright (c) 2022 The Plonky3 Authors
// Use of this source code is governed by a MIT/Apache-2.0 style license that
// can be found in the LICENSE-MIT.plonky3 and the LICENCE-APACHE.plonky3
// file.

#include "tachyon/math/finite_fields/baby_bear/packed_baby_bear_avx512.h"

#include <immintrin.h>

#include "tachyon/math/finite_fields/packed_prime_field32_avx512.h"

namespace tachyon::math {

namespace {

__m512i kP;
__m512i kInv;
__m512i kZero;
__m512i kOne;

__m512i ToVector(const PackedBabyBearAVX512& packed) {
  return _mm512_loadu_si512(
      reinterpret_cast<const __m512i_u*>(packed.values().data()));
}

PackedBabyBearAVX512 FromVector(__m512i vector) {
  PackedBabyBearAVX512 ret;
  _mm512_storeu_si512(reinterpret_cast<__m512i_u*>(ret.values().data()),
                      vector);
  return ret;
}

__m512i Add(__m512i lhs, __m512i rhs) { return AddMod32(lhs, rhs, kP); }

__m512i Sub(__m512i lhs, __m512i rhs) { return SubMod32(lhs, rhs, kP); }

__m512i Negate(__m512i val) { return NegateMod32(val, kP); }

__m512i Mul(__m512i lhs, __m512i rhs) {
  return MontMulMod32(lhs, rhs, kP, kInv);
}

}  // namespace

// static
void PackedBabyBearAVX512::Init() {
  kP = _mm512_set1_epi32(BabyBear::Config::kModulus);
  kInv = _mm512_set1_epi32(BabyBear::Config::kInverse32);
  kZero = _mm512_set1_epi32(0);
  kOne = _mm512_set1_epi32(BabyBear::Config::kOne);
}

// static
PackedBabyBearAVX512 PackedBabyBearAVX512::Zero() { return FromVector(kZero); }

// static
PackedBabyBearAVX512 PackedBabyBearAVX512::One() { return FromVector(kOne); }

// static
PackedBabyBearAVX512 PackedBabyBearAVX512::Broadcast(const PrimeField& value) {
  return FromVector(_mm512_set1_epi32(value.value()));
}

PackedBabyBearAVX512 PackedBabyBearAVX512::Add(
    const PackedBabyBearAVX512& other) const {
  return FromVector(tachyon::math::Add(ToVector(*this), ToVector(other)));
}

PackedBabyBearAVX512 PackedBabyBearAVX512::Sub(
    const PackedBabyBearAVX512& other) const {
  return FromVector(tachyon::math::Sub(ToVector(*this), ToVector(other)));
}

PackedBabyBearAVX512 PackedBabyBearAVX512::Negate() const {
  return FromVector(tachyon::math::Negate(ToVector(*this)));
}

PackedBabyBearAVX512 PackedBabyBearAVX512::Mul(
    const PackedBabyBearAVX512& other) const {
  return FromVector(tachyon::math::Mul(ToVector(*this), ToVector(other)));
}

}  // namespace tachyon::math
