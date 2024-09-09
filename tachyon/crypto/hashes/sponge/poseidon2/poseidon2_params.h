#ifndef TACHYON_CRYPTO_HASHES_SPONGE_POSEIDON2_POSEIDON2_PARAMS_H_
#define TACHYON_CRYPTO_HASHES_SPONGE_POSEIDON2_POSEIDON2_PARAMS_H_

#include <stddef.h>
#include <stdint.h>

#include <type_traits>

#include "tachyon/base/bits.h"
#include "tachyon/math/elliptic_curves/bn/bn254/fr.h"
#include "tachyon/math/finite_fields/finite_field_traits.h"

namespace tachyon::crypto {

template <typename F, size_t Width, uint64_t Alpha>
constexpr size_t GetPoseidon2PartialRounds() {
  using PrimeField = math::MaybeUnpack<F>;

  if constexpr (PrimeField::Config::kModulusBits == 31) {
    if constexpr (Width == 16) {
      switch (Alpha) {
        case 3:
          return 20;
        case 5:
          return 14;
        case 7:
          return 13;
        case 9:
          return 13;
        case 11:
          return 13;
      }
    } else if constexpr (Width == 24) {
      switch (Alpha) {
        case 3:
          return 23;
        case 5:
          return 22;
        case 7:
          return 21;
        case 9:
          return 21;
        case 11:
          return 21;
      }
    }
  } else if constexpr (PrimeField::Config::kModulusBits == 64) {
    if constexpr (Width == 8) {
      switch (Alpha) {
        case 3:
          return 41;
        case 5:
          return 27;
        case 7:
          return 22;
        case 9:
          return 19;
        case 11:
          return 17;
      }
    } else if constexpr (Width == 16) {
      switch (Alpha) {
        case 3:
          return 42;
        case 5:
          return 27;
        case 7:
          return 22;
        case 9:
          return 20;
        case 11:
          return 18;
      }
    } else if constexpr (Width == 24) {
      switch (Alpha) {
        case 3:
          return 47;
        case 5:
          return 27;
        case 7:
          return 22;
        case 9:
          return 20;
        case 11:
          return 18;
      }
    }
  }
  return 56;
}

template <typename _Field, size_t _Rate, uint32_t _Alpha, size_t _Capacity = 1,
          size_t _FullRounds = 8,
          size_t _PartialRounds =
              (GetPoseidon2PartialRounds<_Field, _Rate + _Capacity, _Alpha>())>
struct Poseidon2Params {
  using Field = _Field;
  constexpr static size_t Rate = _Rate;
  constexpr static size_t Capacity = _Capacity;
  constexpr static size_t Width = Rate + Capacity;
  // NOTE(ashjeong): |Alpha| is also referred to as |D|
  constexpr static uint32_t Alpha = _Alpha;
  constexpr static size_t FullRounds = _FullRounds;
  constexpr static size_t PartialRounds = _PartialRounds;
};

}  // namespace tachyon::crypto

#endif  // TACHYON_CRYPTO_HASHES_SPONGE_POSEIDON2_POSEIDON2_PARAMS_H_
