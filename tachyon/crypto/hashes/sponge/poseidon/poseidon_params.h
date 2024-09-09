#ifndef TACHYON_CRYPTO_HASHES_SPONGE_POSEIDON_POSEIDON_PARAMS_H_
#define TACHYON_CRYPTO_HASHES_SPONGE_POSEIDON_POSEIDON_PARAMS_H_

#include <stddef.h>
#include <stdint.h>

#include "tachyon/math/elliptic_curves/bn/bn254/fr.h"
#include "tachyon/math/finite_fields/baby_bear/internal/baby_bear.h"

namespace tachyon::crypto {

template <typename _Field, size_t _Rate, uint32_t _Alpha, size_t _FullRounds,
          size_t _PartialRounds, size_t _Capacity = 1>
struct PoseidonParams {
  using Field = _Field;
  constexpr static size_t Rate = _Rate;
  constexpr static size_t Capacity = _Capacity;
  constexpr static size_t Width = Rate + Capacity;
  // NOTE(ashjeong): |Alpha| is also referred to as |D|
  constexpr static uint32_t Alpha = _Alpha;
  constexpr static size_t FullRounds = _FullRounds;
  constexpr static size_t PartialRounds = _PartialRounds;
};

// NOTE(ashjeong): The variables names' ending number refers to the |Width|;
// however, note that the |PartialRounds| and |FullRounds| depend on both
// |Width| and |Alpha|.
using BabyBearPoseidonParams16 = PoseidonParams<math::BabyBear, 15, 7, 8, 22>;
using BabyBearPoseidonParams24 = PoseidonParams<math::BabyBear, 23, 7, 8, 22>;
using BN254PoseidonParams9 = PoseidonParams<math::bn254::Fr, 8, 5, 8, 63>;
using BN254PoseidonParams5 = PoseidonParams<math::bn254::Fr, 4, 5, 8, 60>;

}  // namespace tachyon::crypto

#endif  // TACHYON_CRYPTO_HASHES_SPONGE_POSEIDON_POSEIDON_PARAMS_H_
