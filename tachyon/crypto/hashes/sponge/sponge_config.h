#ifndef TACHYON_CRYPTO_HASHES_SPONGE_SPONGE_CONFIG_H_
#define TACHYON_CRYPTO_HASHES_SPONGE_SPONGE_CONFIG_H_

#include <stddef.h>

#include "tachyon/base/buffer/copyable.h"
#include "tachyon/export.h"

namespace tachyon::crypto {

template <size_t _Rate, size_t _Capacity>
struct TACHYON_EXPORT SpongeConfig {
  // The rate (in terms of number of field elements).
  // See https://iacr.org/archive/eurocrypt2008/49650180/49650180.pdf
  constexpr static size_t Rate = _Rate;

  // The capacity (in terms of number of field elements).
  constexpr static size_t Capacity = _Capacity;
};

}  // namespace tachyon::crypto

#endif  // TACHYON_CRYPTO_HASHES_SPONGE_SPONGE_CONFIG_H_
