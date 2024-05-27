#ifndef VENDORS_CIRCOM_CIRCOMLIB_BASE_PRIME_FIELD_H_
#define VENDORS_CIRCOM_CIRCOMLIB_BASE_PRIME_FIELD_H_

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "tachyon/base/buffer/endian_auto_reset.h"
#include "tachyon/base/logging.h"
#include "tachyon/math/base/big_int.h"

namespace tachyon::circom {

struct PrimeField {
  std::vector<uint8_t> bytes;

  bool operator==(const PrimeField& other) const {
    return bytes == other.bytes;
  }
  bool operator!=(const PrimeField& other) const {
    return bytes != other.bytes;
  }

  template <size_t N>
  math::BigInt<N> ToBigInt() const {
    CHECK_EQ(bytes.size() / 8, N);
    return math::BigInt<N>::FromBytesLE(bytes);
  }

  template <size_t N>
  static PrimeField FromBigInt(const math::BigInt<N>& big_int) {
    std::array<uint8_t, N * 8> bytes = big_int.ToBytesLE();
    return {{bytes.begin(), bytes.end()}};
  }

  template <bool IsMontgomery, typename F>
  F ToNative() const {
    if constexpr (IsMontgomery) {
      return F::FromMontgomery(ToBigInt<F::kLimbNums>());
    } else {
      return F(ToBigInt<F::kLimbNums>());
    }
  }

  template <bool IsMontgomery, typename F>
  static PrimeField FromNative(const F& prime_field) {
    if constexpr (IsMontgomery) {
      return FromBigInt(prime_field.value());
    } else {
      return FromBigInt(prime_field.ToBigInt());
    }
  }

  template <typename F>
  void Normalize() {
    std::array<uint8_t, F::kLimbNums * 8> bytes_to_copy =
        ToBigInt<F::kLimbNums>().Mod(F::Config::kModulus).ToBytesLE();
    memcpy(bytes.data(), bytes_to_copy.data(), bytes_to_copy.size());
  }

  bool Read(const base::ReadOnlyBuffer& buffer, uint32_t num_bytes = 0) {
    base::EndianAutoReset reset(buffer, base::Endian::kLittle);
    if (num_bytes == 0) {
      if (!buffer.Read(&num_bytes)) return false;
      if (num_bytes % 8 != 0) {
        LOG(ERROR) << "field size is not a multiple of 8";
        return false;
      }
    }
    bytes.resize(num_bytes);
    return buffer.Read(bytes.data(), bytes.size());
  }

  std::string ToString() const;
};

}  // namespace tachyon::circom

#endif  // VENDORS_CIRCOM_CIRCOMLIB_BASE_PRIME_FIELD_H_
