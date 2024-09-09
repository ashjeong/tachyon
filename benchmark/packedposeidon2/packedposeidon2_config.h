#ifndef BENCHMARK_PACKEDPOSEIDON2_PACKEDPOSEIDON2_CONFIG_H_
#define BENCHMARK_PACKEDPOSEIDON2_PACKEDPOSEIDON2_CONFIG_H_

#include <stddef.h>

#include <vector>

// clang-format off
#include "benchmark/config.h"
#include "benchmark/field_type.h"
// clang-format on

namespace tachyon::benchmark {

class PackedPoseidon2Config : public Config {
 public:
  PackedPoseidon2Config();
  PackedPoseidon2Config(const PackedPoseidon2Config& other) = delete;
  PackedPoseidon2Config& operator=(const PackedPoseidon2Config& other) = delete;

  size_t repeating_num() const { return repeating_num_; }
  FieldType prime_field() const { return prime_field_; }

 private:
  // Config methods
  bool Validate() const override;

  size_t repeating_num_;
  FieldType prime_field_;
};

}  // namespace tachyon::benchmark

#endif  // BENCHMARK_PACKEDPOSEIDON2_PACKEDPOSEIDON2_CONFIG_H_
