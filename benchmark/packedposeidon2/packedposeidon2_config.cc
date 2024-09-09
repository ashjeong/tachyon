#include "benchmark/packedposeidon2/packedposeidon2_config.h"

#include <set>

#include "tachyon/base/console/iostream.h"

namespace tachyon::benchmark {

PackedPoseidon2Config::PackedPoseidon2Config() {
  parser_.AddFlag<base::Flag<size_t>>(&repeating_num_)
      .set_short_name("-n")
      .set_default_value(10)
      .set_help("Specify the number of repetition 'n'. By default, 10.");
  parser_.AddFlag<base::Flag<FieldType>>(&prime_field_)
      .set_short_name("-p")
      .set_long_name("--prime_field")
      .set_required()
      .set_help(
          "A prime field to be benchmarked with. (supported prime fields: "
          "baby_bear, bn254_fr)");
  parser_.AddFlag<base::Flag<std::set<Vendor>>>(&vendors_)
      .set_long_name("--vendor")
      .set_help(
          "Vendors to be benchmarked with. (supported vendors: horizen, "
          "plonky3)");
}

bool PackedPoseidon2Config::Validate() const {
  for (const Vendor vendor : vendors_) {
    if ((vendor.value() != Vendor::kHorizen) &&
        (vendor.value() != Vendor::kPlonky3)) {
      tachyon_cerr << "Unsupported vendor " << vendor.ToString() << std::endl;
      return false;
    }
  }
  if (prime_field_.value() == FieldType::kBabyBear && vendors_.size() != 1) {
    tachyon_cerr << "Please run one vendor at a time for Baby Bear!"
                 << std::endl;
    return false;
  }
  if ((prime_field_.value() != FieldType::kBabyBear) &&
      // (prime_field_.value() != FieldType::kBn254Fr) &&
      (prime_field_.value() != FieldType::kPackedBabyBear)) {
    tachyon_cerr << "Unsupported prime field " << prime_field_.ToString()
                 << std::endl;
    return false;
  }
  return true;
}

}  // namespace tachyon::benchmark
