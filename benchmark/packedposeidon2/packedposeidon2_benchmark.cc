#include <iostream>

// clang-format off
#include "benchmark/simple_reporter.h"
#include "benchmark/packedposeidon2/packedposeidon2_benchmark_runner.h"
#include "benchmark/packedposeidon2/packedposeidon2_config.h"
// clang-format on
#include "tachyon/base/containers/contains.h"
#include "tachyon/base/logging.h"
#include "tachyon/base/profiler.h"
#include "tachyon/c/math/finite_fields/baby_bear/baby_bear.h"
#include "tachyon/c/math/finite_fields/baby_bear/baby_bear_type_traits.h"
#include "tachyon/crypto/hashes/sponge/poseidon2/poseidon2_config.h"
#include "tachyon/crypto/hashes/sponge/poseidon2/poseidon2_external_matrix.h"
#include "tachyon/crypto/hashes/sponge/poseidon2/poseidon2_external_matrix_traits_forward.h"
#include "tachyon/crypto/hashes/sponge/poseidon2/poseidon2_horizen_external_matrix.h"
#include "tachyon/crypto/hashes/sponge/poseidon2/poseidon2_params.h"
#include "tachyon/math/finite_fields/baby_bear/baby_bear.h"
#include "tachyon/math/finite_fields/baby_bear/internal/packed_baby_bear.h"
#include "tachyon/math/finite_fields/baby_bear/poseidon2.h"

namespace tachyon::benchmark {

using namespace crypto;

extern "C" tachyon_baby_bear* run_poseidon2_plonky3_baby_bear(
    uint64_t* duration);

template <typename PackedF, typename Fn>
void Run(SimpleReporter& reporter, const PackedPoseidon2Config& config,
         Fn plonky3_fn) {
  // using Field = math::MaybeUnpack<PackedF>;
  PackedF::Init();

  Poseidon2BenchmarkRunner<PackedF> runner(reporter, config);

  PackedF result;
  if constexpr (std::is_same_v<PackedF, math::BabyBear>) {
    if (base::Contains(config.vendors(), Vendor::Plonky3())) {
      using Params = Poseidon2Params<math::BabyBear, 15, 7>;

      crypto::Poseidon2Config<Params> poseidon2_config =
          crypto::Poseidon2Config<Params>::CreateCustom(
              math::GetPoseidon2BabyBearInternalShiftArray<15>());
      result = runner.template Run<Params>(poseidon2_config);

    } else {
      using Params = Poseidon2Params<math::BabyBear, 15, 7>;
      crypto::Poseidon2Config<Params> poseidon2_config =
          crypto::Poseidon2Config<Params>::CreateCustom(
              math::GetPoseidon2BabyBearInternalDiagonalArray<16>());
      result = runner.template Run<Params>(poseidon2_config);
    }
  } else {
    using PackedParams = Poseidon2Params<math::PackedBabyBear, 15, 7>;
    crypto::Poseidon2Config<PackedParams> poseidon2_config =
        crypto::Poseidon2Config<PackedParams>::CreateCustom(
            math::GetPoseidon2BabyBearInternalShiftArray<15>());
    result = runner.template Run<PackedParams>(poseidon2_config);
  }

  // for (const Vendor vendor : config.vendors()) {
  //   PackedF result_vendor;
  //   result_vendor = runner.RunExternal(vendor, plonky3_fn);

  // if (config.check_results()) {
  //   if constexpr (Field::Config::kModulusBits < 32 &&
  //                 vendor.value() == Vendor::kPlonky3) {
  //     CHECK_EQ(result, result_vendor)
  //         << "Tachyon and Plonky3 results do not match";
  //   } else {
  //     CHECK_EQ(result, result_vendor) << "Results do not match";
  //   }
  // }
  // }
}

int RealMain(int argc, char** argv) {
  base::FilePath tmp_file;
  CHECK(base::GetTempDir(&tmp_file));
  tmp_file = tmp_file.Append("poseidon2_benchmark.perfetto-trace");
  base::Profiler profiler({tmp_file});

  profiler.Init();
  profiler.Start();

  PackedPoseidon2Config config;
  if (!config.Parse(argc, argv)) {
    return 1;
  }

  SimpleReporter reporter;
  reporter.set_title("Poseidon2 Benchmark");
  reporter.set_x_label("Trial number");
  reporter.set_column_labels(
      base::CreateVector(config.repeating_num(),
                         [](size_t i) { return base::NumberToString(i); }));

  if (config.prime_field().value() == FieldType::kBabyBear) {
    Run<math::BabyBear>(reporter, config, run_poseidon2_plonky3_baby_bear);
  } else if (config.prime_field().value() == FieldType::kPackedBabyBear) {
    Run<math::PackedBabyBear>(reporter, config,
                              run_poseidon2_plonky3_baby_bear);
  } else {
    NOTREACHED();
  }

  reporter.AddAverageAsLastColumn();
  reporter.Show();

  return 0;
}

}  // namespace tachyon::benchmark

int main(int argc, char** argv) {
  return tachyon::benchmark::RealMain(argc, argv);
}
