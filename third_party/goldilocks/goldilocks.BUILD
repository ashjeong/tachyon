load("@kroma_network_tachyon//bazel:tachyon.bzl", "if_has_avx512")
load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "base_field",
    srcs = [
        "src/goldilocks_base_field.cpp",
    ],
    hdrs = [
        "src/goldilocks_base_field.hpp",
        "src/goldilocks_base_field_avx.hpp",
        "src/goldilocks_base_field_avx512.hpp",
        "src/goldilocks_base_field_batch.hpp",
        "src/goldilocks_base_field_scalar.hpp",
        "src/goldilocks_base_field_tools.hpp",
    ],
    copts = if_has_avx512(
        ["-mavx512f"],
        ["-mavx2"],
    ),
    defines = if_has_avx512(["__AVX512__"]),
    include_prefix = "third_party/goldilocks/include",
    includes = ["src"],
    strip_include_prefix = "src",
    deps = ["@local_config_gmp//:gmp"],
)
