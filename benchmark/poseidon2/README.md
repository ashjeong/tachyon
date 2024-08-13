# Poseidon2 Hash Benchmark

```shell
Run on 13th Gen Intel(R) Core(TM) i9-13900K (32 X 5500 MHz CPU s)
CPU Caches:
  L1 Data 48 KiB (x16)
  L1 Instruction 32 KiB (x16)
  L2 Unified 2048 KiB (x16)
  L3 Unified 36864 KiB (x1)

Run on Apple M3 Pro (12 X 4050 MHz)
CPU Caches:
  L1 Data 64 KiB (x12)
  L1 Instruction 128 KiB (x12)
  L2 Unified 4096 KiB (x12)
```

Note that Poseidon2 runs 100x per test due to some time results being too small when running a single iteration.

## BN254

```shell
bazel run -c opt --//:has_openmp --//:has_rtti --//:has_matplotlib //benchmark/poseidon2:poseidon2_benchmark -- -p bn254_fr --vendor horizen --vendor plonky3 --check_results
```

### On Intel i9-13900K

| Trial Number | Tachyon   | Horizen       | Plonky3   |
| :----------- | --------- | ------------- | --------- |
| 0            | 0.00064   | **0.000538**  | 0.000877  |
| 1            | 0.000641  | **0.000527**  | 0.000876  |
| 2            | 0.000625  | **0.000525**  | 0.00087   |
| 3            | 0.000626  | **0.000531**  | 0.000867  |
| 4            | 0.000621  | **0.000527**  | 0.000867  |
| 5            | 0.000617  | **0.000526**  | 0.00087   |
| 6            | 0.000617  | **0.000525**  | 0.00086   |
| 7            | 0.000623  | **0.000525**  | 0.000853  |
| 8            | 0.000618  | **0.000529**  | 0.000853  |
| 9            | 0.000617  | **0.000522**  | 0.000852  |
| avg          | 0.0006245 | **0.0005275** | 0.0008645 |

![image](/benchmark/poseidon2/poseidon2_benchmark_bn254_ubuntu_i9.png)

### On Mac M3 Pro

<!-- TO UPDATE -->
<!-- | Trial Number | Tachyon | Horizen     | Plonky3  |
| :----------- | ------- | ----------- | -------- |
| 0            | 1.3e-05 | **1.2e-05** | 1.5e-05  |
| 1            | 1e-05   | **8e-06**   | 1.1e-05  |
| 2            | 9e-06   | **7e-06**   | 1e-05    |
| 3            | 9e-06   | **7e-06**   | 1e-05    |
| 4            | 9e-06   | **7e-06**   | 1e-05    |
| 5            | 9e-06   | **7e-06**   | 1e-05    |
| 6            | 9e-06   | **7e-06**   | 1e-05    |
| 7            | 9e-06   | **7e-06**   | 1e-05    |
| 8            | 9e-06   | **7e-06**   | 1e-05    |
| 9            | 9e-06   | **7e-06**   | 1e-05    |
| avg          | 9.5e-06 | **7.6e-06** | 1.06e-05 | -->

![image](/benchmark/poseidon2/poseidon2_benchmark_bn254_mac_m3.png)
