# Poseidon2 Hash Benchmark

```
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

## Baby Bear

Note: Plonky3 and Horizen compute values with a different internal matrix, requiring them to be compared with Tachyon separately.

### Plonky3

```shell
bazel run -c opt --//:has_openmp --//:has_rtti --//:has_matplotlib //benchmark/poseidon2:poseidon2_benchmark -- -p baby_bear --vendor plonky3 --check_results
```

### Horizen

```shell
bazel run -c opt --//:has_openmp --//:has_rtti --//:has_matplotlib //benchmark/poseidon2:poseidon2_benchmark -- -p baby_bear --vendor horizen --check_results
```

### On Intel i9-13900K - Plonky3

| Trial Number | Tachyon   | Plonky3      |
| :----------- | --------- | ------------ |
| 0            | 0.000112  | **6.6e-05**  |
| 1            | 0.000111  | **6.5e-05**  |
| 2            | 0.000111  | **6.6e-05**  |
| 3            | 0.000111  | **6.6e-05**  |
| 4            | 0.00011   | **6.6e-05**  |
| 5            | 0.000116  | **6.6e-05**  |
| 6            | 0.00011   | **6.5e-05**  |
| 7            | 0.000109  | **6.6e-05**  |
| 8            | 0.00011   | **6.6e-05**  |
| 9            | 0.000109  | **6.5e-05**  |
| avg          | 0.0001109 | **6.57e-05** |

![image](/benchmark/poseidon2/poseidon2_benchmark_baby_bear_plonky3_ubuntu_i9.png)

### On Intel i9-13900K - Horizen

| Trial Number | Tachyon       | Horizen   |
| :----------- | ------------- | --------- |
| 0            | **0.000127**  | 0.000381  |
| 1            | **0.000126**  | 0.00036   |
| 2            | **0.000125**  | 0.00037   |
| 3            | **0.000125**  | 0.000356  |
| 4            | **0.000125**  | 0.000354  |
| 5            | **0.000125**  | 0.000354  |
| 6            | **0.000125**  | 0.000354  |
| 7            | **0.000125**  | 0.00036   |
| 8            | **0.000125**  | 0.000359  |
| 9            | **0.000125**  | 0.000353  |
| avg          | **0.0001253** | 0.0003601 |

![image](/benchmark/poseidon2/poseidon2_benchmark_baby_bear_horizen_ubuntu_i9.png)

#### On Mac M3 Pro - Plonky3

<!-- TO UPDATE -->
<!-- | Repetition | Tachyon   | Plonky3       |
| :--------- | --------- | ------------- |
| 0          | 0.000147  | **0.000137**  |
| 1          | 0.00019   | **0.000136**  |
| 2          | 0.000199  | **0.000137**  |
| 3          | 0.000197  | **0.000137**  |
| 4          | 0.000199  | **0.00014**   |
| 5          | 0.000198  | **0.000136**  |
| 6          | 0.000198  | **0.000137**  |
| 7          | 0.000199  | **0.000137**  |
| 8          | 0.000208  | **0.000137**  |
| 9          | 0.00021   | **0.000137**  |
| avg        | 0.0001945 | **0.0001371** | -->

![image](/benchmark/poseidon2/poseidon2_benchmark_baby_bear_plonky3_mac_m2.png)

#### On Mac M3 Pro - Horizen

<!-- TO UPDATE -->
<!-- | Repetition | Tachyon      | Horizen      |
| :--------- | ------------ | ------------ |
| 0          | **0.000166** | 0.000193     |
| 1          | **0.000164** | 0.000172     |
| 2          | **0.000164** | 0.000174     |
| 3          | **0.000165** | 0.000173     |
| 4          | 0.000218     | **0.000173** |
| 5          | 0.000177     | **0.000174** |
| 6          | 0.000181     | **0.000173** |
| 7          | **0.000166** | 0.000184     |
| 8          | **0.000169** | 0.000174     |
| 9          | **0.00016**  | 0.0002       |
| avg        | **0.000173** | 0.000179     | -->

![image](/benchmark/poseidon2/poseidon2_benchmark_baby_bear_horizen_mac_m2.png)
