# Tracking Testing Types Across Example Circuits

## Currently Tested Types

| Circuit             | FloorPlanner       | PCS     | LS      | Fully Tested?      |
| ------------------- | ------------------ | ------- | ------- | ------------------ |
| SimpleCircuit       | SimpleFloorPlanner | SHPlonk | Halo2LS | :heavy_check_mark: |
| SimpleCircuit       | V1FloorPlanner     | SHPlonk | Halo2LS | :heavy_check_mark: |
| SimpleLookupCircuit | SimpleFloorPlanner | SHPlonk | Halo2LS | :heavy_check_mark: |
| SimpleLookupCircuit | V1FloorPlanner     | SHPlonk | Halo2LS | :heavy_check_mark: |
| ShuffleCircuit      | SimpleFloorPlanner | SHPlonk | Halo2LS | :heavy_check_mark: |
| ShuffleCircuit      | V1FloorPlanner     | SHPlonk | Halo2LS | :heavy_check_mark: |
| ShuffleCircuit      | SimpleFloorPlanner | GWC     | Halo2LS | :x:                |
| ShuffleCircuit      | V1FloorPlanner     | GWC     | Halo2LS | :x:                |

## Identical Test Data

The following circuits result in the same test data when using `SimpleFloorPlanner` or `V1FloorPlanner`

- SimpleLookupCircuit
- ShuffleCircuit
