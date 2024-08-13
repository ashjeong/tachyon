use ff::PrimeField;
use p3_baby_bear::{BabyBear, DiffusionMatrixBabyBear};
use p3_bn254_fr::{Bn254Fr, DiffusionMatrixBN254, FFBn254Fr};
use p3_field::AbstractField;
use p3_poseidon2::{Poseidon2, Poseidon2ExternalMatrixHL};
use p3_symmetric::Permutation;
use std::time::Instant;
use tachyon_rs::math::{
    elliptic_curves::bn::bn254::Fr as CppBn254Fr, finite_fields::baby_bear::BabyBear as CppBabyBear,
};
use zkhash::ark_ff::{BigInteger, PrimeField as ark_PrimeField};
use zkhash::fields::{babybear::FpBabyBear as ark_FpBabyBear, bn256::FpBN256 as ark_FpBN256};
use zkhash::poseidon2::{
    poseidon2_instance_babybear::RC16 as BabyBearRC16, poseidon2_instance_bn256::RC3 as BN256RC3,
};

fn bn254_from_ark_ff(input: ark_FpBN256) -> Bn254Fr {
    let bytes = input.into_bigint().to_bytes_le();

    let mut res = <FFBn254Fr as PrimeField>::Repr::default();

    for (i, digit) in res.0.as_mut().iter_mut().enumerate() {
        *digit = bytes[i];
    }

    let value = FFBn254Fr::from_repr(res);

    if value.is_some().into() {
        Bn254Fr {
            value: value.unwrap(),
        }
    } else {
        panic!("Invalid field element")
    }
}

fn baby_bear_from_ark_ff(input: ark_FpBabyBear) -> BabyBear {
    BabyBear::from_canonical_u32(input.into_bigint().0[0] as u32)
}

// #[no_mangle]
// pub extern "C" fn run_poseidon_plonky3_baby_bear(duration: *mut u64, width, d, rounds_f, rounds_p) -> *mut CppBabyBear {
//     const width: usize = 16;
//     const d: u64 = 7;
//     const rounds_f: usize = 8;
//     const rounds_p: usize = 13;

//     // Copy over round constants from zkhash.
//     let mut round_constants: Vec<[BabyBear; width]> = BabyBearRC16
//         .iter()
//         .map(|vec| {
//             vec.iter()
//                 .cloned()
//                 .map(baby_bear_from_ark_ff)
//                 .collect::<Vec<_>>()
//                 .try_into()
//                 .unwrap()
//         })
//         .collect();

//     let internal_start = rounds_f / 2;
//     let internal_end = (rounds_f / 2) + rounds_p;
//     let internal_round_constants = round_constants
//         .drain(internal_start..internal_end)
//         .map(|vec| vec[0])
//         .collect::<Vec<_>>();
//     let external_round_constants = round_constants;

//     let poseidon =
//         Poseidon2::<BabyBear, Poseidon2ExternalMatrixHL, DiffusionMatrixBabyBear, width, d>::new(
//             rounds_f,
//             external_round_constants,
//             Poseidon2ExternalMatrixHL,
//             rounds_p,
//             internal_round_constants,
//             DiffusionMatrixBabyBear,
//         );

//     let mut input = (0..width)
//         .map(|_i| BabyBear::zero())
//         .collect::<Vec<_>>()
//         .try_into()
//         .unwrap();

//     let start = Instant::now();
//     for _ in 0..100 {
//         poseidon.permute_mut(&mut input);
//     }
//     unsafe {
//         duration.write(start.elapsed().as_micros() as u64);
//     }

//     Box::into_raw(Box::new(input[1])) as *mut CppBabyBear
// }

// #[no_mangle]
// pub extern "C" fn run_poseidon_plonky3_bn254_fr(duration: *mut u64) -> *mut CppBn254Fr {
//     const width: usize = 3;
//     const d: u64 = 5;
//     const rounds_f: usize = 8;
//     const rounds_p: usize = 56;

//     // Copy over round constants from zkhash.
//     let mut round_constants: Vec<[Bn254Fr; width]> = BN256RC3
//         .iter()
//         .map(|vec| {
//             vec.iter()
//                 .cloned()
//                 .map(bn254_from_ark_ff)
//                 .collect::<Vec<_>>()
//                 .try_into()
//                 .unwrap()
//         })
//         .collect();
//     let internal_start = rounds_f / 2;
//     let internal_end = (rounds_f / 2) + rounds_p;
//     let internal_round_constants = round_constants
//         .drain(internal_start..internal_end)
//         .map(|vec| vec[0])
//         .collect::<Vec<_>>();
//     let external_round_constants = round_constants;

//     let poseidon =
//         Poseidon2::<Bn254Fr, Poseidon2ExternalMatrixHL, DiffusionMatrixBN254, width, d>::new(
//             rounds_f,
//             external_round_constants,
//             Poseidon2ExternalMatrixHL,
//             rounds_p,
//             internal_round_constants,
//             DiffusionMatrixBN254,
//         );

//     let mut input = (0..width)
//         .map(|_i| Bn254Fr::zero())
//         .collect::<Vec<_>>()
//         .try_into()
//         .unwrap();

//     let start = Instant::now();
//     for _ in 0..100 {
//         poseidon.permute_mut(&mut input);
//     }
//     unsafe {
//         duration.write(start.elapsed().as_micros() as u64);
//     }

//     Box::into_raw(Box::new(input[1])) as *mut CppBn254Fr
// }


fn run_poseidon2<F: PrimeField + std::convert::From<i32>, R>(
  duration: *mut u64,
  RC: &BabyBearRC16,
  width: usize,
  d: u64,
  rounds_f : usize,
  rounds_p:usize

) -> *mut R {
  // Copy over round constants from zkhash.
  let mut round_constants: Vec<[F; width]> = RC
  .iter()
  .map(|vec| {
      vec.iter()
          .cloned()
          .map(baby_bear_from_ark_ff)
          .collect::<Vec<_>>()
          .try_into()
          .unwrap()
  })
  .collect();

let internal_start = rounds_f / 2;
let internal_end = (rounds_f / 2) + rounds_p;
let internal_round_constants = round_constants
  .drain(internal_start..internal_end)
  .map(|vec| vec[0])
  .collect::<Vec<_>>();
let external_round_constants = round_constants;

let poseidon =
  Poseidon2::<F, Poseidon2ExternalMatrixHL, DiffusionMatrixBabyBear, width, d>::new(
      rounds_f,
      external_round_constants,
      Poseidon2ExternalMatrixHL,
      rounds_p,
      internal_round_constants,
      DiffusionMatrixBabyBear,
  );

let mut input = (0..width)
  .map(|_i| F::zero())
  .collect::<Vec<_>>()
  .try_into()
  .unwrap();

let start = Instant::now();
for _ in 0..100 {
  poseidon.permute_mut(&mut input);
}
unsafe {
  duration.write(start.elapsed().as_micros() as u64);
}

Box::into_raw(Box::new(input[1])) as *mut R
}

#[no_mangle]
pub extern "C" fn run_poseidon_plonky3_baby_bear(duration: *mut u64) -> *mut CppBabyBear {
  run_poseidon2::<_, CppBabyBear>(duration, &BabyBearRC16, 16, 7, 8, 13)
}

#[no_mangle]
pub extern "C" fn run_poseidon_plonky3_bn254_fr(duration: *mut u64) -> *mut CppBn254Fr {
  run_poseidon2::<_, CppBn254Fr>(duration, &BN256RC3, 3, 5, 8, 56)
}
