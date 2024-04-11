mod my_circuit;

use halo2_proofs::{circuit::Value, dev::MockProver, pasta::Fp};
use my_circuit::MyCircuit;

fn main() {
    // x^3 + x + 5 = 35

    // max number of rows (2^k)
    let k = 4;

    // private inputs:
    // secret value
    let a = Fp::from(2);
    let b = Fp::from(3);
    let constant = Fp::from(5);

    // public input:
    let res = Fp::from(11);

    // create circuit with private inputs:
    let circuit = MyCircuit {
        constant,
        a: Value::known(a),
        b: Value::known(b),
    };

    let public_inputs = vec![res];
    let prover = MockProver::run(k, &circuit, vec![public_inputs]).unwrap();
    print!("advice cells: {:#?}", prover.advice);
    assert_eq!(prover.verify(), Ok(()));
}
