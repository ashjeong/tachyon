use halo2_proofs::{
    arithmetic::Field,
    circuit::{Layouter, SimpleFloorPlanner, Value},
    plonk::{Circuit, ConstraintSystem, Error},
};

mod my_chip;

use crate::my_circuit::my_chip::NumericInstructions;
use my_chip::{FieldConfig, MyChip};

/// The full circuit implementation.
///
/// In this struct we store the private input variables. We use `Option<F>` because
/// they won't have any value during key generation. During proving, if any of these
/// were `None` we would get an error.
#[derive(Default)]
pub struct MyCircuit<F: Field> {
    pub constant: F,
    pub a: Value<F>,
    pub b: Value<F>,
}

impl<F: Field> Circuit<F> for MyCircuit<F> {
    // Since we are using a single chip for everything, we can just reuse its config.
    type Config = FieldConfig;
    type FloorPlanner = SimpleFloorPlanner;

    fn without_witnesses(&self) -> Self {
        Self::default()
    }

    fn configure(meta: &mut ConstraintSystem<F>) -> Self::Config {
        // We create the two advice columns that FieldChip uses for I/O.
        let advice = [meta.advice_column(), meta.advice_column()];

        // We also need an instance column to store public inputs.
        let instance = meta.instance_column();

        // Create a fixed column to load constants.
        let constant = meta.fixed_column();

        MyChip::configure(meta, advice, instance, constant)
    }

    fn synthesize(
        &self,
        config: Self::Config,
        mut layouter: impl Layouter<F>,
    ) -> Result<(), Error> {
        let field_chip = MyChip::<F>::construct(config);

        // Load our private values into the circuit.
        let a = field_chip.load_private(layouter.namespace(|| "load a"), self.a)?;
        let b = field_chip.load_private(layouter.namespace(|| "load b"), self.b)?;

        // Load the constant factor into the circuit.
        let constant =
            field_chip.load_constant(layouter.namespace(|| "load constant"), self.constant)?;

        // We only have access to plain multiplication.
        // We could implement our circuit as:
        //     asq  = a*a
        //     bsq  = b*b
        //     absq = asq*bsq
        //     c    = constant*asq*bsq
        //
        // but it's more efficient to implement it as:
        //     ab   = a*b
        //     absq = ab^2
        //     c    = constant*absq
        // let aa = field_chip.mul(layouter.namespace(|| "a * a"), a.clone(), a.clone())?;
        // let aaa = field_chip.mul(layouter.namespace(|| "aa * a"), aa, a)?;
        // let aaa_b = field_chip.add(layouter.namespace(|| "aaa + b"), aaa, b)?;
        // // let aaa_b_c = field_chip.add(layouter.namespace(|| "aaa_b + 5"), aaa_b, constant)?;
        let ab = field_chip.mul(layouter.namespace(|| "a * b"), a, b)?;
        let ab_c = field_chip.add(layouter.namespace(|| "ab + constant"), ab, constant)?;

        // Expose the result as a public input to the circuit.
        field_chip.expose_public(layouter.namespace(|| "expose ab_c"), ab_c, 0)
    }
}
