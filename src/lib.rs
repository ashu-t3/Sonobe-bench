#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(clippy::upper_case_acronyms)]

use ark_mnt4_298::{
    Fr as Fr4, MNT4_298, G1Projective as G1, 
    g1::{Config as Config4}, Fq as Fq4
};
use ark_mnt6_298::{
    Fr as Fr6, G1Projective as G2,
    g1::{Config as Config6}, Fq as Fq6
};
use ark_ff::PrimeField;
use ark_groth16::Groth16;
use ark_r1cs_std::{
    alloc::AllocVar,
    groups::curves::short_weierstrass::ProjectiveVar,
    fields::fp::FpVar,
    ToConstraintFieldGadget,
    prelude::CurveVar
};
use ark_relations::r1cs::{ConstraintSystemRef, SynthesisError};
use std::marker::PhantomData;
use std::time::Instant;

use folding_schemes::{
    commitment::{kzg::KZG, pedersen::Pedersen},
    folding::nova::{
        decider_eth::{prepare_calldata, Decider as DeciderEth},
        Nova, PreprocessorParam,
    },
    frontend::FCircuit,
    transcript::poseidon::poseidon_canonical_config,
    Decider, Error, FoldingScheme,
};

// Define constraint field variables for MNT4/MNT6
type FqVar4 = FpVar<Fq4>;
type FqVar6 = FpVar<Fq6>;

// Define the curve variable types using G1 configs
type GVar4 = ProjectiveVar<Config4, FqVar4>;
type GVar6 = ProjectiveVar<Config6, FqVar6>;

/// Test circuit to be folded
#[derive(Clone, Copy, Debug)]
pub struct CubicFCircuit<F: PrimeField> {
    _f: PhantomData<F>,
}

impl<F: PrimeField> FCircuit<F> for CubicFCircuit<F> {
    type Params = ();
    
    fn new(_params: Self::Params) -> Result<Self, Error> {
        Ok(Self { _f: PhantomData })
    }
    
    fn state_len(&self) -> usize {
        1
    }
    
    fn external_inputs_len(&self) -> usize {
        0
    }
    
    fn step_native(
        &self,
        _i: usize,
        z_i: Vec<F>,
        _external_inputs: Vec<F>,
    ) -> Result<Vec<F>, Error> {
        Ok(vec![z_i[0] * z_i[0] * z_i[0] + z_i[0] + F::from(5_u32)])
    }
    
    fn generate_step_constraints(
        &self,
        cs: ConstraintSystemRef<F>,
        _i: usize,
        z_i: Vec<FpVar<F>>,
        _external_inputs: Vec<FpVar<F>>,
    ) -> Result<Vec<FpVar<F>>, SynthesisError> {
        let five = FpVar::<F>::new_constant(cs.clone(), F::from(5u32))?;
        let z_i = z_i[0].clone();
        
        Ok(vec![&z_i * &z_i * &z_i + &z_i + &five])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_mnt() {
        let n_steps = 10;
        // set the initial state
        let z_0 = vec![Fr4::from(3_u32)];
        
        let f_circuit = CubicFCircuit::<Fr4>::new(()).unwrap();
        
        // Define Nova type using MNT4-298/MNT6-298 curve cycle
        pub type N = Nova<
            G1,
            GVar4,
            G2,
            GVar6,
            CubicFCircuit<Fr4>,
            KZG<'static, MNT4_298>,
            Pedersen<G2>,
            false
        >;
        
        // Define Decider type
        pub type D = DeciderEth<
            G1,
            GVar4,
            G2,
            GVar6,
            CubicFCircuit<Fr4>,
            KZG<'static, MNT4_298>,
            Pedersen<G2>,
            Groth16<MNT4_298>,
            N,
        >;
        
        let poseidon_config = poseidon_canonical_config::<Fr4>();
        let mut rng = rand::rngs::OsRng;
        
        // prepare the Nova prover & verifier params
        let nova_preprocess_params = PreprocessorParam::new(poseidon_config.clone(), f_circuit);
        let nova_params = N::preprocess(&mut rng, &nova_preprocess_params).unwrap();
        let pp_hash = nova_params.1.pp_hash().unwrap();
        
        // initialize the folding scheme engine
        let mut nova = N::init(&nova_params, f_circuit, z_0).unwrap();
        
        // prepare the Decider prover & verifier params
        let (decider_pp, decider_vp) = D::preprocess(&mut rng, nova_params, nova.clone()).unwrap();
        
        // run n steps of the folding iteration
        for i in 0..n_steps {
            let start = Instant::now();
            nova.prove_step(rng, vec![], None).unwrap();
            println!("Nova::prove_step {}: {:?}", i, start.elapsed());
        }
        
        let start = Instant::now();
        let proof = D::prove(rng, decider_pp, nova.clone()).unwrap();
        println!("generated Decider proof: {:?}", start.elapsed());
        
        let verified = D::verify(
            decider_vp.clone(),
            nova.i,
            nova.z_0.clone(),
            nova.z_i.clone(),
            &nova.U_i,
            &nova.u_i,
            &proof,
        )
        .unwrap();
        assert!(verified);
        println!("Decider proof verification: {}", verified);
    }
}