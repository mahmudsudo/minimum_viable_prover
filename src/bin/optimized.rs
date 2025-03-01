use ark_bls12_381::{Fr as ScalarField, G1Projective as G1};
use ark_ec::{CurveGroup, VariableBaseMSM};
use ark_ff::{BigInteger, PrimeField, UniformRand, Zero, One};
use ark_poly::{ EvaluationDomain, Radix2EvaluationDomain};
use ark_std::test_rng;
use ark_std::ops::Mul;
use rayon::prelude::*;
use sha2::{Digest, Sha256};

const N: usize = 1 << 17; // 2^17

fn main() {
    // Setup
    let rng = &mut test_rng();
    let domain = Radix2EvaluationDomain::<ScalarField>::new(2 * N).unwrap();

    //  Optimized tau powers calculation
    let tau = ScalarField::rand(rng);
    let mut tau_powers = vec![ScalarField::one(); 2 * N];
    for i in 1..2*N {
        tau_powers[i] = tau_powers[i-1] * tau;
    }

    //  Parallel SRS generation
    let g = G1::rand(rng);
    let srs_monomial: Vec<G1> = tau_powers.par_iter()
        .map(|&tau_i| g.mul(tau_i))
        .collect();

    //  In-place FFT optimization for SRS
    let mut srs_lagrange = srs_monomial.clone();
    domain.ifft_in_place(&mut srs_lagrange);

    //  Parallel polynomial generation and FFT
    let poly_coeffs: Vec<ScalarField> = (0..2*N)
        .into_par_iter()
        .map_init(test_rng, |rng, _| ScalarField::rand(rng))
        .collect();
    
    let mut poly_lagrange = poly_coeffs.clone();
    domain.fft_in_place(&mut poly_lagrange);

    // Prover phase
    // 5. Parallel witness generation with batched hashing
    let witness: Vec<ScalarField> = (0..N)
        .into_par_iter()
        .map_init(
            || (test_rng(), Sha256::new()),
            |(rng_local, hasher), _| {
                let el = ScalarField::rand(rng_local).into_bigint().to_bytes_le();
                hasher.update(&el);
                ScalarField::from_le_bytes_mod_order(hasher.finalize_reset().as_slice())
            }
        )
        .collect();

    // 6. Optimized zero-padding with parallel extend
    let mut padded_witness = witness.clone();
    padded_witness.par_extend((0..N).into_par_iter().map(|_| ScalarField::zero()));

    // 7. In-place FFT for witness
    domain.fft_in_place(&mut padded_witness);

    // 8. Parallel Hadamard product
    let hadamard_product: Vec<ScalarField> = padded_witness.par_iter()
        .zip(poly_lagrange.par_iter())
        .map(|(&w, &p)| w * p)
        .collect();

    // 9. Optimized MSM using arkworks' implementation
    let srs_affine: Vec<_> = srs_lagrange.iter().map(|p| p.into_affine()).collect();
    let commitment = G1::msm(&srs_affine, &hadamard_product)
        .expect("MSM failed");

    println!("Commitment: {:?}", commitment);


    rayon::ThreadPoolBuilder::new()
    .num_threads(10)
    .build_global()
    .unwrap();
}