use ark_bls12_381::{Fr as ScalarField, G1Projective as G1};
use ark_ff::{BigInteger, Field, PrimeField, UniformRand};
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::{ EvaluationDomain, Radix2EvaluationDomain};
use ark_std:: test_rng;
use ark_ff::Zero;
use ark_std::ops::Mul;
use sha2::{Digest, Sha256};



#[allow(non_upper_case_globals)]
const n : usize = 2_usize.pow(17);
fn main() {

    //setup 
    let rng = &mut test_rng();
    let domain = Radix2EvaluationDomain::<ScalarField>::new(2*n).unwrap();
    let tau = ScalarField::rand(rng);
    let tau_powers: Vec<ScalarField> = (0..2*n).map(|i| tau.pow([i as u64])).collect();


    let g = G1::rand(rng);
    let srs_monomial: Vec<G1> = tau_powers.iter().map(|&tau_i| g.mul(tau_i)).collect();
    let srs_lagrange = domain.ifft(&srs_monomial);

    let poly_coeffs: Vec<ScalarField> = (0..2*n).map(|_| ScalarField::rand(rng)).collect();

    let poly: DensePolynomial<ScalarField> = DensePolynomial::<ScalarField>::from_coefficients_vec(poly_coeffs);
    let poly_lagrange = domain.fft(&poly.coeffs);

    // Prover gets SRS in Lagrange basis and polynomial in Lagrange basis
    println!("SRS in Lagrange basis: {:?}", srs_lagrange);
    println!("Polynomial in Lagrange basis: {:?}", poly_lagrange);

     // Witness: Random field elements for the polynomial
     let mut  witness: Vec<ScalarField> = (0..n).map(|_|{ let el = ScalarField::rand(rng).into_bigint().to_bytes_le();
        let mut hasher = Sha256::new();
        hasher.update(&el);
        ScalarField::from_le_bytes_mod_order(hasher.finalize().as_slice())
    }).collect();
witness.extend_from_slice(&vec![ScalarField::zero(); n]);

      //  Convert witness to Lagrange basis using FFT
    let witness_lagrange = domain.fft(&witness);
    // order 1 * 2n 
    // Compute the commitment using Hadamard product and multiscalar multiplication
    let hadamard_product: Vec<ScalarField> = witness_lagrange.iter().zip(poly_lagrange.iter())
        .map(|(&w, &p)| w * p)
        .collect();
    // order 2n * 1 * 1 *2n
    let commitment = srs_lagrange.iter().zip(hadamard_product.iter())
        .map(|(&srs_i, &h_i)| srs_i.mul(h_i))
        .fold(G1::zero(), |acc, x| acc + x);

    println!("Commitment: {:?}", commitment);
    

}
