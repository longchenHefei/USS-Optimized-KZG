//! KZG-based homomorphic vector commitment (KHVC) with FK20 all-openings
//! and precomputed vanishing-polynomial openings for KeyUp.

use crate::fft::{domain, fk20_all_proofs, g1_fft, nth_root};
use ark_bls12_381::{Bls12_381, Fr, G1Affine, G1Projective, G2Affine, G2Projective};
use ark_ec::pairing::Pairing;
use ark_ec::{CurveGroup, PrimeGroup, VariableBaseMSM};
use ark_ff::{Field, PrimeField};
use ark_poly::{EvaluationDomain, Polynomial};
use ark_serialize::CanonicalSerialize;
use ark_std::rand::RngCore;
use ark_std::UniformRand;
use sha2::{Digest, Sha256};

pub fn hash_to_fr(bytes: &[u8]) -> Fr {
    let h = Sha256::digest(bytes);
    Fr::from_le_bytes_mod_order(&h)
}

#[derive(Clone)]
pub struct Khvc {
    pub n: usize,
    pub srs_g1: Vec<G1Affine>,      // [τ^0], …, [τ^n]
    pub srs_g1_proj: Vec<G1Projective>,
    pub g2: G2Affine,
    pub tau_g2: G2Affine,
    pub z_g1: G1Projective,         // [τ^n - 1]_1
    pub omega: Fr,
    pub omegas: Vec<Fr>,
    /// [Z(τ)/(τ − ω^j)]_1 for KeyUp / hiding term of all-open.
    pub z_openings: Vec<G1Projective>,
}

#[derive(Clone, Copy, Debug)]
pub struct Commitment(pub G1Projective);

#[derive(Clone, Copy, Debug)]
pub struct Opening(pub G1Projective);

impl Khvc {
    pub fn setup<R: RngCore>(n: usize, rng: &mut R) -> Self {
        assert!(n.is_power_of_two() && n >= 2);
        let tau = Fr::rand(rng);
        let g1 = G1Projective::generator();
        let g2 = G2Projective::generator();
        let mut srs_g1_proj = Vec::with_capacity(n + 1);
        let mut acc = Fr::ONE;
        for _ in 0..=n {
            srs_g1_proj.push(g1 * acc);
            acc *= tau;
        }
        let srs_g1: Vec<G1Affine> = srs_g1_proj.iter().map(|p| p.into_affine()).collect();
        let tau_g2 = (g2 * tau).into_affine();
        let z_g1 = srs_g1_proj[n] - srs_g1_proj[0];
        let omega = nth_root(n);
        let mut omegas = Vec::with_capacity(n);
        let mut w = Fr::ONE;
        for _ in 0..n {
            omegas.push(w);
            w *= omega;
        }
        let z_openings = compute_z_openings(n, &srs_g1_proj, omega);
        Self {
            n,
            srs_g1,
            srs_g1_proj,
            g2: g2.into_affine(),
            tau_g2,
            z_g1,
            omega,
            omegas,
            z_openings,
        }
    }

    pub fn crs_bytes(&self) -> usize {
        // G1 affine compressed 48 bytes, G2 96 bytes
        self.srs_g1.len() * 48 + 96 * 2 + 48
    }

    /// Interpolate evaluations on Ω, hide with r·Z, commit.
    pub fn commit(&self, values: &[Fr], r: Fr) -> Commitment {
        assert_eq!(values.len(), self.n);
        let coeffs = self.coeffs_of(values, r);
        Commitment(self.commit_coeffs(&coeffs))
    }

    pub fn commit_coeffs(&self, coeffs: &[Fr]) -> G1Projective {
        let bases = &self.srs_g1[..coeffs.len()];
        G1Projective::msm(bases, coeffs).unwrap()
    }

    /// `F = f + r (X^n − 1)` with `f` the IFFT of `values`.
    pub fn coeffs_of(&self, values: &[Fr], r: Fr) -> Vec<Fr> {
        let mut evals = values.to_vec();
        domain(self.n).ifft_in_place(&mut evals);
        // F coeffs: (f_0 - r, f_1, …, f_{n-1}, r)
        evals[0] -= r;
        evals.push(r);
        evals
    }

    pub fn open(&self, i: usize, values: &[Fr], r: Fr) -> Opening {
        let coeffs = self.coeffs_of(values, r);
        Opening(self.open_coeffs(i, &coeffs, values[i]))
    }

    pub fn open_coeffs(&self, i: usize, coeffs: &[Fr], y: Fr) -> G1Projective {
        let w = self.omegas[i];
        let q = synthetic_division(coeffs, w, y);
        let bases = &self.srs_g1[..q.len()];
        G1Projective::msm(bases, &q).unwrap()
    }

    /// All openings via FK20 on the degree `< n` part plus `r · z_openings`.
    /// For `n < 1024` the naive path is faster (FFT overhead); see Table 1.
    pub fn open_all(&self, values: &[Fr], r: Fr) -> Vec<Opening> {
        if self.n < 1024 {
            return self.open_all_naive(values, r);
        }
        self.open_all_from_coeffs(&self.coeffs_of(values, r), r)
    }

    pub fn open_all_from_coeffs(&self, coeffs: &[Fr], r: Fr) -> Vec<Opening> {
        // coeffs has length n+1: [f0-r, f1, …, f_{n-1}, r]
        let low = coeffs[..self.n].to_vec();
        // The degree < n polynomial is f − r, already in `low` (low[0] = f0 - r).
        let proofs_low = fk20_all_proofs(&low, &self.srs_g1_proj[..self.n]);
        proofs_low
            .into_iter()
            .zip(self.z_openings.iter())
            .map(|(p, z)| Opening(p + *z * r))
            .collect()
    }

    /// Naive O(n²) all-open, used to check FK20.
    pub fn open_all_naive(&self, values: &[Fr], r: Fr) -> Vec<Opening> {
        let coeffs = self.coeffs_of(values, r);
        (0..self.n)
            .map(|i| Opening(self.open_coeffs(i, &coeffs, values[i])))
            .collect()
    }

    pub fn verify(&self, c: Commitment, value: Fr, i: usize, opening: Opening) -> bool {
        let g1 = G1Projective::generator();
        let lhs = Bls12_381::pairing(c.0 - g1 * value, G2Projective::from(self.g2));
        let w = self.omegas[i];
        let rhs = Bls12_381::pairing(
            opening.0,
            G2Projective::from(self.tau_g2) - G2Projective::from(self.g2) * w,
        );
        lhs == rhs
    }

    pub fn com_hom(a: Commitment, b: Commitment) -> Commitment {
        Commitment(a.0 + b.0)
    }

    pub fn open_hom(a: Opening, b: Opening) -> Opening {
        Opening(a.0 + b.0)
    }

    pub fn rerand<R: RngCore>(&self, c: Commitment, rng: &mut R) -> (Commitment, Fr) {
        let r = Fr::rand(rng);
        (Commitment(c.0 + self.z_g1 * r), r)
    }

    pub fn serialized_g1_size(p: G1Projective) -> usize {
        let mut v = Vec::new();
        p.into_affine().serialize_compressed(&mut v).unwrap();
        v.len()
    }
}

/// Synthetic division of `f(X) − y` by `X − w`. Remainder must be 0.
pub fn synthetic_division(coeffs: &[Fr], w: Fr, y: Fr) -> Vec<Fr> {
    let n = coeffs.len();
    assert!(n >= 2);
    let mut a = coeffs.to_vec();
    a[0] -= y;
    // q has degree n-2, length n-1
    let mut q = vec![Fr::from(0u64); n - 1];
    q[n - 2] = a[n - 1];
    for k in (0..n - 2).rev() {
        q[k] = a[k + 1] + q[k + 1] * w;
    }
    q
}

/// [ (τ^n − 1)/(τ − ω^j) ]_1 = Σ_k ω^{j(n-1-k)} [τ^k].
fn compute_z_openings(n: usize, srs: &[G1Projective], omega: Fr) -> Vec<G1Projective> {
    // π_j = Σ_{k=0}^{n-1} srs[k] * ω^{j(n-1-k)} = ω^{-j} * Σ_k srs[k] (ω^{-j})^k
    // Σ_k srs[k] ω^{-jk} is IFFT-like: fft with ω^{-1}.
    let omega_inv = omega.inverse().unwrap();
    let evals = g1_fft(&srs[..n], omega_inv); // evals[j] = Σ_k srs[k] (ω^{-1})^{jk} = Σ srs[k] ω^{-jk}
    let mut out = Vec::with_capacity(n);
    let mut winv = Fr::ONE;
    for j in 0..n {
        // multiply by ω^{-j}
        out.push(evals[j] * winv);
        winv *= omega_inv;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn commit_open_verify() {
        let mut rng = ark_std::test_rng();
        let n = 8;
        let khvc = Khvc::setup(n, &mut rng);
        let values: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let r = Fr::rand(&mut rng);
        let c = khvc.commit(&values, r);
        for i in 0..n {
            let pi = khvc.open(i, &values, r);
            assert!(khvc.verify(c, values[i], i, pi));
        }
        let all = khvc.open_all(&values, r);
        let naive = khvc.open_all_naive(&values, r);
        for i in 0..n {
            assert_eq!(all[i].0, naive[i].0, "FK20 mismatch at {i}");
            assert!(khvc.verify(c, values[i], i, all[i]));
        }
    }

    #[test]
    fn hiding_and_hom() {
        let mut rng = ark_std::test_rng();
        let n = 8;
        let khvc = Khvc::setup(n, &mut rng);
        let v1: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let v2: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let r1 = Fr::rand(&mut rng);
        let r2 = Fr::rand(&mut rng);
        let c1 = khvc.commit(&v1, r1);
        let c2 = khvc.commit(&v2, r2);
        let c = Khvc::com_hom(c1, c2);
        let sum: Vec<Fr> = v1.iter().zip(v2.iter()).map(|(a, b)| *a + *b).collect();
        let pi = khvc.open(3, &sum, r1 + r2);
        assert!(khvc.verify(c, sum[3], 3, pi));
        let (c3, r3) = khvc.rerand(c1, &mut rng);
        let pi3 = khvc.open(1, &v1, r1 + r3);
        assert!(khvc.verify(c3, v1[1], 1, pi3));
    }
}
