//! Scalar and G1 radix-2 FFTs aligned with `ark_poly::Radix2EvaluationDomain`.

use ark_bls12_381::{Fr, G1Projective};
use ark_ec::PrimeGroup;
use ark_ff::Field;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Polynomial, Radix2EvaluationDomain};
use ark_poly::univariate::DensePolynomial;

pub fn domain(n: usize) -> Radix2EvaluationDomain<Fr> {
    Radix2EvaluationDomain::<Fr>::new(n).expect("n must be a power of two fitting the field")
}

pub fn nth_root(n: usize) -> Fr {
    domain(n).group_gen
}

fn bitreverse(mut x: usize, log_n: usize) -> usize {
    let mut r = 0;
    for _ in 0..log_n {
        r = (r << 1) | (x & 1);
        x >>= 1;
    }
    r
}

fn bitrev_g1(a: &mut [G1Projective]) {
    let n = a.len();
    let log_n = n.trailing_zeros() as usize;
    for i in 0..n {
        let j = bitreverse(i, log_n);
        if i < j {
            a.swap(i, j);
        }
    }
}

/// In-place Cooley–Tukey FFT on G1. After the call, `a[i] = Σ_k a_k ω^{ik}`.
pub fn g1_fft_in_place(a: &mut [G1Projective], omega: Fr) {
    let n = a.len();
    assert!(n.is_power_of_two());
    bitrev_g1(a);
    let mut m = 1usize;
    while m < n {
        let w_m = omega.pow([(n / (2 * m)) as u64]);
        let mut k = 0;
        while k < n {
            let mut w = Fr::ONE;
            for j in 0..m {
                let t = a[k + j + m] * w;
                let u = a[k + j];
                a[k + j] = u + t;
                a[k + j + m] = u - t;
                w *= w_m;
            }
            k += 2 * m;
        }
        m *= 2;
    }
}

pub fn g1_ifft_in_place(a: &mut [G1Projective], omega: Fr) {
    let n = a.len();
    let omega_inv = omega.inverse().unwrap();
    g1_fft_in_place(a, omega_inv);
    let ninv = Fr::from(n as u64).inverse().unwrap();
    for p in a.iter_mut() {
        *p *= ninv;
    }
}

pub fn g1_fft(input: &[G1Projective], omega: Fr) -> Vec<G1Projective> {
    let mut a = input.to_vec();
    g1_fft_in_place(&mut a, omega);
    a
}

pub fn g1_ifft(input: &[G1Projective], omega: Fr) -> Vec<G1Projective> {
    let mut a = input.to_vec();
    g1_ifft_in_place(&mut a, omega);
    a
}

/// Multiply two length-`n` sequences via a 2n circular convolution in G1 (scalars × points).
pub fn g1_pointwise_mul(points: &[G1Projective], scalars: &[Fr]) -> Vec<G1Projective> {
    assert_eq!(points.len(), scalars.len());
    points.iter().zip(scalars.iter()).map(|(p, s)| *p * *s).collect()
}

pub fn scalar_fft(values: &mut [Fr], n: usize) {
    let d = domain(n);
    let mut v = values.to_vec();
    d.fft_in_place(&mut v);
    values.copy_from_slice(&v);
}

pub fn scalar_ifft(values: &mut [Fr], n: usize) {
    let d = domain(n);
    let mut v = values.to_vec();
    d.ifft_in_place(&mut v);
    values.copy_from_slice(&v);
}

/// FK20: all KZG opening proofs of a degree `< n` polynomial at the n-th roots of unity.
///
/// Port of ethereum/research `fk20_single` (Feist–Khovratovich). `srs` is `[1],[τ],…,[τ^{n-1}]`.
pub fn fk20_all_proofs(coeffs: &[Fr], srs: &[G1Projective]) -> Vec<G1Projective> {
    let n = coeffs.len();
    assert!(n.is_power_of_two());
    assert_eq!(srs.len(), n);

    // x = [τ^{n-2}, τ^{n-3}, …, τ^0, 0]  (length n)
    let mut x = Vec::with_capacity(n);
    for i in (0..n - 1).rev() {
        x.push(srs[i]);
    }
    x.push(G1Projective::default()); // identity

    let xext_fft = {
        let mut xext = x;
        xext.resize(2 * n, G1Projective::default());
        let omega2 = nth_root(2 * n);
        g1_fft_in_place(&mut xext, omega2);
        xext
    };

    // [a_{n-1}, 0, …, 0, a_1, …, a_{n-2}] of length 2n
    // toeplitz_coefficients = polynomial[-1:] + [0]*(n+1) + polynomial[1:-1]
    let mut toeplitz = vec![Fr::from(0u64); 2 * n];
    toeplitz[0] = coeffs[n - 1];
    let rest = &coeffs[1..n - 1]; // length n-2
    let start = 1 + (n + 1); // n+2
    for (i, c) in rest.iter().enumerate() {
        toeplitz[start + i] = *c;
    }

    let d2 = domain(2 * n);
    let mut t_fft = toeplitz.clone();
    d2.fft_in_place(&mut t_fft);
    let mut hext = g1_pointwise_mul(&xext_fft, &t_fft);
    let omega2 = nth_root(2 * n);
    g1_ifft_in_place(&mut hext, omega2);
    let h = hext[..n].to_vec();

    let omega = nth_root(n);
    g1_fft(&h, omega)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ec::PrimeGroup;
    use ark_poly::Polynomial;
    use ark_poly::univariate::DensePolynomial;
    use ark_std::UniformRand;

    #[test]
    fn g1_fft_matches_scalar_fft() {
        let n = 8;
        let mut rng = ark_std::test_rng();
        let coeffs: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let g = G1Projective::generator();
        let pts: Vec<_> = coeffs.iter().map(|c| g * *c).collect();
        let omega = nth_root(n);
        let mut sc = coeffs.clone();
        scalar_fft(&mut sc, n);
        let gp = g1_fft(&pts, omega);
        for i in 0..n {
            assert_eq!(gp[i], g * sc[i]);
        }
    }

    #[test]
    fn fk20_matches_naive_quotients() {
        let n = 16;
        let mut rng = ark_std::test_rng();
        let tau = Fr::rand(&mut rng);
        let mut srs = Vec::with_capacity(n);
        let mut acc = Fr::ONE;
        let g = G1Projective::generator();
        for _ in 0..n {
            srs.push(g * acc);
            acc *= tau;
        }
        let coeffs: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let proofs = fk20_all_proofs(&coeffs, &srs);
        let omega = nth_root(n);
        let mut w = Fr::ONE;
        let poly = DensePolynomial::from_coefficients_vec(coeffs.clone());
        for i in 0..n {
            let y = poly.evaluate(&w);
            // naive q(X) = (f(X)-y)/(X-w)
            let mut num = coeffs.clone();
            num[0] -= y;
            // synthetic division by (X - w)
            let mut q = vec![Fr::from(0u64); n - 1];
            q[n - 2] = num[n - 1];
            for k in (0..n - 2).rev() {
                q[k] = num[k + 1] + q[k + 1] * w;
            }
            let mut pi = G1Projective::default();
            for (k, ck) in q.iter().enumerate() {
                pi += srs[k] * *ck;
            }
            let _ = y;
            assert_eq!(proofs[i], pi, "mismatch at i={i}");
            w *= omega;
        }
    }
}
