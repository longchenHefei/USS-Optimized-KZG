//! Catalano–Fiore homomorphic vector commitment (the ASIACRYPT 2023 USS baseline).
//! Type-3 pairing: verify e(C − [m] h_i, ĥ_i) = e(Λ_i, g2).

use ark_bls12_381::{Bls12_381, Fr, G1Affine, G1Projective, G2Affine, G2Projective};
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup, PrimeGroup};
use ark_ff::UniformRand;
use ark_std::rand::RngCore;

#[derive(Clone)]
pub struct CfHvc {
    pub n: usize,
    pub h: Vec<G1Affine>,
    pub h2: Vec<G2Affine>,
    pub hij: Vec<Vec<G1Affine>>,
    pub g2: G2Affine,
    pub g1: G1Affine,
}

#[derive(Clone, Copy, Debug)]
pub struct CfCommitment(pub G1Projective);

#[derive(Clone, Copy, Debug)]
pub struct CfOpening(pub G1Projective);

impl CfHvc {
    pub fn setup<R: RngCore>(n: usize, rng: &mut R) -> Self {
        let g1 = G1Projective::generator();
        let g2p = G2Projective::generator();
        let g2 = g2p.into_affine();
        let z: Vec<Fr> = (0..=n).map(|_| Fr::rand(rng)).collect();
        let h: Vec<G1Affine> = z.iter().map(|zi| (g1 * *zi).into_affine()).collect();
        let h2: Vec<G2Affine> = z.iter().map(|zi| (g2p * *zi).into_affine()).collect();
        let mut hij = vec![vec![G1Affine::zero(); n + 1]; n + 1];
        for i in 0..=n {
            for j in 0..=n {
                if i != j {
                    hij[i][j] = (g1 * (z[i] * z[j])).into_affine();
                }
            }
        }
        Self {
            n,
            h,
            h2,
            hij,
            g2,
            g1: g1.into_affine(),
        }
    }

    pub fn crs_bytes(&self) -> usize {
        let g1s = (self.n + 1) + (self.n + 1) * self.n;
        let g2s = self.n + 1;
        g1s * 48 + g2s * 96
    }

    pub fn commit(&self, values: &[Fr], r: Fr) -> CfCommitment {
        assert_eq!(values.len(), self.n);
        let mut c = G1Projective::from(self.h[self.n]) * r;
        for (i, m) in values.iter().enumerate() {
            c += G1Projective::from(self.h[i]) * *m;
        }
        CfCommitment(c)
    }

    pub fn open(&self, i: usize, values: &[Fr], r: Fr) -> CfOpening {
        let mut acc = G1Projective::from(self.hij[i][self.n]) * r;
        for (j, m) in values.iter().enumerate() {
            if j != i {
                acc += G1Projective::from(self.hij[i][j]) * *m;
            }
        }
        CfOpening(acc)
    }

    pub fn open_all(&self, values: &[Fr], r: Fr) -> Vec<CfOpening> {
        (0..self.n).map(|i| self.open(i, values, r)).collect()
    }

    pub fn verify(&self, c: CfCommitment, value: Fr, i: usize, opening: CfOpening) -> bool {
        let left_g1 = c.0 - G1Projective::from(self.h[i]) * value;
        let lhs = Bls12_381::pairing(left_g1, G2Projective::from(self.h2[i]));
        let rhs = Bls12_381::pairing(opening.0, G2Projective::from(self.g2));
        lhs == rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cf_roundtrip() {
        let mut rng = ark_std::test_rng();
        let n = 4;
        let cf = CfHvc::setup(n, &mut rng);
        let values: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let r = Fr::rand(&mut rng);
        let c = cf.commit(&values, r);
        for i in 0..n {
            let pi = cf.open(i, &values, r);
            assert!(cf.verify(c, values[i], i, pi));
        }
    }
}
