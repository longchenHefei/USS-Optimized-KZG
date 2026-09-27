//! RISE updatable encryption (Lehmann–Tackmann) instantiated in BLS12-381 G1.

use ark_bls12_381::{Fr, G1Projective};
use ark_ec::PrimeGroup;
use ark_ff::fields::Field;
use ark_std::rand::RngCore;
use ark_std::UniformRand;

#[derive(Clone, Copy, Debug)]
pub struct RiseSk(pub Fr);

#[derive(Clone, Copy, Debug)]
pub struct RisePk(pub G1Projective); // Y = sk · G

#[derive(Clone, Copy, Debug)]
pub struct RiseCt {
    pub c1: G1Projective, // r · Y
    pub c2: G1Projective, // r · G + M
}

#[derive(Clone, Copy, Debug)]
pub struct RiseToken {
    pub delta: Fr,        // sk' / sk
    pub y_new: G1Projective,
}

impl RiseSk {
    pub fn keygen<R: RngCore>(rng: &mut R) -> (Self, RisePk) {
        let sk = Fr::rand(rng);
        let y = G1Projective::generator() * sk;
        (Self(sk), RisePk(y))
    }

    pub fn pk(&self) -> RisePk {
        RisePk(G1Projective::generator() * self.0)
    }
}

pub fn rise_enc<R: RngCore>(pk: RisePk, m: G1Projective, rng: &mut R) -> RiseCt {
    let r = Fr::rand(rng);
    RiseCt {
        c1: pk.0 * r,
        c2: G1Projective::generator() * r + m,
    }
}

pub fn rise_dec(sk: RiseSk, ct: RiseCt) -> G1Projective {
    let kinv = sk.0.inverse().unwrap();
    ct.c2 - ct.c1 * kinv
}

pub fn rise_next(old: RiseSk, new: RiseSk) -> RiseToken {
    RiseToken {
        delta: new.0 * old.0.inverse().unwrap(),
        y_new: G1Projective::generator() * new.0,
    }
}

pub fn rise_upd<R: RngCore>(tk: RiseToken, ct: RiseCt, rng: &mut R) -> RiseCt {
    let r = Fr::rand(rng);
    RiseCt {
        c1: ct.c1 * tk.delta + tk.y_new * r,
        c2: ct.c2 + G1Projective::generator() * r,
    }
}

/// Multiplicative (here additive) homomorphism on G1 messages.
pub fn rise_hom(a: RiseCt, b: RiseCt) -> RiseCt {
    RiseCt {
        c1: a.c1 + b.c1,
        c2: a.c2 + b.c2,
    }
}

pub fn rise_ct_bytes() -> usize {
    48 * 2
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rise_enc_dec_upd_hom() {
        let mut rng = ark_std::test_rng();
        let (sk, pk) = RiseSk::keygen(&mut rng);
        let m = G1Projective::generator() * Fr::rand(&mut rng);
        let ct = rise_enc(pk, m, &mut rng);
        assert_eq!(rise_dec(sk, ct), m);

        let (sk2, _) = RiseSk::keygen(&mut rng);
        let tk = rise_next(sk, sk2);
        let ct2 = rise_upd(tk, ct, &mut rng);
        assert_eq!(rise_dec(sk2, ct2), m);

        let m2 = G1Projective::generator() * Fr::rand(&mut rng);
        let ct_b = rise_enc(pk, m2, &mut rng);
        let h = rise_hom(ct, ct_b);
        assert_eq!(rise_dec(sk, h), m + m2);
    }
}
