//! One RLWE ciphertext per file, holding the Ring-LWR mask secret \(\sigma\in R_q\).
//!
//! KeyUp samples an ephemeral \(\delta\), adds \(F_{\Delta+\delta}\) to the body, and adds
//! \(\delta\) into this ciphertext under the new header key. \(\delta\) is not part of the
//! token. A later snapshot that has the token, the old body, and the new body therefore
//! cannot test \(c'-c=F_\Delta\).
//!
//! Ciphertexts live modulo the Goldilocks prime \(2^{64}-2^{32}+1\). The message is a
//! vector in \(\mathbb{Z}_q^N\), placed in the high bits (`SCALE`) so decryption rounds
//! off the RLWE noise. The key-switching key is sent once per epoch.

use rand::RngCore;

use crate::rlwr_ue::{self, N, Q};

/// Goldilocks prime \(2^{64}-2^{32}+1\).
const QH: u128 = 0xFFFF_FFFF_0000_0001;

/// Message \(\sigma\) occupies bits `[SCALE_BITS, SCALE_BITS+30)`.
const SCALE_BITS: u32 = 26;
const SCALE: u128 = 1u128 << SCALE_BITS;

/// Gadget base for key switching. Sixteen 4-bit digits cover a 64-bit coefficient.
const BASE_BITS: u32 = 4;
const DIGITS: usize = 16;
const BASE: u64 = 1u64 << BASE_BITS;

#[derive(Clone)]
pub struct HdrPk {
    pub a: Vec<u64>,
    pub b: Vec<u64>,
}

#[derive(Clone)]
pub struct HdrSk {
    pub s: Vec<u64>,
}

#[derive(Clone)]
pub struct HdrCt {
    pub c1: Vec<u64>,
    pub c2: Vec<u64>,
}

#[derive(Clone)]
pub struct SwitchKey {
    pub parts: Vec<HdrCt>,
}

fn addh(a: u64, b: u64) -> u64 {
    let s = a as u128 + b as u128;
    if s >= QH { (s - QH) as u64 } else { s as u64 }
}
fn subh(a: u64, b: u64) -> u64 {
    if a >= b { a - b } else { (a as u128 + QH - b as u128) as u64 }
}
fn mulh(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % QH) as u64
}
fn powh(mut b: u64, mut e: u128) -> u64 {
    let mut r = 1u64;
    while e > 0 {
        if e & 1 == 1 {
            r = mulh(r, b);
        }
        b = mulh(b, b);
        e >>= 1;
    }
    r
}
fn neg_one() -> u64 {
    (QH - 1) as u64
}

struct Tables {
    psi: Vec<u64>,
    psi_inv: Vec<u64>,
    omega: u64,
    omega_inv: u64,
    ninv: u64,
}

fn tables() -> &'static Tables {
    use std::sync::OnceLock;
    static T: OnceLock<Tables> = OnceLock::new();
    T.get_or_init(|| {
        let qh_minus_1 = QH - 1;
        let gen = 7u64;
        let psi1 = powh(gen, qh_minus_1 / (2 * N as u128));
        let psi_inv1 = powh(psi1, QH - 2);
        let mut psi = vec![0u64; N];
        let mut psi_inv = vec![0u64; N];
        let mut p = 1u64;
        let mut ip = 1u64;
        for i in 0..N {
            psi[i] = p;
            psi_inv[i] = ip;
            p = mulh(p, psi1);
            ip = mulh(ip, psi_inv1);
        }
        let omega = mulh(psi1, psi1);
        Tables {
            psi,
            psi_inv,
            omega,
            omega_inv: powh(omega, QH - 2),
            ninv: powh(N as u64, QH - 2),
        }
    })
}

fn bitrev(a: &mut [u64]) {
    let n = a.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            a.swap(i, j);
        }
    }
}

fn ntt_cyclic(a: &mut [u64], omega: u64) {
    let n = a.len();
    bitrev(a);
    let mut len = 2;
    while len <= n {
        let wlen = powh(omega, (n / len) as u128);
        for i in (0..n).step_by(len) {
            let mut w = 1u64;
            for j in 0..len / 2 {
                let u = a[i + j];
                let v = mulh(a[i + j + len / 2], w);
                a[i + j] = addh(u, v);
                a[i + j + len / 2] = subh(u, v);
                w = mulh(w, wlen);
            }
        }
        len <<= 1;
    }
}

fn twist_fwd(a: &[u64]) -> Vec<u64> {
    let t = tables();
    let mut b = vec![0u64; N];
    for i in 0..N {
        b[i] = mulh(a[i], t.psi[i]);
    }
    ntt_cyclic(&mut b, t.omega);
    b
}

fn twist_inv(hat: &[u64]) -> Vec<u64> {
    let t = tables();
    let mut b = hat.to_vec();
    ntt_cyclic(&mut b, t.omega_inv);
    for i in 0..N {
        b[i] = mulh(mulh(b[i], t.ninv), t.psi_inv[i]);
    }
    b
}

fn pointwise(a: &[u64], s: &[u64]) -> Vec<u64> {
    let mut c = vec![0u64; N];
    for i in 0..N {
        c[i] = mulh(a[i], s[i]);
    }
    c
}

fn mul_poly(a: &[u64], b: &[u64]) -> Vec<u64> {
    twist_inv(&pointwise(&twist_fwd(a), &twist_fwd(b)))
}

fn add_poly(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut c = vec![0u64; N];
    for i in 0..N {
        c[i] = addh(a[i], b[i]);
    }
    c
}

fn uniform_poly(rng: &mut impl RngCore) -> Vec<u64> {
    let mut a = vec![0u64; N];
    for x in &mut a {
        let mut buf = [0u8; 16];
        rng.fill_bytes(&mut buf);
        *x = (u128::from_le_bytes(buf) % QH) as u64;
    }
    a
}

fn ternary(rng: &mut impl RngCore) -> Vec<u64> {
    let mut a = vec![0u64; N];
    let mut buf = [0u8; 1];
    for x in &mut a {
        rng.fill_bytes(&mut buf);
        *x = match buf[0] % 3 {
            0 => 0,
            1 => 1,
            _ => neg_one(),
        };
    }
    a
}

/// Uniform mask secret in \(R_q\), the same ring as the Ring-LWR key.
pub fn sample_sigma(rng: &mut impl RngCore) -> Vec<u64> {
    let mut a = vec![0u64; N];
    for x in &mut a {
        let mut buf = [0u8; 8];
        rng.fill_bytes(&mut buf);
        *x = u64::from_le_bytes(buf) % Q;
    }
    a
}

pub fn keygen(rng: &mut impl RngCore) -> (HdrSk, HdrPk) {
    let s = ternary(rng);
    let a = uniform_poly(rng);
    let e = ternary(rng);
    let b = add_poly(&mul_poly(&a, &s), &e);
    (HdrSk { s }, HdrPk { a, b })
}

fn scale_message(msg: &[u64]) -> Vec<u64> {
    let mut out = vec![0u64; N];
    for i in 0..N {
        out[i] = ((msg[i] as u128 * SCALE) % QH) as u64;
    }
    out
}

fn encrypt_raw(pk: &HdrPk, raw: &[u64], rng: &mut impl RngCore) -> HdrCt {
    let r = ternary(rng);
    let e1 = ternary(rng);
    let e2 = ternary(rng);
    let c1 = add_poly(&mul_poly(&pk.a, &r), &e1);
    let c2 = add_poly(&add_poly(&mul_poly(&pk.b, &r), &e2), raw);
    HdrCt { c1, c2 }
}

pub fn encrypt(pk: &HdrPk, sigma: &[u64], rng: &mut impl RngCore) -> HdrCt {
    encrypt_raw(pk, &scale_message(sigma), rng)
}

fn center(x: u64) -> i128 {
    let half = (QH / 2) as u64;
    if x >= half {
        x as i128 - QH as i128
    } else {
        x as i128
    }
}

pub fn decrypt(sk: &HdrSk, ct: &HdrCt) -> Vec<u64> {
    let prod = mul_poly(&ct.c1, &sk.s);
    let mut msg = vec![0u64; N];
    for i in 0..N {
        let d = center(subh(ct.c2[i], prod[i]));
        let rounded = (d + (SCALE as i128 / 2)).div_euclid(SCALE as i128);
        let mut m = rounded.rem_euclid(Q as i128) as u64;
        m %= Q;
        msg[i] = m;
    }
    msg
}

pub fn add(x: &HdrCt, y: &HdrCt) -> HdrCt {
    HdrCt {
        c1: add_poly(&x.c1, &y.c1),
        c2: add_poly(&x.c2, &y.c2),
    }
}

pub fn switch_key(old: &HdrSk, new_pk: &HdrPk, rng: &mut impl RngCore) -> SwitchKey {
    let mut parts = Vec::with_capacity(DIGITS);
    let mut place = vec![0u64; N];
    // `old.s` is ternary in \(\{0,1,-1\}\). Digit \(i\) encrypts \(s\cdot\mathrm{BASE}^i\).
    let mut pow = 1u64;
    for _ in 0..DIGITS {
        for j in 0..N {
            place[j] = mulh(old.s[j], pow);
        }
        parts.push(encrypt_raw(new_pk, &place, rng));
        pow = mulh(pow, BASE);
    }
    SwitchKey { parts }
}

fn digits(x: u64) -> [u64; DIGITS] {
    let mut d = [0u64; DIGITS];
    let mut v = x;
    for i in 0..DIGITS {
        d[i] = v % BASE;
        v /= BASE;
    }
    d
}

fn sub_poly(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut c = vec![0u64; N];
    for i in 0..N {
        c[i] = subh(a[i], b[i]);
    }
    c
}

pub fn keyswitch(ksk: &SwitchKey, ct: &HdrCt) -> HdrCt {
    // c2 - c1·s_old = m. Digits of c1 times a key-switching key that encrypts
    // s_old·BASE^k produce s_old·c1, so the new ciphertext is
    // (Σ d_k a_k, c2 - Σ d_k b_k).
    let mut acc1 = vec![0u64; N];
    let mut acc2 = vec![0u64; N];
    for k in 0..DIGITS {
        let mut dig_poly = vec![0u64; N];
        for i in 0..N {
            dig_poly[i] = digits(ct.c1[i])[k];
        }
        acc1 = add_poly(&acc1, &mul_poly(&ksk.parts[k].c1, &dig_poly));
        acc2 = add_poly(&acc2, &mul_poly(&ksk.parts[k].c2, &dig_poly));
    }
    // Negate the switched `c1` so that `c2 - c1·s_new` recovers the message.
    let mut c1 = vec![0u64; N];
    for i in 0..N {
        c1[i] = subh(0, acc1[i]);
    }
    HdrCt {
        c1,
        c2: sub_poly(&ct.c2, &acc2),
    }
}

pub fn ct_bytes() -> usize {
    2 * N * 8
}

pub fn ksk_bytes() -> usize {
    DIGITS * ct_bytes()
}

/// Body key \(s+\sigma\) in the Ring-LWR NTT domain.
pub fn combined_dem(dem: &rlwr_ue::Sk, sigma: &[u64]) -> rlwr_ue::Sk {
    rlwr_ue::add_sk(dem, &rlwr_ue::sk_from_poly(sigma))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn ntt_roundtrip_and_product() {
        let mut rng = StdRng::seed_from_u64(1);
        let a = uniform_poly(&mut rng);
        let back = twist_inv(&twist_fwd(&a));
        assert_eq!(a, back);
        let b = ternary(&mut rng);
        let mut exp = vec![0u64; N];
        for i in 0..N {
            for j in 0..N {
                let prod = mulh(a[i], b[j]);
                if i + j < N {
                    exp[i + j] = addh(exp[i + j], prod);
                } else {
                    exp[i + j - N] = subh(exp[i + j - N], prod);
                }
            }
        }
        assert_eq!(mul_poly(&a, &b), exp);
    }

    #[test]
    fn encrypt_decrypt_and_rotate() {
        let mut rng = StdRng::seed_from_u64(7);
        let (sk, pk) = keygen(&mut rng);
        let sigma = sample_sigma(&mut rng);
        let mut ct = encrypt(&pk, &sigma, &mut rng);
        assert_eq!(decrypt(&sk, &ct), sigma);

        let mut cur_sk = sk;
        let mut acc = sigma.clone();
        // A 127-step run of this loop stays inside the rounding budget; the
        // short loop below is the regression check.
        for step in 0..8 {
            let (nsk, npk) = keygen(&mut rng);
            let ksk = switch_key(&cur_sk, &npk, &mut rng);
            let delta = sample_sigma(&mut rng);
            ct = add(&keyswitch(&ksk, &ct), &encrypt(&npk, &delta, &mut rng));
            for i in 0..N {
                acc[i] = (acc[i] + delta[i]) % Q;
            }
            assert_eq!(decrypt(&nsk, &ct), acc, "rotation {step}");
            cur_sk = nsk;
        }
    }
}
