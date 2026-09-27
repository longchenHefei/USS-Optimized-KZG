//! CPA Ring-LWR updatable encryption from an almost key-homomorphic PRF.
//!
//! Random-oracle shape used by BLMR / BPR and recalled in Boneh–Eskandarian–Kim–Fisch
//! (ASIACRYPT 2020, ePrint 2020/222, §6): for each block index \(i\),
//! \(F_s(i)=\lfloor (p/q)\,(a_i\cdot s)\rfloor\) in \(R_p\), with \(a_i=H(i)\) public.
//! One negacyclic product per block. The update token is \(\Delta=s'-s\in R_q\),
//! independent of the message length. Addition is in \(\mathbb{Z}_p\); each update
//! adds a rounding error in \(\{-1,0,1\}\) per coefficient, so the low bits of
//! each coefficient are padding, not payload.
//!
//! This is a single-threaded reference NTT, not an AVX implementation.
//! Concrete cost of these parameters is estimated in the paper, not in this file.

use rand::RngCore;
use sha2::{Digest, Sha256};

pub const N: usize = 1024;
pub const Q: u64 = 998_244_353;
pub const P: u32 = 65_536;
pub const MSG_BITS: u32 = 8;
pub const NOISE_BITS: u32 = 8;
const MSG_MASK: u32 = (1 << MSG_BITS) - 1;
const NOISE_MASK: u32 = (1 << NOISE_BITS) - 1;

/// Plaintext bytes carried by one ring element (8 bits × N coefficients).
pub const BLOCK: usize = N;

#[derive(Clone)]
pub struct Public {
    /// Twisted NTT of \(a_i\), one per block index.
    pub a_ntt: Vec<Vec<u64>>,
}

#[derive(Clone)]
pub struct Sk {
    pub s_ntt: Vec<u64>,
}

#[derive(Clone)]
pub struct Ct {
    pub len: usize,
    /// Rounded ring elements, `P`-domain, `nblocks * N` coefficients.
    pub y: Vec<u32>,
}

fn addq(a: u64, b: u64) -> u64 {
    let s = a + b;
    if s >= Q { s - Q } else { s }
}
fn subq(a: u64, b: u64) -> u64 {
    if a >= b { a - b } else { a + Q - b }
}
fn mulq(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % Q as u128) as u64
}
fn powq(mut b: u64, mut e: u64) -> u64 {
    let mut r = 1u64;
    while e > 0 {
        if e & 1 == 1 {
            r = mulq(r, b);
        }
        b = mulq(b, b);
        e >>= 1;
    }
    r
}

struct NttTables {
    psi: Vec<u64>,
    psi_inv: Vec<u64>,
    omega: u64,
    omega_inv: u64,
    ninv: u64,
}

fn tables() -> &'static NttTables {
    use std::sync::OnceLock;
    static T: OnceLock<NttTables> = OnceLock::new();
    T.get_or_init(|| {
        let psi1 = powq(3, (Q - 1) / (2 * N as u64));
        let mut psi = vec![0u64; N];
        let mut psi_inv = vec![0u64; N];
        let inv = powq(psi1, Q - 2);
        let mut p = 1u64;
        let mut ip = 1u64;
        for i in 0..N {
            psi[i] = p;
            psi_inv[i] = ip;
            p = mulq(p, psi1);
            ip = mulq(ip, inv);
        }
        NttTables {
            psi,
            psi_inv,
            omega: mulq(psi1, psi1),
            omega_inv: powq(mulq(psi1, psi1), Q - 2),
            ninv: powq(N as u64, Q - 2),
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
        let wlen = powq(omega, (n / len) as u64);
        for i in (0..n).step_by(len) {
            let mut w = 1u64;
            for j in 0..len / 2 {
                let u = a[i + j];
                let v = mulq(a[i + j + len / 2], w);
                a[i + j] = addq(u, v);
                a[i + j + len / 2] = subq(u, v);
                w = mulq(w, wlen);
            }
        }
        len <<= 1;
    }
}

fn twist_fwd(a: &[u64]) -> Vec<u64> {
    let t = tables();
    let mut b = vec![0u64; N];
    for i in 0..N {
        b[i] = mulq(a[i], t.psi[i]);
    }
    ntt_cyclic(&mut b, t.omega);
    b
}

fn twist_inv(hat: &[u64]) -> Vec<u64> {
    let t = tables();
    let mut b = hat.to_vec();
    ntt_cyclic(&mut b, t.omega_inv);
    for i in 0..N {
        b[i] = mulq(mulq(b[i], t.ninv), t.psi_inv[i]);
    }
    b
}

fn pointwise(a: &[u64], s: &[u64]) -> Vec<u64> {
    let mut c = vec![0u64; N];
    for i in 0..N {
        c[i] = mulq(a[i], s[i]);
    }
    c
}

/// \(\lfloor p\cdot x / q\rfloor\) for \(x\in\mathbb{Z}_q\).
fn round_p(x: u64) -> u32 {
    ((x as u128 * P as u128) / Q as u128) as u32
}

fn prf(a_ntt: &[u64], s_ntt: &[u64]) -> Vec<u32> {
    twist_inv(&pointwise(a_ntt, s_ntt))
        .into_iter()
        .map(round_p)
        .collect()
}

fn encode_coeff(m: u8) -> u32 {
    (m as u32) << NOISE_BITS
}

fn decode_coeff(v: u32) -> u8 {
    let low = v & NOISE_MASK;
    let mut high = v >> NOISE_BITS;
    if low >= (1 << (NOISE_BITS - 1)) {
        high = high.wrapping_add(1);
    }
    (high & MSG_MASK) as u8
}

fn random_poly(rng: &mut impl RngCore) -> Vec<u64> {
    let mut a = vec![0u64; N];
    for x in &mut a {
        let mut buf = [0u8; 8];
        rng.fill_bytes(&mut buf);
        *x = u64::from_le_bytes(buf) % Q;
    }
    a
}

fn poly_from_tweak(tweak: &[u8], block: usize) -> Vec<u64> {
    let mut poly = vec![0u64; N];
    let mut filled = 0usize;
    let mut counter = 0u64;
    while filled < N {
        let mut h = Sha256::new();
        h.update(b"rlwr-a");
        h.update(tweak);
        h.update(&(block as u64).to_le_bytes());
        h.update(&counter.to_le_bytes());
        let dig = h.finalize();
        for chunk in dig.chunks(8) {
            if filled >= N {
                break;
            }
            let mut buf = [0u8; 8];
            buf.copy_from_slice(chunk);
            poly[filled] = u64::from_le_bytes(buf) % Q;
            filled += 1;
        }
        counter += 1;
    }
    poly
}

/// One public ring element \(a_b=H(\mathrm{tweak}\Vert b)\), already in the NTT domain.
pub fn a_ntt_block(tweak: &[u8], block: usize) -> Vec<u64> {
    twist_fwd(&poly_from_tweak(tweak, block))
}

/// Public ring elements \(a_b=H(\mathrm{tweak}\Vert b\Vert \mathrm{ctr})\), in the NTT domain.
/// Same tweak and block count always rebuild the same table, which is what KeyUp needs.
pub fn public_from_tweak(tweak: &[u8], nblocks: usize) -> Public {
    let mut a_ntt = Vec::with_capacity(nblocks);
    for b in 0..nblocks {
        a_ntt.push(a_ntt_block(tweak, b));
    }
    Public { a_ntt }
}

pub fn setup(nblocks: usize, rng: &mut impl RngCore) -> Public {
    let mut a_ntt = Vec::with_capacity(nblocks);
    for _ in 0..nblocks {
        a_ntt.push(twist_fwd(&random_poly(rng)));
    }
    Public { a_ntt }
}

pub fn keygen(rng: &mut impl RngCore) -> Sk {
    Sk { s_ntt: twist_fwd(&random_poly(rng)) }
}

/// \(\Delta\) such that updating under \(\Delta\) moves ciphertexts from `old` to `new`.
/// Sent on the wire as \(N\) coefficients in \(\mathbb{Z}_q\) (see [`token_bytes`]).
pub fn token(old: &Sk, new: &Sk) -> Sk {
    let old_c = twist_inv(&old.s_ntt);
    let new_c = twist_inv(&new.s_ntt);
    let mut d = vec![0u64; N];
    for i in 0..N {
        d[i] = subq(new_c[i], old_c[i]);
    }
    Sk { s_ntt: twist_fwd(&d) }
}

pub fn token_bytes() -> usize {
    N * 4
}

/// Encrypt one block into `P`-domain coefficients. `msg` is at most [`BLOCK`] bytes.
pub fn enc_block(a_ntt: &[u64], sk: &Sk, msg: &[u8]) -> Vec<u16> {
    let f = prf(a_ntt, &sk.s_ntt);
    let mut y = Vec::with_capacity(N);
    for j in 0..N {
        let m = if j < msg.len() { msg[j] } else { 0 };
        y.push((encode_coeff(m).wrapping_add(f[j]) % P) as u16);
    }
    y
}

pub fn dec_block(a_ntt: &[u64], sk: &Sk, y: &[u16], nbytes: usize) -> Vec<u8> {
    let f = prf(a_ntt, &sk.s_ntt);
    let mut msg = Vec::with_capacity(nbytes);
    for j in 0..nbytes {
        let v = (y[j] as u32).wrapping_sub(f[j]) % P;
        msg.push(decode_coeff(v));
    }
    msg
}

pub fn upd_block(a_ntt: &[u64], delta: &Sk, y: &mut [u16]) {
    let f = prf(a_ntt, &delta.s_ntt);
    for j in 0..N {
        y[j] = y[j].wrapping_add(f[j] as u16);
    }
}

pub fn enc(pp: &Public, sk: &Sk, msg: &[u8]) -> Ct {
    let nblocks = msg.len().div_ceil(BLOCK);
    assert!(pp.a_ntt.len() >= nblocks, "public a_i shorter than the message");
    let mut y = vec![0u32; nblocks * N];
    for b in 0..nblocks {
        let f = prf(&pp.a_ntt[b], &sk.s_ntt);
        for j in 0..N {
            let m = if b * BLOCK + j < msg.len() {
                msg[b * BLOCK + j]
            } else {
                0
            };
            y[b * N + j] = encode_coeff(m).wrapping_add(f[j]) % P;
        }
    }
    Ct { len: msg.len(), y }
}

pub fn dec(pp: &Public, sk: &Sk, ct: &Ct) -> Vec<u8> {
    let nblocks = ct.y.len() / N;
    let mut msg = vec![0u8; ct.len];
    for b in 0..nblocks {
        let f = prf(&pp.a_ntt[b], &sk.s_ntt);
        for j in 0..N {
            let idx = b * BLOCK + j;
            if idx >= ct.len {
                break;
            }
            let v = ct.y[b * N + j].wrapping_sub(f[j]) % P;
            msg[idx] = decode_coeff(v);
        }
    }
    msg
}

pub fn upd(pp: &Public, delta: &Sk, ct: &Ct) -> Ct {
    let nblocks = ct.y.len() / N;
    let mut y = vec![0u32; ct.y.len()];
    for b in 0..nblocks {
        let f = prf(&pp.a_ntt[b], &delta.s_ntt);
        for j in 0..N {
            y[b * N + j] = ct.y[b * N + j].wrapping_add(f[j]) % P;
        }
    }
    Ct { len: ct.len, y }
}

pub fn ct_bytes(ct: &Ct) -> usize {
    ct.y.len() * 2
}

/// Largest absolute accumulated rounding error after `steps` homomorphic updates,
/// over every coefficient of `nblocks` fresh PRF samples. Error is centered in \((-P/2,P/2]\).
pub fn max_abs_rounding_error(nblocks: usize, steps: usize, rng: &mut impl RngCore) -> u32 {
    let pp = setup(nblocks, rng);
    let mut sk = keygen(rng);
    let mut acc = vec![vec![0u32; N]; nblocks];
    for b in 0..nblocks {
        acc[b] = prf(&pp.a_ntt[b], &sk.s_ntt);
    }
    for _ in 0..steps {
        let nxt = keygen(rng);
        let d = token(&sk, &nxt);
        for b in 0..nblocks {
            let fd = prf(&pp.a_ntt[b], &d.s_ntt);
            for j in 0..N {
                acc[b][j] = acc[b][j].wrapping_add(fd[j]) % P;
            }
        }
        sk = nxt;
    }
    let mut max_e = 0u32;
    for b in 0..nblocks {
        let fresh = prf(&pp.a_ntt[b], &sk.s_ntt);
        for j in 0..N {
            let mut diff = acc[b][j].wrapping_sub(fresh[j]) % P;
            if diff > P / 2 {
                diff = P - diff;
            }
            max_e = max_e.max(diff);
        }
    }
    max_e
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    fn schoolbook(a: &[u64], b: &[u64]) -> Vec<u64> {
        let mut c = vec![0u64; N];
        for i in 0..N {
            for j in 0..N {
                let prod = mulq(a[i], b[j]);
                if i + j < N {
                    c[i + j] = addq(c[i + j], prod);
                } else {
                    c[i + j - N] = subq(c[i + j - N], prod);
                }
            }
        }
        c
    }

    #[test]
    fn ntt_matches_schoolbook() {
        let mut rng = StdRng::seed_from_u64(1);
        for _ in 0..3 {
            let a = random_poly(&mut rng);
            let b = random_poly(&mut rng);
            let got = twist_inv(&pointwise(&twist_fwd(&a), &twist_fwd(&b)));
            let exp = schoolbook(&a, &b);
            assert_eq!(got, exp);
        }
    }

    #[test]
    fn one_update_error_is_at_most_one() {
        let mut rng = StdRng::seed_from_u64(2);
        let e = max_abs_rounding_error(2, 1, &mut rng);
        assert!(e <= 1, "one-update rounding error {e}");
    }

    #[test]
    fn tweak_public_roundtrip() {
        let mut rng = StdRng::seed_from_u64(4);
        let msg = b"tweak-domain-file";
        let pp = public_from_tweak(b"file-tweak-01", 1);
        let pp2 = public_from_tweak(b"file-tweak-01", 1);
        assert_eq!(pp.a_ntt, pp2.a_ntt);
        let sk = keygen(&mut rng);
        let ct = enc(&pp, &sk, msg);
        assert_eq!(dec(&pp2, &sk, &ct), msg);
    }

    #[test]
    fn roundtrip_survives_padding_budget() {
        let mut rng = StdRng::seed_from_u64(3);
        let msg = vec![0xA5u8; 4 * 1024];
        let pp = setup(4, &mut rng);
        let mut sk = keygen(&mut rng);
        let mut ct = enc(&pp, &sk, &msg);
        assert_eq!(dec(&pp, &sk, &ct), msg);
        for _ in 0..32 {
            let nxt = keygen(&mut rng);
            let d = token(&sk, &nxt);
            ct = upd(&pp, &d, &ct);
            sk = nxt;
            assert_eq!(dec(&pp, &sk, &ct), msg);
        }
    }
}
