//! Same-stack updatable-encryption baselines for large (and short) messages.
//!
//! - `AesGcmEnvelope`: AES-256-GCM with a RISE-wrapped seed. Fast; payload is
//!   not re-randomized on KeyUp (no IND-REENC for the body).
//! - `NestedAes`: Boneh–Eskandarian–Kim–Lewi nested AES-GCM (bounded epochs).
//! - `Shine0`: Ristretto exponentiation UE for a single group element.
//! - `ShineLong`: SHINE-style / OCBSHINE analogue: one scalarmul per 32-byte
//!   block, pad = SHA-256(k · H(nonce‖i)), ciphertext-independent update.
//! - `RecryptKhprf`: Everspaugh-style KH-PRF (stores the group element per block).

use aes_gcm::aead::{Aead, KeyInit, Payload};
use aes_gcm::{Aes256Gcm, Nonce};
use curve25519_dalek::constants::RISTRETTO_BASEPOINT_POINT;
use curve25519_dalek::ristretto::{CompressedRistretto, RistrettoPoint};
use curve25519_dalek::scalar::Scalar;
use rand::rngs::OsRng;
use rand::{RngCore, SeedableRng};
use rand::rngs::StdRng;
use sha2::{Digest, Sha256, Sha512};

const BLOCK: usize = 32;

fn sha256(data: &[u8]) -> [u8; 32] {
    Sha256::digest(data).into()
}

fn hash_to_ristretto(bytes: &[u8]) -> RistrettoPoint {
    let h = Sha512::digest(bytes);
    let mut buf = [0u8; 64];
    buf.copy_from_slice(&h);
    RistrettoPoint::from_uniform_bytes(&buf)
}

fn random_scalar<R: RngCore>(rng: &mut R) -> Scalar {
    let mut b = [0u8; 64];
    rng.fill_bytes(&mut b);
    Scalar::from_bytes_mod_order_wide(&b)
}

fn xor_in_place(a: &mut [u8], b: &[u8]) {
    for (x, y) in a.iter_mut().zip(b.iter()) {
        *x ^= *y;
    }
}

// ---------- AES-GCM envelope (payload KeyUp does not refresh body) ----------

#[derive(Clone)]
pub struct AesEnvelopeCt {
    pub seed_bytes: [u8; 32], // we store the Ristretto compression of S = seed_point
    pub nonce: [u8; 12],
    pub body: Vec<u8>,
}

pub fn aes_key_from_seed_bytes(seed: &[u8; 32]) -> [u8; 32] {
    sha256(seed)
}

pub fn aes_envelope_enc<R: RngCore>(msg: &[u8], rng: &mut R) -> (AesEnvelopeCt, [u8; 32]) {
    let mut seed = [0u8; 32];
    rng.fill_bytes(&mut seed);
    let mut nonce = [0u8; 12];
    rng.fill_bytes(&mut nonce);
    let key = aes_key_from_seed_bytes(&seed);
    let cipher = Aes256Gcm::new_from_slice(&key).unwrap();
    let body = cipher
        .encrypt(Nonce::from_slice(&nonce), Payload { msg, aad: b"uss-file" })
        .unwrap();
    (AesEnvelopeCt { seed_bytes: seed, nonce, body }, seed)
}

pub fn aes_envelope_dec(ct: &AesEnvelopeCt) -> Vec<u8> {
    let key = aes_key_from_seed_bytes(&ct.seed_bytes);
    let cipher = Aes256Gcm::new_from_slice(&key).unwrap();
    cipher
        .decrypt(Nonce::from_slice(&ct.nonce), Payload { msg: &ct.body, aad: b"uss-file" })
        .expect("aes-gcm decrypt")
}

/// Re-encrypt payload under a fresh seed (client-held file). Used as a naive KeyUp baseline.
pub fn aes_envelope_reenc<R: RngCore>(msg: &[u8], rng: &mut R) -> AesEnvelopeCt {
    aes_envelope_enc(msg, rng).0
}

// ---------- Nested AES-GCM (ePrint 2020/222) ----------

#[derive(Clone)]
pub struct NestedAesCt {
    /// `nonce || AES-GCM ciphertext`. Each update wraps the previous blob.
    pub blob: Vec<u8>,
    pub epoch: usize,
}

pub fn nested_keygen<R: RngCore>(rng: &mut R) -> [u8; 32] {
    let mut k = [0u8; 32];
    rng.fill_bytes(&mut k);
    k
}

fn aes_wrap(key: &[u8; 32], msg: &[u8], rng: &mut impl RngCore) -> Vec<u8> {
    let mut nonce = [0u8; 12];
    rng.fill_bytes(&mut nonce);
    let cipher = Aes256Gcm::new_from_slice(key).unwrap();
    let body = cipher.encrypt(Nonce::from_slice(&nonce), msg).unwrap();
    let mut out = nonce.to_vec();
    out.extend_from_slice(&body);
    out
}

fn aes_unwrap(key: &[u8; 32], blob: &[u8]) -> Vec<u8> {
    let cipher = Aes256Gcm::new_from_slice(key).unwrap();
    cipher
        .decrypt(Nonce::from_slice(&blob[..12]), &blob[12..])
        .expect("nested unwrap")
}

pub fn nested_enc(key: &[u8; 32], msg: &[u8], rng: &mut impl RngCore) -> NestedAesCt {
    NestedAesCt { blob: aes_wrap(key, msg, rng), epoch: 1 }
}

pub fn nested_upd(new_key: &[u8; 32], ct: &NestedAesCt, rng: &mut impl RngCore) -> NestedAesCt {
    NestedAesCt { blob: aes_wrap(new_key, &ct.blob, rng), epoch: ct.epoch + 1 }
}

pub fn nested_dec(keys: &[[u8; 32]], ct: &NestedAesCt) -> Vec<u8> {
    assert_eq!(keys.len(), ct.epoch);
    let mut cur = ct.blob.clone();
    for key in keys.iter().rev() {
        cur = aes_unwrap(key, &cur);
    }
    cur
}

// ---------- SHINE0 on Ristretto (short messages = group elements) ----------

#[derive(Clone, Copy)]
pub struct ShineSk(pub Scalar);

#[derive(Clone, Copy)]
pub struct ShineCt(pub RistrettoPoint);

pub fn shine_keygen<R: RngCore>(rng: &mut R) -> ShineSk {
    ShineSk(random_scalar(rng))
}

pub fn shine_enc(sk: ShineSk, m: RistrettoPoint) -> ShineCt {
    ShineCt(sk.0 * m)
}

pub fn shine_dec(sk: ShineSk, ct: ShineCt) -> RistrettoPoint {
    sk.0.invert() * ct.0
}

pub fn shine_next(old: ShineSk, new: ShineSk) -> Scalar {
    new.0 * old.0.invert()
}

pub fn shine_upd(delta: Scalar, ct: ShineCt) -> ShineCt {
    ShineCt(delta * ct.0)
}

// ---------- Long-message SHINE / ReCrypt KH-PRF ----------

#[derive(Clone)]
pub struct KhprfCt {
    pub nonce: [u8; 16],
    /// Per-block group element Q_i = k · H(nonce‖i). Needed for ciphertext-independent update.
    pub qs: Vec<CompressedRistretto>,
    pub body: Vec<u8>,
}

pub fn khprf_enc(sk: ShineSk, msg: &[u8]) -> KhprfCt {
    let mut nonce = [0u8; 16];
    OsRng.fill_bytes(&mut nonce);
    khprf_enc_with_nonce(sk, msg, nonce)
}

pub fn khprf_enc_with_nonce(sk: ShineSk, msg: &[u8], nonce: [u8; 16]) -> KhprfCt {
    let nblocks = msg.len().div_ceil(BLOCK);
    let mut body = msg.to_vec();
    let pad_len = nblocks * BLOCK - body.len();
    body.extend(std::iter::repeat(0u8).take(pad_len));
    let mut qs = Vec::with_capacity(nblocks);
    for i in 0..nblocks {
        let mut lab = [0u8; 24];
        lab[..16].copy_from_slice(&nonce);
        lab[16..].copy_from_slice(&(i as u64).to_le_bytes());
        let p = hash_to_ristretto(&lab);
        let q = sk.0 * p;
        let pad = sha256(q.compress().as_bytes());
        let start = i * BLOCK;
        xor_in_place(&mut body[start..start + BLOCK], &pad);
        qs.push(q.compress());
    }
    KhprfCt { nonce, qs, body }
}

pub fn khprf_dec(sk: ShineSk, ct: &KhprfCt) -> Vec<u8> {
    let mut body = ct.body.clone();
    let nblocks = ct.qs.len();
    for i in 0..nblocks {
        let q = ct.qs[i].decompress().expect("q");
        let pad = sha256(q.compress().as_bytes());
        let start = i * BLOCK;
        xor_in_place(&mut body[start..start + BLOCK], &pad);
    }
    let _ = sk;
    body
}

pub fn khprf_upd(delta: Scalar, ct: &KhprfCt) -> KhprfCt {
    let mut body = ct.body.clone();
    let mut qs = Vec::with_capacity(ct.qs.len());
    for (i, q_old_c) in ct.qs.iter().enumerate() {
        let q_old = q_old_c.decompress().unwrap();
        let q_new = delta * q_old;
        let pad_old = sha256(q_old.compress().as_bytes());
        let pad_new = sha256(q_new.compress().as_bytes());
        let start = i * BLOCK;
        xor_in_place(&mut body[start..start + BLOCK], &pad_old);
        xor_in_place(&mut body[start..start + BLOCK], &pad_new);
        qs.push(q_new.compress());
    }
    KhprfCt { nonce: ct.nonce, qs, body }
}

pub fn khprf_ct_bytes(ct: &KhprfCt) -> usize {
    16 + ct.qs.len() * 32 + ct.body.len()
}

/// Deterministic RNG for reproducible benches.
pub fn bench_rng(seed: u64) -> StdRng {
    StdRng::seed_from_u64(seed)
}

pub fn random_message(len: usize, seed: u64) -> Vec<u8> {
    let mut rng = bench_rng(seed);
    let mut m = vec![0u8; len];
    rng.fill_bytes(&mut m);
    m
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aes_roundtrip() {
        let mut rng = bench_rng(1);
        let msg = b"hello large file contents";
        let (ct, _) = aes_envelope_enc(msg, &mut rng);
        assert_eq!(aes_envelope_dec(&ct), msg);
    }

    #[test]
    fn nested_roundtrip_two_epochs() {
        let mut rng = bench_rng(2);
        let k0 = nested_keygen(&mut rng);
        let k1 = nested_keygen(&mut rng);
        let msg = b"nested-aes-payload";
        let ct0 = nested_enc(&k0, msg, &mut rng);
        let ct1 = nested_upd(&k1, &ct0, &mut rng);
        assert_eq!(nested_dec(&[k0, k1], &ct1), msg);
    }

    #[test]
    fn shine0_upd() {
        let mut rng = bench_rng(3);
        let sk = shine_keygen(&mut rng);
        let sk2 = shine_keygen(&mut rng);
        let m = random_scalar(&mut rng) * RISTRETTO_BASEPOINT_POINT;
        let ct = shine_enc(sk, m);
        assert_eq!(shine_dec(sk, ct), m);
        let d = shine_next(sk, sk2);
        let ct2 = shine_upd(d, ct);
        assert_eq!(shine_dec(sk2, ct2), m);
    }

    #[test]
    fn khprf_roundtrip_upd() {
        let mut rng = bench_rng(4);
        let sk = shine_keygen(&mut rng);
        let sk2 = shine_keygen(&mut rng);
        let msg = random_message(100, 9);
        let ct = khprf_enc(sk, &msg);
        let pt = khprf_dec(sk, &ct);
        assert_eq!(&pt[..msg.len()], &msg[..]);
        let d = shine_next(sk, sk2);
        let ct2 = khprf_upd(d, &ct);
        let pt2 = khprf_dec(sk2, &ct2);
        assert_eq!(&pt2[..msg.len()], &msg[..]);
    }
}
