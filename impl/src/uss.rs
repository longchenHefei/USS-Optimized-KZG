//! Updatable secure storage: KHVC stub + RISE-wrapped openings + Ring-LWR file bodies.
//!
//! The DEM key is one Ring-LWR secret for the epoch, plus a per-file mask \(\sigma\)
//! stored in an RLWE header. KeyUp sends \(\Delta=s'-s\) and a header key-switching
//! key. The server samples an ephemeral \(\delta\), adds \(F_{\Delta+\delta}\) to every
//! body, and folds \(\delta\) into the header. \(\delta\) is erased, so the public
//! token does not determine the new body. Openings stay under RISE. Plaintext
//! integrity is the client-side stub, not a MAC on the DEM.
//!
//! FileUp openings use the one-hot closed form (`Khvc::open_one_hot`). KeyUp reuses
//! cached public `a_b` NTTs stored with each file (not secret; rebuildable from tweak).

use crate::hdr::{self, HdrCt, HdrPk, HdrSk};
use crate::khvc::{hash_to_fr, Commitment, Khvc, Opening};
use crate::rise::{rise_ct_bytes, rise_dec, rise_enc, rise_hom, rise_next, rise_upd, RiseCt, RisePk, RiseSk};
use crate::rlwr_ue::{self, a_ntt_block, dec_block, enc_block, upd_block, Sk};
use ark_bls12_381::Fr;
use ark_ff::UniformRand;
use ark_std::rand::RngCore;

#[derive(Clone)]
pub struct UssPp {
    pub khvc: Khvc,
}

#[derive(Clone)]
pub struct UssSk {
    pub rise: RiseSk,
    pub rise_pk: RisePk,
    pub dem: Sk,
    pub hdr: HdrSk,
    pub hdr_pk: HdrPk,
}

#[derive(Clone)]
pub struct FileRecord {
    /// Domain separator for a_b = H(tweak ‖ b). Stored in the clear; it is not a key.
    pub tweak: [u8; 16],
    /// Ring-LWR ciphertext of ν ‖ file, one `u16` per coefficient (value in `0..p`).
    pub body: Vec<u16>,
    pub body_len: usize,
    /// Cached twisted NTT of each public `a_b`. Rebuildable from `tweak`; speeds KeyUp/Rev.
    pub a_ntt: Vec<Vec<u64>>,
    pub proof_ct: RiseCt,
    /// RLWE encryption of the per-file mask secret \(\sigma\).
    pub mask: HdrCt,
}

#[derive(Clone)]
pub struct Repository {
    pub files: Vec<FileRecord>,
    pub stub: Commitment,
}

pub fn pargen<R: RngCore>(n: usize, rng: &mut R) -> UssPp {
    UssPp {
        khvc: Khvc::setup(n, rng),
    }
}

pub fn keygen<R: RngCore>(rng: &mut R) -> UssSk {
    let (rise, rise_pk) = RiseSk::keygen(rng);
    let (hdr, hdr_pk) = hdr::keygen(rng);
    UssSk {
        rise,
        rise_pk,
        dem: rlwr_ue::keygen(rng),
        hdr,
        hdr_pk,
    }
}

fn pack(nu: &[u8; 32], file: &[u8]) -> Vec<u8> {
    let mut inner = Vec::with_capacity(32 + file.len());
    inner.extend_from_slice(nu);
    inner.extend_from_slice(file);
    inner
}

fn body_of(tweak: &[u8; 16], dem: &Sk, msg: &[u8]) -> (Vec<u16>, usize, Vec<Vec<u64>>) {
    let nblocks = msg.len().div_ceil(rlwr_ue::BLOCK);
    let mut y = Vec::with_capacity(nblocks * rlwr_ue::N);
    let mut a_ntt = Vec::with_capacity(nblocks);
    for b in 0..nblocks {
        let a = a_ntt_block(tweak, b);
        let start = b * rlwr_ue::BLOCK;
        let end = (start + rlwr_ue::BLOCK).min(msg.len());
        y.extend(enc_block(&a, dem, &msg[start..end]));
        a_ntt.push(a);
    }
    (y, msg.len(), a_ntt)
}

fn open_body(a_ntt: &[Vec<u64>], dem: &Sk, y: &[u16], len: usize) -> Vec<u8> {
    let mut msg = Vec::with_capacity(len);
    let nblocks = y.len() / rlwr_ue::N;
    for b in 0..nblocks {
        if msg.len() >= len {
            break;
        }
        let need = (len - msg.len()).min(rlwr_ue::BLOCK);
        msg.extend(dec_block(&a_ntt[b], dem, &y[b * rlwr_ue::N..], need));
    }
    msg
}

pub struct StagedFile {
    tweak: [u8; 16],
    body: Vec<u16>,
    body_len: usize,
    a_ntt: Vec<Vec<u64>>,
    h: Fr,
    mask: HdrCt,
}

pub fn stage_file<R: RngCore>(sk: &UssSk, file: &[u8], rng: &mut R) -> StagedFile {
    let mut nu = [0u8; 32];
    let mut tweak = [0u8; 16];
    rng.fill_bytes(&mut nu);
    rng.fill_bytes(&mut tweak);
    let mut pre = file.to_vec();
    pre.extend_from_slice(&nu);
    let h = hash_to_fr(&pre);
    let packed = pack(&nu, file);
    let sigma = hdr::sample_sigma(rng);
    let dem = hdr::combined_dem(&sk.dem, &sigma);
    let (body, body_len, a_ntt) = body_of(&tweak, &dem, &packed);
    StagedFile {
        tweak,
        body,
        body_len,
        a_ntt,
        h,
        mask: hdr::encrypt(&sk.hdr_pk, &sigma, rng),
    }
}

pub fn seal<R: RngCore>(pp: &UssPp, sk: &UssSk, staged: Vec<StagedFile>, rng: &mut R) -> Repository {
    let n = pp.khvc.n;
    assert_eq!(staged.len(), n);
    let hs: Vec<Fr> = staged.iter().map(|s| s.h).collect();
    let r = Fr::rand(rng);
    let stub = pp.khvc.commit(&hs, r);
    let openings = pp.khvc.open_all(&hs, r);
    let files = staged
        .into_iter()
        .enumerate()
        .map(|(i, s)| FileRecord {
            tweak: s.tweak,
            body: s.body,
            body_len: s.body_len,
            a_ntt: s.a_ntt,
            proof_ct: rise_enc(sk.rise_pk, openings[i].0, rng),
            mask: s.mask,
        })
        .collect();
    Repository { files, stub }
}

pub fn store<R: RngCore>(pp: &UssPp, sk: &UssSk, files: &[Vec<u8>], rng: &mut R) -> Repository {
    let n = pp.khvc.n;
    assert_eq!(files.len(), n);
    let staged: Vec<StagedFile> = files.iter().map(|f| stage_file(sk, f, rng)).collect();
    seal(pp, sk, staged, rng)
}

fn file_dem(sk: &UssSk, rec: &FileRecord) -> Sk {
    let sigma = hdr::decrypt(&sk.hdr, &rec.mask);
    hdr::combined_dem(&sk.dem, &sigma)
}

pub fn retrieve(pp: &UssPp, sk: &UssSk, repo: &Repository, i: usize) -> Option<Vec<u8>> {
    let rec = &repo.files[i];
    let dem = file_dem(sk, rec);
    let inner = open_body(&rec.a_ntt, &dem, &rec.body, rec.body_len);
    if inner.len() < 32 {
        return None;
    }
    let mut ri = [0u8; 32];
    ri.copy_from_slice(&inner[..32]);
    let file = inner[32..].to_vec();
    let mut pre = file.clone();
    pre.extend_from_slice(&ri);
    let hi = hash_to_fr(&pre);
    let lam = Opening(rise_dec(sk.rise, rec.proof_ct));
    if !pp.khvc.verify(repo.stub, hi, i, lam) {
        return None;
    }
    Some(file)
}

/// Update file `i`.
/// Wire token: `(i, Δh, ρ_δ, C_δ)` plus a fresh UE body (`O(1)+|m'|`, independent of `n`).
/// The server recomputes one-hot openings from the public SRS and RISE-encrypts under `Y`.
/// This co-located harness applies the opening refresh locally; the returned size matches the wire token.
pub fn file_up<R: RngCore>(
    pp: &UssPp,
    sk: &UssSk,
    repo: &mut Repository,
    i: usize,
    new_file: &[u8],
    rng: &mut R,
) -> usize {
    let dem_i = file_dem(sk, &repo.files[i]);
    let inner = open_body(&repo.files[i].a_ntt, &dem_i, &repo.files[i].body, repo.files[i].body_len);
    let old_nu = &inner[..32];
    let old_file = &inner[32..];
    let mut pre_old = old_file.to_vec();
    pre_old.extend_from_slice(old_nu);
    let h_old = hash_to_fr(&pre_old);

    let mut new_ri = [0u8; 32];
    rng.fill_bytes(&mut new_ri);
    let mut pre_new = new_file.to_vec();
    pre_new.extend_from_slice(&new_ri);
    let h_new = hash_to_fr(&pre_new);

    let delta_val = h_new - h_old;
    let r_delta = Fr::rand(rng);
    let mut delta_vec = vec![Fr::from(0u64); pp.khvc.n];
    delta_vec[i] = delta_val;
    let c_delta = pp.khvc.commit(&delta_vec, r_delta);
    repo.stub = Khvc::com_hom(repo.stub, c_delta);
    let openings_delta = pp.khvc.open_one_hot(i, delta_val, r_delta);
    for j in 0..pp.khvc.n {
        let enc_delta = rise_enc(sk.rise_pk, openings_delta[j].0, rng);
        repo.files[j].proof_ct = rise_hom(repo.files[j].proof_ct, enc_delta);
    }

    let mut tweak = [0u8; 16];
    rng.fill_bytes(&mut tweak);
    let packed = pack(&new_ri, new_file);
    let sigma = hdr::sample_sigma(rng);
    let dem = hdr::combined_dem(&sk.dem, &sigma);
    let (body, body_len, a_ntt) = body_of(&tweak, &dem, &packed);
    repo.files[i].tweak = tweak;
    repo.files[i].body = body;
    repo.files[i].body_len = body_len;
    repo.files[i].a_ntt = a_ntt;
    repo.files[i].mask = hdr::encrypt(&sk.hdr_pk, &sigma, rng);

    // Wire size: |m'| + ν (32) + tweak (16) + C_δ (48) + Δh||ρ_δ (64)
    // + one RLWE mask header. Opening envelopes are produced by the server.
    new_file.len() + 32 + 16 + 48 + 64 + hdr::ct_bytes()
}

/// Rotate keys. Header token is the RISE fields; the DEM token is one ring element.
/// Returns the sum, in bytes.
pub fn key_up<R: RngCore>(
    pp: &UssPp,
    old: &UssSk,
    new: &UssSk,
    repo: &mut Repository,
    rng: &mut R,
) -> usize {
    let tk = rise_next(old.rise, new.rise);
    let r0 = Fr::rand(rng);
    let zeros = vec![Fr::from(0u64); pp.khvc.n];
    let c0 = pp.khvc.commit(&zeros, r0);
    repo.stub = Khvc::com_hom(repo.stub, c0);
    let openings0: Vec<Opening> = pp
        .khvc
        .z_openings
        .iter()
        .map(|z| Opening(*z * r0))
        .collect();
    let delta = rlwr_ue::token(&old.dem, &new.dem);
    // One key-switching key per epoch, not per file. The stub update below is the
    // client-local half of this co-located harness; the server never stores `stub`.
    let ksk = hdr::switch_key(&old.hdr, &new.hdr_pk, rng);
    for j in 0..pp.khvc.n {
        repo.files[j].proof_ct = rise_upd(tk, repo.files[j].proof_ct, rng);
        let enc0 = rise_enc(new.rise_pk, openings0[j].0, rng);
        repo.files[j].proof_ct = rise_hom(repo.files[j].proof_ct, enc0);
        let delta_mask = hdr::sample_sigma(rng);
        let step = rlwr_ue::add_sk(&delta, &rlwr_ue::sk_from_poly(&delta_mask));
        let file = &mut repo.files[j];
        let nblocks = file.body.len() / rlwr_ue::N;
        debug_assert_eq!(file.a_ntt.len(), nblocks);
        for b in 0..nblocks {
            let start = b * rlwr_ue::N;
            upd_block(&file.a_ntt[b], &step, &mut file.body[start..start + rlwr_ue::N]);
        }
        let switched = hdr::keyswitch(&ksk, &file.mask);
        let fresh = hdr::encrypt(&new.hdr_pk, &delta_mask, rng);
        file.mask = hdr::add(&switched, &fresh);
    }
    // RISE: Δ_R (32) + Y' (48) + ρ0 (32) + epoch word (32).
    // DEM: Δ in R_q (4096). Header: one RLWE key-switching key.
    32 + 48 + 32 + 32 + rlwr_ue::token_bytes() + hdr::ksk_bytes()
}

pub fn repo_bytes(repo: &Repository) -> usize {
    repo.files
        .iter()
        .map(|f| f.body.len() * 2 + 16 + rise_ct_bytes())
        .sum::<usize>()
        + 48
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::test_rng;

    #[test]
    fn store_rev_fileup_keyup() {
        let mut rng = test_rng();
        let n = 4;
        let pp = pargen(n, &mut rng);
        let sk = keygen(&mut rng);
        let files: Vec<Vec<u8>> = (0..n)
            .map(|i| format!("file-{}-contents", i).into_bytes())
            .collect();
        let mut repo = store(&pp, &sk, &files, &mut rng);
        for i in 0..n {
            assert_eq!(retrieve(&pp, &sk, &repo, i).unwrap(), files[i]);
            assert!(!repo.files[i].a_ntt.is_empty());
        }
        let newf = b"updated-file-1".to_vec();
        file_up(&pp, &sk, &mut repo, 1, &newf, &mut rng);
        assert_eq!(retrieve(&pp, &sk, &repo, 1).unwrap(), newf);
        assert_eq!(retrieve(&pp, &sk, &repo, 0).unwrap(), files[0]);

        let before = repo.files[1].body.clone();
        let sk2 = keygen(&mut rng);
        let tok = key_up(&pp, &sk, &sk2, &mut repo, &mut rng);
        assert!(tok > rlwr_ue::token_bytes());
        assert_ne!(repo.files[1].body, before);
        assert!(retrieve(&pp, &sk, &repo, 1).is_none());
        assert_eq!(retrieve(&pp, &sk2, &repo, 1).unwrap(), newf);
        assert_eq!(retrieve(&pp, &sk2, &repo, 2).unwrap(), files[2]);
        let sk3 = keygen(&mut rng);
        key_up(&pp, &sk2, &sk3, &mut repo, &mut rng);
        assert_eq!(retrieve(&pp, &sk3, &repo, 1).unwrap(), newf);
        assert!(retrieve(&pp, &sk2, &repo, 1).is_none());
    }
}
