//! USS with a client-local KHVC stub and a Ring-LWR body with RLWE mask headers.
//! All files share the public block-polynomial sequence. Each file has an
//! independent encrypted mask, and KeyUp adds a fresh mask to each file.

use crate::hdr::{self, HdrCt, HdrPk, HdrSk, SwitchKey};
use crate::khvc::{hash_to_fr, Commitment, Khvc, Opening};
use crate::rise::{
    rise_ct_bytes, rise_dec, rise_enc, rise_hom, rise_next, rise_upd, RiseCt, RisePk, RiseSk,
    RiseToken,
};
use crate::rlwr_ue::{self, a_ntt_block, dec_block, enc_block, upd_block, Sk};
use ark_bls12_381::Fr;
use ark_ff::UniformRand;
use ark_std::rand::RngCore;
use std::sync::{Arc, RwLock, RwLockReadGuard};

const NONCE_BYTES: usize = 32;
const SHARED_BODY_DOMAIN: &[u8] = b"uss-shared-body-v1";

#[derive(Clone)]
pub struct UssPp {
    pub khvc: Khvc,
    block_cache: Arc<RwLock<Vec<Vec<u64>>>>,
}

impl UssPp {
    fn blocks(&self, n: usize) -> RwLockReadGuard<'_, Vec<Vec<u64>>> {
        let blocks = self.block_cache.read().unwrap();
        if blocks.len() >= n {
            return blocks;
        }
        drop(blocks);
        {
            let mut blocks = self.block_cache.write().unwrap();
            let start = blocks.len();
            blocks.extend((start..n).map(|i| a_ntt_block(SHARED_BODY_DOMAIN, i)));
        }
        self.block_cache.read().unwrap()
    }
}

#[derive(Clone)]
pub struct UssSk {
    pub rise: RiseSk,
    pub rise_pk: RisePk,
    pub dem: Sk,
    pub hdr: HdrSk,
    pub hdr_pk: HdrPk,
}

/// The server may see the commitment, but only the client owns this state.
#[derive(Clone)]
pub struct ClientState {
    stub: Commitment,
    blind: Fr,
}

impl ClientState {
    pub fn stub(&self) -> Commitment {
        self.stub
    }

    pub fn rerandomize<R: RngCore>(&mut self, pp: &UssPp, proof_delta: Fr, rng: &mut R) {
        let blind = Fr::rand(rng);
        self.stub = Khvc::com_hom(
            self.stub,
            Commitment(pp.khvc.z_g1 * (proof_delta + blind - self.blind)),
        );
        self.blind = blind;
    }
}

#[derive(Clone)]
pub struct FileRecord {
    pub body: Vec<u16>,
    pub body_len: usize,
    pub proof_ct: RiseCt,
    pub mask: HdrCt,
}

#[derive(Clone)]
pub struct Repository {
    pub files: Vec<FileRecord>,
}

/// The same epoch token can update any repository using the old epoch key.
pub struct KeyUpToken {
    rise: RiseToken,
    body_delta: Sk,
    switch_key: SwitchKey,
    new_header_pk: HdrPk,
}

impl KeyUpToken {
    pub fn bytes(&self) -> usize {
        32 + 48 + rlwr_ue::token_bytes() + hdr::ksk_bytes() + hdr::ct_bytes()
    }
}

pub fn pargen<R: RngCore>(n: usize, rng: &mut R) -> UssPp {
    UssPp {
        khvc: Khvc::setup(n, rng),
        block_cache: Arc::default(),
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

pub struct StagedFile {
    body: Vec<u16>,
    body_len: usize,
    digest: Fr,
    mask: HdrCt,
}

pub fn stage_file<R: RngCore>(pp: &UssPp, sk: &UssSk, file: &[u8], rng: &mut R) -> StagedFile {
    let mut nonce = [0u8; NONCE_BYTES];
    rng.fill_bytes(&mut nonce);
    let mut packed = Vec::with_capacity(file.len() + NONCE_BYTES);
    packed.extend_from_slice(file);
    packed.extend_from_slice(&nonce);
    let digest = hash_to_fr(&packed);

    let sigma = hdr::sample_sigma(rng);
    let dem = hdr::combined_dem(&sk.dem, &sigma);
    let n = packed.len().div_ceil(rlwr_ue::BLOCK);
    let blocks = pp.blocks(n);
    let mut body = Vec::with_capacity(n * rlwr_ue::N);
    for (chunk, a) in packed.chunks(rlwr_ue::BLOCK).zip(blocks.iter()) {
        body.extend(enc_block(a, &dem, chunk));
    }
    StagedFile {
        body,
        body_len: packed.len(),
        digest,
        mask: hdr::encrypt(&sk.hdr_pk, &sigma, rng),
    }
}

pub fn seal<R: RngCore>(
    pp: &UssPp,
    sk: &UssSk,
    staged: Vec<StagedFile>,
    rng: &mut R,
) -> (ClientState, Repository) {
    assert_eq!(staged.len(), pp.khvc.n);
    let digests: Vec<Fr> = staged.iter().map(|s| s.digest).collect();
    let blind = Fr::rand(rng);
    let proof_blind = Fr::rand(rng);
    let stub = pp.khvc.commit(&digests, proof_blind + blind);
    // Proofs have their own blinding. The client adds its private offset on verification.
    let openings = pp.khvc.open_all(&digests, proof_blind);
    let files = staged
        .into_iter()
        .zip(openings)
        .map(|(s, opening)| FileRecord {
            body: s.body,
            body_len: s.body_len,
            proof_ct: rise_enc(sk.rise_pk, opening.0, rng),
            mask: s.mask,
        })
        .collect();
    (ClientState { stub, blind }, Repository { files })
}

pub fn store<R: RngCore>(
    pp: &UssPp,
    sk: &UssSk,
    files: &[Vec<u8>],
    rng: &mut R,
) -> (ClientState, Repository) {
    assert_eq!(files.len(), pp.khvc.n);
    let staged = files.iter().map(|f| stage_file(pp, sk, f, rng)).collect();
    seal(pp, sk, staged, rng)
}

fn verified_file(
    pp: &UssPp,
    sk: &UssSk,
    client: &ClientState,
    repo: &Repository,
    i: usize,
) -> Option<(Vec<u8>, Fr)> {
    let blind_opening = *pp.khvc.z_openings.get(i)? * client.blind;
    let record = repo.files.get(i)?;
    let n = record.body_len.div_ceil(rlwr_ue::BLOCK);
    if record.body_len < NONCE_BYTES
        || n.checked_mul(rlwr_ue::N)? != record.body.len()
        || !hdr::valid_ciphertext(&record.mask)
    {
        return None;
    }

    let sigma = hdr::decrypt(&sk.hdr, &record.mask);
    let dem = hdr::combined_dem(&sk.dem, &sigma);
    let blocks = pp.blocks(n);
    let mut packed = Vec::with_capacity(record.body_len);
    for (a, chunk) in blocks.iter().zip(record.body.chunks_exact(rlwr_ue::N)) {
        let need = (record.body_len - packed.len()).min(rlwr_ue::BLOCK);
        packed.extend(dec_block(a, &dem, chunk, need));
    }

    let digest = hash_to_fr(&packed);
    let opening = Opening(rise_dec(sk.rise, record.proof_ct) + blind_opening);
    if !pp.khvc.verify(client.stub, digest, i, opening) {
        return None;
    }
    packed.truncate(record.body_len - NONCE_BYTES);
    Some((packed, digest))
}

pub fn retrieve(
    pp: &UssPp,
    sk: &UssSk,
    client: &ClientState,
    repo: &Repository,
    i: usize,
) -> Option<Vec<u8>> {
    verified_file(pp, sk, client, repo, i).map(|(file, _)| file)
}

/// Update only after verifying the old record against the client's own stub.
/// Each update blinds its opening delta independently. The server receives that
/// scalar with the update message; it must be hidden when the transcript is withheld.
/// Returns the client-to-server byte count, or None without mutating either state.
pub fn file_up<R: RngCore>(
    pp: &UssPp,
    sk: &UssSk,
    client: &mut ClientState,
    repo: &mut Repository,
    i: usize,
    new_file: &[u8],
    rng: &mut R,
) -> Option<usize> {
    if repo.files.len() != pp.khvc.n {
        return None;
    }
    let (_, old_digest) = verified_file(pp, sk, client, repo, i)?;
    let staged = stage_file(pp, sk, new_file, rng);
    let delta = staged.digest - old_digest;
    let opening_blind = Fr::rand(rng);
    let new_blind = Fr::rand(rng);
    let mut delta_vec = vec![Fr::from(0u64); pp.khvc.n];
    delta_vec[i] = delta;
    let new_stub = Khvc::com_hom(
        client.stub,
        pp.khvc
            .commit(&delta_vec, opening_blind + new_blind - client.blind),
    );
    let openings = pp.khvc.open_one_hot(i, delta, opening_blind);
    for (record, opening) in repo.files.iter_mut().zip(openings) {
        record.proof_ct = rise_hom(record.proof_ct, rise_enc(sk.rise_pk, opening.0, rng));
    }

    let record = &mut repo.files[i];
    record.body = staged.body;
    record.body_len = staged.body_len;
    record.mask = staged.mask;
    client.stub = new_stub;
    client.blind = new_blind;

    // Slot (u64), delta digest and opening blind (Fr), RISE public key, and payload.
    Some(8 + 32 + 32 + 48 + payload_bytes(record))
}

pub fn key_up_token<R: RngCore>(old: &UssSk, new: &UssSk, rng: &mut R) -> KeyUpToken {
    KeyUpToken {
        rise: rise_next(old.rise, new.rise),
        body_delta: rlwr_ue::token(&old.dem, &new.dem),
        switch_key: hdr::switch_key(&old.hdr, &new.hdr_pk, rng),
        new_header_pk: new.hdr_pk.clone(),
    }
}

pub fn rotate_repository<R: RngCore>(
    pp: &UssPp,
    token: &KeyUpToken,
    proof_delta: Fr,
    repo: &mut Repository,
    rng: &mut R,
) {
    assert_eq!(repo.files.len(), pp.khvc.n);
    for (record, opening_shift) in repo.files.iter_mut().zip(&pp.khvc.z_openings) {
        record.proof_ct = rise_upd(token.rise, record.proof_ct, rng);
        record.proof_ct.c2 += *opening_shift * proof_delta;
        let mask_delta = hdr::sample_sigma(rng);
        let step = rlwr_ue::add_sk(&token.body_delta, &rlwr_ue::sk_from_poly(&mask_delta));
        let blocks = pp.blocks(record.body.len() / rlwr_ue::N);
        for (a, chunk) in blocks.iter().zip(record.body.chunks_exact_mut(rlwr_ue::N)) {
            upd_block(a, &step, chunk);
        }
        let switched = hdr::keyswitch(&token.switch_key, &record.mask);
        let fresh = hdr::encrypt(&token.new_header_pk, &mask_delta, rng);
        record.mask = hdr::add(&switched, &fresh);
    }
}

pub fn key_up<R: RngCore>(
    pp: &UssPp,
    old: &UssSk,
    new: &UssSk,
    client: &mut ClientState,
    repo: &mut Repository,
    rng: &mut R,
) -> usize {
    let token = key_up_token(old, new, rng);
    let proof_delta = Fr::rand(rng);
    rotate_repository(pp, &token, proof_delta, repo, rng);
    client.rerandomize(pp, proof_delta, rng);
    token.bytes() + 32
}

fn payload_bytes(record: &FileRecord) -> usize {
    8 + record.body.len() * 2 + hdr::ct_bytes()
}

pub fn repo_bytes(repo: &Repository) -> usize {
    repo.files
        .iter()
        .map(|f| payload_bytes(f) + rise_ct_bytes())
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::test_rng;

    #[test]
    fn fileup_masks_the_two_update_candidate_source_comparison() {
        let mut rng = test_rng();
        let pp = pargen(2, &mut rng);
        let first_key = keygen(&mut rng);
        let candidates = [
            vec![b"source-a".to_vec(), b"other-aa".to_vec()],
            vec![b"source-b".to_vec(), b"other-bb".to_vec()],
        ];
        let stored: Vec<_> = candidates
            .iter()
            .map(|files| store(&pp, &first_key, files, &mut rng))
            .collect();
        let candidate_digests: Vec<_> = stored
            .iter()
            .map(|(client, repo)| verified_file(&pp, &first_key, client, repo, 0).unwrap().1)
            .collect();
        assert_ne!(candidate_digests[0], candidate_digests[1]);
        let second_key = keygen(&mut rng);
        let leaked_key = keygen(&mut rng);
        let unit_opening = pp
            .khvc
            .open(1, &[Fr::from(1u64), Fr::from(0u64)], Fr::from(0u64))
            .0;

        for source in 0..2 {
            let (mut client, mut repo) = stored[source].clone();
            key_up(
                &pp,
                &first_key,
                &second_key,
                &mut client,
                &mut repo,
                &mut rng,
            );
            key_up(
                &pp,
                &second_key,
                &leaked_key,
                &mut client,
                &mut repo,
                &mut rng,
            );

            // The first overwrite makes slot 1 visible to the observer, who
            // has the later epoch key but neither connecting KeyUp transcript.
            file_up(
                &pp,
                &leaked_key,
                &mut client,
                &mut repo,
                1,
                b"visible",
                &mut rng,
            )
            .unwrap();
            let before = rise_dec(leaked_key.rise, repo.files[1].proof_ct);
            let unchanged = repo.files[1].clone();
            let uploaded = file_up(
                &pp,
                &leaked_key,
                &mut client,
                &mut repo,
                0,
                b"replaced",
                &mut rng,
            )
            .unwrap();
            let after = rise_dec(leaked_key.rise, repo.files[1].proof_ct);
            let new_digest = verified_file(&pp, &leaked_key, &client, &repo, 0)
                .unwrap()
                .1;

            // The new digest is available by decrypting the replacement with
            // the leaked epoch key. Neither candidate satisfies the former
            // equality B'_1 - B_1 = (h'_0 - h_{0,source}) A_{0,1}.
            for candidate in &candidate_digests {
                assert_ne!(after - before, unit_opening * (new_digest - candidate));
            }
            assert_eq!(uploaded, 8 + 32 + 32 + 48 + payload_bytes(&repo.files[0]));
            assert_eq!(repo.files[1].body, unchanged.body);
            assert_eq!(repo.files[1].body_len, unchanged.body_len);
            assert_eq!(repo.files[1].mask.c1, unchanged.mask.c1);
            assert_eq!(repo.files[1].mask.c2, unchanged.mask.c2);
            assert_eq!(
                retrieve(&pp, &leaked_key, &client, &repo, 0).unwrap(),
                b"replaced"
            );
            assert_eq!(
                retrieve(&pp, &leaked_key, &client, &repo, 1).unwrap(),
                b"visible"
            );
        }
    }

    #[test]
    fn fileup_samples_a_fresh_opening_blind_for_every_update() {
        let mut rng = test_rng();
        let pp = pargen(4, &mut rng);
        let sk = keygen(&mut rng);
        let mut files = vec![b"initial".to_vec(); 4];
        let (mut client, mut repo) = store(&pp, &sk, &files, &mut rng);
        let mut previous_mask = None;

        for replacement in [b"first".as_slice(), b"second", b"third"] {
            let old_digest = verified_file(&pp, &sk, &client, &repo, 0).unwrap().1;
            let before: Vec<_> = repo
                .files
                .iter()
                .map(|r| rise_dec(sk.rise, r.proof_ct))
                .collect();
            file_up(&pp, &sk, &mut client, &mut repo, 0, replacement, &mut rng).unwrap();
            let new_digest = verified_file(&pp, &sk, &client, &repo, 0).unwrap().1;
            let mut delta = vec![Fr::from(0u64); 4];
            delta[0] = new_digest - old_digest;

            for (i, record) in repo.files.iter().enumerate() {
                let difference = rise_dec(sk.rise, record.proof_ct) - before[i];
                let unblinded = pp.khvc.open(i, &delta, Fr::from(0u64)).0;
                assert_ne!(
                    difference, unblinded,
                    "slot {i} must include the opening blind"
                );
                if i == 1 {
                    let mask = difference - unblinded;
                    if let Some(previous) = previous_mask {
                        assert_ne!(mask, previous, "each FileUp needs independent randomness");
                    }
                    previous_mask = Some(mask);
                }
            }
            files[0] = replacement.to_vec();
            for (i, expected) in files.iter().enumerate() {
                assert_eq!(
                    retrieve(&pp, &sk, &client, &repo, i).as_ref(),
                    Some(expected)
                );
            }
        }
    }

    #[test]
    fn fileup_and_rotation_preserve_current_versions() {
        let mut rng = test_rng();
        let pp = pargen(4, &mut rng);
        let mut sk = keygen(&mut rng);
        let files = vec![
            vec![],
            vec![0x42; 1025],
            b"unchanged".to_vec(),
            vec![0xff; 991],
        ];
        let (mut client, mut repo) = store(&pp, &sk, &files, &mut rng);
        for (i, file) in files.iter().enumerate() {
            assert_eq!(retrieve(&pp, &sk, &client, &repo, i).as_ref(), Some(file));
        }

        let new_file = b"updated";
        assert!(file_up(&pp, &sk, &mut client, &mut repo, 1, new_file, &mut rng).is_some());
        for _ in 0..8 {
            let next = keygen(&mut rng);
            key_up(&pp, &sk, &next, &mut client, &mut repo, &mut rng);
            for (i, file) in files.iter().enumerate() {
                let expected = if i == 1 {
                    new_file.as_slice()
                } else {
                    file.as_slice()
                };
                assert_eq!(retrieve(&pp, &next, &client, &repo, i).unwrap(), expected);
            }
            assert!(retrieve(&pp, &sk, &client, &repo, 1).is_none());
            sk = next;
        }
        assert!(file_up(
            &pp,
            &sk,
            &mut client,
            &mut repo,
            1,
            b"next version",
            &mut rng
        )
        .is_some());
        assert_eq!(
            retrieve(&pp, &sk, &client, &repo, 1).unwrap(),
            b"next version"
        );
    }

    #[test]
    fn fileup_rejects_stale_and_malformed_records_before_mutating_state() {
        let mut rng = test_rng();
        let pp = pargen(2, &mut rng);
        let sk = keygen(&mut rng);
        let (mut client, mut repo) =
            store(&pp, &sk, &[b"old".to_vec(), b"other".to_vec()], &mut rng);
        let mut stale = repo.clone();
        assert!(file_up(&pp, &sk, &mut client, &mut repo, 0, b"current", &mut rng).is_some());
        let state_before = client.clone();
        let untouched_proof = stale.files[1].proof_ct.c1;
        assert!(file_up(&pp, &sk, &mut client, &mut stale, 0, b"rejected", &mut rng).is_none());
        assert_eq!(stale.files[1].proof_ct.c1, untouched_proof);
        assert_eq!(client.stub.0, state_before.stub.0);
        assert_eq!(client.blind, state_before.blind);

        let mut malformed = repo.clone();
        malformed.files[0].body_len = usize::MAX;
        assert!(file_up(
            &pp,
            &sk,
            &mut client,
            &mut malformed,
            0,
            b"rejected",
            &mut rng
        )
        .is_none());
        malformed.files[0] = repo.files[0].clone();
        malformed.files[0].mask.c1.clear();
        assert!(file_up(
            &pp,
            &sk,
            &mut client,
            &mut malformed,
            0,
            b"rejected",
            &mut rng
        )
        .is_none());
        assert!(file_up(
            &pp,
            &sk,
            &mut client,
            &mut malformed,
            2,
            b"rejected",
            &mut rng
        )
        .is_none());
        assert_eq!(client.stub.0, state_before.stub.0);
        assert_eq!(client.blind, state_before.blind);
        assert_eq!(retrieve(&pp, &sk, &client, &repo, 0).unwrap(), b"current");
    }

    #[test]
    fn server_can_rotate_two_repositories_with_one_token_without_client_state() {
        let mut rng = test_rng();
        let pp = pargen(2, &mut rng);
        let sk = keygen(&mut rng);
        let next = keygen(&mut rng);
        let files = vec![vec![0u8; 12], vec![1u8; 12]];
        let (mut first, mut left) = store(&pp, &sk, &files, &mut rng);
        let (mut second, mut right) = store(&pp, &sk, &files, &mut rng);
        let before = left.files[0].clone();
        let old_first = first.stub.0;
        let old_second = second.stub.0;
        let token = key_up_token(&sk, &next, &mut rng);
        let first_delta = Fr::rand(&mut rng);
        let second_delta = Fr::rand(&mut rng);
        rotate_repository(&pp, &token, first_delta, &mut left, &mut rng);
        rotate_repository(&pp, &token, second_delta, &mut right, &mut rng);
        first.rerandomize(&pp, first_delta, &mut rng);
        second.rerandomize(&pp, second_delta, &mut rng);

        assert_eq!(token.bytes(), 282_704);
        assert_eq!(left.files[0].body_len, before.body_len);
        assert_ne!(left.files[0].body, before.body);
        assert_ne!(left.files[0].mask.c1, before.mask.c1);
        assert_ne!(left.files[0].proof_ct.c1, before.proof_ct.c1);
        assert_ne!(first.stub.0 - old_first, pp.khvc.z_g1 * first_delta);
        assert_ne!(second.stub.0 - old_second, pp.khvc.z_g1 * second_delta);
        assert_eq!(retrieve(&pp, &next, &first, &left, 0).unwrap(), files[0]);
        assert_eq!(retrieve(&pp, &next, &second, &right, 0).unwrap(), files[0]);
    }
}
