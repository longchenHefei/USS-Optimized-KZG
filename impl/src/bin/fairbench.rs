//! Fair cost accounting for the same-stack UE comparison.
//!
//! Problem with `bench --suite ue`: RISE and SHINE0 are single-group-element
//! schemes. They are measured on ONE element (48 B / 32 B), once, and printed
//! in the same "ms" columns as the 4 KB / 1 MB / 100 MB symmetric baselines.
//! A single-element latency placed next to a bulk throughput is not a
//! comparison: the table implies RISE (0.29 ms) is ~1.8x slower than
//! AES-GCM (0.16 ms), while per plaintext byte the gap is ~150x at 4 KB and
//! ~1400x at bulk sizes.
//!
//! This harness fixes the accounting:
//!   1. every scheme is measured at the SAME plaintext sizes;
//!   2. long-message algebraic UE is measured block-wise (RISE: 48 B blocks,
//!      SHINE0: 32 B blocks) so the bulk cost is real, not extrapolated;
//!   3. every cell is repeated and reported as median with a p10/p90 spread;
//!   4. the CSV carries normalised columns (ns per plaintext byte, MB/s)
//!      and the ciphertext size, so expansion and communication are visible.
//!
//! Usage: cargo run --release --bin fairbench -- [--quick]

use ark_bls12_381::{Fr, G1Projective};
use ark_ec::PrimeGroup;
use curve25519_dalek::constants::RISTRETTO_BASEPOINT_POINT;
use curve25519_dalek::scalar::Scalar;
use std::time::Instant;
use uss_khvc::rise::{rise_dec, rise_enc, rise_next, rise_upd, RiseCt, RiseSk};
use uss_khvc::ue::{
    aes_envelope_dec, aes_envelope_enc, bench_rng, khprf_ct_bytes, khprf_dec, khprf_enc,
    khprf_upd, nested_dec, nested_enc, nested_keygen, nested_upd, random_message, shine_dec,
    shine_enc, shine_keygen, shine_next, shine_upd, ShineCt,
};

const RISE_BLOCK: usize = 48; // one BLS12-381 G1 element
const SHINE_BLOCK: usize = 32; // one Ristretto point

fn ms(d: std::time::Duration) -> f64 {
    d.as_secs_f64() * 1000.0
}

fn record(
    rows: &mut Vec<String>,
    scheme: &str,
    size: usize,
    op: &str,
    mut samples: Vec<f64>,
    ct_bytes: usize,
) {
    samples.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = samples.len();
    let med = samples[n / 2];
    let p10 = samples[((n as f64 * 0.1) as usize).min(n - 1)];
    let p90 = samples[((n as f64 * 0.9) as usize).min(n - 1)];
    let (ns_per_byte, mbps) = if size > 0 && med > 0.0 {
        (med * 1e6 / size as f64, (size as f64 / 1e6) / (med / 1e3))
    } else {
        (-1.0, -1.0)
    };
    rows.push(format!(
        "{scheme},{size},{op},{n},{med:.4},{p10:.4},{p90:.4},{ct_bytes},{ns_per_byte:.2},{mbps:.4}"
    ));
}

fn main() {
    let quick = std::env::args().any(|a| a == "--quick");
    let mut rows: Vec<String> = vec![
        "scheme,size,op,reps,median_ms,p10_ms,p90_ms,ct_bytes,ns_per_byte,MBps".to_string(),
    ];

    let big = 1024 * 1024;
    let sym_sizes: Vec<usize> = if quick {
        vec![4 * 1024, 64 * 1024]
    } else {
        vec![4 * 1024, 64 * 1024, big]
    };

    // ---------- Single-element latency (Table 2a: seeds / openings) ----------
    {
        let mut rng = bench_rng(200);
        let (sk, pk) = RiseSk::keygen(&mut rng);
        let (sk2, _) = RiseSk::keygen(&mut rng);
        let tk = rise_next(sk, sk2);
        let m = G1Projective::generator() * Fr::from(7u64);
        for _ in 0..3 {
            let ct = rise_enc(pk, m, &mut rng);
            let _ = rise_dec(sk, ct);
            let _ = rise_upd(tk, ct, &mut rng);
        }
        let reps = 21;
        let mut encs = Vec::new();
        let mut upds = Vec::new();
        let mut decs = Vec::new();
        for _ in 0..reps {
            let t = Instant::now();
            let ct = rise_enc(pk, m, &mut rng);
            encs.push(ms(t.elapsed()));
            let t = Instant::now();
            let ct2 = rise_upd(tk, ct, &mut rng);
            upds.push(ms(t.elapsed()));
            let t = Instant::now();
            let _ = rise_dec(sk2, ct2);
            decs.push(ms(t.elapsed()));
        }
        record(&mut rows, "RISE-G1", 48, "enc", encs, 96);
        record(&mut rows, "RISE-G1", 48, "upd", upds, 96);
        record(&mut rows, "RISE-G1", 48, "dec", decs, 96);

        let mut rng = bench_rng(2001);
        let sks = shine_keygen(&mut rng);
        let sks2 = shine_keygen(&mut rng);
        let d = shine_next(sks, sks2);
        let mp = Scalar::from(9u64) * RISTRETTO_BASEPOINT_POINT;
        for _ in 0..3 {
            let ct = shine_enc(sks, mp);
            let _ = shine_dec(sks, ct);
            let _ = shine_upd(d, ct);
        }
        let mut encs = Vec::new();
        let mut upds = Vec::new();
        let mut decs = Vec::new();
        for _ in 0..reps {
            let t = Instant::now();
            let ct = shine_enc(sks, mp);
            encs.push(ms(t.elapsed()));
            let t = Instant::now();
            let ct2 = shine_upd(d, ct);
            upds.push(ms(t.elapsed()));
            let t = Instant::now();
            let _ = shine_dec(sks2, ct2);
            decs.push(ms(t.elapsed()));
        }
        record(&mut rows, "SHINE0", 32, "enc", encs, 32);
        record(&mut rows, "SHINE0", 32, "upd", upds, 32);
        record(&mut rows, "SHINE0", 32, "dec", decs, 32);
    }

    // ---------- RISE (one G1 element per 48 B block) ----------
    let mut rng = bench_rng(201);
    let (sk, pk) = RiseSk::keygen(&mut rng);
    let (sk2, _) = RiseSk::keygen(&mut rng);
    let tk = rise_next(sk, sk2);
    let m = G1Projective::generator() * Fr::from(7u64);
    let rise_sizes: Vec<usize> = if quick {
        vec![4 * 1024, 64 * 1024]
    } else {
        vec![4 * 1024, 64 * 1024, big]
    };
    for &size in &rise_sizes {
        let blocks = size.div_ceil(RISE_BLOCK);
        let reps = if size >= big { 3 } else if size >= 64 * 1024 { 5 } else { 15 };
        eprintln!("RISE size={size} blocks={blocks} reps={reps}");
        {
            let mut cts: Vec<RiseCt> = Vec::with_capacity(blocks);
            for _ in 0..blocks {
                cts.push(rise_enc(pk, m, &mut rng));
            }
            let nxt: Vec<RiseCt> = cts.iter().map(|ct| rise_upd(tk, *ct, &mut rng)).collect();
            for ct in &nxt {
                let _ = rise_dec(sk2, *ct);
            }
        }
        let mut encs = Vec::new();
        let mut upds = Vec::new();
        let mut decs = Vec::new();
        for _ in 0..reps {
            let t = Instant::now();
            let mut cts: Vec<RiseCt> = Vec::with_capacity(blocks);
            for _ in 0..blocks {
                cts.push(rise_enc(pk, m, &mut rng));
            }
            encs.push(ms(t.elapsed()));
            let t = Instant::now();
            let mut nxt: Vec<RiseCt> = Vec::with_capacity(blocks);
            for ct in &cts {
                nxt.push(rise_upd(tk, *ct, &mut rng));
            }
            upds.push(ms(t.elapsed()));
            let t = Instant::now();
            for ct in &nxt {
                let _ = rise_dec(sk2, *ct);
            }
            decs.push(ms(t.elapsed()));
        }
        record(&mut rows, "RISE-G1", size, "enc", encs, 2 * 48 * blocks);
        record(&mut rows, "RISE-G1", size, "upd", upds, 2 * 48 * blocks);
        record(&mut rows, "RISE-G1", size, "dec", decs, 2 * 48 * blocks);
    }

    // ---------- SHINE0 (one Ristretto point per 32 B block) ----------
    let mut rng = bench_rng(202);
    let sks = shine_keygen(&mut rng);
    let sks2 = shine_keygen(&mut rng);
    let d = shine_next(sks, sks2);
    let mp = Scalar::from(9u64) * RISTRETTO_BASEPOINT_POINT;
    let shine_sizes: Vec<usize> = if quick {
        vec![4 * 1024, 64 * 1024]
    } else {
        vec![4 * 1024, 64 * 1024, big]
    };
    for &size in &shine_sizes {
        let blocks = size.div_ceil(SHINE_BLOCK);
        let reps = if size >= big { 3 } else if size >= 64 * 1024 { 5 } else { 15 };
        eprintln!("SHINE0 size={size} blocks={blocks} reps={reps}");
        {
            let mut cts: Vec<ShineCt> = Vec::with_capacity(blocks);
            for _ in 0..blocks {
                cts.push(shine_enc(sks, mp));
            }
            let nxt: Vec<ShineCt> = cts.iter().map(|ct| shine_upd(d, *ct)).collect();
            for ct in &nxt {
                let _ = shine_dec(sks2, *ct);
            }
        }
        let mut encs = Vec::new();
        let mut upds = Vec::new();
        let mut decs = Vec::new();
        for _ in 0..reps {
            let t = Instant::now();
            let mut cts: Vec<ShineCt> = Vec::with_capacity(blocks);
            for _ in 0..blocks {
                cts.push(shine_enc(sks, mp));
            }
            encs.push(ms(t.elapsed()));
            let t = Instant::now();
            let mut nxt: Vec<ShineCt> = Vec::with_capacity(blocks);
            for ct in &cts {
                nxt.push(shine_upd(d, *ct));
            }
            upds.push(ms(t.elapsed()));
            let t = Instant::now();
            for ct in &nxt {
                let _ = shine_dec(sks2, *ct);
            }
            decs.push(ms(t.elapsed()));
        }
        record(&mut rows, "SHINE0", size, "enc", encs, 32 * blocks);
        record(&mut rows, "SHINE0", size, "upd", upds, 32 * blocks);
        record(&mut rows, "SHINE0", size, "dec", decs, 32 * blocks);
    }

    // ---------- Symmetric baselines at floating-point sizes ----------
    for &size in &sym_sizes {
        let msg = random_message(size, 7);
        let reps = if size >= big { 7 } else { 21 };
        eprintln!("AES/Nested/KHPRF size={size} reps={reps}");
        {
            let mut rng = bench_rng(99);
            let (ct, _) = aes_envelope_enc(&msg, &mut rng);
            let _ = aes_envelope_dec(&ct);
        }

        let mut rng = bench_rng(101);
        let mut encs = Vec::new();
        let mut decs = Vec::new();
        let mut reenc = Vec::new();
        let mut ctb = 0usize;
        for _ in 0..reps {
            let t = Instant::now();
            let (ct, _) = aes_envelope_enc(&msg, &mut rng);
            encs.push(ms(t.elapsed()));
            ctb = ct.body.len() + 12 + 32;
            let t = Instant::now();
            let _ = aes_envelope_dec(&ct);
            decs.push(ms(t.elapsed()));
            let t = Instant::now();
            let (ct2, _) = aes_envelope_enc(&msg, &mut rng);
            reenc.push(ms(t.elapsed()));
            ctb = ct2.body.len() + 12 + 32;
        }
        record(&mut rows, "AES-GCM", size, "enc", encs, ctb);
        record(&mut rows, "AES-GCM", size, "dec", decs, size);
        record(&mut rows, "AES-GCM", size, "client_reenc", reenc, ctb);

        let mut rng = bench_rng(102);
        let k0 = nested_keygen(&mut rng);
        let k1 = nested_keygen(&mut rng);
        let mut encs = Vec::new();
        let mut upds = Vec::new();
        let mut decs = Vec::new();
        let mut nctb = 0usize;
        for _ in 0..reps {
            let t = Instant::now();
            let nct = nested_enc(&k0, &msg, &mut rng);
            encs.push(ms(t.elapsed()));
            nctb = nct.blob.len();
            let t = Instant::now();
            let nct1 = nested_upd(&k1, &nct, &mut rng);
            upds.push(ms(t.elapsed()));
            let t = Instant::now();
            let pt = nested_dec(&[k0, k1], &nct1);
            decs.push(ms(t.elapsed()));
            assert_eq!(pt, msg);
        }
        record(&mut rows, "Nested-AES", size, "enc", encs, nctb);
        record(&mut rows, "Nested-AES", size, "upd", upds, nctb);
        record(&mut rows, "Nested-AES", size, "dec", decs, nctb);

        let mut rng = bench_rng(103);
        let sk = shine_keygen(&mut rng);
        let sk2 = shine_keygen(&mut rng);
        let dd = shine_next(sk, sk2);
        let mut encs = Vec::new();
        let mut upds = Vec::new();
        let mut decs = Vec::new();
        let mut kctb = 0usize;
        for _ in 0..reps {
            let t = Instant::now();
            let ct = khprf_enc(sk, &msg);
            encs.push(ms(t.elapsed()));
            kctb = khprf_ct_bytes(&ct);
            let t = Instant::now();
            let ct2 = khprf_upd(dd, &ct);
            upds.push(ms(t.elapsed()));
            let t = Instant::now();
            let pt = khprf_dec(sk2, &ct2);
            decs.push(ms(t.elapsed()));
            assert_eq!(&pt[..msg.len()], &msg[..]);
        }
        record(&mut rows, "ReCrypt-KHPRF", size, "enc", encs, kctb);
        record(&mut rows, "ReCrypt-KHPRF", size, "upd", upds, kctb);
        record(&mut rows, "ReCrypt-KHPRF", size, "dec", decs, kctb);
    }

    println!("{}", rows.join("\n"));
}
