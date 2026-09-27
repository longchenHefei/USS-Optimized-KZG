//! Micro-benchmarks and end-to-end USS timings. Writes CSV to stdout or a file.

use clap::Parser;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;
use uss_khvc::cf_hvc::CfHvc;
use uss_khvc::khvc::Khvc;
use uss_khvc::rise::{rise_dec, rise_enc, rise_next, rise_upd, RiseSk};
use uss_khvc::ue::{
    aes_envelope_dec, aes_envelope_enc, bench_rng, khprf_ct_bytes, khprf_dec, khprf_enc, khprf_upd,
    nested_dec, nested_enc, nested_keygen, nested_upd, random_message, shine_dec, shine_enc,
    shine_keygen, shine_next, shine_upd,
};
use uss_khvc::uss::{file_up, key_up, keygen, pargen, repo_bytes, retrieve, seal, stage_file, store};
use ark_bls12_381::{Fr, G1Projective};
use ark_ec::PrimeGroup;
use ark_ff::UniformRand;
use curve25519_dalek::constants::RISTRETTO_BASEPOINT_POINT;

#[derive(Parser, Debug)]
struct Args {
    #[arg(long, default_value = "all")]
    suite: String,
    #[arg(long)]
    out: Option<PathBuf>,
    /// Skip 100MB+ algebraic UE and 1GB USS (faster smoke run).
    #[arg(long)]
    smoke: bool,
}

fn ms(d: std::time::Duration) -> f64 {
    d.as_secs_f64() * 1000.0
}

fn main() {
    let args = Args::parse();
    let mut rows: Vec<String> = vec![
        "suite,scheme,n,k_or_size,op,ms,bytes,notes".to_string(),
    ];
    match args.suite.as_str() {
        "hvc" => bench_hvc(&mut rows, args.smoke),
        "ue" => bench_ue(&mut rows, args.smoke),
        "uss" => bench_uss(&mut rows, args.smoke),
        "uss1g" => run_uss_lowmem(&mut rows, 4, 1024 * 1024 * 1024),
        "ussmid" => {
            run_uss(&mut rows, 64, 1024 * 1024, "small");
            run_uss(&mut rows, 8, 100 * 1024 * 1024, "large");
        }
        _ => {
            bench_hvc(&mut rows, args.smoke);
            bench_ue(&mut rows, args.smoke);
            bench_uss(&mut rows, args.smoke);
        }
    }
    let text = rows.join("\n") + "\n";
    if let Some(p) = args.out {
        if let Some(dir) = p.parent() {
            fs::create_dir_all(dir).ok();
        }
        File::create(&p).unwrap().write_all(text.as_bytes()).unwrap();
        eprintln!("wrote {}", p.display());
    } else {
        print!("{}", text);
    }
}

fn bench_hvc(rows: &mut Vec<String>, smoke: bool) {
    let mut ns = vec![64usize, 256, 1024];
    if !smoke {
        ns.extend_from_slice(&[4096, 16384]);
    }
    for n in ns {
        eprintln!("HVC n={n}");
        let mut rng = ark_std::rand::rngs::StdRng::seed_from_u64(0xC0FFEE);
        use rand::SeedableRng;
        let t0 = Instant::now();
        let khvc = Khvc::setup(n, &mut rng);
        let setup_ms = ms(t0.elapsed());
        rows.push(format!(
            "hvc,KHVC,{n},0,setup,{setup_ms:.4},{},",
            khvc.crs_bytes()
        ));

        let values: Vec<Fr> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let r = Fr::rand(&mut rng);
        let t0 = Instant::now();
        let c = khvc.commit(&values, r);
        rows.push(format!("hvc,KHVC,{n},0,commit,{:.4},48,", ms(t0.elapsed())));

        let t0 = Instant::now();
        let pi = khvc.open(0, &values, r);
        let open_ms = ms(t0.elapsed());
        let t0 = Instant::now();
        assert!(khvc.verify(c, values[0], 0, pi));
        rows.push(format!("hvc,KHVC,{n},0,open,{open_ms:.4},48,"));
        rows.push(format!("hvc,KHVC,{n},0,verify,{:.4},0,", ms(t0.elapsed())));

        if n <= 1024 {
            let t0 = Instant::now();
            let _ = khvc.open_all_naive(&values, r);
            rows.push(format!(
                "hvc,KHVC,{n},0,open_all_naive,{:.4},{},",
                ms(t0.elapsed()),
                48 * n
            ));
        }
        let t0 = Instant::now();
        let all = khvc.open_all(&values, r);
        rows.push(format!(
            "hvc,KHVC,{n},0,open_all_fk20,{:.4},{},",
            ms(t0.elapsed()),
            48 * n
        ));
        assert!(khvc.verify(c, values[1], 1, all[1]));

        // sparse FileUp-style: k changed positions, still all-open of a sparse poly
        for k in [1usize, 8, 64] {
            if k > n {
                continue;
            }
            let mut delta = vec![Fr::from(0u64); n];
            for i in 0..k {
                delta[i] = Fr::rand(&mut rng);
            }
            let rd = Fr::rand(&mut rng);
            let t0 = Instant::now();
            let _ = khvc.open_all(&delta, rd);
            rows.push(format!(
                "hvc,KHVC,{n},{k},fileup_all_open,{:.4},{},k={k}",
                ms(t0.elapsed()),
                48 * n
            ));
        }
        let t0 = Instant::now();
        let r0 = Fr::rand(&mut rng);
        let _pre = khvc
            .z_openings
            .iter()
            .map(|z| *z * r0)
            .collect::<Vec<_>>();
        rows.push(format!(
            "hvc,KHVC,{n},0,keyup_precomputed,{:.4},{},",
            ms(t0.elapsed()),
            48 * n
        ));

        if n <= 1024 {
            let t0 = Instant::now();
            let cf = CfHvc::setup(n, &mut rng);
            rows.push(format!(
                "hvc,CF-HVC,{n},0,setup,{:.4},{},",
                ms(t0.elapsed()),
                cf.crs_bytes()
            ));
            let t0 = Instant::now();
            let cc = cf.commit(&values, r);
            rows.push(format!("hvc,CF-HVC,{n},0,commit,{:.4},48,", ms(t0.elapsed())));
            let t0 = Instant::now();
            let op = cf.open(0, &values, r);
            rows.push(format!("hvc,CF-HVC,{n},0,open,{:.4},48,", ms(t0.elapsed())));
            let t0 = Instant::now();
            assert!(cf.verify(cc, values[0], 0, op));
            rows.push(format!("hvc,CF-HVC,{n},0,verify,{:.4},0,", ms(t0.elapsed())));
            let t0 = Instant::now();
            let _ = cf.open_all(&values, r);
            rows.push(format!(
                "hvc,CF-HVC,{n},0,open_all,{:.4},{},",
                ms(t0.elapsed()),
                48 * n
            ));
        }
    }
}

fn bench_ue(rows: &mut Vec<String>, smoke: bool) {
    let mut sizes = vec![4 * 1024usize, 1024 * 1024];
    if !smoke {
        sizes.push(100 * 1024 * 1024);
    }
    for size in sizes {
        eprintln!("UE size={size}");
        let msg = random_message(size, 7);
        let mut rng = bench_rng(11);

        // AES-GCM envelope
        let t0 = Instant::now();
        let (ct, _) = aes_envelope_enc(&msg, &mut rng);
        let enc_ms = ms(t0.elapsed());
        let t0 = Instant::now();
        let _ = aes_envelope_dec(&ct);
        rows.push(format!(
            "ue,AES-GCM,{0},{size},enc,{enc_ms:.4},{1},",
            1,
            ct.body.len() + 12 + 32
        ));
        rows.push(format!(
            "ue,AES-GCM,1,{size},dec,{:.4},{},",
            ms(t0.elapsed()),
            size
        ));
        let t0 = Instant::now();
        let (ct2, _) = aes_envelope_enc(&msg, &mut rng);
        rows.push(format!(
            "ue,AES-GCM,1,{size},reenc_naive,{:.4},{},client reencrypt",
            ms(t0.elapsed()),
            ct2.body.len()
        ));

        // Nested AES
        let mut rng = bench_rng(12);
        let k0 = nested_keygen(&mut rng);
        let k1 = nested_keygen(&mut rng);
        let t0 = Instant::now();
        let nct = nested_enc(&k0, &msg, &mut rng);
        rows.push(format!(
            "ue,Nested-AES,1,{size},enc,{:.4},{},",
            ms(t0.elapsed()),
            nct.blob.len()
        ));
        let t0 = Instant::now();
        let nct1 = nested_upd(&k1, &nct, &mut rng);
        rows.push(format!(
            "ue,Nested-AES,1,{size},upd,{:.4},{},epoch=2",
            ms(t0.elapsed()),
            nct1.blob.len()
        ));
        let t0 = Instant::now();
        let pt = nested_dec(&[k0, k1], &nct1);
        assert_eq!(pt, msg);
        rows.push(format!(
            "ue,Nested-AES,1,{size},dec,{:.4},{},epochs=2",
            ms(t0.elapsed()),
            size
        ));

        // RISE on a single G1 message (short) — report per-op, independent of file size
        if size == 4 * 1024 {
            let mut rng = ark_std::rand::rngs::StdRng::seed_from_u64(13);
            use rand::SeedableRng;
            let (sk, pk) = RiseSk::keygen(&mut rng);
            let m = G1Projective::generator() * Fr::rand(&mut rng);
            let t0 = Instant::now();
            let ct = rise_enc(pk, m, &mut rng);
            rows.push(format!("ue,RISE-G1,1,48,enc,{:.4},96,", ms(t0.elapsed())));
            let t0 = Instant::now();
            assert_eq!(rise_dec(sk, ct), m);
            rows.push(format!("ue,RISE-G1,1,48,dec,{:.4},48,", ms(t0.elapsed())));
            let (sk2, _) = RiseSk::keygen(&mut rng);
            let tk = rise_next(sk, sk2);
            let t0 = Instant::now();
            let ct2 = rise_upd(tk, ct, &mut rng);
            rows.push(format!("ue,RISE-G1,1,48,upd,{:.4},96,", ms(t0.elapsed())));
            assert_eq!(rise_dec(sk2, ct2), m);

            let sks = shine_keygen(&mut rng);
            let sks2 = shine_keygen(&mut rng);
            let mp = shine_keygen(&mut rng).0 * RISTRETTO_BASEPOINT_POINT;
            let t0 = Instant::now();
            let sct = shine_enc(sks, mp);
            rows.push(format!("ue,SHINE0,1,32,enc,{:.4},32,", ms(t0.elapsed())));
            let t0 = Instant::now();
            let _ = shine_dec(sks, sct);
            rows.push(format!("ue,SHINE0,1,32,dec,{:.4},32,", ms(t0.elapsed())));
            let d = shine_next(sks, sks2);
            let t0 = Instant::now();
            let _ = shine_upd(d, sct);
            rows.push(format!("ue,SHINE0,1,32,upd,{:.4},32,", ms(t0.elapsed())));
        }

        // KH-PRF / ReCrypt-style: skip 100MB in smoke; at 100MB still run once
        let run_khprf = size <= 1024 * 1024 || !smoke;
        if run_khprf && size <= 100 * 1024 * 1024 {
            let mut rng = bench_rng(14);
            let sk = shine_keygen(&mut rng);
            let sk2 = shine_keygen(&mut rng);
            // For 100MB this is slow (one scalarmul per 32B). Still measure Enc; Upd similar.
            if size >= 100 * 1024 * 1024 {
                // subsample: encrypt 1MB and scale comment
                let msg_s = random_message(1024 * 1024, 15);
                let t0 = Instant::now();
                let ct = khprf_enc(sk, &msg_s);
                let enc_ms = ms(t0.elapsed());
                rows.push(format!(
                    "ue,ReCrypt-KHPRF,1,{size},enc_est,{:.4},{},measured 1MB then x100",
                    enc_ms * 100.0,
                    khprf_ct_bytes(&ct) * 100
                ));
                let d = shine_next(sk, sk2);
                let t0 = Instant::now();
                let _ = khprf_upd(d, &ct);
                rows.push(format!(
                    "ue,ReCrypt-KHPRF,1,{size},upd_est,{:.4},0,measured 1MB then x100",
                    ms(t0.elapsed()) * 100.0
                ));
            } else {
                let t0 = Instant::now();
                let ct = khprf_enc(sk, &msg);
                rows.push(format!(
                    "ue,ReCrypt-KHPRF,1,{size},enc,{:.4},{},",
                    ms(t0.elapsed()),
                    khprf_ct_bytes(&ct)
                ));
                let t0 = Instant::now();
                let pt = khprf_dec(sk, &ct);
                assert_eq!(&pt[..msg.len()], &msg[..]);
                rows.push(format!(
                    "ue,ReCrypt-KHPRF,1,{size},dec,{:.4},{},",
                    ms(t0.elapsed()),
                    size
                ));
                let d = shine_next(sk, sk2);
                let t0 = Instant::now();
                let ct2 = khprf_upd(d, &ct);
                rows.push(format!(
                    "ue,ReCrypt-KHPRF,1,{size},upd,{:.4},{},",
                    ms(t0.elapsed()),
                    khprf_ct_bytes(&ct2)
                ));
                let _ = khprf_dec(sk2, &ct2);
            }
        }
    }
}

fn bench_uss(rows: &mut Vec<String>, smoke: bool) {
    // small files
    for &(n, fsize) in &[(8usize, 4 * 1024), (32, 4 * 1024), (64, 1024 * 1024)] {
        if smoke && fsize > 4 * 1024 {
            continue;
        }
        if smoke && n > 32 {
            continue;
        }
        run_uss(rows, n, fsize, "small");
    }
    if !smoke {
        for &(n, fsize) in &[
            (8usize, 100 * 1024 * 1024),
            (4, 1024 * 1024 * 1024),
        ] {
            run_uss(rows, n, fsize, "large");
        }
    }
}

fn run_uss_lowmem(rows: &mut Vec<String>, n: usize, fsize: usize) {
    eprintln!("USS lowmem n={n} file={fsize}");
    use rand::SeedableRng;
    use sha2::{Digest, Sha256};
    let mut rng = ark_std::rand::rngs::StdRng::seed_from_u64(99);
    let t0 = Instant::now();
    let pp = pargen(n, &mut rng);
    rows.push(format!(
        "uss,KHVC-USS,{n},{fsize},pargen,{:.4},{},lowmem",
        ms(t0.elapsed()),
        pp.khvc.crs_bytes()
    ));
    let sk = keygen(&mut rng);
    let mut digests = Vec::with_capacity(n);
    let mut staged = Vec::with_capacity(n);
    let t0 = Instant::now();
    for i in 0..n {
        let file = random_message(fsize, 1000 + i as u64);
        digests.push(Sha256::digest(&file));
        staged.push(stage_file(&sk.dem, &file, &mut rng));
        eprintln!("  staged file {i} in {:.1}s", t0.elapsed().as_secs_f64());
    }
    let mut repo = seal(&pp, &sk, staged, &mut rng);
    rows.push(format!(
        "uss,KHVC-USS,{n},{fsize},store,{:.4},{},lowmem",
        ms(t0.elapsed()),
        repo_bytes(&repo)
    ));
    let t0 = Instant::now();
    let got = retrieve(&pp, &sk, &repo, 0).expect("rev");
    assert_eq!(Sha256::digest(&got), digests[0]);
    drop(got);
    rows.push(format!(
        "uss,KHVC-USS,{n},{fsize},rev,{:.4},{},lowmem",
        ms(t0.elapsed()),
        fsize
    ));
    let newf = random_message(fsize, 42);
    let new_digest = Sha256::digest(&newf);
    let t0 = Instant::now();
    let tok = file_up(&pp, &sk, &mut repo, 0, &newf, &mut rng);
    drop(newf);
    rows.push(format!(
        "uss,KHVC-USS,{n},{fsize},fileup,{:.4},{tok},lowmem",
        ms(t0.elapsed())
    ));
    let got = retrieve(&pp, &sk, &repo, 0).expect("rev after fileup");
    assert_eq!(Sha256::digest(&got), new_digest);
    drop(got);
    let sk2 = keygen(&mut rng);
    let t0 = Instant::now();
    let tok = key_up(&pp, &sk, &sk2, &mut repo, &mut rng);
    rows.push(format!(
        "uss,KHVC-USS,{n},{fsize},keyup,{:.4},{tok},lowmem",
        ms(t0.elapsed())
    ));
    let got = retrieve(&pp, &sk2, &repo, 0).expect("rev after keyup");
    assert_eq!(Sha256::digest(&got), new_digest);
    assert!(retrieve(&pp, &sk, &repo, 2).is_none());
}

fn run_uss(rows: &mut Vec<String>, n: usize, fsize: usize, note: &str) {
    eprintln!("USS n={n} file={fsize} ({note})");
    use rand::SeedableRng;
    let mut rng = ark_std::rand::rngs::StdRng::seed_from_u64(99);
    let t0 = Instant::now();
    let pp = pargen(n, &mut rng);
    rows.push(format!(
        "uss,KHVC-USS,{n},{fsize},pargen,{:.4},{},{note}",
        ms(t0.elapsed()),
        pp.khvc.crs_bytes()
    ));
    let sk = keygen(&mut rng);
    let files: Vec<Vec<u8>> = (0..n)
        .map(|i| random_message(fsize, 1000 + i as u64))
        .collect();
    let t0 = Instant::now();
    let mut repo = store(&pp, &sk, &files, &mut rng);
    rows.push(format!(
        "uss,KHVC-USS,{n},{fsize},store,{:.4},{},{note}",
        ms(t0.elapsed()),
        repo_bytes(&repo)
    ));
    let t0 = Instant::now();
    let got = retrieve(&pp, &sk, &repo, 0).unwrap();
    assert_eq!(got, files[0]);
    rows.push(format!(
        "uss,KHVC-USS,{n},{fsize},rev,{:.4},{},{note}",
        ms(t0.elapsed()),
        fsize
    ));
    let newf = random_message(fsize, 42);
    let t0 = Instant::now();
    let tok = file_up(&pp, &sk, &mut repo, 0, &newf, &mut rng);
    rows.push(format!(
        "uss,KHVC-USS,{n},{fsize},fileup,{:.4},{tok},{note}",
        ms(t0.elapsed())
    ));
    assert_eq!(retrieve(&pp, &sk, &repo, 0).unwrap(), newf);
    let sk2 = keygen(&mut rng);
    let t0 = Instant::now();
    let tok = key_up(&pp, &sk, &sk2, &mut repo, &mut rng);
    rows.push(format!(
        "uss,KHVC-USS,{n},{fsize},keyup,{:.4},{tok},{note}",
        ms(t0.elapsed())
    ));
    assert_eq!(retrieve(&pp, &sk2, &repo, 0).unwrap(), newf);
}
