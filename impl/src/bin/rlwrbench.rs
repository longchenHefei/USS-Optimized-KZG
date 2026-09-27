//! Single-thread Ring-LWR KH-PRF UE, same accounting as fairbench.
//!
//! `a_i` is treated as a public parameter already in the NTT domain, so Enc/Upd/Dec
//! time is the ring product plus rounding. Setup time and the byte size of `{a_i}`
//! are reported separately: they grow with the number of blocks.
//!
//! Usage: cargo run --release --bin rlwrbench -- [--quick]

use rand::rngs::StdRng;
use rand::SeedableRng;
use std::time::Instant;
use uss_khvc::rlwr_ue::{
    ct_bytes, dec, enc, keygen, max_abs_rounding_error, setup, token, token_bytes, upd, BLOCK,
    N, P,
};

fn ms(d: std::time::Duration) -> f64 {
    d.as_secs_f64() * 1000.0
}

fn median(mut xs: Vec<f64>) -> (f64, f64, f64) {
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = xs.len();
    let med = xs[n / 2];
    let p10 = xs[((n as f64 * 0.1) as usize).min(n - 1)];
    let p90 = xs[((n as f64 * 0.9) as usize).min(n - 1)];
    (med, p10, p90)
}

fn main() {
    let quick = std::env::args().any(|a| a == "--quick");
    let sizes: Vec<usize> = if quick {
        vec![4 * 1024]
    } else {
        vec![4 * 1024, 64 * 1024, 1024 * 1024]
    };
    let mut rows = vec![
        "scheme,size,op,reps,median_ms,p10_ms,p90_ms,ct_bytes,ns_per_byte,MBps".to_string(),
    ];

    println!(
        "Ring-LWR UE  N={N}  q=998244353  p={P}  payload={BLOCK} B/block  noise_bits=8  msg_bits=8  token={} B",
        token_bytes()
    );
    println!("single-threaded reference NTT; a_i NTT excluded from Enc/Upd/Dec");

    let mut rng = StdRng::seed_from_u64(20260324);
    println!("\nrounding error |e|_inf (8-bit padding fails once |e| >= 128)");
    for (blocks, t) in [(2usize, 1usize), (2, 32), (2, 128), (2, 255), (64, 64), (64, 128)] {
        let e = max_abs_rounding_error(blocks, t, &mut rng);
        println!("  blocks={blocks:<3} t={t:<3} max |e|={e}");
    }

    for size in sizes {
        let nblocks = size.div_ceil(BLOCK);
        let reps = if size <= 4 * 1024 {
            21
        } else if size <= 64 * 1024 {
            7
        } else {
            5
        };
        let msg = vec![0x5Au8; size];

        let t0 = Instant::now();
        let pp = setup(nblocks, &mut rng);
        let setup_ms = ms(t0.elapsed());
        let a_bytes = nblocks * N * 4;
        println!(
            "\nsize={size} B  blocks={nblocks}  setup={setup_ms:.1} ms  |a_i|={a_bytes} B ({:.2}× plaintext)",
            a_bytes as f64 / size as f64
        );

        let sk = keygen(&mut rng);
        let sk2 = keygen(&mut rng);
        let delta = token(&sk, &sk2);
        let ct0 = enc(&pp, &sk, &msg);
        let ct_len = ct_bytes(&ct0);
        let ct1 = upd(&pp, &delta, &ct0);
        let pt = dec(&pp, &sk2, &ct1);
        assert_eq!(pt, msg, "roundtrip failed at {size}");

        let mut enc_s = Vec::new();
        let mut upd_s = Vec::new();
        let mut dec_s = Vec::new();
        // warm-up
        let _ = enc(&pp, &sk, &msg);
        let warm = enc(&pp, &sk, &msg);
        let _ = upd(&pp, &delta, &warm);
        let _ = dec(&pp, &sk, &warm);

        for _ in 0..reps {
            let t = Instant::now();
            let ct = enc(&pp, &sk, &msg);
            enc_s.push(ms(t.elapsed()));
            let t = Instant::now();
            let ct2 = upd(&pp, &delta, &ct);
            upd_s.push(ms(t.elapsed()));
            let t = Instant::now();
            let _ = dec(&pp, &sk2, &ct2);
            dec_s.push(ms(t.elapsed()));
        }

        for (op, samples) in [("enc", enc_s), ("upd", upd_s), ("dec", dec_s)] {
            let (med, p10, p90) = median(samples.clone());
            let ns = med * 1e6 / size as f64;
            let mbps = (size as f64 / 1e6) / (med / 1e3);
            println!(
                "  {op:3}  median {med:.3} ms  p10 {p10:.3}  p90 {p90:.3}  {ns:.1} ns/B  {mbps:.2} MB/s  |ct|={ct_len} ({:.2}×)",
                ct_len as f64 / size as f64
            );
            rows.push(format!(
                "Ring-LWR,{size},{op},{reps},{med:.4},{p10:.4},{p90:.4},{ct_len},{ns:.2},{mbps:.4}"
            ));
        }
        rows.push(format!(
            "Ring-LWR,{size},setup,1,{setup_ms:.4},{setup_ms:.4},{setup_ms:.4},{a_bytes},-1,-1"
        ));
    }

    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../experiments/results/rlwr_ue.csv");
    std::fs::write(&path, rows.join("\n") + "\n").expect("write csv");
    println!("\nwrote {}", path.display());
}
