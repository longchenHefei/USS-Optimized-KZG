//! Metric 3.4 smoke bench: 100 attribute commitments.
use ark_bls12_381::Fr;
use ark_ec::CurveGroup;
use ark_ff::UniformRand;
use ark_serialize::CanonicalSerialize;
use ark_std::rand::{rngs::StdRng, SeedableRng};
use std::time::Instant;
use uss_khvc::khvc::Khvc;

fn main() {
    let mut rng = StdRng::seed_from_u64(0x34_00_03_04);
    // n=128 >= 100 attributes (power-of-two evaluation domain)
    let n = 128usize;
    let khvc = Khvc::setup(n, &mut rng);

    // Exactly 100 attribute values; pad remaining slots with zero
    let mut values: Vec<Fr> = (0..100).map(|_| Fr::rand(&mut rng)).collect();
    values.resize(n, Fr::from(0u64));
    let r = Fr::rand(&mut rng);

    // Warm-up
    let _ = khvc.commit(&values, r);

    let rounds = 50usize;
    let t0 = Instant::now();
    let mut c = khvc.commit(&values, r);
    for _ in 1..rounds {
        c = khvc.commit(&values, r);
    }
    let elapsed = t0.elapsed();
    let avg_ms = elapsed.as_secs_f64() * 1000.0 / rounds as f64;

    let mut bytes = Vec::new();
    c.0.into_affine().serialize_compressed(&mut bytes).unwrap();

    println!("metric=3.4");
    println!("n_domain={}", n);
    println!("num_attributes=100");
    println!("commit_size_bytes={}", bytes.len());
    println!("commit_size_kb={:.4}", bytes.len() as f64 / 1024.0);
    println!("commit_avg_ms={:.4}", avg_ms);
    println!("commit_avg_s={:.6}", avg_ms / 1000.0);
    println!("rounds={}", rounds);
    println!("limit_size_kb=2");
    println!("limit_time_s=5");
    println!("pass_size={}", (bytes.len() as f64 / 1024.0) <= 2.0);
    println!("pass_time={}", (avg_ms / 1000.0) <= 5.0);

    let pi = khvc.open(0, &values, r);
    assert!(khvc.verify(c, values[0], 0, pi));
    println!("open_verify_ok=true");
}
