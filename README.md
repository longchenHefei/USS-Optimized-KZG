# Scalable USS from KHVC

Implementation and English paper for an optimized Updatable Secure Storage
scheme built from a KZG-based homomorphic vector commitment (KHVC), with a
Ring-LWR key-homomorphic PRF as the file-body DEM.

Repository: <https://github.com/longchenHefei/USS-Optimized-KZG>

```
uss-khvc/
  impl/           Rust crate (KHVC, CF-HVC, RISE, UE baselines, Ring-LWR DEM, USS)
  paper/          English article (main.tex / main.pdf)
  experiments/    bench CSVs and machine notes
```

## Build and test

```bash
cargo test --release -p uss-khvc
cargo run --release --bin fairbench
cargo run --release --bin rlwrbench
cargo run --release --bin bench -- --suite uss --smoke
```

Large end-to-end runs (memory-aware for 1 GB bodies):

```bash
cargo run --release --bin bench -- --suite hvc --out experiments/results/hvc_sparse.csv
cargo run --release --bin bench -- --suite ussmid --out experiments/results/uss_mid.csv
cargo run --release --bin bench -- --suite uss1g  --out experiments/results/uss_e2e.csv
```

CSV output: `experiments/results/{hvc,ue,uss,fair_ue,rlwr_ue,...}.csv`.

## Headline numbers (Apple M5 Pro)

- KHVC verification CRS at $n=1024$: 48.3 KB vs CF-HVC 50.5 MB
- One-hot FileUp openings at $n=1024$, $k=1$: **194 ms** (was 1.62 s with dense FK20)
- Ring-LWR DEM at 1 MB: ≈32.6 MB/s Enc/Upd, $2\times$ ciphertext, 4096-byte token
- End-to-end USS KeyUp of $4\times 1$ GB: **254 s**, **4240-byte** token (cached $a_b$ NTTs; was 404 s with per-block rebuild)
- End-to-end USS KeyUp of $8\times 100$ MB: **58.9 s** / 4240 B
- End-to-end USS KeyUp of $64\times 1$ MB: **2.22 s** / 4240 B
