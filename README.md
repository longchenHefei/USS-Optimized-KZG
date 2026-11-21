# Scalable USS from KHVC

The Rust prototype separates the client's current KHVC stub and blinding offset
from the server repository. FileUp verifies the old record against that state
before updating it. Each accepted FileUp independently blinds its opening
increment and includes the same blind in the stub update. KeyUp shares a public block-polynomial sequence across
repositories and refreshes independent, RLWE-encrypted file masks. The current
record format has no per-file public tweak. The public KeyUp token does not
include the client's blinding offset.

The paper uses the four games from IACR ePrint 2023/1421 and replaces the
single-transition exposure rule with its transitive closure. It proves OF-PTXT
and gives an adaptive confidentiality reduction for an isolated suffix of
epochs. The reduction permits keys from earlier epochs and transcript exposure
within the suffix, assumes negligible honest verification failure, and does
not cover key exposures within the suffix. The earlier zero-blind FileUp allowed a two-update source test. The current
opening blind masks that difference, with regression coverage. The paper proves
joint off-slot masking and an adaptive continuation for aligned opening states;
the full adaptive payload reduction remains open.

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

For measurements of the shared-polynomial format:

```bash
cargo run --release --bin bench -- --suite uss --smoke --out experiments/results/uss_shared_smoke.csv
```

The 4 KiB smoke results are in `experiments/results/uss_shared_smoke.csv`.
The encoded size is 282,704 bytes for the reusable epoch token plus a 32-byte
opening shift per repository. Earlier USS CSVs used different record formats
and are historical measurements. The revised format has not been measured at
gigabyte scale.
