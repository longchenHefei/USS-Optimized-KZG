#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/experiments/results"
mkdir -p "$OUT"
cd "$ROOT"
echo "machine:" | tee "$OUT/run.log"
cat "$ROOT/experiments/machine.txt" | tee -a "$OUT/run.log"
echo "==== HVC ====" | tee -a "$OUT/run.log"
cargo run --release -p uss-khvc --bin bench -- --suite hvc --out "$OUT/hvc.csv" 2>>"$OUT/run.log"
echo "==== UE ====" | tee -a "$OUT/run.log"
cargo run --release -p uss-khvc --bin bench -- --suite ue --out "$OUT/ue.csv" 2>>"$OUT/run.log"
echo "==== USS ====" | tee -a "$OUT/run.log"
cargo run --release -p uss-khvc --bin bench -- --suite uss --out "$OUT/uss.csv" 2>>"$OUT/run.log"
python3 - <<'PY'
from pathlib import Path
root = Path(__file__).resolve().parent / "results" if False else Path.cwd()
PY
# merge CSVs
{
  head -n 1 "$OUT/hvc.csv"
  tail -n +2 "$OUT/hvc.csv"
  tail -n +2 "$OUT/ue.csv"
  tail -n +2 "$OUT/uss.csv"
} > "$OUT/all.csv"
echo "wrote $OUT/all.csv"
