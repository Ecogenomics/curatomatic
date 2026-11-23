#!/usr/bin/env bash
# Benchmark the current Rust curatomatic against the stored Python output and the curated releases.
#
# For each dataset (default: r226_ar53 r226_bac120 r232_ar53 r232_bac120) this:
#   1. converts the Python inputs to the Rust format     (scripts/convert_python_inputs.py)
#   2. runs the Rust version, with --min-bs 95.0          (target/release/curatomatic)
#   3. compares it with the stored Python output          (scripts/compare_outputs.py)
#   4. limits the Rust report to taxa placed differently  (scripts/report_placed_differently.py)
#   5. scores the input, Rust and Python against the truth (scripts/compare_truth.py)
#
# Usage (from the repository root):
#   scripts/benchmark.sh <out_dir> [dataset ...]
#
# <out_dir> must not exist. The inputs are read from test_data/ (see docs/testing.md), which is
# never written to. Step 3 needs dendropy: set PY to a Python that has it (default: the conda
# environment `general`). The other steps only need the standard library (python3).
set -euo pipefail

PY=${PY:-/home/aaron/.conda/envs/general/bin/python}

if [[ $# -lt 1 ]]; then
    sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi
out=$1
shift
datasets=("$@")
[[ ${#datasets[@]} -gt 0 ]] || datasets=(r226_ar53 r226_bac120 r232_ar53 r232_bac120)

repo=$(cd "$(dirname "$0")/.." && pwd)
py_in=$repo/test_data/python/input
py_out=$repo/test_data/python/output
truth=$repo/test_data/curated

if [[ -e $out ]]; then
    echo "error: $out already exists" >&2
    exit 1
fi
"$PY" -c 'import dendropy' || { echo "error: $PY can't import dendropy (set PY)" >&2; exit 1; }

(cd "$repo" && cargo build --release)
mkdir -p "$out"
out=$(cd "$out" && pwd)
commit=$(git -C "$repo" rev-parse --short HEAD)
git -C "$repo" diff --quiet HEAD || commit="$commit-dirty"
echo "$commit" > "$out/commit.txt"

for d in "${datasets[@]}"; do
    echo "== $d"
    w=$out/$d
    stem=$d.decorated.scaled.red

    # 1. The Python inputs of this dataset only, converted to the Rust format.
    mkdir -p "$w/python_input" "$w/input" "$w/python"
    for f in decorated.dict decorated.scaled.red.tree metadata.tsv; do
        ln -s "$py_in/$d.$f" "$w/python_input/$d.$f"
    done
    python3 -I "$repo/scripts/convert_python_inputs.py" "$w/python_input" "$w/input"
    # The names scripts/compare_outputs.py looks for in --rs-input.
    ln -s "$stem.tree" "$w/input/$d.nwk"
    ln -s "$d.taxonomy.tsv" "$w/input/$d.decorated_taxonomy.tsv"
    # The names it (and compare_truth.py) looks for in --py-out.
    ln -s "$py_out/$stem.curatomatic.tree" "$w/python/$stem.curatomatic.tree"
    ln -s "$py_out/$d.log" "$w/python/curatomatic.log"

    # 2. Rust.
    "$repo/target/release/curatomatic" --min-bs 95.0 \
        --taxonomy "$w/input/$d.taxonomy.tsv" \
        --novel-name-file "$w/input/$d.metadata.tsv" --novel-name-column ncbi_wgs_formatted \
        "$w/input/$stem.tree" "$w/input/$d.median_red.tsv" "$w/rust" > /dev/null 2>&1 \
        || { echo "error: curatomatic failed, see $w/rust/curatomatic.log" >&2; exit 1; }

    # 3. and 4. Rust vs Python.
    "$PY" -I "$repo/scripts/compare_outputs.py" --name "$d" \
        --input-tree "$py_in/$stem.tree" --rs-input "$w/input" \
        --py-out "$w/python" --rs-out "$w/rust" --out "$w/python_comparison" > /dev/null
    python3 -I "$repo/scripts/report_placed_differently.py" \
        "$w/rust/$stem.curatomatic.report.tsv" "$w/python_comparison" "$w/python_comparison"

    # 5. Input, Rust and Python vs the truth.
    python3 -I "$repo/scripts/compare_truth.py" --name "$d" --truth "$truth/$d.tsv" \
        --baseline "$w/input/$d.taxonomy.tsv" --rs-out "$w/rust" --py-out "$w/python" \
        --out "$w/truth_comparison"
done

# All the scores in one file.
first=1
for d in "${datasets[@]}"; do
    if [[ $first == 1 ]]; then head -1 "$out/$d/truth_comparison/scores.tsv"; first=0; fi
    tail -n +2 "$out/$d/truth_comparison/scores.tsv"
done > "$out/scores.tsv"
echo "Done: $out (commit $commit)"
