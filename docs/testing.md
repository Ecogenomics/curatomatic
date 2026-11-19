# Testing curatomatic against curated releases

This file is the full brief for a testing session. Everything needed is here; read `AGENTS.md`
first for the project rules (commit after each task, never push, no new crates without approval,
tests in separate files).

## The tasks

1. **Rust vs truth.** Run the current Rust version on each `test_data/<release>_<domain>/rs_input`
   dataset and compare its output taxonomy with the curated taxonomy of that release (the source
   of truth, see [Datasets](#datasets)).
2. **Python vs truth.** Run (or reuse the stored output of) the Python version on the same
   datasets, compare it with the same source of truth, and compare how well each version does.

The point is to measure how close each version's *pre-curation* gets to what the curators finally
published, and where the Rust version is better or worse than the Python one.

## Background

- Curatomatic pre-curates a GTDB tree: it moves existing taxon labels to better nodes (by RED and
  bootstrap) and creates novel taxa for genomes missing a rank. A curator then fixes the rest by
  hand. The curated taxonomy is what the curators published after that.
- **RED**: relative evolutionary divergence, 0 at the root and 1 at the leaves. Each rank has a
  median RED (phylum to genus) that defines where that rank should be on the tree.
- Ranks: `d__ p__ c__ o__ f__ g__ s__`. A placeholder is the prefix alone (`g__`), not a name.
  `s__{unresolved} NAME` species are never changed by curatomatic.
- The rules the Rust version follows are in `docs/refinement.md`. The intended differences from
  Python are in `docs/differences.md` (L1 to L7). Upcoming work is in `TODO.md`.

## Datasets

Four datasets have inputs for both versions and a source of truth. All paths are relative to the
repository root (`/home/aaron/gtdb/curatomatic`).

| Dataset | Source of truth | Truth rows | Tree genomes |
| --- | --- | ---: | ---: |
| `r226_ar53` | `test_data/r226_ar53/ar53_taxonomy_r226_curated.tsv` | 17,245 | 6,968 |
| `r226_bac120` | `test_data/r226_bac120/bac120_taxonomy_r226_curated.tsv` | 715,230 | 136,646 |
| `r232_ar53` | `test_data/r232_ar53/ar53_taxonomy_r232.tsv` | 22,343 | 10,122 |
| `r232_bac120` | `test_data/r232_bac120/bac120_taxonomy_r232.tsv` | 878,998 | 189,801 |

Note the r232 truth files are **not** named `_curated`.

A fifth set, `test_data/r237/`, has no Python input (no old-format tree or `.dict`), so it can
only be used for task 1. Its truth is `test_data/r237/final_taxonomy.tsv`, and its Rust inputs are
`test_data/r237/median_red_marked.nwk`, `median_red_rank_distribution.tsv`,
`gtdb_r237_bac120_metadata.tsv` (or `example/gtdb_r237_bac120_metadata_slim.tsv`), and the
taxonomy `example/gtdb_r237_decorated_taxonomy.tsv`. It is optional.

### Source of truth files

- Tab separated, no header, two columns: genome ID and a `;`-joined taxonomy (no spaces).
- r226 IDs are canonical (`G945873965`). r232 IDs are accessions (`RS_GCF_964307495.1`,
  `GB_GCA_048462655.1`). **Make IDs canonical before matching**: `RS_GCF_000008885.1` and
  `GB_GCA_000008885.1` are both `G000008885` (regex `^(?:RS_|GB_)?GC[AF]_(\d{9})\.\d+$`, then
  `G` + the 9 digits). The r226 `_curated.tsv` files may lack a final newline.
- The truth has **every** genome of the release (representatives and non-representatives). Only
  the genomes on the tree matter: every tree genome is in its truth file (checked, all four).
- For every tree genome the truth has all 7 ranks, **no placeholders and no
  `s__{unresolved}` species** (checked, all four).

### Rust inputs: `test_data/<dataset>/rs_input/`

Converted from the old (Python) format in an earlier session. Each has four files, for dataset
`<d>` (e.g. `r226_ar53`) and the original tree stem `<stem>`:

| File | Contents |
| --- | --- |
| `<stem>.nwk` | The tree. Every dummy pair `('G1|RED=1.0':0,'D-G1|RED=1.0':0)'<taxa>|RED=1.0':L` was collapsed to the leaf `G1:L`. Every other label, bootstrap, RED and branch length is byte for byte the original. Internal labels keep the old style, e.g. `'100.0:f__X; g__Y|RED=0.957017'`, which the Rust parser reads. |
| `gtdb_<d>.decorated_taxonomy.tsv` | Each genome's taxonomy read **literally** off the original tree, joined with `; `: every taxon on the path including the dummy node, so gaps and repeated ranks are kept as on the tree. |
| `gtdb_<d>.median_red_rank_distribution.tsv` | `rank` and `median` columns, Phylum to Genus, from the old `.dict`. |
| `gtdb_<d>_metadata_slim.tsv` | `accession` and `ncbi_wgs_formatted` columns of the old metadata file. |

The stems are `gtdb_r226_ar53.decorated.scaled.red`, `gtdb_r226_bac120.decorated.scaled.red`,
`gtdb_r232_ar53.decorated.red_decorated` and `gtdb_r232_bac120.decorated.red_decorated`.

### Python inputs: `test_data/<dataset>/`

The Python version cannot read `rs_input` (it needs the dummy nodes and the JSON RED dict), so it
uses the original files the `rs_input` files were made from. "The same input data" means these:

| Dataset | Tree | RED dict | Metadata (`--meta`) |
| --- | --- | --- | --- |
| `r226_ar53` | `gtdb_r226_ar53.decorated.scaled.red.tree` | `gtdb_r226_ar53.decorated.dict` | `metadata.tsv` |
| `r226_bac120` | `gtdb_r226_bac120.decorated.scaled.red.tree` | `gtdb_r226_bac120.decorated.dict` | `metadata.tsv` |
| `r232_ar53` | `gtdb_r232_ar53.decorated.red_decorated.tree` | `gtdb_r232_ar53.decorated.dict` | `curatomatic_metadata.tsv` |
| `r232_bac120` | `gtdb_r232_bac120.decorated.red_decorated.tree` | `gtdb_r232_bac120.decorated.dict` | `curatomatic_metadata.tsv` |

- The r226 `metadata.tsv` columns are `accession, taxonomy, ncbi_wgs_formatted`. Python reads a
  prior taxonomy only from a column named `gtdb_taxonomy`, so **the r226 Python runs had no prior
  taxonomy and did no downward refinement**.
- The r232 `curatomatic_metadata.tsv` columns are `accession, ncbi_wgs_formatted, gtdb_taxonomy`
  (last release's taxonomy), so **the r232 Python runs did refine downward**.

### Stored outputs

- **Python:** `test_data/<dataset>/curatomatic_bs95.0/`, run with `--min-bs 95.0` by Python
  v0.0.2. Files: `curatomatic.log`, `<stem>.curatomatic.tree` and
  `<stem>.curatomatic.tree.taxonomy.tsv`.
  - The r226 runs were made from this repository. The `r226_bac120` output was reproduced byte
    for byte from the `main` branch with the driver below. The other three have not been checked
    yet.
  - The r232 runs were made on a server, from `/srv/projects/gtdb/release232/...` paths (see the
    first line of each log). They are assumed to be the same files as the local ones, which has
    not been checked.
- **Rust:** `test_output/<dataset>/`, made by an **older** Rust build (before commit `a3aee1e`).
  The refinement logic hasn't changed since, but task 1 asks for the current version, so **rerun
  Rust** (see [Running Rust](#running-the-rust-version)).

## Running the Rust version

Build, then run each dataset. The output directory must not exist (unless `--log debug`).

```shell
cargo build --release
d=r226_ar53; stem=gtdb_r226_ar53.decorated.scaled.red    # repeat for each dataset
./target/release/curatomatic --min-bs 95.0 \
    --taxonomy test_data/$d/rs_input/gtdb_$d.decorated_taxonomy.tsv \
    --novel-name-file test_data/$d/rs_input/gtdb_${d}_metadata_slim.tsv \
    --novel-name-column ncbi_wgs_formatted \
    test_data/$d/rs_input/$stem.nwk \
    test_data/$d/rs_input/gtdb_$d.median_red_rank_distribution.tsv \
    <out_dir>/$d
```

Each run takes a few seconds. Outputs in `<out_dir>/<d>/`:

- `<stem>.curatomatic.taxonomy.tsv`: one row per genome, `G…` ID, all 7 ranks joined with `;`,
  placeholders for empty ranks (`f__;g__`). Each rank is the taxon nearest the leaf. **Use this
  file** for the comparison with the truth.
- `<stem>.curatomatic.tree`: leaves are bare IDs. Taxa held by leaves (e.g. species) are only in
  the taxonomy file, not on the tree.
- `<stem>.curatomatic.report.tsv`: one row per moved or created taxon (`taxon, rank, action, flag,
  n_genomes_before/after, bootstrap_before/after, red_before/after, prior_taxon, nested_in,
  sister_taxa`). `action` is `moved_up`, `moved_down` or `novel`.
- `curatomatic.log`. Expect many `Shared name` warnings (from repeated ranks in the literal
  taxonomy) and `UNTESTED CODE` warnings. Both are normal.

## Running the Python version

The Python implementation is on the `main` branch (`curatomatic/`). Its CLI needs `typer`, which
isn't installed. Rather than install it, export `main` and run its steps with a small driver; this
reproduced the stored `r226_bac120` output byte for byte.

**Environment:** the conda environment `general`, at `/home/aaron/.conda/envs/general/bin/python`.
It has dendropy 5.1.3 and no typer. The shell's default `python` is the system one and has no
dendropy. The sandbox can't read `~/.conda`, so these commands need the sandbox disabled. The user
allows installing curatomatic's Python dependencies **only** if you are certain you are in the
`general` environment, but the driver makes that unnecessary.

```shell
PY=/home/aaron/.conda/envs/general/bin/python
PYDIR=<scratch dir>/pymain                 # a fresh, empty directory
git archive main curatomatic | tar -x -C $PYDIR
# write the driver below to $PYDIR/run_py.py, then, from inside $PYDIR:
cd $PYDIR && $PY run_py.py <tree> <dict> <out_dir> <metadata>
```

`run_py.py`, the body of `run()` in `main:curatomatic/__main__.py` without typer, with
`--min-bs 95.0`:

```python
"""Run the Python curatomatic (main) without typer: the body of __main__.run."""
import sys
from pathlib import Path
from curatomatic.model.log import LogLevel
from curatomatic.model.metadata import Metadata
from curatomatic.model.rank import RANKS
from curatomatic.model.red import RED
from curatomatic.model.tree import Tree
from curatomatic.util.log import log_setup
from curatomatic.util.tree import fill_in_missing_ranks

tree_p, red_p, out_dir, meta_p = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
min_bs = 95.0
out_dir.mkdir(parents=True, exist_ok=True)
logger = log_setup(out_dir=out_dir, level=LogLevel.INFO)
tree = Tree.from_file(tree_p)
metadata = Metadata.from_file(meta_p)
metadata.subset_to_ids(set(tree.gid_to_node.keys()))
red_dict = RED.from_file(red_p)
for rank in RANKS:
    if rank in {'d', 's'}:
        continue
    logger.info(f'Refining existing rank: {rank}')
    ref = tree.refine_existing_labels(rank, min_bs, red_dict, metadata)
    tree.process_refinement(ref, rank, min_bs, red_dict, metadata)
    tree.reload()
missing = tree.get_leaf_nodes_missing_rank()
fill_in_missing_ranks(missing, tree.node_to_descs, min_bs, red_dict, metadata, set(tree.taxon_labels_to_node.keys()))
tree.remove_red_labels()
tree.to_file(out_dir / f'{tree_p.stem}.curatomatic.tree')
```

- Runtime: about 70 s for `r226_bac120`. `r232_bac120` is larger.
- The driver writes the tree and `curatomatic.log` only, like the original. The stored
  `.taxonomy.tsv` files were made separately, and are the **raw path** (every taxon on the path
  joined with `; `, gaps and repeats kept), so they are not comparable with the truth as is.
- To get one taxon per rank from a Python tree, read it with dendropy and take, for each rank, the
  taxon nearest the leaf. The dummy pair's parent holds the singleton taxa (species, and any
  rank the tree couldn't hold). `scripts/compare_outputs.py` already does this:
  `Tree(path).genome_taxonomy()` returns `{gid: {rank_letter: taxon}}` with placeholders dropped,
  and handles both Python trees (dummy pairs) and Rust trees (`Tree.add_leaf_taxa(<taxonomy.tsv>)`
  puts the leaf taxa back).
- Python's novel taxa without a WGS name are `<rank>__NOVEL-<gid>`. Rust's are `<rank>__<gid>`
  (same genome ID in every case seen so far).

## Comparing with the truth: things to get right

- **Restrict to tree genomes.** The truth has many more genomes than the tree.
- **Canonical IDs** for r232 truth (above).
- **Species.** Curatomatic never names species: they come from the input, and `s__{unresolved}`
  species stay unresolved, while the truth names every species. Species agreement measures the
  input, not curatomatic, so report it separately or leave it out.
- **Placeholders.** The truth has none, so a placeholder in an output (Rust writes `f__`, `g__`
  etc.; Python leaves the rank out) is always a miss. Count "left empty" separately from "wrong
  name".
- **Novel names.** Curators may rename a novel taxon (or keep its WGS name). Compare in two ways:
  - **Exact name** per genome and rank.
  - **Grouping** per rank, ignoring names: do the genomes that share a taxon in the output share
    one in the truth? For example, for each output taxon the best matching truth taxon by overlap,
    or pairwise precision and recall of "same taxon", or the adjusted Rand index. This shows
    whether a taxon is placed right even when its name changed.
- **Baseline.** Also score the input (`rs_input/gtdb_<d>.decorated_taxonomy.tsv`, taking the
  taxon nearest the leaf for each rank) against the truth. This is "no curatomatic at all", and
  shows how much each version helps or hurts.
- **Bias (unconfirmed).** The r232 truth was probably curated starting from the stored Python
  output, and r226 from an earlier Python run. If so, the truth favours Python wherever a curator
  accepted its suggestion without checking. Ask the user before drawing conclusions from small
  differences.
- **Per rank, per dataset.** Report d to g separately for all four datasets, plus the counts behind
  each percentage. For the bigger differences, list example taxa and genomes, and explain them by
  the causes below.

## Known differences between the two versions

From an earlier session (`scripts/compare_outputs.py`, see its docstring). Every difference
between the stored Python and Rust outputs was explained by one of these:

| Cause (`docs/differences.md`) | Effect on these datasets |
| --- | --- |
| **L1** Nested taxa never move up in Rust | 1 / 13 / 0 / 22 taxa (r226_ar53 / r226_bac120 / r232_ar53 / r232_bac120). E.g. Python moved `c__Bacilli_A` up and took 3,528 genomes from `c__Bacilli`. |
| **L2** Downward refinement | Rust made **no** downward moves anywhere: its prior taxonomy is `--taxonomy`, which here is read off the same tree, so rule D2 never fires. Python moved 146 (r232_ar53) and 1,451 (r232_bac120) taxa down, using `gtdb_taxonomy`. None in r226. |
| **L4** Novel names without WGS | 3 / 181 / 3 / 205: `NOVEL-G…` vs `G…`, same genome ID. |
| **L6** Leaf taxa from `--taxonomy` | 12 / 167 / 13 / 287 singleton taxa are not on the Rust tree: a leaf only gets taxa ranked below the lowest rank on its path, so nested singletons (e.g. `g__Phormidium` inside `g__PCC-6304`) are dropped. Side effects: the genome gets the outer taxon (G017313335, *s__Phormidium pseudopriestleyi*, is reported in `g__PCC-6304`), and the dropped taxon no longer blocks upward moves (`g__PCC-6304` and `f__Sporomusaceae` move up in Rust only). **The user considers this unacceptable: Python is more correct here. Not fixed yet.** |
| **L7** No novel taxa above named taxa | 15 / 578 / 187 / 2,811 novel taxa made only by Python. Mostly a genus above a named species: e.g. `s__JAHITG01 sp965373035` with no `g__JAHITG01` on its path, where Python made `g__CBDVEP01`. Also, once Rust skips a rank for a genome it fills none of its lower missing ranks: G041740975 (r226_bac120) gets `p__;…;f__;g__` in Rust, but `p__CACXFL01` and `f__JBGSKV01; g__JBGSKV01` (its WGS prefix) in Python. |

Expect these to show up in the comparison with the truth. They are the first places to look when
one version beats the other.

Per-genome counts that differ between the two versions (r226_ar53 / r226_bac120 / r232_ar53 /
r232_bac120, after mapping `NOVEL-` names): p 0 / 628 / 53 / 5,110, c 0 / 3,620 / 26 / 3,981,
o 1 / 429 / 16 / 1,924, f 4 / 329 / 72 / 899, g 35 / 753 / 254 / 3,776, s 0 everywhere.

### Open issues (state at the time of writing)

- `TODO.md` Job 10: an upward move can put a taxon above a higher rank. Rust puts
  `f__Sporomusaceae` above `o__Sporomusales` in `r226_bac120`. Both versions have this gap.
- L6 (nested singletons dropped): to be fixed. Python's behaviour is the reference.
- L2: whether Rust should take a separate prior taxonomy so downward refinement can run is still
  being discussed.
- The end-of-run check (`check_taxa_unchanged`, `src/main.rs`) covers the taxa on the tree after
  the leaf assignment. It doesn't cover taxa in `--taxonomy` that were never assigned.

## Tools already in the repository

- `scripts/compare_outputs.py`: Python vs Rust, node by node (see its docstring for arguments;
  run with the conda Python above, sandbox disabled). It writes `summary.md`,
  `existing_taxa.tsv`, `novel_taxa.tsv`, `genome_diff.tsv` and `log_moves.tsv`. Its `Tree` class
  (genome-set keyed nodes, `genome_taxonomy()`, `add_leaf_taxa()`), `canonical()` and
  `read_taxonomy()` can be reused for the comparison with the truth. A new script for the truth
  comparison belongs in `scripts/` too.
- `scripts/extract_taxonomy.py` (untracked, the user's): reads the raw path taxonomy off a tree.
  It doesn't strip `|RED=…`, so don't use it on RED-decorated trees.

## Practical notes

- Write scratch output (runs, reports) to the session's scratchpad or `$TMPDIR`, not the
  repository. Give the user full absolute paths to anything they should inspect.
- `test_data/` and `test_output/` are the user's. Don't modify or delete them.
- Files for the user to review go in a temporary directory, with the path given to the user.
- The repository has untracked dotfiles in its root (`.bashrc`, `.zshrc`, `.claude/*`, ...) and a
  modified `.gitignore` that aren't part of the project. Leave them out of commits.
- After any code change: `cargo fmt`, `cargo clippy --all-targets` (lints are `deny`) and
  `cargo nextest run`, then commit with a semantic message.
