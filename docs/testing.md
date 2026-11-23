# Testing curatomatic against curated releases

This file is the full brief for a testing session. Everything needed is here; read `AGENTS.md`
first for the project rules (commit after each task, never push, no new crates without approval,
tests in separate files).

## The tasks

1. **Rust vs truth.** Run the current Rust version on each dataset and compare its output taxonomy
   with the curated taxonomy of that release (the source of truth, see [Datasets](#datasets)).
2. **Python vs truth.** Compare the stored output of the Python version with the same source of
   truth, and compare how well each version does.
3. **Rust vs Python.** Compare the two outputs node by node, and list the taxa they place
   differently.

`scripts/benchmark.sh` does all three (see [Running the benchmark](#running-the-benchmark)).

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
  Python are in `docs/differences.md` (L1 to L9). Upcoming work is in `TODO.md`.

## Datasets

Four datasets (`<d>` below) have inputs for both versions, a stored Python output and a source of
truth. All paths are relative to the repository root. `test_data/` is not in git (it is the
user's, and is never written to): its layout is

```
test_data/
├── curated/<d>.tsv                                     the source of truth
├── python/input/<d>.decorated.scaled.red.tree          Python input: the tree, with dummy nodes
├── python/input/<d>.decorated.dict                     Python input: median RED of each rank (JSON)
├── python/input/<d>.metadata.tsv                       Python input: metadata (--meta)
├── python/output/<d>.decorated.scaled.red.curatomatic.tree   stored Python output
├── python/output/<d>.log                               its log
└── rust/input/gtdb_r237_*                              the r237 Rust inputs (no Python, no truth)
```

| Dataset `<d>` | Truth rows | Tree genomes |
| --- | ---: | ---: |
| `r226_ar53` | 17,245 | 6,968 |
| `r226_bac120` | 715,230 | 136,646 |
| `r232_ar53` | 22,343 | 10,122 |
| `r232_bac120` | 878,998 | 189,801 |

There are no Rust inputs for these four in `test_data/`: they are made from the Python inputs by
`scripts/convert_python_inputs.py` on every run (it takes a few seconds), see
[Rust inputs](#rust-inputs).

### Source of truth files

- Tab separated, no header, two columns: genome ID and a `;`-joined taxonomy (no spaces).
- r226 IDs are canonical (`G945873965`). r232 IDs are accessions (`RS_GCF_964307495.1`,
  `GB_GCA_048462655.1`). **Make IDs canonical before matching**: `RS_GCF_000008885.1` and
  `GB_GCA_000008885.1` are both `G000008885` (regex `^(?:RS_|GB_)?GC[AF]_(\d{9})\.\d+$`, then
  `G` + the 9 digits). The r226 files may lack a final newline.
- The truth has **every** genome of the release (representatives and non-representatives). Only
  the genomes on the tree matter: every tree genome is in its truth file (checked, all four).
- For every tree genome the truth has all 7 ranks, **no placeholders and no
  `s__{unresolved}` species** (checked, all four).

### Python inputs

- The r226 metadata columns are `accession, taxonomy, ncbi_wgs_formatted`. Python reads a prior
  taxonomy only from a column named `gtdb_taxonomy`, so **the r226 Python runs had no prior
  taxonomy and did no downward refinement**.
- The r232 metadata columns are `accession, ncbi_wgs_formatted, gtdb_taxonomy` (last release's
  taxonomy), so **the r232 Python runs did refine downward**.

### Rust inputs

`python scripts/convert_python_inputs.py <dir> <out_dir>` converts every `<d>` in `<dir>` (to
convert one dataset, give it a directory with links to that dataset's three files only). It
writes:

| File | Contents |
| --- | --- |
| `<d>.decorated.scaled.red.tree` | The tree. Every dummy pair `('G1|RED=1.0':0,'D-G1|RED=1.0':0)'<taxa>|RED=1.0':L` is collapsed to the leaf `G1:L`. Internal labels are unchanged apart from `; ` becoming `;` and `|RED=` becoming `|red=`. |
| `<d>.taxonomy.tsv` | `--taxonomy`: each genome's taxa off its path in the original tree, dummy node included, grouped by rank (placeholders for missing ranks, repeated ranks kept). |
| `<d>.median_red.tsv` | `rank` and `median` columns, Phylum to Genus, from the `.dict`. |
| `<d>.metadata.tsv` | `--novel-name-file`: the `accession` and `ncbi_wgs_formatted` columns. |

These inputs round trip with 0 differences (see the README and `scripts/check_round_trip.py`).

### Prior taxonomy (`--prior-taxonomy`)

Moving taxa down needs the taxonomy of a previous release (`docs/refinement.md`, rule D2). Make it
from the Python metadata (both domains of a release share one metadata file; genomes that aren't on
the tree are ignored):

- **r226**: column 2 `taxonomy`, 596,626 named rows. It is probably r220 (not checked). The `none`
  rows must be dropped:
  `awk -F'\t' 'NR>1 && $2!="none"{print $1"\t"$2}' test_data/python/input/$d.metadata.tsv`
- **r232**: column 3 `gtdb_taxonomy`, 731,926 named rows. For every genome named in both, it is
  identical to the curated r226 taxonomy. Its 320,439 placeholder-only rows
  (`d__Bacteria;p__;…;s__`, genomes new in r232) must be dropped, or those genomes count as having
  a prior and rule D2a never applies to them:
  `awk -F'\t' 'NR>1{n=split($3,a,";"); for(i=2;i<=n;i++) if(length(a[i])>3){print $1"\t"$3; next}}' test_data/python/input/$d.metadata.tsv`
  It keeps a row with any named rank below the domain: a partly placeholder row would now stop the run, as every
  rank of a prior row must be named (none of r232's rows is one).

All four of these, `test_data/curated/*.tsv` and `example/gtdb_r232_bac120_taxonomy.tsv` load as priors (checked on
2026-11-23). `test_data/r237/final_taxonomy.tsv` does not (5,070 rows with `g__;s__{unresolved} …`, and 605 whose
species is in another genus), but it is the curated truth of r237, not a prior. Requiring the full
prior row (`TODO.md` Job 23) left the output with these four priors byte-identical, and the curator
workload below at 8,547.

These give the prior files that the [curator workload](#curator-workload-at-be526b9) was measured
with.

`example/gtdb_r232_bac120_taxonomy.tsv` (the prior for the r237 example) is cut down so that it can
be pushed to GitHub (the full file is 126 MB, over its 100 MB limit). It keeps every row whose
genome is on `example/gtdb_r237_median_red_marked.nwk`, plus one row for each taxon that only
off-tree genomes have: the report's set of prior taxa is read before the prior is subset to the
tree (`read_metadata`, `src/main.rs`). The output is byte-identical to that of the full file.

### Stored Python outputs

- Run with `--min-bs 95.0` by Python v0.0.2. Each is a tree and a log only.
- The r226 runs were made from this repository. The `r226_bac120` output was reproduced byte for
  byte from the `main` branch with the driver below. The other three have not been checked.
- The r232 runs were made on a server, from `/srv/projects/gtdb/release232/...` paths (see the
  first line of each log). They are assumed to be the same files as the local ones, which has not
  been checked.

## Running the benchmark

```shell
PY=/home/aaron/.conda/envs/general/bin/python scripts/benchmark.sh <out_dir> [<d> ...]
```

With no datasets, all four are run. `<out_dir>` must not exist. It takes about two minutes for
all four (most of it in `compare_outputs.py`). `PY` must have dendropy (only
`compare_outputs.py` needs it, see [Running the Python version](#running-the-python-version) for
the environment); the sandbox can't read `~/.conda`, so run it with the sandbox disabled. The
script builds the release binary first, and records the commit it ran (`-dirty` if the working
tree has changes) in `<out_dir>/commit.txt`.

For each dataset, `<out_dir>/<d>/` holds:

| Path | Contents |
| --- | --- |
| `python_input/` | Links to the dataset's Python inputs. |
| `input/` | The converted Rust inputs (above), plus the links `<d>.nwk` and `<d>.decorated_taxonomy.tsv` that `compare_outputs.py --rs-input` looks for. |
| `python/` | Links to the stored Python tree and log, as `compare_outputs.py` and `compare_truth.py --py-out` look for them. |
| `rust/` | The Rust output (see [below](#rust-output)). |
| `python_comparison/` | `scripts/compare_outputs.py` (`summary.md`, `existing_taxa.tsv`, `novel_taxa.tsv`, `genome_diff.tsv`, `log_moves.tsv`), and `scripts/report_placed_differently.py`: `report.placed_differently.tsv` (the rows of the Rust report that Python places differently) and `report.placed_differently.all.tsv` (every taxon placed differently, including those Rust left alone). |
| `truth_comparison/` | `scripts/compare_truth.py`: `scores.tsv`, `genomes.tsv`, `taxa.tsv`, `changes.tsv`, `py_vs_rs.tsv` (see its docstring). The `baseline` version is the input taxonomy, i.e. no curatomatic at all. |

`<out_dir>/scores.tsv` has the scores of every dataset in one file.

The first sanity check in `python_comparison/summary.md` must say that the topology of all three
trees and the taxa of the converted input match the input.

The script doesn't pass `--prior-taxonomy` yet, so it measures curatomatic without downward
refinement, and it doesn't score the curator workload (`TODO.md` Job 18).

### Results at `b44e6b6`

Per cent of tree genomes whose taxon is the truth's (`pct_correct` in `scores.tsv`):

| Dataset | Version | p | c | o | f | g |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| r226_ar53 | input | 83.0 | 84.2 | 99.1 | 96.0 | 90.1 |
| | Rust | 83.0 | 84.2 | 99.2 | 96.5 | 93.8 |
| | Python | 83.0 | 84.2 | 99.2 | 96.5 | 93.8 |
| r226_bac120 | input | 78.3 | 92.5 | 97.4 | 97.1 | 91.4 |
| | Rust | 78.3 | 92.5 | 97.5 | 97.4 | 94.6 |
| | Python | 78.3 | 90.0 | 97.6 | 97.1 | 94.3 |
| r232_ar53 | input | 99.5 | 99.3 | 94.6 | 95.8 | 85.0 |
| | Rust | 99.5 | 99.4 | 94.9 | 96.8 | 91.3 |
| | Python | 100.0 | 99.6 | 94.9 | 97.2 | 92.8 |
| r232_bac120 | input | 94.6 | 96.5 | 95.0 | 93.2 | 89.6 |
| | Rust | 94.6 | 96.5 | 95.2 | 93.7 | 93.0 |
| | Python | 94.6 | 96.6 | 94.3 | 93.7 | 94.1 |

- Rust never scores below the input. Most of its gain is at genus and family, from filling ranks
  the input left empty. Its remaining misses are mostly taxa that were already wrong in the input.
- Python's ambiguous genome-ranks (a node holding two taxa of the rank) count as neither correct
  nor wrong, e.g. 1,993 orders in r232_bac120.
- Rust vs Python (`report.placed_differently*.tsv`): Rust never moved a taxon that Python left in
  place, and never made a novel taxon that Python didn't. Only 3 rows of the Rust reports are
  placed differently, all in r232_bac120 (`f__JBAZCA01`, `g__JAFGQJ01`, `g__JALSRK01`): Python
  moved them down (L2), Rust moved them up, and the truth agrees with Python on all 9 genome-ranks
  they change. The other differences are taxa only Python moved (1 / 13 / 146 / 1,470) and novel
  taxa only Python made (15 / 578 / 245 / 3,120), for r226_ar53 / r226_bac120 / r232_ar53 /
  r232_bac120.

### Curator workload at `be526b9`

The user's metric of success is the curator's manual work, not the genome-weighted scores above.
It counts, per rank p to g and with each taxon counted once however many genomes it has, the
operations a curator needs to turn an output into the truth: `move` (a truth taxon whose name is in
the output, on a different genome set), `create` (a truth taxon with no counterpart, including
ranks left empty) and `delete` (an output taxon with no counterpart in the truth). Renames are left
out: curatomatic never edits a name. The scorer is `wip/scripts/workload.py`, which is not
committed yet (`TODO.md` Job 18).

Operations (move / create / delete), with the prior files [above](#prior-taxonomy---prior-taxonomy):

| Dataset | Input | Python | `32d3350` | Jobs 10 to 13 (`2622c2f`) | `be526b9` (Jobs 10 to 13 and 15) |
| --- | ---: | ---: | ---: | ---: | ---: |
| r226_ar53 | 469 (108 / 337 / 24) | 248 (105 / 104 / 39) | 263 (105 / 119 / 39) | 173 | 179 |
| r226_bac120 | 7,621 (1,726 / 5,449 / 446) | 3,860 (1,570 / 1,431 / 851) | 4,103 (1,573 / 1,837 / 693) | 2,656 | 2,758 |
| r232_ar53 | 1,237 (296 / 822 / 119) | 470 (168 / 66 / 236) | 693 (258 / 255 / 180) | 474 | 482 |
| r232_bac120 | 13,580 (3,060 / 9,400 / 1,120) | 5,030 (1,844 / 1,062 / 2,107) | 7,918 (2,835 / 3,570 / 1,513) | 5,088 | 5,128 |
| Total | 22,907 | 9,608 | 12,977 | 8,391 | 8,547 |

`32d3350` is the code before Jobs 10 to 15. Each job's effect on the total, against `32d3350`:
Job 10 +3, Job 11 -659, Job 12 -236, Jobs 11 and 12 -903, Jobs 10 to 12 -900, Job 13 alone -2,073,
Jobs 10 to 13 -4,586. Job 15 then adds +156. Applied to every genus, it added +3,003 (11,394 in
total): it stopped 1,992 moves down, and for most of them the truth was closer to the move, as new
genomes' species echo the decoration. The user chose the narrowed rule (Latin names inside the RED
window, `docs/differences.md` L9), accepting +156, to keep genera such as `g__Noviluteimonas`
(r237) whole. Without `--prior-taxonomy`, Jobs 10 to 13 give the same output as Jobs
10 to 12.

Job 14 (an unresolvable novel taxon is reported, never written) changes only the genomes that got
`<rank>__UNRESOLVABLE`: G020350785 (`g__`) and G030406045 (`o__`, `f__`, `g__`) in r226_bac120,
and G020350445 (`g__`) in r232_bac120 (none in the ar53 datasets). Each is now a placeholder and
one `unresolvable` report row, and every other output is unchanged. With the `wip/bench_v3`
priors, the total falls from 8,547 to 8,543 (r226_bac120 2,758 to 2,755, r232_bac120 5,128 to
5,127): the sentinels were taxa to delete.

Caveat: the r232 prior is the curated r226, and the r232 truth was probably curated from the Python
output, which made the same moves down (146 / 1,451 taxa, equal to Rust's). So r232 favours Job 13.
r226 (prior probably r220; Python made no moves down there) is the less biased evidence: Jobs 10 to
13 save 90 operations in r226_ar53 and 1,447 in r226_bac120.

## Running the Rust version

To run one dataset by hand, from inputs converted as above (`<in>`):

```shell
cargo build --release
d=r226_ar53
./target/release/curatomatic --min-bs 95.0 \
    --taxonomy <in>/$d.taxonomy.tsv \
    --novel-name-file <in>/$d.metadata.tsv --novel-name-column ncbi_wgs_formatted \
    <in>/$d.decorated.scaled.red.tree <in>/$d.median_red.tsv <out_dir>
```

Each run takes a few seconds. The output directory must not exist (unless `--log debug`).

### Rust output

- `<stem>.curatomatic.taxonomy.tsv`: one row per genome, `G…` ID, all 7 ranks joined with `;`,
  placeholders for empty ranks (`f__;g__`). Each rank is the taxon nearest the leaf. **Use this
  file** for the comparison with the truth.
- `<stem>.curatomatic.tree`: leaves are bare IDs. Taxa held by leaves (e.g. species) are only in
  the taxonomy and leaf taxonomy files, not on the tree.
- `<stem>.curatomatic.leaf-taxonomy.tsv`: one row per genome, the taxa held by its leaf only (no
  placeholders, no taxa of its path). With the tree, it gives back the full taxonomy of every genome
  (see the README). On r237 (the README command, with Job 24) all 229,307 genomes have a row,
  G053840885's is `g__Neosphaerochaeta;s__Neosphaerochaeta fermentans`, the tree plus this file give
  back `<stem>.curatomatic.taxonomy.tsv` for every genome, and the taxa they hold are the input taxa
  plus the 5,447 novel taxa of the report. The tree, report and taxonomy files were unchanged by Job 24.
  Since Job 14 the report has 5,446 novel taxa and one `unresolvable` row (G965247515, namesake
  `g__G965247515`), which is no longer written as `g__UNRESOLVABLE`.
- `<stem>.curatomatic.report.tsv`: one row per taxon that was moved or created, one per novel taxon
  that could not be named, and nothing else (`taxon, rank, action, flag, n_genomes_before/after,
  bootstrap_before/after, red_before/after, prior_taxon, nested_in, sister_taxa, namesakes`).
  `action` is `moved_up`, `moved_down`, `novel` or `unresolvable`.
- `curatomatic.log`.

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
- **Baseline.** Also score the input (the converted `<d>.taxonomy.tsv`, taking the
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

From an earlier session (`scripts/compare_outputs.py`, see its docstring), before Jobs 10 to 15
(L8 and L9 did not exist yet, and the L2 and L7 rows have changed since, see below
the table). Every difference between the stored Python and Rust outputs was explained by one of
these:

| Cause (`docs/differences.md`) | Effect on these datasets |
| --- | --- |
| **L1** Nested taxa never move up in Rust | 1 / 13 / 0 / 22 taxa (r226_ar53 / r226_bac120 / r232_ar53 / r232_bac120). E.g. Python moved `c__Bacilli_A` up and took 3,528 genomes from `c__Bacilli`. |
| **L2** Downward refinement | Rust made **no** downward moves anywhere: its prior taxonomy is `--taxonomy`, which here is read off the same tree, so rule D2 never fires. Python moved 146 (r232_ar53) and 1,451 (r232_bac120) taxa down, using `gtdb_taxonomy`. None in r226. |
| **L4** Novel names without WGS | 3 / 181 / 3 / 205: `NOVEL-G…` vs `G…`, same genome ID. |
| **L6** Leaf taxa from `--taxonomy` | 12 / 167 / 13 / 287 singleton taxa are not on the Rust tree: a leaf only gets taxa ranked below the lowest rank on its path, so nested singletons (e.g. `g__Phormidium` inside `g__PCC-6304`) are dropped. Side effects: the genome gets the outer taxon (G017313335, *s__Phormidium pseudopriestleyi*, is reported in `g__PCC-6304`), and the dropped taxon no longer blocks upward moves (`g__PCC-6304` and `f__Sporomusaceae` move up in Rust only). **Fixed in `cd09ece`** (leaf taxa are matched by taxon, not rank). Rescored after the fix: all 12 / 169 / 13 / 288 genome-ranks it changed now agree with Python. |
| **L7** No novel taxa above named taxa | 15 / 578 / 187 / 2,811 novel taxa made only by Python. Mostly a genus above a named species: e.g. `s__JAHITG01 sp965373035` with no `g__JAHITG01` on its path, where Python made `g__CBDVEP01`. Also, once Rust skips a rank for a genome it fills none of its lower missing ranks: G041740975 (r226_bac120) gets `p__;…;f__;g__` in Rust, but `p__CACXFL01` and `f__JBGSKV01; g__JBGSKV01` (its WGS prefix) in Python. |

Expect these to show up in the comparison with the truth. They are the first places to look when
one version beats the other.

Per-genome counts that differ between the two versions (r226_ar53 / r226_bac120 / r232_ar53 /
r232_bac120, after mapping `NOVEL-` names): p 0 / 628 / 53 / 5,110, c 0 / 3,620 / 26 / 3,981,
o 1 / 429 / 16 / 1,924, f 4 / 329 / 72 / 899, g 35 / 753 / 254 / 3,776, s 0 everywhere.

Changes since the table was made:

- **L2**: with `--prior-taxonomy` (Job 13, `2622c2f`), Rust moves taxa down: 37 / 816 / 146 /
  1,451 taxa with Jobs 10 to 13. The r232 counts equal Python's. Without it, Rust still makes no
  moves down.
- **L7**: Jobs 11 and 12 (`fd77429`, `7513c7d`) removed both causes in the row: a novel genus is
  now created above a named species, and skipping a rank no longer stops the lower ranks being
  filled. Rust still makes no novel taxon above a named taxon of a lower rank, other than a species.
- **L8** (Job 10, `9a6125a`): Rust no longer moves `f__Sporomusaceae` above `o__Sporomusales` in
  r226_bac120. Python has this gap.

### Open issues

- Novel taxa can break rank order (`TODO.md` Job 16, needs approval).
- Rewritten labels round the RED and bootstrap (`TODO.md` Job 17).
- `scripts/benchmark.sh` doesn't use a prior or score the curator workload (`TODO.md` Job 18).

Settled: L6 was fixed in `cd09ece`. Since then every named taxon of a `--taxonomy` row is on the
genome's path or assigned to its leaf, and anything else fails the run, so the end-of-run check
(`check_taxa_unchanged`) covers every input taxon.

## Tools already in the repository

- `scripts/benchmark.sh`: runs everything below for each dataset (see
  [Running the benchmark](#running-the-benchmark)).
- `scripts/convert_python_inputs.py`: the Python inputs to the Rust inputs.
- `scripts/compare_outputs.py`: Python vs Rust, node by node (see its docstring for arguments;
  run with the conda Python above, sandbox disabled). It writes `summary.md`,
  `existing_taxa.tsv`, `novel_taxa.tsv`, `genome_diff.tsv` and `log_moves.tsv`. Its `Tree` class
  (genome-set keyed nodes, `genome_taxonomy()`, `add_leaf_taxa()`), `canonical()` and
  `read_taxonomy()` can be reused.
- `scripts/report_placed_differently.py`: the Rust report, limited to the taxa that
  `compare_outputs.py` found placed differently.
- `scripts/compare_truth.py`: the input, Rust and Python vs the truth.
- `scripts/check_round_trip.py`: checks that the converted inputs load without losing a taxon.

## Practical notes

- Write scratch output (runs, reports) to the session's scratchpad or `$TMPDIR`, not the
  repository. Give the user full absolute paths to anything they should inspect.
- `test_data/` and `test_output/` are the user's. Don't modify or delete them.
- Files for the user to review go in a temporary directory, with the path given to the user.
- The dotfiles that the Claude Code sandbox creates in the root (`.bashrc`, `.zshrc`, ...) are
  ignored by `.gitignore`. So is `wip/`, for work-in-progress notes and runs: never commit it.
- After any code change: `cargo fmt`, `cargo clippy --all-targets` (lints are `deny`) and
  `cargo nextest run`, then commit with a semantic message.
