# TODO

Finished jobs are removed from this file: their rules are in `docs/refinement.md` and `docs/differences.md`, and
their measurements in `docs/testing.md`. Jobs that need the user's approval say so: don't start them before it is
given.

## Job 16: Don't create a novel taxon that breaks rank order

**Status: needs the user's approval** (proposed as P6 on 2026-11-20).

### Problem

`fill_in_missing_ranks` (`src/tree.rs`) can create a novel taxon that is out of rank order with the named taxa around
it. There are two forms:

- A novel taxon **below** a named taxon of a lower rank on its path. In r237 bac120, 10 novel orders are created
  inside `g__MED-G20`, and `o__G000829155` inside `g__Karelsulcia` (RED 0.891).
- A novel taxon **above** a named taxon of a higher rank. In r237 bac120, `g__CAZDDW01` is created above
  `f__DNKR01` and the novel `f__G965231015`.

With Jobs 10 to 13 there are 7 such taxa in r237 bac120 (Job 15 left that unchanged), and 13 / 15 in r226_bac120 /
r232_bac120 (6 / 10 at `32d3350`, before Jobs 11 and 12; not measured after Job 15). Every one checked is placed
elsewhere by the truth (`o__Flavobacteriales_D`, `o__Bacillales`, `f__CAZDDW01;g__G965231015`).

### Rule (proposed)

A novel taxon of rank `R` is skipped if its candidate node is below a named taxon of a rank lower than `R` on its
path, or has a named taxon of a rank higher than `R` below it. Novel taxa of higher ranks made earlier in the same
run count. Skipping only skips that rank (as for the existing check, `docs/refinement.md` "Novel taxa").

### Tasks

1. Add the check to `fill_in_missing_ranks`, next to the `lowest_valid_rank_below` check.
2. Tests in `src/tree/tests.rs`, one for each form, and one where a novel taxon of a higher rank made earlier in the
   run blocks a lower one.
3. Docs: `docs/refinement.md` "Novel taxa", `docs/differences.md` (new entry: Python has no such check).
4. Measure the curator workload on the four datasets (see `docs/testing.md`), and count the rank inversions left in
   r237 bac120.

## Job 17: Keep RED and bootstrap values as they were read when a label is rewritten

**Status: needs the user's decision.**

### Problem

`NodeLabel::to_label` (`src/node_label.rs`) writes the RED with 3 decimal places and the bootstrap with 1
(`100.0`), so a rewritten label differs from the input in format and precision (1,574 nodes in the r237 bac120
output with Jobs 10 to 13). The later refinement passes read the rounded RED back. This copies the Python behaviour,
but `AGENTS.md` says curatomatic never changes RED or bootstrap values.

### Question for the user

Should a rewritten label keep the RED and bootstrap exactly as read (no rounding), and should the later passes use
the unrounded RED? Either answer changes some outputs, so measure it on the four datasets.

## Job 18: Benchmark with a prior taxonomy and score the curator workload

**Status: not started. Confirm with the user before starting.**

### Problem

- `scripts/benchmark.sh` runs curatomatic without `--prior-taxonomy`, so it measures the behaviour without downward
  refinement (Jobs 13 and 15 are not measured).
- It compares against the stored Python outputs, whose r226 runs had no prior, so Python is not given the same input
  as Rust.
- The curator workload, the user's metric of success, is scored by `wip/scripts/workload.py` and
  `wip/scripts/summary.py`, which are not committed (`wip/` is git-ignored). The Python metadata with the prior is
  made by `wip/scripts/make_py_meta.py`, also not committed.

The benchmark was run by hand on 2026-11-23 (`docs/testing.md`, "Curator workload at `3dc483a`", steps in "Running
it"). This job scripts it.

### Tasks

1. Take the prior files as arguments (the official GTDB releases joined across ar53 and bac120: r220 for the r226
   datasets, r226 for the r232 datasets, see `docs/testing.md`, "Prior taxonomy") and pass them with
   `--prior-taxonomy`.
2. Move `make_py_meta.py` into `scripts/`, and re-run Python (`main`) with the same prior instead of linking the
   stored outputs (r232_bac120 takes about 40 min).
3. Move the workload scorer into `scripts/`, run it in `benchmark.sh`, and document its metrics in
   `docs/testing.md`.
4. Check that `benchmark.sh` gives the numbers in `docs/testing.md`, "Curator workload at `3dc483a`".

## Job 26: Rule D2a for non-Latin genera on a supported node

**Status: requested by the user on 2026-11-23, after the measurement below (+284 operations) was shown to them.**

### Problem

In r237 bac120 (the `example/` data, with `--prior-taxonomy example/gtdb_r232_bac120_taxonomy.tsv`), the new genome
G057028815 is split from its genus and given a novel genus `g__JBXLYT01`. Its species is `s__JAHELV01 sp057028815`,
which comes from an earlier species clustering and puts it in `g__JAHELV01`:

```
100:f__JAHELV01;g__JAHELV01|red=0.8758   (9 genomes)
├── G057028815                       new (not in r232)
└── 100:red=0.9095   (8 genomes)     all 5 genomes in r232 are g__JAHELV01, plus 3 new ones
```

Rule D2 (`docs/refinement.md`) follows the prior genomes into the 8-genome child, which is supported and closer to
the genus median (0.9189), so the genus is moved there (`Moving g__JAHELV01 from 100:f__JAHELV01;g__JAHELV01|red=0.8758
to 100:red=0.9095`). G057028815 is left with no genus, and `fill_in_missing_ranks` creates `g__JBXLYT01` on its leaf.
The curators would keep `g__JAHELV01` on the family node, with all 9 genomes.

Rule D2a (a new genome counts as being in the genus its species names) would stop this, but it only applies to a
Latin genus name (`Taxon::is_latin_name`), and `g__JAHELV01` is not one.

### Rule

When refining the genus, a new genome (one without a prior taxonomy) whose species is a binomial name counts as being
in the genus that the first word gives (`s__JAHELV01 sp057028815` gives `g__JAHELV01`; `s__{unresolved} NAME` gives
none). This applies to a genus only when its original node:

1. is inside the genus RED window, for every genus; and
2. for a **non-Latin** genus name (`g__JAHELV01`, `g__CAG-698`), is also **supported** (bootstrap `≥ --min-bs`).

A Latin genus name (`g__Noviluteimonas`, `g__Luteimonas_D`) needs the window only, as now. A genome in the prior is
always given its prior genus, never the genus of its species. When the rule applies, the descent works as before
(Rule D2): if more than one child contains the genus, it stays on its original node.

| Genus name | Original node | D2a now | D2a after this job |
| --- | --- | --- | --- |
| Latin | in window | applies | applies |
| Latin | outside window | no | no |
| non-Latin | in window, supported | no | **applies** |
| non-Latin | in window, unsupported | no | no |
| non-Latin | outside window | no | no |

The support condition is not applied to Latin genera: `g__Noviluteimonas` (r237, `91:g__Noviluteimonas|red=0.9579`,
bootstrap below 95) would lose its new genome G057139745 again, which is the case D2a was written for.

### Measurement (2026-11-23, made before the job, in `wip/`)

A variant of `3dc483a` with this rule (`wip/bin/jahelv01_opt1`, outputs in `wip/runs/jahelv01_opt1`, log at the end
of `wip/task.md`). Curator workload (`ops_no_rename`, `docs/testing.md` "Curator workload at `3dc483a`"), with the
previous release as the prior (r220 for r226, r226 for r232):

| Version | r226_ar53 | r226_bac120 | r232_ar53 | r232_bac120 | Total |
| --- | ---: | ---: | ---: | ---: | ---: |
| `3dc483a` | 179 | 2,755 | 482 | 5,127 | 8,543 |
| This rule | 181 | 2,800 | 480 | 5,366 | 8,827 (+284) |
| D2a for every genus, window only (rejected) | 210 | 3,349 | 570 | 6,448 | 10,577 (+2,034) |

All of the change is at genus. Of the genera it changes, the truth equals `3dc483a` for 744 and this rule for 460
(r226_bac120: 233 and 189). The r232 truth was probably curated from the Python output, which made these moves down,
so it favours `3dc483a`. In r237 bac120, G057028815 stays in `g__JAHELV01`, genera moved down go from 802 to 597, and
novel genera from 4,677 to 4,456.

### Tasks

1. `Tree::refine_existing_labels` and `Tree::species_genera_of_new_genomes` (`src/tree.rs`): drop the
   `is_latin_name` filter from the species genera, and choose them per taxon instead of per node: use them when the
   node is in the window and (`taxon.is_latin_name()` or `original_bs >= min_bs`). Update the doc comments.
2. Tests in `src/tree/tests.rs`. The fixture `genus_above_a_new_genome_at` hard-codes the bootstrap `91` of the
   genus node: add a bootstrap parameter. Then:
   - a non-Latin genus on a supported node in the window is not moved down away from a new genome its species names
     (the `g__JAHELV01` case);
   - a non-Latin genus on an unsupported node is still moved down
     (`a_genus_with_a_placeholder_name_is_moved_down_away_from_a_new_genome_its_species_names`, bootstrap 91, should
     pass unchanged: rename it to say "unsupported");
   - a Latin genus on an unsupported node in the window is still kept
     (`a_genus_is_not_moved_down_away_from_a_new_genome_that_its_species_names`, bootstrap 91, unchanged);
   - a non-Latin genus on a supported node outside the window is moved down.
3. Docs: `docs/refinement.md` Rule D2a (conditions, and the workload paragraph), `docs/differences.md` L9 (table row
   and section, with the `g__JAHELV01` example), and `docs/testing.md` (the new workload numbers, as for Job 15).
4. Measure with the steps in `docs/testing.md` "Running it" (Rust only: the Python outputs don't change). The total
   should be 8,827 (181 / 2,800 / 480 / 5,366); if not, find out why before committing. Check r237 bac120:
   G057028815 in `g__JAHELV01`, and no novel `g__JBXLYT01`.
5. `cargo fmt`, `cargo clippy --all-targets`, `cargo nextest run`, then commit (`feat(tree): ...`).

## Open questions (ask the user before starting)

- Moving down (Job 13) adds work as well as saving it: the operations added against the input went from 764 at
  `32d3350` to 1,897 with Jobs 10 to 13 (Python adds 1,814), and the taxon specificity in r226_bac120 dropped from
  99.95 to 99.56. Should the added operations be looked at?
- The effect of Job 10 was only measured without a prior (+3 operations, r226_bac120 only). It was not measured
  separately with Job 13.
