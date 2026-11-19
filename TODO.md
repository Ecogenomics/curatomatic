# TODO


## Job 1 (done): Write a review table for curators

Done in `src/report.rs`. Notes on how the spec was read:
- Several flags can apply to one taxon: `flag` holds all of them, comma separated. Three extra columns:
  `prior_taxon` (`displaced_from_prior`), `nested_in` (`nested_in_same_rank`) and `sister_taxa`
  (`sister_join_past_child_thresh`, the highest taxa of the rank under the other children of the first join).
- `sister_join_past_child_thresh` does not require the join to be supported. A nested taxon has no sister.
- `displaced_from_prior` is set when the members had any prior taxon of the rank, not only exactly one.
- `new_in_decoration` uses every `gtdb_taxonomy` in the metadata file (before it is subset to the tree), and is
  skipped when the metadata has no taxonomies.
- Flags are for ranks p to g only. Rows are sorted by rank, then taxon. `absorbed` waits for an absorb rule.


### Problem
Some curator decisions can't be automated from RED and bootstrap, but curatomatic can point at them. The log is the
only record of what curatomatic changed, and it isn't machine readable.

### What to do
Write `<stem>.curatomatic.report.tsv` next to the output tree, one row per taxon that was changed or flagged:
`taxon, rank, action, flag, n_genomes_before, n_genomes_after, bootstrap_before, bootstrap_after, red_before,
red_after`. Actions: `moved_up`, `moved_down`, `absorbed` (Job 5), `novel`. Flags, all computed from data already in
memory:
- `nested_in_same_rank`: the taxon sits inside another taxon of the same rank (Job 2 skips these). Curators merged 8
  of the 13 found in bac120. See Job 9 for the pairs that are already nested in the input.
- `novel_low_red`: a novel taxon with two or more genomes whose node has `q < 0.5`. Curators split 38% of these
  genera in R226 bac120 (50/130) and 57% in R232 bac120 (130/227), against 2% at `q >= 0.5` in both.
- `sister_join_past_child_thresh`: two taxa of rank `r` whose common ancestor has RED above `child_thresh(r)`.
  For genera with a supported join, curators merged one into the other in 27% (R226 bac120) and 16% (R232 bac120)
  of cases, against about 1% for joins below `q = 0.5`. For phylum to order in R232 bac120 it was 56% (23/41).
- `unsupported_label`: the label stays on a node with bootstrap below `--min-bs`.
- `large_gain`: a move that multiplies the number of genomes by 5 or more. Curators kept 0 of 7 (R226 bac120) and 1
  of 12 (R232 bac120).

Flags added after the R232 review:
- `split_candidate`: an unmoved taxon on a supported node with `q < 0.5`, whose members from the previous release
  are all under one child while the other children hold only new genomes. Curators split 36–48% of these genera in
  R232 bac120 (Job 6, step 2).
- `displaced_from_prior`: a novel taxon whose members already had a taxon of that rank in the previous release
  (Job 8). Include the previous name in an extra column.
- `new_in_decoration`: a taxon that is on the input tree but not in any `gtdb_taxonomy` value in the metadata
  (names added before curatomatic runs, mostly Latin names from NCBI). Only 33–38% of these were kept unchanged in
  bac120 (42 of 129 in R232, 44 of 117 in R226), against 92–93% for names carried over from the previous release.

The report must not change the tree. Python and Rust must write identical files.


## Job 2: Name fallback novel taxa `<rank>__<accession>`, not `<rank>__NOVEL-<accession>`

### Problem
When no unused WGS prefix is available, `generate_unique_rank` names the taxon `g__NOVEL-G963674675`. Curators use the
bare accession: `g__G963674675`. In bac120, 298 novel taxa got a `NOVEL-` name; for 187 of them the curated name is
exactly the name with `NOVEL-` removed, on the same genomes (the other 111 were merged into a neighbour or given a
Latin name). ar53 has only 4 `NOVEL-` names, none kept.

### R232 check
Still valid. R232 bac120 has 365 `NOVEL-` names: 206 were curated to the bare accession on the same genomes, 72 were
merged into a neighbour and 87 got another name. ar53: 2 of 9. Many of the "other name" cases are WGS-style names
(`g__NOVEL-G020072925` became `g__JAIQGA01`), which means the genome has a WGS prefix that the metadata file didn't
supply. 13,146 of the 189,801 R232 bac120 genomes have an empty `ncbi_wgs_formatted`.

### Where
- Rust: `src/taxon.rs` → `generate_unique_rank` (last line).
- Python: `curatomatic/util/taxon.py` (the same function).

### What to do
1. Emit `{prefix}__{gid}`. If that name is already in `known_taxa`, fall back to the next lowest genome ID in the
   group, and only then to a suffixed name. Keep the result deterministic.
2. Update `src/taxon/tests.rs`.

### Done when
No `NOVEL-` appears in any output tree, and the bac120 exact count rises by about 187.


## Job 3: Replace the hand-rolled `thousands` with the `thousands` crate

### Problem
`src/main.rs` → `thousands` formats counts with comma separators to match Python's `f'{n:,}'`. The standard library
can't group digits, so the function was written by hand. The `thousands` crate does the same thing, so there's no
need to maintain our own.

### Where
- Rust: `src/main.rs` → `thousands` and its callers in `curate` and `taxonomy_coverage_messages`.
- Tests: `src/main/tests.rs` → `thousands_matches_python_comma_formatting`.

### What to do
1. `cargo add thousands`.
2. Delete `thousands` and call `n.separate_with_commas()` (`use thousands::Separable;`) at each call site.
3. Keep the test cases from `thousands_matches_python_comma_formatting`, but have them call `separate_with_commas`,
   so they still pin the output to Python's (always a comma, never locale-dependent).

### Done when
`thousands` is gone from `src/main.rs`, `cargo test` passes, and the log lines for every test_data run are
unchanged.


## Job 4: Read the median RED values from a TSV instead of the JSON `.dict`

### Problem
`Red::from_file` reads the median RED for each rank from a JSON object, e.g.
`test_data/r226_ar53/gtdb_r226_ar53.decorated.dict`:
`{"phylum":0.231,"class":0.356,"order":0.529,"family":0.723,"genus":0.907,"species":1.0}`. The RED tooling now writes
these values as a TSV, so curatomatic should read that instead.

### Input format
`test_data/r237/median_red_rank_distribution.tsv`: tab separated, with a header and one row per rank:
```
rank	num_taxa	num_taxa_pass_filter	median	mean	std_dev	p1	p5	p10	p90	p95	p99
Phylum	161	46	0.3095	0.3314	0.0692	0.2261	0.2368	0.2549	0.4243	0.4684	0.4987
...
Genus	20066	19408	0.9189	0.9175	0.0394	0.8208	0.8508	0.8657	0.9694	0.9804	0.9918
```
- Look up the `rank` and `median` columns by header name, not position. Ignore the other columns.
- Rank names are capitalised (`Phylum`); `Rank::of_taxon` already uses only the first character, but match the name
  case-insensitively and reject names that aren't a rank.
- There is no `Domain` or `Species` row. Domain stays `0.0` (as now). Species must be `1.0` (`RedValue::LEAF`):
  `fill_in_missing_ranks`/upward refinement for genera call `red_dict.get(Species)` via `rank.child()`, and the old
  `.dict` supplied `species: 1.0`.

### Where
- Rust: `src/red.rs` → `Red::from_file` / `Red::from_json` (replace with `from_tsv`), `src/main.rs` (`red_dict` arg doc
  and the "Reading RED dictionary" log line), `src/error.rs` (`Red`/`MissingRed` messages say "dictionary").
- Tests: `src/red/tests.rs`, and the `from_json` call in `src/tree/tests.rs`.

### What to do
1. Parse the TSV into `Red` with the rules above. Error (with the path) on a missing `rank`/`median` column, an
   unparsable median, an unknown rank, or a duplicate rank.
2. Drop the JSON parser and `serde_json` if nothing else uses it.
3. Update the tests to build `Red` from TSV text, including the missing-column, bad-value and duplicate-rank errors.
4. Convert the `.dict` files under `test_data/r226_*` to this TSV format (rank + median columns are enough) so the
   existing runs still work, and update the README/usage if they mention the `.dict`.

### Done when
`curatomatic` runs on `test_data/r237` with `median_red_rank_distribution.tsv`, the `debug` line for the RED values
lists every rank including `s=1.000`, and the r226 runs give the same trees as with the old `.dict` files.
