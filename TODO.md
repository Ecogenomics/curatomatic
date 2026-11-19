# TODO

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
