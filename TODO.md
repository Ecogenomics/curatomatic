# TODO

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


## Job 5: Replace `--meta` with `--novel-name-file PATH --novel-name-column NAME`

### Problem
`--meta` needs a three column file (`accession`, `ncbi_wgs_formatted`, `gtdb_taxonomy`) just to get the WGS prefix used
to name novel taxa. Building it means joining the GTDB metadata with a taxonomy file first
(`scripts/create_input_file.py`). The GTDB metadata file already has the WGS column, so curatomatic should read it
directly and let the caller choose the column.

### Input format
`test_data/r237/gtdb_r237_bac120_metadata.tsv` with `--novel-name-column ncbi_wgs_formatted`: tab separated, with a
header, many columns. Only `accession` and the named column are read; the rest are ignored.
- Find both columns by header name. Error (with the path) when `accession` or the named column is missing.
- Values `none`, `n/a`, `null` or empty (after trimming, case-insensitive) mean "no name" and are ignored. They are not
  an error (empty is an error now). In the r237 bac120 file 15,313 rows are `none`.

### Where
- Rust: `src/main.rs` (`Cli::meta` → `novel_name_file` + `novel_name_column`, the "Reading metadata" block),
  `src/metadata.rs` (`from_file`/`from_reader` take the column name, `not_none` handles the new placeholder set),
  `src/error.rs` if messages mention the metadata file.
- Tests: `src/metadata/tests.rs`, `src/main/tests.rs`.
- Docs: `README.md` (`--meta` usage and the "Metadata file (--meta)" section), `example/meta.tsv` and the example
  command.

### What to do
1. Add the two flags. `--novel-name-column` is required when `--novel-name-file` is given (clap `requires`); remove
   `--meta`.
2. Read the named column in place of the hard-coded `ncbi_wgs_formatted`, with the placeholder rules above.
3. Tests: custom column name, missing column error, and each placeholder (`none`, `N/A`, `null`, empty) ignored.

### Prior taxonomy
`--meta` also supplies `gtdb_taxonomy`, which drives downward refinement, the report's prior taxa
(`displaced_from_prior`, `new_in_decoration`) and the taxonomy coverage log lines. This moves to `--taxonomy` (Job 6):
the novel name file supplies only the WGS names.

### Done when
`curatomatic` runs on `test_data/r237` with `--novel-name-file test_data/r237/gtdb_r237_bac120_metadata.tsv
--novel-name-column ncbi_wgs_formatted`, novel taxa get WGS names, and no `--meta` remains in the code or README.


## Job 6: Add `--taxonomy` to name the ranks the tree can't hold

### Problem
Input trees used to put a dummy leaf next to each genome (e.g. `('G028276145|RED=1.0','D-G028276145|RED=1.0')`), so
that the species could be written on their parent node. The new trees (`test_data/r237/median_red_marked.nwk`) have no
dummy leaves and no `s__` labels, so there is nowhere in the tree for the species, or for any rank below the lowest
one decorated above a genome. These names must come from a taxonomy file, and be stored on the leaf.

The r237 bac120 tree has 229,307 leaves, none with a species. The lowest rank above each leaf is: genus 209,791,
family 17,507, order 1,463, class 449, phylum 71, domain 26.

### Input format
`test_data/r237/final_taxonomy.tsv`: tab separated, no header, two columns: the genome ID and a full seven rank
taxonomy string, e.g.
```
G003697165	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
```
- Every species name is unique to one genome.
- Some ranks are empty placeholders (`g__` 5,067, `f__` 668, `o__` 139, `c__` 31, `p__` 3). A placeholder is not a
  name: it is never assigned.
- Match genome IDs the same way as the metadata file (`gid_to_canonical`).
- A tree leaf (not a zombie) missing from the file is an error. A file genome missing from the tree is a warning.

### Rules
- **The tree always wins.** Never change or remove a taxon already in the tree, even when the file disagrees
  (r237 bac120: 1,148 genera, 304 families, 69 orders, 13 classes and 7 phyla differ).
- For each leaf, go up the tree to the lowest rank assigned on the path. Assign every rank **below** it from the
  file to the leaf. E.g. if the lowest rank above the leaf is `p__`, the leaf gets `c__` to `s__` from the file;
  if it is `g__`, the leaf gets only `s__`. Don't check the names for consistency with the tree (e.g. a species of
  `g__Neosphaerochaeta` under `g__JAAYFO01` is kept as is).
- A leaf with a lowest rank of `s__` (an old tree with dummy leaves) gets nothing.
- This runs first, straight after reading the tree and before any refinement.

### Output
A new file, `<out_dir>/<stem>.curatomatic.taxonomy.tsv`, written after all refinement steps, using the same format as
the input: tab separated, no header, the accession then the taxonomy string, e.g.
```
G003697165	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
```
- One row per leaf in the final tree (zombies skipped), in tree order.
- The taxonomy is read from the final tree: the taxa on the leaf and on every ancestor, ordered domain to species,
  joined with `;` (no spaces).
- A rank that appears on more than one node of the path (nested taxa) uses the one nearest the leaf.

### Where
- Rust:
  - `src/node_label.rs`: leaves must hold a taxonomy like internal nodes do (`from_node`/`to_label`; `to_label`
    currently rejects it). Pick a leaf label format that reads back the same, e.g. `'<gid>|<taxonomy>|RED=<red>'`.
  - `src/node_label.rs`, `src/tree.rs` (`RED_SEPARATOR`, `split("|RED")`), `parse_node_label`: accept both `|RED=` and
    `|red=`. The r237 leaves have no RED at all: treat a missing leaf RED as `RedValue::LEAF`.
  - `src/tree.rs`: a new step that assigns the file ranks to the leaves. Remove the `s__` placeholder hack in
    `Tree::from_newick` if nothing needs it. Check every place that assumes a leaf has no taxonomy:
    `get_leaf_nodes_missing_rank` (starts at the parent), `taxon_nodes`, the refinement checks (e.g. "cannot be moved
    to the leaf node"), `fill_in_missing_ranks`, and the report (`src/report.rs`).
  - A new module (e.g. `src/taxonomy_file.rs`) to read the file, using `Taxonomy::from_string`.
  - `src/main.rs`: the `--taxonomy PATH` argument. It replaces `gtdb_taxonomy` from `--meta` as the prior taxonomy,
    which drives downward refinement, the report's prior taxa and the taxonomy coverage log lines
    (`taxonomy_coverage_messages`). See Job 5. Write the output taxonomy file in `curate`, next to
    `path_tree_out`/`path_report_out`.
- Tests: `src/node_label/tests.rs`, `src/tree/tests.rs`, `src/main/tests.rs`, and the tests for the new module.
- Docs: `README.md` (usage, input files).

### What to do
1. Let a leaf hold a taxonomy, and accept `|red=` and a leaf without RED.
2. Read the taxonomy file, with the ID and placeholder rules above.
3. Add the leaf assignment step and run it first.
4. Make the rest of the pipeline treat a taxon on a leaf the same as one on an internal node.
5. Point the prior taxonomy at `--taxonomy`.
6. Write the output taxonomy file after the last refinement step.
7. Tests: a leaf gets the ranks below its lowest ancestor rank and nothing above it, the tree's taxa are never
   changed, placeholders are skipped, an old style tree gets nothing, a missing leaf is an error, an extra genome is a
   warning, and leaf labels round trip with a taxonomy and with both RED spellings. For the output file: ranks come from the leaf and its
   ancestors in rank order, the nearest of two nested taxa wins, zombies are skipped, and the format matches the
   input file.

### Open question
The rule above fills only the ranks below the lowest one in the tree. A leaf can still have a gap higher up: in
r237 bac120, 5,155 leaves have no phylum on their path but only 26 have nothing below the domain, so about 5,100 of
them miss a phylum above a decorated class or lower (likewise class, order and family). These stay missing after this
step, and `fill_in_missing_ranks` will create novel taxa for them. Decide whether that's right, or whether the file
should fill these gaps too.

### Done when
`curatomatic` runs on `test_data/r237` with `--taxonomy test_data/r237/final_taxonomy.tsv`. Before refinement, every
leaf has every named rank from the file below the lowest rank in the tree, including its species. After the full run, every leaf has one taxon at each of the seven ranks on its path,
and `<stem>.curatomatic.taxonomy.tsv` has one seven rank row per genome that matches the final tree.
`cargo test` passes, and none of the tree's taxa, bootstrap values or structure were changed.
