# TODO

Finished jobs are removed from this file: their rules are in `docs/refinement.md` and `docs/differences.md`, and
their measurements in `docs/testing.md`. Jobs that need the user's approval say so: don't start them before it is
given.

## Job 14: Report a namesake genome outside its taxon

**Status: approved.**

### Problem

`generate_unique_rank` (`src/taxon.rs`) names a novel taxon without a WGS name after its lowest genome ID whose name
isn't already a taxon. If every genome's name is taken, the taxon is `<rank>__UNRESOLVABLE`. This happens when the
input tree already has a taxon named after one of the genomes, but that genome isn't inside it. For example,
`G020350785` in r226_bac120 gets `g__UNRESOLVABLE`, because `g__G020350785` is already on the tree elsewhere (the
curators put the genome into `g__G020350785`). In r237 bac120, G965247515 gets `g__UNRESOLVABLE` the same way. Refusing
the name is right, but the report doesn't say why.

### Rule

When a novel taxon is `UNRESOLVABLE`, its report row gets a new flag, `namesake_outside_taxon`, and names the
existing taxa that are named after its genomes (e.g. `g__G020350785`). Nothing else changes.

### Tasks

1. Make `generate_unique_rank` (or its caller in `fill_in_missing_ranks`) return the taxa that blocked each name, and
   pass them to the report (`src/report.rs`, new `Flag` variant with a doc comment).
2. Tests in `src/taxon/tests.rs` and `src/report/tests.rs`.
3. Docs: `docs/refinement.md` "Novel taxa", the report columns in `README.md`.
4. Check the r226_bac120 cases (`G020350785`, `G030406045`), the r232_bac120 case (`G020350445`) and the r237 bac120
   case (`G965247515`).

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
- The curator workload, the user's metric of success, is scored by `wip/scripts/workload.py` and
  `wip/scripts/summary.py`, which are not committed (`wip/` is git-ignored).

### Tasks

1. Make the prior files in `benchmark.sh` from the Python metadata (recipes in `docs/testing.md`, "Prior taxonomy")
   and pass them with `--prior-taxonomy`.
2. Move the workload scorer into `scripts/`, run it in `benchmark.sh`, and document its metrics in
   `docs/testing.md`.
3. Check that `benchmark.sh` gives the numbers in `docs/testing.md`, "Curator workload".

## Job 22: Should the output taxonomy contain placeholders?

**Status: needs the user's decision. Discuss before starting.**

### Problem

The user does not want placeholders written to the output taxonomy (`<stem>.curatomatic.taxonomy.tsv`). Today
`Tree::leaf_taxonomies` (`src/tree.rs`) builds each output row one rank at a time, and writes the bare prefix
(`o__`) for a rank with no named taxon on the genome's path. `Tree::full_taxonomies`, which writes the `--debug`
taxonomy, does the same.

In a normal run every such placeholder comes from a novel taxon that `fill_in_missing_ranks` skipped, because a
valid name of a lower rank is at or below its node (`docs/differences.md` L7, `docs/refinement.md` "Novel taxa":
"the genomes of a skipped taxon keep the placeholder for that rank only"). A missing domain stops the run, and every
genome has its species.

### Options to discuss

1. Leave the rank out of the row (e.g. no `o__`). Rows then have varying lengths, which may break downstream tools
   that expect 7 ranks.
2. Always create the novel taxon: drop the L7 skip, as Python did. This brings back the curator work L7 removed.
3. Write a marked name for a curator to resolve, e.g. `o__UNRESOLVABLE`, the existing convention for a novel taxon
   that can't be named.
4. Keep placeholders in the output only (they are banned on the tree and stripped from the inputs, see `AGENTS.md`).

Also decide the `--debug` taxonomy: `scripts/check_round_trip.py` compares it with the input file, which has
placeholders, so the script must change if the debug output does.

### Tasks (once decided)

1. Count how many genomes and ranks get a placeholder in the r237, r226 and r232 outputs. The log prints
   `Skipped N novel <rank>__ taxa that have a valid name below them.` for each rank.
2. Implement the chosen option, with tests, and update `README.md` ("Output") and `docs/refinement.md`.

## Job 23: Require a full seven rank taxonomy in every row of `--prior-taxonomy`

**Status: approved** (requested by the user on 2026-11-23).

### Problem

`--prior-taxonomy` is read by the same parser as `--taxonomy` (`TaxonomyFile::from_reader` and `parse_row`,
`src/taxonomy_file.rs`), so it is as lenient. Placeholders are stripped (`d__Bacteria;p__;c__;…;s__X y` loads as
`d__Bacteria;s__X y`), a rank can be missing (`d__B;p__P;s__X y`), and a rank can hold several taxa
(`g__A;g__B`). It is only required to be in rank order, give each taxon once, and have at most one domain and one
species (commit `da1a5bf`). A prior taxonomy is a previous release, where every genome is fully classified, so a
row that is not is a broken input that should stop the run, not a gap downward refinement silently works around.

### Rule

Every row of `--prior-taxonomy` (including genomes that are not on the tree) must be exactly seven taxa, one per
rank, in order: `d__`, `p__`, `c__`, `o__`, `f__`, `g__`, `s__`, each with a name. Anything else is a
`CuratomaticError::TaxonomyFile` error on that line, naming the problem, e.g.:

- a placeholder: `d__B;p__;c__C;…` → `the phylum rank is the placeholder p__, but a prior taxonomy names every rank`
- a missing rank: `d__B;p__P;c__C;o__O;f__F;s__X y` → `the genus rank is missing, …`
- two taxa of a rank: `…;g__A;g__B;s__X y` → `the genus rank has g__A and g__B, …`
- out of order or too many taxa: already caught by the shared checks, keep their messages.

The species must be a valid binomial name (decided by the user on 2026-11-23): exactly two words separated by one
space, where the first word is the name of the row's genus (`g__Escherichia;s__Escherichia coli`,
`g__CAG-698;s__CAG-698 sp000123`). So `s__{unresolved} NAME` is an error in the prior, though it is still allowed in
`--taxonomy`. Errors, e.g.:

- `s__{unresolved} sp1` → `the species s__{unresolved} sp1 is not a binomial name, but a prior taxonomy names every
  species`
- `s__Escherichia` (one word) or `s__Escherichia coli K12` (three words) → the same message.
- `g__Escherichia;s__Salmonella enterica` → `the species s__Salmonella enterica is not in the genus g__Escherichia`

`--taxonomy` keeps its current rules: placeholders stripped, missing ranks and nested taxa allowed, and unresolved
species allowed.

### Tasks

1. Give the prior its own reader, checking the raw row before placeholders are stripped. Prefer a newtype (e.g.
   `PriorTaxonomyFile` wrapping `TaxonomyFile`, or a row format parameter that `from_reader` takes) so that
   `main.rs` (`let prior = match &cli.prior_taxonomy`, around line 294) can't read the prior with the lenient reader.
   `read_metadata` takes `Option<&TaxonomyFile>`: update its signature to match.
2. Tests in `src/taxonomy_file/tests.rs`: a full row loads, and one test per error above (placeholder in each
   position including `d__` and `s__`, a missing rank, two taxa of a rank, a row of placeholders only, a domain-only
   row, an unresolved species, a one word and a three word species, a double space, and a species whose first word
   is not the row's genus). Check that `--taxonomy` still accepts all of these, apart from the ones the shared checks
   reject.
3. Docs: `AGENTS.md` ("Taxonomy files"), `README.md` (the `--prior-taxonomy` option and "Taxonomy file"), and
   `docs/differences.md` (the placeholder table at the "In a row of `--taxonomy` / `--prior-taxonomy`" row).
4. Check the real prior files still load (all passed this rule on 2026-11-23, 0 rows failing, the binomial and
   genus checks included):
   `example/gtdb_r232_bac120_taxonomy.tsv`, `test_data/curated/*.tsv`, and the four files that the recipe in
   `docs/testing.md` ("Prior taxonomy") makes from `test_data/python/input/*.metadata.tsv`. The r232 recipe drops
   placeholder-only rows; if it ever kept a partly placeholder row, the run would now stop, so say so there.
   `test_data/r237/final_taxonomy.tsv` fails (5,070 rows with `g__;s__{unresolved} …`), but it is the curated truth
   of r237, not a prior.
5. Run the benchmark with a prior (`docs/testing.md`, "Curator workload") and check that the numbers do not change.

## Open questions (ask the user before starting)

- Moving down (Job 13) adds work as well as saving it: the operations added against the input went from 764 at
  `32d3350` to 1,897 with Jobs 10 to 13 (Python adds 1,814), and the taxon specificity in r226_bac120 dropped from
  99.95 to 99.56. Should the added operations be looked at?
- The effect of Job 10 was only measured without a prior (+3 operations, r226_bac120 only). It was not measured
  separately with Job 13.
