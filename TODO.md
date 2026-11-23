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

## Job 20: No placeholders on the tree, and placeholders stripped from the input files

**Status: approved** (decided with the user on 2026-11-23). Job 19, which it needed, is done.

### Problem

Placeholders (`p__`, `g__`, ...) are not names, but the code carries them around:

- A placeholder in a tree label is accepted. None of the trees in `example/` or `test_data/` has one (checked
  2026-11-23), but if one did: two nodes with the same placeholder fail to load with `NonUniqueTaxon`
  (`Tree::get_taxon_labels` only exempts `UNRESOLVABLE`). Refinement would move a placeholder like a name. A
  placeholder on a genome's path stops `Tree::get_leaf_nodes_missing_rank` from reporting that rank as missing. And
  `check_taxa_unchanged` would not notice if one disappeared.
- The rows of `--taxonomy` and `--prior-taxonomy` keep their placeholders. Every reader then filters them out again
  with `Taxon::is_named`.

### Rule

- **Tree:** a placeholder in any node label is an error that stops the run (new error,
  `CuratomaticError::PlaceholderOnTree`, naming the label). `d__Bacteria|RED=0.000` is fine; `100:g__|RED=0.9`,
  `100:f__F;g__|RED=0.9` and `d__|RED=0.0` are not.
- **Taxonomy files** (`--taxonomy` and `--prior-taxonomy`, both read by `TaxonomyFile::from_reader`): each row is
  checked, then its placeholders are stripped:
  1. A rank that holds both a placeholder and a name (`g__;g__X`) fails with a `TaxonomyFile` error on that line.
     This must be checked before stripping.
  2. The placeholders are removed, so the row's `Taxonomy` holds named taxa only. A row whose domain or species is a
     placeholder is caught later, by the leaf checks (Job 21 and `Tree::check_leaves_have_species`).
- Placeholders are not assigned to leaves. They already aren't, and with stripping there are none to assign.
- Suggested approach: make it impossible for a `Taxonomy` to hold a placeholder (`Taxonomy::add_taxon` refuses one).
  `TaxonomyFile` strips before building the row, and `NodeLabel` turns the refusal into `PlaceholderOnTree`. The
  output rows built by `Tree::leaf_taxonomies` and `Tree::full_taxonomies` are `Vec<Taxon>`, not `Taxonomy`, so
  they are not affected (see Job 22).
- The `is_named` filters on tree taxa (`would_cause_conflicts`, `nested_inside` and others) become redundant. Keep
  them unless removing one is trivial; this job must not change behaviour on valid input.

### Tasks

1. `src/node_label.rs`, `src/taxonomy_file.rs`, `src/taxonomy.rs`, `src/error.rs`.
2. Tests:
   - `src/node_label/tests.rs`: `a_placeholder_on_a_tree_label_is_an_error`, with the three labels above.
   - `src/taxonomy_file/tests.rs`:
     - `placeholders_are_stripped_from_a_row`: `G1\td__B;p__;c__;o__;f__;g__;s__X y` holds only `d__B` and
       `s__X y`.
     - `a_rank_with_a_placeholder_and_a_name_is_an_error`: `g__;g__X` gives `TaxonomyFile { line, .. }`.
   - Update the existing tests that expect placeholders in a `TaxonomyFile` row (e.g. `src/taxonomy_file/tests.rs`
     line 29).
3. Docs: replace the AGENTS.md text saying the tree may hold placeholders (Taxonomy section) with this rule. Also
   update "Input preparation" in `docs/refinement.md`, the taxonomy file format in `README.md`, and add an entry to
   `docs/differences.md`.
4. Check the r237 run and the round trip of the r226 and r232 Python trees (README, "Output") still pass with 0
   differences.

## Job 21: Each genome's taxonomy on the tree must match its row, and include a named domain

**Status: approved** (decided with the user on 2026-11-23, option A). Needs Job 20.

### Problem

`Tree::assign_leaf_taxonomy` (`src/tree.rs`) only checks the **anchor**: the nearest node above a leaf that holds
a named taxon. Its taxa must be in the genome's row. Taxa further up the path are never compared with the row, so
these pass silently today. Both use the tree `(((G1,G2)g__X,G3)f__F,G4)d__Bacteria`:

- **A name in the row that disagrees with a name on the tree.** G1's row is `d__Bacteria;f__Z;g__X;s__X one`. The
  anchor `g__X` is in the row. `f__Z` is not on G1's path, so G1's leaf is given `f__Z`, and G1 ends up in both
  `f__F` and `f__Z`.
- **A placeholder in the row where the tree has a name.** G3's row is `d__Bacteria;f__;g__;s__Y three`. After
  stripping (Job 20) the row has no family, but the tree puts G3 in `f__F`. The anchor check only rejects a
  *named* anchor taxon that is missing from the row.

Nothing checks that a genome has a named domain, though `AGENTS.md` says every genome has one.

### Rule

- After the leaves are assigned, for every leaf: the named taxa on its path plus those on its leaf must **equal**
  the taxa of its row (placeholders stripped). They are matched by taxon, as sets: never by rank, and order is
  ignored. Any difference fails with `TaxonomyFileConflict`, naming the genome and the taxa that are on the tree but
  not in the row (and the reverse, though assignment should leave none). Both examples above fail.
- Keep the existing anchor check, which gives a clearer message for the case it covers.
- Rows for genomes that are not on the tree are not checked; they stay a warning.
- **Domain:** every leaf's path plus its leaf must hold a named `d__` taxon. This is checked at the same two points
  as the species, after assignment and before output, and fails with a new `CuratomaticError::LeafWithoutDomain
  { n, first }`, like `LeafWithoutSpecies`. `--prior-taxonomy` rows are not checked.
- This makes the AGENTS.md rule "removing the dummy nodes and assigning the file must give back exactly the taxa
  they held" a check that runs every time, not only a test and a script.

### Tasks

1. `src/tree.rs` (`assign_leaf_taxonomy`, a domain check next to `check_leaves_have_species`), `src/main.rs`,
   `src/error.rs`.
2. Tests in `src/tree/tests.rs`, written first and seen to fail:

   ```rust
   #[test]
   fn a_tree_taxon_above_the_anchor_that_is_not_in_the_row_is_an_error() {
       // The tree puts G1 in p__A, its row says p__B. The anchor (g__X) agrees, so only a check of
       // the whole path catches it. Before this job G1 silently got p__B on its leaf.
       let newick = "((G1:0.1,G2:0.1)'100:g__X|red=0.9':0.1,G3:0.1)'d__D;p__A|red=0.0';";
       let taxonomy = "\
   G1\td__D;p__B;c__;o__;f__;g__X;s__X one
   G2\td__D;p__A;c__;o__;f__;g__X;s__X two
   G3\td__D;p__A;c__;o__;f__;g__;s__Y three
   ";
       assert!(matches!(
           conflict(newick, taxonomy),
           CuratomaticError::TaxonomyFileConflict { gid, message }
               if gid.as_str() == "G1" && message.contains("p__A")
       ));
   }

   #[test]
   fn a_row_placeholder_where_the_tree_has_a_name_is_an_error() {
       // The tree puts G3 in f__F, its row has no family. The anchor is f__F, but the anchor check
       // only rejects a named taxon missing from the row.
       let newick = "(((G1:0.1,G2:0.1)'100:g__X|red=0.9':0.1,G3:0.1)'100:f__F|red=0.7':0.1,G4:0.1)'d__D|red=0.0';";
       let taxonomy = "\
   G1\td__D;p__;c__;o__;f__F;g__X;s__X one
   G2\td__D;p__;c__;o__;f__F;g__X;s__X two
   G3\td__D;p__;c__;o__;f__;g__;s__Y three
   G4\td__D;p__;c__;o__;f__;g__;s__Z four
   ";
       assert!(matches!(
           conflict(newick, taxonomy),
           CuratomaticError::TaxonomyFileConflict { gid, message }
               if gid.as_str() == "G3" && message.contains("f__F")
       ));
   }
   ```

   Also:
   - `every_genomes_full_taxonomy_matches_its_row_after_assignment`: for every existing assignment fixture,
     including the nested order and out-of-rank-order r226 cases, the named taxa of `Tree::full_taxonomies` equal
     the row. This guards the correct nested cases.
   - `a_leaf_without_a_named_domain_is_an_error`: a row with `d__` and no domain on the tree gives
     `LeafWithoutDomain`.
   - `src/main/tests.rs`: a run whose taxonomy file disagrees above the anchor fails with `TaxonomyFileConflict`
     and writes no tree, taxonomy or report.
   - Give a domain to the existing fixtures that have none, e.g. the rows with `d__` in `src/tree/tests.rs`
     (around lines 463 and 506) and their trees.
3. Docs: `AGENTS.md` (the taxonomy file and species rules), "Input preparation" and "Invariants" in
   `docs/refinement.md`, and `README.md`.
4. Check that the r237 run and the round trip of the r226 and r232 Python trees still pass. They round trip with 0
   differences today, so the new check should not fire.

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
4. Keep placeholders in the output only (they are banned on the tree and stripped from the inputs by Job 20).

Also decide the `--debug` taxonomy: `scripts/check_round_trip.py` compares it with the input file, which has
placeholders, so the script must change if the debug output does.

### Tasks (once decided)

1. Count how many genomes and ranks get a placeholder in the r237, r226 and r232 outputs. The log prints
   `Skipped N novel <rank>__ taxa that have a valid name below them.` for each rank.
2. Implement the chosen option, with tests, and update `README.md` ("Output") and `docs/refinement.md`.

## Open questions (ask the user before starting)

- Moving down (Job 13) adds work as well as saving it: the operations added against the input went from 764 at
  `32d3350` to 1,897 with Jobs 10 to 13 (Python adds 1,814), and the taxon specificity in r226_bac120 dropped from
  99.95 to 99.56. Should the added operations be looked at?
- The effect of Job 10 was only measured without a prior (+3 operations, r226_bac120 only). It was not measured
  separately with Job 13.
