# Differences between the Python (`main`) and Rust (`dev`) implementations

This compares the Python package on the `main` branch (`curatomatic/`) with the Rust rewrite on
`dev` (`src/`). It lists every difference in behaviour, with **logic changes first**. Each item
names the code on both sides.

See [refinement.md](refinement.md) for the full rules that the Rust version follows.

## Summary of logic changes

| # | Area | Python (`main`) | Rust (`dev`) | Changes output? |
| --- | --- | --- | --- | --- |
| L1 | Refining up | A nested taxon can move up | A nested taxon never moves up | **Yes** |
| L2 | Refining down | Prior taxonomy from the `gtdb_taxonomy` column of `--meta` | Prior taxonomy from `--prior-taxonomy` | **Yes**, if a prior taxonomy was given |
| L3 | Species placeholder | Replaces the **whole** label of the node with `s__{unresolved} <gid>` | Not rewritten: species come from `--taxonomy`, onto the leaves | **Yes**, for trees with an `s__` placeholder |
| L4 | Novel names without WGS | `<rank>__NOVEL-<lowest gid>`, may duplicate an existing taxon | `<rank>__<lowest gid not already used>`, else nothing is written and the taxon is reported as `unresolvable` | **Yes** |
| L5 | Novel domain | Uses the species median as the parent | Error | Stops the run instead of a wrong placement |
| L6 | Leaf taxa | Leaves can't hold a taxon, the tree needs a species node above each leaf | `--taxonomy` assigns the ranks the tree can't hold to the leaves | **Yes**, a new input |
| L7 | Novel taxa above names | A novel taxon is created for every missing rank, even above named taxa | Only created where no valid name of a lower rank is below it | **Yes** |
| L8 | Refining up | A taxon can be moved up above a taxon of a higher rank | The upward walk stops at a node with a higher rank below it | **Yes** |
| L9 | Refining down | Genomes without a prior taxonomy are in no genus | A new genome is in the genus its binomial species names, so a genus inside the RED window (and, if its name is not Latin, on a supported node) is not moved down away from it | **Yes**, if a prior taxonomy was given |
| L10 | Novel taxa out of rank order | A novel taxon can be created below a named taxon of a lower rank, or above one of a higher rank | Not created | **Yes** |

All other refinement and novel-rank rules are the same; see
[Unchanged logic](#unchanged-logic).

## Logic changes in detail

### L1. Nested taxa are not moved up

- **Python** (`Tree.refine_existing_labels`, `model/tree.py`): every taxon on a node with a
  bootstrap is walked up the tree. The only stop is a conflict below the candidate ancestor.
- **Rust** (`Tree::refine_existing_labels` / `Tree::nested_inside`, `src/tree.rs`): before
  walking up, it checks whether any ancestor holds a different **named** taxon of the same rank.
  If so the taxon is nested inside it, and the upward search is skipped entirely. The downward
  search still runs.

Why: the conflict check only looks *below* the candidate ancestor. A taxon nested under another
of the same rank could climb onto, or past, the node of the outer taxon and take its genomes. The
Python version did have a nested-rank detector (`Tree.get_nested_ranks`), but only in the old,
unreachable code path (see [Python code that was not ported](#python-code-that-was-not-ported)).

### L2. Moving down reads the prior taxonomy from `--prior-taxonomy`

The downward rules themselves are identical (see [Unchanged logic](#unchanged-logic)). What
changed is where the prior taxonomy comes from.

- **Python** (`Metadata.from_file`, `model/metadata.py`): `--meta` reads the `gtdb_taxonomy`
  column, if present, as each genome's prior taxonomy. That enables downward refinement.
- **Rust** (`main.rs`, `read_metadata`): the prior taxonomy is the `--prior-taxonomy` file (a
  previous release, in the format of `--taxonomy`, read with `src/taxonomy_file.rs`).
  `--novel-name-file` reads only `accession` and the `--novel-name-column`. Without
  `--prior-taxonomy` the log warns *"No prior taxonomy was given, downward refinement is
  disabled."* `--taxonomy` is never the prior: it is the decoration of the same tree, so no taxon
  could ever move down.

The r226 test data's `metadata.tsv` has a `taxonomy` column (r220), not `gtdb_taxonomy`, so the
stored r226 Python reference outputs did no downward refinement. The r232 metadata has
`gtdb_taxonomy` (r226), so the r232 Python outputs did move taxa down. To compare the two versions,
give both the same prior from the previous release (`docs/testing.md`, "Prior taxonomy").

### L3. The species placeholder is not rewritten

- **Python** (`Tree.from_file`, `model/tree.py`): a leaf's parent labelled with `s__` has its
  whole label set to `s__{unresolved} <gid>|RED=…`. Any other taxon on that node (e.g. `g__X; s__`)
  is **dropped**.
- **Rust**: the tree is read as is. No input tree has a species placeholder any more: the species
  of each genome comes from `--taxonomy`, and is put on the leaf (see L6). A placeholder of any
  rank on the tree is an error (see [Placeholders](#placeholders)).

### L4. Naming novel taxa without a WGS prefix

`generate_unique_rank` (`util/taxon.py` vs `src/taxon.rs`):

| | Python | Rust |
| --- | --- | --- |
| A WGS prefix is available | Most readable, then alphabetical, excluding known taxa | Same |
| No WGS prefix | `<rank>__NOVEL-<lowest gid>` | `<rank>__<gid>`: the lowest genome ID whose name isn't already a known taxon |
| Every genome ID already taken | `<rank>__NOVEL-<gid>` (not checked; may duplicate) | Not created: its genomes keep the placeholder (`g__`), and the report has a row `<rank>__UNRESOLVABLE` with action `unresolvable` and its `namesakes` |

The `NOVEL-` prefix was removed in commit `598ac97` ("remove __NOVEL"). Python doesn't check the
fallback name against the known taxa. If two novel taxa of the same rank have the same lowest
genome, they get the same name.

### L5. A novel domain cannot be placed

- **Python** (`fill_in_missing_ranks`, `util/tree.py`): the parent rank is
  `RANKS[RANKS.index(rank) - 1]`. For `d` this wraps to `RANKS[-1]`, i.e. the species (RED `1.0`),
  so the RED threshold is wrong.
- **Rust** (`fill_in_missing_ranks`, `src/tree.rs`): `parent_rank` returns `NoParentRank` for the
  domain, which stops the run. Refining the domain is also an error, but it's never attempted.

### L6. The ranks the tree can't hold are put on the leaves

- **Python**: a leaf is labelled with its genome ID only. The species (and any rank below the
  lowest one decorated above a genome) had to be on a dummy species node above each leaf.
- **Rust** (`Tree::assign_leaf_taxonomy`, `src/tree.rs`): a leaf can hold taxa, which are kept
  in the `Taxonomy` of its `Leaf`, never in its label. `--taxonomy` gives each leaf the taxa that
  were on its dummy node: every taxon of its row that is not on its path, before refinement (see
  [refinement.md](refinement.md#input-preparation)). Taxa can be moved down onto a leaf, and
  novel taxa placed on one. The curated taxonomy of every genome is written to
  `<stem>.curatomatic.taxonomy.tsv`, and the taxa held by each leaf to
  `<stem>.curatomatic.leaf-taxonomy.tsv` (Python has no such file: its leaf taxa were on the dummy
  nodes of its tree). Leaf labels are never changed, so the output tree labels
  each leaf with its genome ID, as in Python.

### L7. Novel taxa are not created above named taxa

- **Python** (`fill_in_missing_ranks`, `util/tree.py`): a novel taxon is created on the highest
  candidate node for every rank a genome is missing, even when named taxa of lower ranks are
  below it. In R237, low bootstrap support above `c__Chlorobiia` split it into six novel phyla.
- **Rust** (`fill_in_missing_ranks`, `src/tree.rs`): the novel taxon is skipped if a valid name
  (`Taxon::is_valid_name`: not a placeholder or `s__{unresolved}` species) of a lower rank, other
  than a species, is on the candidate node or below it. So the two versions agree on a novel genus
  above a named species: both create it. Its genomes keep the placeholder for that rank, and each of their
  lower missing ranks is decided on its own candidate node. See
  [refinement.md](refinement.md#novel-taxa).

### L8. A taxon is never moved up above a taxon of a higher rank

- **Python** (`Tree.refine_existing_labels`, `model/tree.py`): the upward walk only checks for
  conflicts of the same rank, so a family can be moved above its own order.
- **Rust** (`Tree::refine_existing_labels` / `Tree::highest_named_rank_below`, `src/tree.rs`): the
  walk stops at the first candidate ancestor that has a named taxon of a higher rank on a node
  strictly below it (the original node and leaves included). In r226 bac120 this stops
  `f__Sporomusaceae` being moved above `o__Sporomusales`. See
  [refinement.md](refinement.md#rule-u2-walk-the-ancestors), step 5.

### L9. A genus is not moved down away from new genomes that its species names

- **Python** (`Tree.refine_existing_labels`, `model/tree.py`): moving down follows only the genomes
  that the prior taxonomy puts in the taxon. A genome without a prior taxonomy (new in this
  release) counts as being in no taxon, so a genus can be moved down onto the part of its node
  that holds its old genomes, leaving behind a new genome whose species names the genus.
- **Rust** (`Tree::species_genera_of_new_genomes` / `Tree::child_rank_taxa`, `src/tree.rs`): when
  refining the genus, a new genome with a binomial species (`s__X epithet`, not
  `s__{unresolved}`) counts as being in `g__X`, if the original node of the genus is inside the
  genus RED window and, when `g__X` is not a Latin name (`Taxon::is_latin_name`), also supported
  (bootstrap at least `--min-bs`). If the new genome is under a different child from the genus'
  old genomes, the descent stops, so the genus is not split. In r237 bac120 (with the r232 prior):
  - `g__Noviluteimonas` (Latin, so the window is enough) stays on its node (bootstrap 91, RED
    0.958) with the new genome G057139745 (`s__Noviluteimonas sp057139745`), rather than moving to
    its 8 genome child (bootstrap 100) and leaving G057139745 to a novel genus.
  - `g__JAHELV01` (not Latin, on a supported node) stays on
    `100:f__JAHELV01;g__JAHELV01|red=0.8758` with all 9 genomes, as the curators would keep it,
    rather than moving to its 8 genome child (`100:red=0.9095`, which holds the 5 prior genomes)
    and leaving the new genome G057028815 (`s__JAHELV01 sp057028815`) to a novel `g__JBXLYT01`.

  See
  [refinement.md](refinement.md#rule-d2-follow-the-genomes-of-the-taxon), rule D2a.

### L10. Novel taxa are not created out of rank order

- **Python** (`fill_in_missing_ranks`, `util/tree.py`): the only checks on the highest candidate
  node are the RED, the bootstrap and that no taxon of the same rank is below it, so a novel taxon
  can be created inside a named taxon of a lower rank, or above one of a higher rank.
- **Rust** (`fill_in_missing_ranks`, `src/tree.rs`): the novel taxon is skipped if a named taxon of
  a lower rank is on a node above its candidate node, or a named taxon of a higher rank (a novel
  one made earlier in the run included) is on a node below it. As for L7, only that rank is
  skipped. With the r220 / r226 priors this skips 12 novel taxa in r226 bac120 and 13 in r232
  bac120 (none in ar53), such as 10 novel orders inside `g__MED-G20` (r232 bac120) and 7 novel
  orders inside families or `g__Karelsulcia` (r237 bac120). See
  [refinement.md](refinement.md#novel-taxa) and
  [testing.md](testing.md#job-16-novel-taxa-out-of-rank-order).

## Failures that are now reported errors

These are the same logical outcome (the run stops), but Python failed with an `assert`,
`TypeError` or `IndexError`, while Rust returns a named error with a message.

| Situation | Python | Rust |
| --- | --- | --- |
| An ancestor inside the RED window has no bootstrap | `TypeError` comparing `None >= float` | `Invariant("expected a node with bootstrap support")` |
| The best move is onto a leaf | `AssertionError` in `NodeLabel.to_string` | `Invariant("… cannot be moved to the leaf node …")` |
| A move would leave a node with neither a bootstrap nor a taxonomy | `AssertionError` in `NodeLabel.to_string` | `Invariant("… neither a bootstrap nor a taxonomy")` |
| A leaf missing a rank has no candidate node | `IndexError` on `candidates[-1]` | `NoCandidateNode { rank, gid }` |
| A leaf's parent has RED ≠ 1.0, or has a bootstrap | `assert` | `InvalidNodeLabel` |
| A label with more than one `:` or `\|RED=` | `ValueError` on unpacking | `InvalidNodeLabel` |
| A taxon to remove isn't on the node | `ValueError` from `list.remove` | `TaxonNotFound` |
| A taxon is on more than one node | `logger.error` + `sys.exit(1)` | `NonUniqueTaxon` (same outcome) |
| The end-of-run taxa check fails | `assert` | `Invariant(…)` |

The upward branch that moves a taxon further from the median to escape an unsupported node
`print`s `UNTESTED CODE` to stdout in Python. Rust is silent, as the branch is tested. Both record
the candidate.

## Input changes

### Median RED file

| | Python | Rust |
| --- | --- | --- |
| Format | JSON object, e.g. `{"phylum": 0.34, …}` (`*.decorated.dict`) | TSV with a header, using the `rank` and `median` columns (e.g. `gtdb_r237_median_red_rank_distribution.tsv`) |
| Rank keys | First letter of the key | Full rank name, any case; unknown names are errors |
| Domain | Forced to `0.0` | Fixed at `0.0`; a row with any other value is an error |
| Species | Taken from the file | Fixed at `1.0`; a row with any other value is an error |
| Duplicates | Last one wins | Error |

The reference logs in `test_data` were made with the JSON `.dict` file, so the Rust version needs
a TSV converted from it.

### Metadata / novel names

| | Python `--meta PATH` | Rust `--novel-name-file PATH --novel-name-column NAME` |
| --- | --- | --- |
| Name column | Hard-coded `ncbi_wgs_formatted` (optional) | Given by `--novel-name-column` (required with the file) |
| Prior taxonomy | `gtdb_taxonomy` column (optional) | Not read (see L2) |
| "No name" values | Exactly `none` | `""`, `none`, `n/a`, `null`; trimmed, any case |
| A short row | Lines with ≤1 column get no name; a missing later column is an `IndexError` | `MissingValue { line, column }` error |
| Trailing empty columns | Lost, because the whole line is `strip()`ped before splitting | Kept |
| Missing `accession` column | `KeyError` | `MissingColumn` error |

The genome ID canonicalisation (`gid_to_canonical`) is the same in both versions.

### Placeholders

A placeholder is a rank prefix with no name (`g__`).

| | Python | Rust |
| --- | --- | --- |
| In a tree label (`100:g__\|RED=0.9`) | Read as a taxon like any other, and refined like a name | `PlaceholderOnTree` error, naming the label |
| In a row of `--taxonomy` | (no such file) | Stripped as the row is read |
| In a row of `--prior-taxonomy` | (no such file) | `TaxonomyFile` error on that line: every rank of a prior row must be named |
| A rank with both a placeholder and a name (`g__;g__X`) | (no such file) | `TaxonomyFile` error on that line |

None of the trees in `example/` or `test_data/` has a placeholder on it, so this does not change
their output.

### Command line

- `--meta` was replaced by `--novel-name-file` and `--novel-name-column` (see above).
- The second positional argument is now the median RED TSV (`MEDIAN_RED`), not the JSON dict.
- Any input whose name ends in `.gz` is decompressed as it is read (`src/input.rs`); Python read every
  input as plain text. The outputs are not compressed, and are named after the tree without `.gz`.
- `--log`, `--min-bs` (default `95.0`) and the "output directory exists" check (skipped at debug
  level) are unchanged.

## Output changes

| | Python | Rust |
| --- | --- | --- |
| Tree | `<stem>.curatomatic.tree`, written by dendropy | Same name, written by `NewickTree::to_file` (`src/newick.rs`) |
| Report | Not written (the code is commented out) | `<stem>.curatomatic.report.tsv`: one row per moved or novel taxon, with flags (`src/report.rs`) |
| Leaf taxonomy | Not written: the taxa of a leaf are on its dummy node in the tree | `<stem>.curatomatic.leaf-taxonomy.tsv`: the taxa held by each leaf, which the tree can't hold (`Tree::leaf_own_taxonomies`) |
| Moves | `process_refinement` returns nothing | Returns each `Move`, which feeds the report |
| Log | — | Also logs how many tree genomes have a prior taxonomy, and warns when downward refinement is disabled |
| Rewritten labels | The bootstrap is rounded to `.1f` and the RED to `.3f` (`100` becomes `100.0`, `0.8894` becomes `0.889`), and the later passes read the rounded RED back | The bootstrap and RED are written exactly as read, so every pass (and the report) uses the full value, and the report's columns write them as read. Only the RED key is written `\|RED=` |

The Rust Newick writer quotes labels that contain `_` or protected characters, and turns spaces
into `_` in unquoted labels, to match dendropy with `preserve_underscores=True`. Branch lengths are
written back exactly as they were read. This document doesn't claim the two writers are
byte-for-byte identical.

## Unchanged logic

These rules match line for line. The Rust version is
documented in [refinement.md](refinement.md).

- **Rank order**: phylum → genus, skipping the domain and species. Reload after each rank.
- **Considered nodes**: nodes with a taxonomy *and* a bootstrap, in post-order.
- **RED window**: upper = midpoint(parent median, rank median), lower = midpoint(rank median,
  child median). Walking up stops when `RED < upper`, and a node is considered only if
  `RED ≤ lower`.
- **Upward bootstrap rule**: `bs ≥ min_bs` or `bs > original bs`.
- **Upward conflict rule**: stop if a descendant holds a named taxon of the rank, and not this
  taxon. "Named" means longer than 3 characters in both versions.
- **Upward record rule**: record if the RED delta is strictly smaller, or (untested branch) if the
  original node is unsupported.
- **Downward rules**: skip when supported and `RED > lower`; stop on polyphyly (more than one
  child containing the taxon in the prior taxonomy); the record/continue table for nodes without
  a bootstrap, supported nodes and improved nodes.
- **Choosing the move**: the candidates plus the original node, sorted by supported first,
  then RED delta ascending, then bootstrap descending. A missing bootstrap counts as `100`. The
  sort is stable, with ties in insertion order. Nothing happens if the original node is best.
- **Transferring a taxon**: remove it from the old node, add it to the new one, rewrite both
  labels. Python rounds the bootstrap and RED of both, Rust doesn't (see
  [Output changes](#output-changes)).
- **Leaves missing ranks**: walk up from each leaf, collecting taxa until every rank is
  assigned.
- **Novel rank placement**: stop at the first ancestor that has any taxon of the rank below it
  (placeholders included). Candidates need bootstrap ≥ `min_bs` (or none), and RED ≥
  `parent + (rank − parent) / 2`. The Rust code uses `RedValue::halfway_to` to keep Python's
  floating-point form. The highest candidate is used, leaves that share it share a taxon, and
  ranks are processed highest first. (Rust then skips it if a name is below it, see L7.)
- **Unique taxa**: a taxon on more than one node stops the run.
- **End-of-run check**: the taxa after, minus the novel taxa, must equal the taxa before.
- **Label removal**: everything from the first `|` is stripped from every label.

## Python code that was not ported

This code is on `main` but is never reached from `run()`. Most of it is after the `sys.exit(0)`
under `"""OLD CODE BELOW, IGNORE"""` in `__main__.py`:

- `Tree.find_taxa_for_refinement`: the earlier "move up only" refinement.
- `Tree.get_nested_ranks`: nested-rank detection (its idea survives as L1).
- `check_for_candidate_sharing`, `refine_nodes` (`util/tree.py`): applying the old refinement.
- `Tree.find_low_bootstrap_support_nodes_with_taxon_placement`: an unfinished draft that returns
  nothing.
- `create_novel_ranks`, `mrca_faster`, `get_tax_strings_for_leaf_nodes`, `add_taxon_to_string`,
  `create_node_label`: unused or stubs.
- In `refine_existing_labels`: `d_taxon_to_gids`, `singleton_taxa` and `new_gids` are computed
  but never used.
- `model/report.py`: never called.
