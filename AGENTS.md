# AGENTS.md

## About

Curatomatic is a command-line tool that performs pre-curation of a domain-specific prokaryotic phylogenetic tree.
It is designed to automate the manual steps that curators of the tree perform, with the goal of reducing the amount of
manual work required.

The tree uses information based on Relative Evolutionary Divergence (RED) and a pre-decoration of the tree using a
different software package.

Curatomatic is written in Rust. It replaces an earlier Python implementation, which is no longer maintained.

## Documentation

- `README.md` = how to build and run curatomatic, and the format of its input files.
- `docs/refinement.md` = the rules curatomatic follows when refining the taxa already on the tree.
- `docs/differences.md` = the differences in behaviour between the Python and Rust implementations.
- `docs/testing.md` = the test datasets, and how to benchmark against the Python version and the curated releases
  (`scripts/benchmark.sh`).
- `TODO.md` = specifications for upcoming work.
- `docs/architecture.md` = diagrams of the pipeline, data structures and rules, with an interactive version in
  `docs/architecture.html`. Both are a once-off snapshot, not kept up to date: update them only when asked, following
  the "Regenerating" section of `docs/architecture.md`.

## Relative evolutionary divergence (RED)

RED is a measure of the relative evolutionary distance of a node in a phylogenetic tree.
The root of the tree has a RED value of 0, and the leaves have a RED value of 1.
Median RED values are given for Phylum, Class, Order, Family, Genus. These define the bands for where those ranks
should be found in the tree.

## Taxonomy

The taxonomy follows a seven rank taxonomic structure: Domain, Phylum, Class, Order, Family, Genus, Species.

A taxon is always `<prefix>__<name>`: one of the rank prefixes `d`, `p`, `c`, `o`, `f`, `g`, `s`, then `__`, then its
name (`g__Escherichia`, `s__Escherichia coli`, `g__CAG-698`, `s__{unresolved} sp1`). The name must not start with `_`
or whitespace, and must not end with whitespace. Anything else stops the run (`CuratomaticError::InvalidTaxon`, or
`UnknownRank` for an unknown prefix such as `x__Foo`). This is checked by `Rank::split_taxon` (`src/rank.rs`), which
every taxon passes through on its way into a `Taxonomy`, from the tree, the taxonomy files and novel naming.

A placeholder is a rank prefix with no name (`p__`, `c__`, `o__`, `f__`, `g__`). It is not a name, and is never held
in a `Taxonomy`: `Taxonomy::add_taxon` (`src/taxonomy.rs`) refuses one with `CuratomaticError::Placeholder`.

- **Tree:** a placeholder in any node label stops the run (`CuratomaticError::PlaceholderOnTree`, naming the label).
  `d__Bacteria|RED=0.000` is fine; `100:g__|RED=0.9`, `100:f__F;g__|RED=0.9` and `d__|RED=0.0` are not. A rank
  with no named taxon on a genome's path is simply missing, and can be filled with a novel taxon.
- **Taxonomy files** (`--taxonomy` and `--prior-taxonomy`, read by `TaxonomyFile::from_reader`): a rank of a row that
  holds both a placeholder and a name (`g__;g__X`) is a `CuratomaticError::TaxonomyFile` error on that line. Otherwise
  the placeholders are stripped, so each row holds named taxa only.
- The output taxonomy files still write a placeholder for a rank with no named taxon (see TODO Job 22).

The domain and species never have a placeholder: every genome always has a named domain (`d__`), and its species
(`s__`) is always given in the taxonomy file.

## Tree format

Previously, trees were given with "Dummy" nodes. These were used so that taxon label(s) could be placed for
leaf nodes or singletons in the tree, where no bifurcation was present to place those labels from.

To supplement the tree file, an additional taxonomy file is provided to store the taxon label(s) for those nodes.

The tree will be annotated with the RED values for each internal node.

## Development

After implementing a task, the following steps must be performed:

```shell
cargo fmt
cargo clippy --all-targets
cargo nextest run
```

`Cargo.toml` sets the lints to `deny`, so any clippy warning fails the check.

After finishing a task, commit the changes using semantic commit messages. Never push.

### Notes for agents

- Development happens on several machines and accounts. Anything an agent needs to remember about this project goes in
  this file (or the docs it links to) and is committed. Never keep project notes in an agent's own memory outside the
  repository.
- Before recommending a change in behaviour, check it against the intent in this file and `docs/refinement.md`, not
  only against scores. The curated test truth was made from the Python output, so it favours Python's behaviour. A
  score gain is not evidence that current behaviour is a bug: present such ideas as questions.

## Test data

Test data are available here:

- `example/gtdb_r237_decorated_taxonomy.tsv` = taxonomy file
- `example/gtdb_r237_median_red_marked.nwk` = tree file
- `example/gtdb_r237_bac120_metadata_slim.tsv` = novel rank creation candidate file
- `example/gtdb_r237_median_red_rank_distribution.tsv` = red distribution per-rank file

## Rules

- Curatomatic must never change the structure of the tree, its bootstrap values, or RED values: the only changes to
  the final tree are the taxon labels at each node (which form the taxonomy).
- Curatomatic can only ever move taxa or create novel taxa. It never removes a taxon, and never edits the text of one.
  This is checked at the end of every run, before any output
  is written, by `check_taxa_unchanged` (`src/main.rs`), which fails the run with `CuratomaticError::TaxaChanged`.
- Tests live in a separate file, never inline: declare `#[cfg(test)] mod tests;` in `foo.rs` and put the tests in
  `foo/tests.rs` (see `src/taxonomy.rs`).
- When possible, use newtypes.
- The taxon labels on the input (decorated) tree are placed there for a specific statistical reason and are
  authoritative. Never propose overriding the tree's labels with the file.
- The inputs (the tree and the taxonomy file) are gospel until curation starts. Loading them must never drop, change or
  second-guess a taxon: anything that does not line up is an error (`CuratomaticError::TaxonomyFileConflict`), never
  a silent skip.
- The taxonomy file (`--taxonomy`) replaces the dummy nodes of the old trees. Each row is the full taxonomy of a genome:
  the taxa on its path in the tree, plus the taxa that were on its dummy node. A leaf gets every named taxon of its row
  that is not on its path, matched **by taxon, never by rank** (`Tree::assign_leaf_taxonomy`). The row is in rank
  order, not path order, so a taxon's position or rank in the row says nothing about whether it is the leaf's: a
  single-genome genus inside a tree-labelled genus (r237 G053840885, `g__Neosphaerochaeta` inside `g__JAAYFO01`), or a
  single-genome order inside a tree-labelled family (r226 bac120 G017303975), are the leaf's own. Removing the dummy
  nodes and assigning the file must give back exactly the taxa they held. This is checked on every run: once the
  leaves are assigned, the named taxa on each genome's path and leaf must equal its row as sets
  (`Tree::check_leaves_match_rows`, `TaxonomyFileConflict`), and it is tested by
  `removing_the_dummy_nodes_and_assigning_the_file_restores_their_taxa` (`src/tree/tests.rs`).
  `scripts/convert_python_inputs.py` converts the Python inputs to this format. On real data, check it with
  `--debug` and `scripts/check_round_trip.py` (see the README): the r226 and r232 Python trees and the r237 inputs all
  round trip with 0 differences, and any change to how the inputs are loaded must keep it that way.
- Every leaf holds at least one taxon, and one of them is always its named species (never a placeholder, never on an
  internal node). This is checked by `Tree::check_leaves_have_species` (`src/tree.rs`) once the taxonomy file is
  assigned and again before any output is written, failing the run with `CuratomaticError::LeafWithoutSpecies`. This
  is why `--taxonomy` is required.
- Every leaf has a named domain, on itself or on a node above it. This is checked by `Tree::check_leaves_have_domain`
  (`src/tree.rs`) at the same two points as the species, failing the run with `CuratomaticError::LeafWithoutDomain`.
  `--prior-taxonomy` rows are not checked.
- Leaf node labels are immutable: a leaf is always labelled with its genome accession, exactly as it was read. The
  taxa that a leaf holds are stored in the `Taxonomy` of its `Leaf` (`src/tree.rs`) and are written to the taxonomy
  file, never to the tree. `NewickTree::set_internal_label` refuses to change a leaf.
- Species names can be given in the format `s__{unresolved} NAME`. These are not to be changed by curatomatic.
- An internal node within the tree can have multiple taxa in the same rank. These cases can still be refined by
  curatomatic, but deciding which of them should be kept must be done manually.
- If an external crate is required, then approval must be obtained from a human. Present the case, and include the
  number of downloads.
- Intent and rules must be enforced as lints or tests so that the program does not diverge from the intended
  specification. Anything that can not be represented as rules or lints, must be recorded here.
