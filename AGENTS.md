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
- `TODO.md` = specifications for upcoming work.

## Relative evolutionary divergence (RED)

RED is a measure of the relative evolutionary distance of a node in a phylogenetic tree.
The root of the tree has a RED value of 0, and the leaves have a RED value of 1.
Median RED values are given for Phylum, Class, Order, Family, Genus. These define the bands for where those ranks
should be found in the tree.

## Taxonomy

The taxonomy follows a seven rank taxonomic structure: Domain, Phylum, Class, Order, Family, Genus, Species.

Some taxa within the tree may have no name assigned to them, they are given in the placeholder format (prefix only, i.e.
`p__`, `c__`, `o__`, `f__`, `g__`).
These can be refined by curatomatic at will.

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
  Placeholders (e.g. `g__`) are not names, and are exempt. This is checked at the end of every run, before any output
  is written, by `check_taxa_unchanged` (`src/main.rs`), which fails the run with `CuratomaticError::TaxaChanged`.
- Tests live in a separate file, never inline: declare `#[cfg(test)] mod tests;` in `foo.rs` and put the tests in
  `foo/tests.rs` (see `src/taxonomy.rs`).
- When possible, use newtypes.
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
