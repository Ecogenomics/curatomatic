# TODO



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
- **`s__{unresolved}` names are intentional: never touch them.** 5,066 species in the file are named like
  `s__{unresolved} sp053939565` (e.g. G053939565). They are real names, not placeholders. Assign them to the leaf
  like any other species, and write them to the output file exactly as read, `{unresolved}` included. Never rename,
  refine, replace or drop them.
- A leaf with a lowest rank of `s__` (an old tree with dummy leaves) gets nothing.
- Gaps above the lowest rank are never filled from the file (see below).
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
   changed, placeholders are skipped, an old style tree gets nothing, a gap above the lowest rank is not filled, an `s__{unresolved} sp…` name reaches the leaf and the output file unchanged, a missing leaf is an error, an extra genome is a
   warning, and leaf labels round trip with a taxonomy and with both RED spellings. For the output file: ranks come from the leaf and its
   ancestors in rank order, the nearest of two nested taxa wins, zombies are skipped, and the format matches the
   input file.

### Gaps above the lowest rank
The file fills only the ranks below the lowest one in the tree, never a gap higher up. In r237 bac120, 5,129 leaves
have a class but no phylum on their path (e.g. G003554345: its parent is
`c__T1SED10-198M;o__T1SED10-198M;f__T1SED10-198M`, then only unlabelled nodes up to `d__Bacteria`). That leaf gets
`g__T1SED10-198M` and `s__T1SED10-198M sp003554345` from the file, but not `p__Sumerlaeota`: the missing phylum is
left for `fill_in_missing_ranks`, which creates a novel taxon for it as now.

### Done when
`curatomatic` runs on `test_data/r237` with `--taxonomy test_data/r237/final_taxonomy.tsv`. Before refinement, every
leaf has every named rank from the file below the lowest rank in the tree, including its species. After the full run, every leaf has one taxon at each of the seven ranks on its path,
and `<stem>.curatomatic.taxonomy.tsv` has one seven rank row per genome that matches the final tree.
`cargo test` passes, and none of the tree's taxa, bootstrap values or structure were changed.
