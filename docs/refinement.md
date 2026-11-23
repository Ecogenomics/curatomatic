# Refinement of existing ranks

This document lists the rules that curatomatic follows when it *refines* the taxa already on the
tree, i.e. moves an existing label to a better node. It covers `Tree::refine_existing_labels`
(finding candidate nodes) and `Tree::process_refinement` (choosing and applying a move), both in
`src/tree.rs`, and the loop that drives them in `curate` (`src/main.rs`).

Creating novel taxa for genomes that are missing ranks (`fill_in_missing_ranks`) runs after
refinement and is not covered here.

Refinement only ever moves a label between nodes. It never changes the structure of the tree, its
bootstrap values or its RED values, and it never creates, renames or deletes a taxon.

## Terms

| Term | Meaning |
| --- | --- |
| RED | Relative evolutionary divergence of a node: `0.0` at the root, `1.0` at the leaves. |
| Median RED | The median RED of a rank, read from the `MEDIAN_RED` file. The domain is fixed at `0.0` and the species at `1.0`. |
| RED delta | `|median RED of the rank − RED of the node|`. Smaller is better. |
| `min_bs` | The `--min-bs` argument (default `95.0`). A node is *supported* if its bootstrap is `≥ min_bs`. |
| Original node | The node the taxon is on when the rank's pass starts. |
| Named taxon | A taxon with a name after its prefix (`g__Escherichia`). A bare prefix (`g__`) is a placeholder. Placeholders are banned on the tree and stripped from the input files (see [Input preparation](#input-preparation)), so every taxon on the tree is named. |
| Prior taxonomy | The taxonomy of each genome in the `--prior-taxonomy` file (a previous release). Only used for moving labels down (and by the report). |

## Order of processing

1. The ranks are refined one at a time, from highest to lowest: **phylum, class, order, family,
   genus**. The domain and species are never refined.
2. Within a rank, every candidate move is found first, on the tree as it was at the start of the
   pass. The moves are then applied one taxon at a time, in the order the taxa were found
   (post-order of their original nodes). Candidates are not recalculated between moves in the same
   pass.
3. After each rank the node labels are reloaded, so the next (lower) rank sees the moves just made.

## RED window of a rank

Each rank has a window of RED values that are acceptable when moving a label **up**:

- **Upper bound** (towards the root): the midpoint of the parent rank's median and this rank's
  median.
- **Lower bound** (towards the leaves): the midpoint of this rank's median and the child rank's
  median.

A node is inside the window when `upper bound ≤ RED ≤ lower bound`.

Using the r226 bac120 medians as an example:

| Rank | Median | Window |
| --- | --- | --- |
| phylum | 0.341 | 0.171 – 0.401 |
| class | 0.461 | 0.401 – 0.539 |
| order | 0.617 | 0.539 – 0.689 |
| family | 0.760 | 0.689 – 0.840 |
| genus | 0.920 | 0.840 – 0.960 |

## Which taxa are considered

For the rank being refined, every taxon of that rank on every node is considered, **as long as the
node has both a taxonomy and a bootstrap value**. Nodes without a bootstrap (leaves, and in old
style trees the species nodes directly above them) are skipped: a taxon on a leaf is never moved
up, but a taxon can be moved down onto a leaf.

Each taxon is searched in two directions, up and down. Candidates from both are pooled and the
best is chosen (see [Choosing the move](#choosing-the-move)).

## Moving up

Moving up merges the original node with its sister lineages under a higher ancestor.

### Rule U1: nested taxa do not move up

If any ancestor of the original node holds a different **named** taxon of the same rank, the taxon
is nested inside it and is not moved up at all. (Moving up would climb onto or past the other
taxon's node and take its genomes.) It can still move down.

### Rule U2: walk the ancestors

Otherwise each ancestor is visited in turn, from the parent of the original node towards the root.
For each ancestor:

1. **Stop** if its RED is above the window (`RED < upper bound`). No higher ancestor is visited.
2. **Skip** (but keep walking) if its RED is below the window (`RED > lower bound`).
3. **Skip** (but keep walking) unless the bootstrap would be acceptable: the ancestor must be
   supported (`bs ≥ min_bs`) **or** have a higher bootstrap than the original node.
   An ancestor without a bootstrap is a malformed tree and is an error.
4. **Stop** if moving here would cause a conflict: any node below the ancestor holds a named taxon
   of this rank, and none of the named taxa of this rank on that node is the taxon being moved.
5. **Stop** if a named taxon of a **higher** rank is on any node strictly below the ancestor,
   leaves and the original node included: the taxon would be placed above it (e.g. a family above
   its own order). A higher rank on the ancestor itself does not count. Every higher ancestor has
   the same taxon below it, so the walk ends here.
6. Otherwise, the ancestor is **recorded as a candidate** if either:
   - its RED delta is strictly smaller than the original node's (it brings the taxon closer to the
     median), or
   - the original node is unsupported (`bs < min_bs`). Here the move may take the taxon further
     from the median, as the trade is for support (tested by
     `an_unsupported_taxon_is_moved_up_to_a_supported_node_further_from_the_median`).

   If neither holds, the ancestor is not recorded, and the walk continues.

Note that the checks in steps 4 and 5 are only reached for ancestors that pass steps 1–3. An
ancestor skipped in step 2 or 3 is not checked, but every higher ancestor contains the same
descendants, so a conflict below it stops the walk at the next ancestor that is checked.

## Moving down

Moving down narrows the taxon to the part of the tree that its genomes actually belong to,
according to the prior taxonomy.

> Moving down needs a prior taxonomy, which is read from `--prior-taxonomy` (a previous release).
> Without it the log warns "No prior taxonomy was given, downward refinement is disabled." and
> nothing is moved down. It is never read from `--taxonomy`: that is the decoration of the same
> tree, so every genome below a taxon has it, more than one child always contains it (Rule D2),
> and the descent would stop at once.

### Rule D1: well supported, deep taxa stay

A taxon is not moved down if the original node is supported **and** its RED is already below the
window (`RED > lower bound`): going lower would only move further from the median.

### Rule D2: follow the genomes of the taxon

For a node, each child is assigned the set of taxa of this rank that the prior taxonomy gives to
the genomes at or below it (a child may be a leaf itself). Genomes without a prior taxonomy
contribute nothing, with one exception for the genus:

- **Rule D2a: new genomes keep the genus their species names.** When refining the genus,
  a genome without a prior taxonomy (new in this release) whose species (from `--taxonomy`) is a
  binomial name, e.g. `s__Noviluteimonas sp057139745`, is taken to be in the genus that the first
  part of the name gives (`g__Noviluteimonas`, see `Taxon::genus_of_species`). So a genus is not
  moved down away from a new genome that its species puts in the genus: the child holding it also
  contains the genus, and the descent stops as below. This only applies to a genus whose
  original node:
  - is inside the genus [RED window](#red-window-of-a-rank), for every genus. A genus decorated
    at a node too deep or too shallow for a genus is moved down as before; and
  - for a genus that is not a Latin name (`Taxon::is_latin_name`), a placeholder style name such
    as `g__JAHELV01` or `g__CAG-698`, is also supported (bootstrap at least `--min-bs`). A Latin
    name (`g__Noviluteimonas`, `g__Luteimonas_D`) needs the window only.

  | Genus name | Original node | Rule D2a |
  | --- | --- | --- |
  | Latin | in window | applies |
  | Latin | outside window | no |
  | not Latin | in window, supported | applies |
  | not Latin | in window, unsupported | no |
  | not Latin | outside window | no |

  Unresolved species (`s__{unresolved} NAME`) give no genus. A genome that is in the prior
  taxonomy is never given the genus of its species, even if the prior has a different genus (or a
  placeholder) for it. Tested by `a_genus_is_not_moved_down_away_from_a_new_genome_that_its_species_names`
  and the tests after it (`src/tree/tests.rs`).

  The decoration names each new species after the genus it decorated, so the species of a new
  genome usually echoes the decoration rather than being independent evidence. Applied to every
  genus, the rule raised the curator workload of r226 and r232 by 36% (8,391 to 11,394
  operations); limited to Latin names inside the window, by 2% (to 8,547), in exchange for keeping
  genera such as `g__Noviluteimonas` in r237 whole (`91:g__Noviluteimonas|red=0.9579`, an
  unsupported node, which is why the support condition is not applied to Latin names). Extending
  it to non-Latin names on a supported node raised the workload by a further 3% (8,543 to 8,827,
  `docs/testing.md`), in exchange for keeping genera such as `g__JAHELV01` in r237 whole: on
  `100:f__JAHELV01;g__JAHELV01|red=0.8758`, it was moved down to the child holding its prior
  genomes, leaving the new genome G057028815 (`s__JAHELV01 sp057028815`) to a novel genus. The
  curators keep such a genus on its node.

- If **more than one** child contains the taxon, moving down would make it polyphyletic, so the
  descent stops here.
- If **no** child contains it, the descent stops here.
- Otherwise the descent continues into the one child that contains it.

This is first applied to the original node, then to each node the descent reaches.

### Rule D3: what is recorded on the way down

At each node reached (the RED window is **not** applied on the way down):

| Node | Recorded as a candidate? | Keep descending? |
| --- | --- | --- |
| No bootstrap (a leaf, or the species node above one in an old style tree) | Yes | No |
| Supported (`bs ≥ min_bs`) | Yes | Only if its RED delta is not larger than the original node's |
| Unsupported, but higher bootstrap than the original node | Yes | Yes |
| Anything else | No | Yes |

"Keep descending" means Rule D2 is applied again from that node.

## Choosing the move

All candidates for a taxon (from both directions) are ranked together with the original node. For
a candidate without a bootstrap, the bootstrap is taken to be `100`. They are sorted by:

1. **Supported first**: `bs ≥ min_bs` before `bs < min_bs`.
2. **Smallest RED delta** to the rank's median.
3. **Highest bootstrap**.

Ties keep the order the nodes were found (the first candidate found comes before the original
node).

- If the best is the original node, nothing is moved.
- Otherwise the taxon is removed from the original node and added to the best node. The move is
  logged (`Moving <taxon> from <label> to <label>`) and recorded for the report.

So the original node is kept whenever it is at least as good: a supported original node is only
beaten by a supported candidate closer to the median (or equally close with a higher bootstrap).

## Moves that are rejected

A chosen move fails with an error, stopping the run, if removing the taxon would leave an internal
original node with **neither a bootstrap nor a taxonomy**, so it would have nothing to be labelled
with. A leaf can hold a taxon (on its `Leaf`, not its label), so its last taxon can always be
moved away.

Neither node is changed in that case.

## Invariants

These are checked during or after the refinement, and break the run if violated:

- Every taxon is on **exactly one** node (checked whenever the labels are loaded), with no
  exception: an input tree with `g__UNRESOLVABLE` on two nodes is a `NonUniqueTaxon` error like any
  other. A novel taxon on more than one node fails `check_taxa_unchanged`.
- Assigning the taxonomy file to the leaves only adds taxa, it never removes one, and afterwards
  each genome's taxa on the tree are exactly those of its row (`TaxonomyFileConflict`).
- Every leaf holds its named species and has a named domain, after the leaf assignment and again
  before any output is written (`LeafWithoutSpecies`, `LeafWithoutDomain`).
- Curatomatic only moves taxa or creates novel taxa, it never removes a taxon or edits its text
  (`check_taxa_unchanged`). After the novel taxa are created, and before anything is written, the
  named taxa on the tree are compared with those after the leaf assignment, counting the nodes
  that hold each one. Every taxon must still be on as many nodes as before (moved, not removed or
  copied), and with the same text. Every other taxon must be a novel taxon that does not reuse an
  existing name. A violation fails the run with an
  error that lists each taxon that was removed, copied, reused or gained.
- The label of a leaf is never changed (`NewickTree::set_internal_label` refuses a leaf): the
  output tree labels each leaf exactly as it was read. The taxa held by a leaf are kept on its
  `Leaf`, and are written to the taxonomy file and the leaf taxonomy file only
  (`<stem>.curatomatic.leaf-taxonomy.tsv`, `Tree::leaf_own_taxonomies`). The RED is removed from every internal node
  (`Tree::remove_red_labels`).

## Input preparation

Before refinement, when the tree is read:

- A leaf is labelled `<gid>[|RED=<red>]`. A leaf without a RED has a RED of `1.0`, and any other
  leaf RED is invalid. A leaf label that holds a taxonomy (e.g. `<gid>|g__A`) is an error
  (`TaxonomyOnLeaf`): the taxa of a leaf come only from `--taxonomy`. An internal node is labelled
  `[<bootstrap>][:][<taxonomy>]|RED=<red>`. The RED key may be written `RED=` or `red=`, after a
  `|`, or after the `:` of a bootstrap with no taxonomy (`100:red=0.889`).
- Every taxon, in a tree label, `--taxonomy`, `--prior-taxonomy` or a novel name, is a rank
  prefix, then `__`, then its name (`g__Escherichia`, `s__{unresolved} sp1`). The name must not
  start with `_` or whitespace, nor end with whitespace, and is otherwise kept exactly as read. A
  taxon with an unknown prefix (`x__Foo`, `G__Foo`) is `UnknownRank`; any other malformed taxon
  (`gFoo`, `g_Foo`, `g___Foo`, `g__ Foo`) is `InvalidTaxon`. Both stop the run. This is checked by
  `Rank::split_taxon` (`src/rank.rs`), which every taxon passes through on its way into a
  `Taxonomy`. A bare prefix (`g__`) is a placeholder: `Taxon::is_named` is false for it.
- A placeholder is never held by a `Taxonomy` (`Taxonomy::add_taxon` refuses one with
  `Placeholder`). A placeholder in a node label of the tree stops the run (`PlaceholderOnTree`,
  naming the label): `100:g__|RED=0.9`, `100:f__F;g__|RED=0.9` and `d__|RED=0.0` are errors. In
  a row of `--taxonomy` or `--prior-taxonomy`, a rank that holds both a placeholder and a name
  (`g__;g__X`) is a `TaxonomyFile` error on that line; otherwise the placeholders are stripped,
  so the row holds named taxa only. A row whose species is a placeholder has no species, and is
  caught by the leaf checks below.
- If `--taxonomy` is given, every leaf is given the taxa the tree can't hold, before anything else
  runs (`Tree::assign_leaf_taxonomy`). The file holds what the dummy nodes of the old trees held:
  each row is the genome's full taxonomy, the taxa on its path plus its own.
  - The nearest node above the leaf with a named taxon is found, and its taxa must all be in the
    row. If not, the tree and the file disagree, and the run fails (`TaxonomyFileConflict`).
  - The leaf gets every named taxon of its row that is **not on its path**, matched by taxon and
    never by rank. A single-genome genus nested in a genus on the tree, or a single-genome order
    nested in a family on the tree, is the leaf's. The row is in rank order, so the position of a
    taxon in it says nothing. If every taxon is on the path (an old style tree, with the dummy
    nodes still on it), the leaf gets nothing.
  - The taxa already on the tree are never changed.
  - `s__{unresolved} <name>` species are names, and are assigned as read. There are no
    placeholders to assign: they were stripped from the rows.
  - A taxon to assign that is already elsewhere on the tree, or that more than one leaf would get,
    is an error (`TaxonomyFileConflict`). Nothing is skipped silently.
  - A leaf whose genome is not in the file is an error. A genome in the file that is not on the
    tree is a warning.
  - Once every leaf is assigned, the named taxa on each leaf and its path must **equal** the taxa
    of its row, matched by taxon as sets (never by rank, in any order). A taxon on the path above
    the nearest one that is not in the row (the tree puts the genome in `p__A`, its row says
    `p__B`), or a name on the tree where the row has a placeholder, is a `TaxonomyFileConflict`
    naming the genome and the taxa. This makes the round trip of the old trees (removing the
    dummy nodes and assigning the file gives back exactly the taxa they held) a check on every
    run.
- Every leaf must then hold its own named species (`Tree::check_leaves_have_species`,
  `LeafWithoutSpecies`), and have a named domain on itself or a node above it
  (`Tree::check_leaves_have_domain`, `LeafWithoutDomain`). Both are checked again before any
  output is written. `--prior-taxonomy` rows are not checked for a domain.

## Novel taxa

Novel taxa are placed by `fill_in_missing_ranks` after refinement. The leaf itself is a candidate
node, as is every node above it. A novel taxon without a WGS name is named after the lowest genome
ID under it whose name is not already a taxon (`g__G963674675`).

If every name is taken (every WGS name and every genome ID of its genomes is already a taxon, e.g.
the input tree named a taxon after a genome that is not in it), the novel taxon is
**unresolvable**, and is **never written**:

- Nothing is placed on its candidate node. Its genomes keep the placeholder of that rank (e.g.
  `g__`) in `<stem>.curatomatic.taxonomy.tsv`, exactly as for a skipped taxon, and nothing is added
  to the curated tree or the leaf taxonomy file.
- Only that rank is affected. The candidates and the names below each node are worked out before
  any novel taxon is created, so every other novel taxon is placed as before, including the lower
  ranks of its genomes.
- It is a row of the report: taxon `<rank>__UNRESOLVABLE`, action `unresolvable`, flag
  `namesake_outside_taxon`, the `*_after` columns of the candidate node it would have been placed
  on, and the `namesakes` column listing the taxa that refused each name (e.g. `g__G020350785`),
  which also says which genomes it is. A curator names it. The run logs the number of unresolvable
  taxa for each rank, beside the number skipped.
- `<rank>__UNRESOLVABLE` is only a label of the report, never a taxon on the tree.

A novel taxon is only created where **no valid name is below it**. A curator never places a novel
rank above existing names:

- For a missing rank `R` on its highest candidate node `N`, the taxon is created only if no node at
  or below `N` (leaves included) holds a valid name of a rank below `R`. Species do not count:
  every genome has a named species, so they would stop every novel genus. A lower rank on `N`
  itself counts (e.g. `c__Chlorobiia` on the node of a would-be phylum).
- A valid name is any taxon that is not a placeholder (`c__`, `g__`, ...) and not an unresolved
  species (`s__{unresolved} NAME`), see `Taxon::is_valid_name`. Named species, species named after
  a placeholder genus (`s__CAYYKD01 sp964471825`) and WGS-style taxa (`g__CAYYKD01`) are all valid.
- The names are read from the tree after refinement, before any novel taxon is created.
- The genomes of a skipped taxon keep the placeholder (e.g. `p__`) for that rank only. Each lower
  missing rank is decided on its own highest candidate node, by the same check. No report row is
  added; the run logs the number skipped for each rank.
