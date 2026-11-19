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
| Named taxon | A taxon with a name after its prefix (`g__Escherichia`). A bare prefix (`g__`) is a placeholder, and the rules below that say "named" ignore placeholders. |
| Prior taxonomy | The taxonomy of each genome in the `--taxonomy` file. Only used for moving labels down (and by the report). |

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
5. Otherwise, the ancestor is **recorded as a candidate** if either:
   - its RED delta is strictly smaller than the original node's (it brings the taxon closer to the
     median), or
   - the original node is unsupported (`bs < min_bs`). Here the move may take the taxon further
     from the median, as the trade is for support. This branch logs `UNTESTED CODE`.

   If neither holds, the ancestor is not recorded, and the walk continues.

Note that the conflict check in step 4 is only reached for ancestors that pass steps 1–3. An
ancestor skipped in step 2 or 3 is not checked, but every higher ancestor contains the same
descendants, so a conflict below it stops the walk at the next ancestor that is checked.

## Moving down

Moving down narrows the taxon to the part of the tree that its genomes actually belong to,
according to the prior taxonomy.

> Moving down needs a prior taxonomy, which is read from `--taxonomy`. Without it the log warns
> "No prior taxonomy was given, downward refinement is disabled." and nothing is moved down.

### Rule D1: well supported, deep taxa stay

A taxon is not moved down if the original node is supported **and** its RED is already below the
window (`RED > lower bound`): going lower would only move further from the median.

### Rule D2: follow the genomes of the taxon

For a node, each child is assigned the set of taxa of this rank that the prior taxonomy gives to
the genomes at or below it (a child may be a leaf itself). Genomes without a prior taxonomy
contribute nothing.

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
with. A leaf can hold a taxon, and keeps its genome ID when its last taxon is moved away.

Neither node is changed in that case.

## Invariants

These are checked during or after the refinement, and break the run if violated:

- Every taxon is on **exactly one** node (checked whenever the labels are loaded). The only
  exception is the sentinel `<rank>__UNRESOLVABLE` (see [Novel taxa](#novel-taxa)).
- Assigning the taxonomy file to the leaves only adds taxa, it never removes one.
- At the end of the run, the set of taxa on the tree equals the set after the leaf assignment plus
  the novel taxa that were created: refinement moves taxa, it never adds or loses them.
- The output tree labels each leaf with its genome ID only, and no node with a RED value
  (`Tree::set_output_labels`). The taxa held by the leaves are written to the taxonomy file.

## Input preparation

Before refinement, when the tree is read:

- A leaf is labelled `<gid>[|<taxonomy>][|RED=<red>]`. A leaf without a RED has a RED of `1.0`, and
  any other leaf RED is invalid. An internal node is labelled
  `[<bootstrap>][:][<taxonomy>]|RED=<red>`. The RED key may be written `RED=` or `red=`, after a
  `|`, or after the `:` of a bootstrap with no taxonomy (`100:red=0.889`).
- If `--taxonomy` is given, every leaf is given the ranks the tree can't hold,
  before anything else runs (`Tree::assign_leaf_taxonomy`):
  - The leaf gets every named taxon of its row whose rank is **below** the lowest rank of any taxon
    on its path to the root. If the lowest is a genus, it gets only its species. If it is a
    species (an old style tree, with a species node above each leaf), it gets nothing.
  - A gap above the lowest rank is never filled (e.g. a missing phylum above a class): it is left
    for the novel taxa.
  - The taxa already on the tree are never changed, even when the file disagrees.
  - Placeholders (`g__`) are never assigned. `s__{unresolved} <name>` species are names, and are
    assigned as read.
  - A taxon that is already on the tree, or that more than one leaf would get, is not assigned (it
    would be on more than one node). The run logs a warning with the number skipped.
  - A leaf whose genome is not in the file is an error. A genome in the file that is not on the
    tree is a warning.

## Novel taxa

Novel taxa are placed by `fill_in_missing_ranks` after refinement. The leaf itself is a candidate
node, as is every node above it. A novel taxon without a WGS name is named after the lowest genome
ID under it whose name is not already a taxon (`g__G963674675`). If every genome ID is taken (e.g.
the input tree named a taxon after a genome that is not in it), it is the sentinel
`<rank>__UNRESOLVABLE`, which may be on more than one node, for a curator to name.
