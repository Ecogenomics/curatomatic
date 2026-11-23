# Architecture

How curatomatic runs, the data structures it uses, and the rules it applies, as diagrams. The same diagrams,
with interactive explorers, are in [`architecture.html`](architecture.html) (see [Regenerating](#regenerating)).

This is a snapshot of the code at commit `27f27d5` (2026-11-23). The **code is the source of truth**: if this
document and the code disagree, the code is what runs. The detailed rules are in [refinement.md](refinement.md).

## Pipeline

`curate` (`src/main.rs`) runs these steps in order. Any error stops the run before output is written, except the
output directory check, which runs before the log is opened.

```mermaid
flowchart TD
    A["Parse the command line (Cli)"] --> B{"OUT_DIR exists and --log is not debug?"}
    B -- yes --> X["Exit: output directory already exists"]
    B -- no --> C["Create OUT_DIR, start the log (curatomatic.log)"]
    C --> D["Tree::from_file: parse the Newick tree and every node label"]
    D --> E["read_taxonomy_file: every genome on the tree must be in --taxonomy"]
    E --> F["Tree::assign_leaf_taxonomy: leaves get the taxa the tree can't hold"]
    F --> G["Tree::check_leaves_have_species"]
    G --> H{"--debug?"}
    H -- yes --> H1["write_debug_inputs: zombie tree and full taxonomy"]
    H -- no --> I
    H1 --> I["Snapshot the taxa: Tree::taxon_counts"]
    I --> J["read_metadata: --prior-taxonomy and --novel-name-file"]
    J --> K["Red::from_file: the median RED of each rank"]
    K --> L["For phylum, class, order, family, genus: refine_existing_labels, process_refinement, reload"]
    L --> M["Tree::get_leaf_nodes_missing_rank"]
    M --> N["fill_in_missing_ranks: create novel taxa"]
    N --> O["check_taxa_unchanged: only moves and novel taxa"]
    O --> P["Tree::check_leaves_have_species (again)"]
    P --> Q["build_report"]
    Q --> R["Write the taxonomy TSV"]
    R --> S["Tree::remove_red_labels, write the tree"]
    S --> T["Write the report TSV"]
```

| Output | Written by |
| --- | --- |
| `<stem>.curatomatic.taxonomy.tsv` | `taxonomy_file::write` from `Tree::leaf_taxonomies` (nearest taxon of each rank on the path) |
| `<stem>.curatomatic.tree` | `Tree::to_file`, after `Tree::remove_red_labels` (leaf labels unchanged) |
| `<stem>.curatomatic.report.tsv` | `write_report` from `build_report` |
| `curatomatic.log` | `log::init` |
| `<stem>.debug.zombie.tree`, `<stem>.debug.taxonomy.tsv` | `write_debug_inputs` (`--debug` only) |

## Data structures

### The tree

```mermaid
classDiagram
    class Tree {
        newick: NewickTree
        leaves: map NodeId to Leaf
        node_labels: Vec of NodeLabel, by NodeId
        gids: set of GenomeId
        taxon_labels: set of Taxon
        reload()
    }
    class NewickTree {
        nodes: Vec of Node, post-order
        root: NodeId
        set_internal_label()
    }
    class Node {
        label: Option String
        length: Option String
        parent: Option NodeId
        children: Vec of NodeId
        first: NodeId
    }
    class Leaf {
        gid: GenomeId
        taxonomy: Taxonomy
    }
    class NodeLabel {
        red: RedValue
        bootstrap: Option BootstrapValue
        taxonomy: Option Taxonomy
        gid: Option GenomeId
    }
    class Taxonomy {
        d, p, c, o, f, g, s: Vec of Taxon
    }
    class Taxon {
        String, e.g. g__Escherichia
        rank()
        is_named()
        is_valid_name()
        is_latin_name()
    }
    class Rank {
        Domain to Species
        parent()
        child()
    }
    Tree *-- NewickTree
    Tree *-- Leaf
    Tree *-- NodeLabel
    NewickTree *-- Node
    Leaf *-- Taxonomy
    NodeLabel *-- Taxonomy
    Taxonomy *-- Taxon
    Taxon ..> Rank
```

- `NodeId` is the position of a node in post-order. Every node comes after its descendants, and the descendants of
  a node are the contiguous ids just before it (`Node.first` to the node).
- `node_labels` is derived from the Newick labels, plus the `Leaf` taxonomy for leaves. It is only up to date after
  `Tree::reload`.
- The taxa of an internal node live in its Newick label. The taxa of a leaf live on its `Leaf`, never in its label.
- `GenomeId`, `Taxon`, `RedValue` and `BootstrapValue` are newtypes over `String` and `f64`.

### Inputs, results and the report

```mermaid
classDiagram
    class TaxonomyFile {
        rows: map GenomeId to Taxonomy
    }
    class Metadata {
        rows: map GenomeId to MetadataRow
        has_prior_taxonomy: bool
    }
    class MetadataRow {
        accession: GenomeId
        taxonomy: Option Taxonomy, from --prior-taxonomy
        wgs: Option String, from --novel-name-file
    }
    class Red {
        data: map Rank to RedValue
        get()
    }
    class Refinement {
        entries: Vec of Taxon and its moves, to and from
    }
    class Move {
        taxon: Taxon
        from: NodeId
        to: NodeId
    }
    class TaxonCounts {
        map Taxon to number of nodes
    }
    class ReportRow {
        taxon, rank, action, flags
        before: Option NodeStats
        after: NodeStats
        prior_taxa, nested_in, sister_taxa
    }
    class Action {
        MovedUp
        MovedDown
        Novel
    }
    class Flag {
        NestedInSameRank
        NovelLowRed
        SisterJoinPastChildThresh
        UnsupportedLabel
        LargeGain
        DisplacedFromPrior
        NewInDecoration
    }
    Metadata *-- MetadataRow
    Refinement ..> Move : process_refinement makes
    ReportRow --> Action
    ReportRow --> Flag
```

### Node labels

| Node | Label | Example |
| --- | --- | --- |
| Leaf | `<gid>[\|RED=1.0]`. Anything else after the genome ID is `TaxonomyOnLeaf`. | `G964471825` |
| Internal | `[<bootstrap>][:][<taxonomy>]\|RED=<red>`, with a bootstrap, a taxonomy or both. The key may be `RED=` or `red=`, after `\|`, or after the `:` of a bootstrap with no taxonomy. | `100:f__CAYYKD01;g__CAYYKD01\|red=0.8386`, `100:red=0.8894` |

## RED windows

For a rank `R` with parent `P` and child `C` (the domain is fixed at 0.0 and the species at 1.0):

- **upper bound** = midpoint(median `P`, median `R`)
- **lower bound** = midpoint(median `R`, median `C`)
- a node is **inside the window** when `upper ≤ RED ≤ lower`
- **RED delta** = `|median R − RED of node|`
- a novel taxon of `R` may only go on a node with `RED ≥ halfway(median P, median R)`, which is the upper bound

## Assigning the taxonomy file to the leaves

`Tree::assign_leaf_taxonomy` (`src/tree.rs`) matches by taxon, never by rank.

```mermaid
flowchart TD
    L["Each leaf, genome G"] --> R{"G in --taxonomy?"}
    R -- no --> E1["Error: MissingFromTaxonomyFile"]
    R -- yes --> U["Walk up and collect the named taxa on the path. The first node with a named taxon is the anchor."]
    U --> A{"Every anchor taxon in the row of G?"}
    A -- no --> E2["Error: TaxonomyFileConflict"]
    A -- yes --> D["The leaf gets every named taxon of the row that is not on the path"]
    D --> X{"Taxon already on the tree, or wanted by another leaf?"}
    X -- yes --> E3["Error: TaxonomyFileConflict"]
    X -- no --> P["Add it to Leaf.taxonomy"]
```

## Refinement

### One pass per rank

```mermaid
flowchart TD
    S["Pass for rank R: phylum, class, order, family, then genus"] --> W["Window of R from the median RED"]
    W --> N["Every node in post-order with a taxonomy and a bootstrap"]
    N --> T["Every taxon of rank R on it"]
    T --> U["Up candidates"]
    T --> D["Down candidates"]
    U --> F["Refinement: the candidates of each taxon, in the order found"]
    D --> F
    F --> P["process_refinement: move each taxon to its best node"]
    P --> RL["Tree::reload, then the next rank"]
```

Candidates are found on the tree as it was at the start of the pass, and are not recalculated between moves.

### Moving up (`Tree::refine_existing_labels`)

```mermaid
flowchart TD
    S(["Taxon of rank R on its original node"]) --> U1{"U1: a different named taxon of R on an ancestor?"}
    U1 -- yes --> NO(["No up candidates"])
    U1 -- no --> NX["Go to the parent"]
    NX --> TOP{"Past the root?"}
    TOP -- yes --> DONE(["Done"])
    TOP -- no --> S1{"RED above the window?"}
    S1 -- yes --> DONE
    S1 -- no --> S2{"RED below the window?"}
    S2 -- yes --> NX
    S2 -- no --> S3{"Supported, or a higher bootstrap than the original? No bootstrap is an error."}
    S3 -- no --> NX
    S3 -- yes --> S4{"A node below holds a different named taxon of R?"}
    S4 -- yes --> DONE
    S4 -- no --> S5{"A named taxon of a higher rank strictly below?"}
    S5 -- yes --> DONE
    S5 -- no --> S6{"Closer to the median than the original, or the original is unsupported?"}
    S6 -- yes --> REC["Record a candidate"] --> NX
    S6 -- no --> NX
```

### Moving down (`Tree::refine_existing_labels`, `Tree::child_rank_taxa`)

```mermaid
flowchart TD
    S(["Taxon of rank R on its original node"]) --> D1{"D1: original supported and RED below the window?"}
    D1 -- yes --> NO(["No down candidates"])
    D1 -- no --> CT["D2: the taxa of R that the prior taxonomy gives the genomes under each child. Genus only: a new genome adds the Latin genus of its species, if the original node is inside the window."]
    CT --> K{"How many children contain the taxon?"}
    K -- "none, or more than one" --> STOP(["Stop"])
    K -- "exactly one" --> V["Visit that child"]
    V --> KIND{"D3: what is the child?"}
    KIND -- "a leaf" --> RL(["Record, stop"])
    KIND -- "supported" --> RS["Record"]
    RS --> FAR{"RED delta larger than the original?"}
    FAR -- yes --> STOP
    FAR -- no --> CT
    KIND -- "unsupported, higher bootstrap than the original" --> RU["Record"] --> CT
    KIND -- "anything else" --> CT
```

Without `--prior-taxonomy` no genome has a prior taxonomy, so no child contains the taxon and nothing moves down.

### Choosing the move (`Tree::process_refinement`)

```mermaid
flowchart LR
    C["The candidates in the order found, with the original node after the first"] --> SO["Stable sort: supported first, then smallest RED delta, then highest bootstrap (a leaf counts as 100)"]
    SO --> B{"Is the best the original node?"}
    B -- yes --> K(["No move"])
    B -- no --> M["transfer_taxon: remove from the original, add to the best"]
    M --> R(["Record the Move for the report"])
```

## Novel taxa (`fill_in_missing_ranks`)

```mermaid
flowchart TD
    M["Tree::get_leaf_nodes_missing_rank: leaves whose path lacks a rank"] --> E["Each leaf and missing rank R, from phylum down"]
    E --> W["Walk up from the leaf, the leaf included"]
    W --> A{"A node below holds a taxon of R?"}
    A -- yes --> H["Take the highest candidate. None is NoCandidateNode."]
    A -- no --> OK{"No bootstrap or supported, and RED at least halfway from the parent median to the R median?"}
    OK -- yes --> C["Candidate"] --> W
    OK -- no --> W
    H --> G["Leaves with the same highest node share one taxon"]
    G --> V{"A valid name of a rank lower than R, species excluded, at or below the node?"}
    V -- yes --> SK["Skip: the leaves keep the placeholder of R"]
    V -- no --> NM["generate_unique_rank: the most readable WGS name, else the lowest unused genome ID, else R__UNRESOLVABLE"]
    NM --> P["Add the novel taxon to the node"]
```

The candidates and the valid names below each node are read before any novel taxon is created. A missing domain
has no parent rank, and is a `NoParentRank` error.

## Checks

| Check | Where | When | Error |
| --- | --- | --- | --- |
| A leaf label is its genome ID, with a RED of 1.0 if any | `parse_leaf_label` (`src/node_label.rs`) | Reading the tree | `TaxonomyOnLeaf`, `InvalidNodeLabel` |
| Every taxon is on one node, except `<rank>__UNRESOLVABLE` | `Tree::get_taxon_labels` | Every load and `reload` | `NonUniqueTaxon` |
| Every genome on the tree is in `--taxonomy` | `read_taxonomy_file` (`src/main.rs`) | After reading the tree | `MissingFromTaxonomyFile` |
| The tree and the file agree, and no taxon is given twice | `Tree::assign_leaf_taxonomy` | Leaf assignment | `TaxonomyFileConflict` |
| Leaf assignment only adds taxa | `assign_leaf_taxonomy` (`src/main.rs`) | Leaf assignment | `Invariant` |
| Every leaf holds its named species | `Tree::check_leaves_have_species` | After assignment, and before output | `LeafWithoutSpecies` |
| An ancestor inside the window has a bootstrap | `required_bootstrap` | Moving up | `Invariant` |
| A move leaves the original node with a bootstrap or a taxon | `transfer_taxon` | Applying a move | `Invariant` |
| A leaf label is never changed | `NewickTree::set_internal_label`, `NodeLabel::to_label` | Any label write | `Invariant` |
| Taxa are only moved or created | `check_taxa_unchanged` (`src/main.rs`) | Before output | `TaxaChanged` |
| Every move and novel taxon has a report row | `check_every_change_has_a_row` (`src/report.rs`) | Building the report | `Invariant` |

## Code map

| File | Holds |
| --- | --- |
| `src/main.rs` | `curate`: the pipeline, input reading, `check_taxa_unchanged`, `--debug` output |
| `src/cli.rs` | `Cli`: the command line |
| `src/tree.rs` | `Tree`, `Leaf`, leaf assignment, refinement, novel taxa, taxonomy output |
| `src/newick.rs` | `NewickTree`, `NodeId`: parsing and writing Newick, post-order arena |
| `src/node_label.rs` | `NodeLabel`: reading and writing node labels |
| `src/taxon.rs` | `Taxon`, `TaxonCounts`, `generate_unique_rank` |
| `src/taxonomy.rs` | `Taxonomy`: the taxa of each rank on one node |
| `src/taxonomy_file.rs` | `TaxonomyFile`: reading and writing the taxonomy TSV |
| `src/metadata.rs` | `Metadata`: novel names and the prior taxonomy |
| `src/red.rs` | `Red`, `RedValue`: the median RED file |
| `src/rank.rs` | `Rank` |
| `src/bootstrap.rs`, `src/genome.rs` | `BootstrapValue`, `GenomeId`, `gid_to_canonical` |
| `src/report.rs` | `build_report`, `ReportRow`, `Action`, `Flag` |
| `src/error.rs` | `CuratomaticError` |
| `src/log.rs` | The logger |

## Regenerating

This document and `architecture.html` were made once, and are not kept up to date automatically. To update them
when asked:

1. Read the code, not only the docs: `curate` in `src/main.rs`, then `src/tree.rs`, `src/node_label.rs`,
   `src/taxon.rs`, `src/report.rs` and the other files in the [code map](#code-map).
2. Compare what the code does with `AGENTS.md`, `README.md`, `docs/refinement.md` and `docs/differences.md`. Where
   they disagree, the diagrams follow the code. List each disagreement for the user, who decides whether the docs or
   the code are wrong. Do not fix either without being asked.
3. Update the diagrams here, then copy them unchanged into the `<pre class="mermaid">` blocks of
   `architecture.html`.
4. `architecture.html` has three explorers written in JavaScript that copy rules from the code. Check each against
   the code and update it:
   - the label reader copies `split_red`, `parse_leaf_label` and `parse_node_label` (`src/node_label.rs`);
   - the RED window explorer copies the bounds in `Tree::refine_existing_labels` and `fill_in_missing_ranks`, and
     holds the medians of `example/gtdb_r237_median_red_rank_distribution.tsv` and of `docs/refinement.md`;
   - the upward walk explorer copies the upward walk of `Tree::refine_existing_labels` and the sort in
     `Tree::process_refinement`.
5. Update the commit and date at the top of this document and in the footer of `architecture.html`.
6. Publish `architecture.html` to its artifact link (below) with the Artifact tool, updating it in place.
