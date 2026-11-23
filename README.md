# curatomatic

Semi-automated curation of GTDB reference trees for new releases.

_This tool was developed for internal use, it is not supported for public use._

See [docs/refinement.md](docs/refinement.md) for the rules curatomatic follows when refining the taxa on the tree.

## Example

```shell
cargo run --release -- --min-bs 95.0 \
    --taxonomy example/gtdb_r237_decorated_taxonomy.tsv \
    --prior-taxonomy example/gtdb_r232_bac120_taxonomy.tsv \
    --novel-name-file example/gtdb_r237_bac120_metadata_slim.tsv \
    --novel-accession db_name --novel-name-column ncbi_wgs_formatted \
    example/gtdb_r237_median_red_marked.nwk example/gtdb_r237_median_red_rank_distribution.tsv /tmp/curatomatic
```

Run `curatomatic --help` for every option. `OUT_DIR` must not already exist (unless `--log debug` is given).

## Installation

### Pre-built binary

A static Linux binary (`x86_64-unknown-linux-musl`) is attached to every
[release](https://github.com/Ecogenomics/curatomatic/releases), with its changes in [CHANGELOG.md](CHANGELOG.md):

### Building from source

```shell
git clone https://github.com/Ecogenomics/curatomatic.git
cd curatomatic
cargo build --release   # ./target/release/curatomatic
```

### Testing

```shell
cargo nextest run
scripts/benchmark.sh /tmp/benchmark   # against Python and the curated releases, needs test_data/ (docs/testing.md)
```

## Input

All files are tab separated. Genome IDs are matched in their canonical form (`GB_GCA_003697165.1` is `G003697165`).
A taxon is a rank prefix (`d`, `p`, `c`, `o`, `f`, `g`, `s`), `__`, then its name; a bare prefix (`g__`) is a
placeholder. Anything malformed or inconsistent stops the run, naming the line (the full rules are in
[AGENTS.md](AGENTS.md), "Taxonomy").

**`TREE`**: the decorated tree from gDecorate (`decorated.red.tree`), with bootstrap, taxa and RED on internal nodes,
and leaves labelled with genome IDs. A placeholder on the tree is an error.

```
((G964471825:0.161,G964608985:0.161)'100:f__CAYYKD01;g__CAYYKD01|red=0.8386':0.194, ...
```

**`MEDIAN_RED`**: the median RED rank distribution from gDecorate. Only the `rank` and `median` columns are read
(Phylum to Genus; the domain is always 0 and the species 1).

```
rank	num_taxa	num_taxa_pass_filter	median	mean	...
Phylum	161	46	0.3095	0.3314	...
```

**`--taxonomy`** (required): the full taxonomy of each genome, from gDecorate (`decorate_classifications.tsv`). No
header. It holds the taxa the tree can't: each leaf gets every taxon of its row that is not on its path (always its
species). The tree and file must agree on every other taxon.

```
G964471825	d__Bacteria;p__Spirochaetota;c__GWE2-31-10;o__GWE2-31-10;f__CAYYKD01;g__CAYYKD01;s__CAYYKD01 sp964471825
G982720455	d__Bacteria;p__Spirochaetota;c__GWE2-31-10;o__GWE2-31-10;f__;g__;s__{unresolved} sp982720455
```

**`--prior-taxonomy`**: the taxonomy of the previous GTDB release, in the same format. Without it no taxon is moved
down. Every row must have seven named taxa, with a binomial species in the row's genus. Genomes not on the tree are
ignored.

```
RS_GCF_003697165.2	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
```

**`--novel-name-file`**: the GTDB metadata file, giving the name of novel taxa created for each genome. Columns are
found by header: `--novel-accession` (default `accession`) and `--novel-name-column` (e.g. `ncbi_wgs_formatted`). A
genome named `none`, `n/a`, `null` or empty gets a novel taxon named after its genome ID.

```
db_name	ncbi_wgs_formatted
G000008885	none
```

## Output

Written to `OUT_DIR`, where `<stem>` is the file name of `TREE` without its last extension:

| File | Contents |
| --- | --- |
| `<stem>.curatomatic.tree` | The curated tree. Only the internal taxon labels change; leaves keep their genome ID. |
| `<stem>.curatomatic.taxonomy.tsv` | The curated taxonomy of each genome, in the `--taxonomy` format, with all 7 ranks. A placeholder (`f__`) is a rank left for a curator to fill (see [docs/refinement.md](docs/refinement.md), "Novel taxa"). |
| `<stem>.curatomatic.leaf-taxonomy.tsv` | The taxa held by each leaf itself (never on its path, no placeholders), e.g. `G053840885	g__Neosphaerochaeta;s__Neosphaerochaeta fermentans`. The curated tree plus this file give every genome's full taxonomy, nested taxa included. |
| `<stem>.curatomatic.report.tsv` | One row per taxon moved or created, and per novel taxon that could not be named (below). |
| `curatomatic.log` | The log of the run. |
| `<stem>.debug.zombie.tree`, `<stem>.debug.taxonomy.tsv` | With `--debug` only: the inputs as loaded, before curation, as an old style tree with dummy nodes and the full taxonomy of each genome (see [docs/testing.md](docs/testing.md), "Round trip"). |

Report columns (several values are joined by `,`):

| Column | Contents |
| --- | --- |
| `taxon`, `rank` | The taxon and its rank prefix. An unresolvable novel taxon is `<rank>__UNRESOLVABLE`. |
| `action` | `moved_up`, `moved_down`, `novel`, or `unresolvable` (not created: every name it could have is taken). |
| `flag` | To review: `nested_in_same_rank`, `novel_low_red`, `sister_join_past_child_thresh`, `unsupported_label`, `large_gain`, `displaced_from_prior`, `new_in_decoration`, `namesake_outside_taxon`. |
| `n_genomes_*`, `bootstrap_*`, `red_*` | The taxon's node `_before` (input tree, empty if novel) and `_after` (curated tree, or where an unresolvable taxon would go). The bootstrap and RED are written exactly as read from the tree. |
| `prior_taxon` | The taxa of this rank that a novel taxon's genomes had in `--prior-taxonomy`. |
| `nested_in`, `sister_taxa` | The closest taxon of the same rank it is nested in, and those joining it past the child threshold. |
| `namesakes` | For an unresolvable taxon, the taxa named after its genomes that took each name. |
