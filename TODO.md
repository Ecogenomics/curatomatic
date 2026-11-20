# TODO

## Job 10: Never move a taxon up above a taxon of a higher rank

**Status: done** (`docs/differences.md` L8).

### Problem

The upward search in `Tree::refine_existing_labels` (`src/tree.rs`) only checks for conflicts of the **same** rank
(`Tree::would_cause_conflicts`). Nothing stops a taxon from being moved onto a node that has a taxon of a **higher**
rank below it, which places e.g. a family above its own order.

Found when comparing the Python and Rust outputs on `test_data/r226_bac120` (`scripts/compare_outputs.py`):

```
G000445445/19  bs=71  RED=0.749                                   <- Rust moves f__Sporomusaceae here
└─ G000445445/18  bs=54  RED=0.789  o__Sporomusales; f__Sporomusaceae
```

Rust logs `Moving f__Sporomusaceae from 54.0:o__Sporomusales; f__Sporomusaceae|RED=0.789394 to 71.0|RED=0.748746`
(report flag `unsupported_label`). The output tree has `f__Sporomusaceae` on `G000445445/19` and `o__Sporomusales`
on its child `G000445445/18`: the family is above its order.

The Python version has the same gap. It did not make this move only because a singleton `f__JAJDHX01` on a dummy
node below happened to be a same-rank conflict.

### Rule

When moving a taxon of rank `R` **up**, the walk stops at the first candidate node `C` that has a named taxon (not a
placeholder) of a rank **higher** than `R` on any node strictly below `C`, leaves included.

- Taxa on `C` itself do not count: a family on the same node as its order is fine.
- The walk **stops** (it does not skip `C` and carry on), because every node above `C` also has that taxon below it.
- The original node counts as "below `C`": in the example, `o__Sporomusales` on `G000445445/18` stops the move of
  `f__Sporomusaceae` to `G000445445/19`.
- Unchanged: the same-rank conflict check, the nested-taxon rule (L1), the RED window, the bootstrap rules, and the
  downward search.

### Tasks

1. Add the check to the upward walk in `Tree::refine_existing_labels`, next to `would_cause_conflicts`.
2. Tests (in `src/tree/tests.rs`):
   - A family on a node with its order is not moved up to the parent (the Sporomusaceae case).
   - A genus is not moved up to a node whose other child holds a family.
   - A higher rank on the candidate node itself does not stop the move.
   - The walk stops at the first such node, it does not continue to a higher candidate.
3. Document the rule in `docs/refinement.md` (moving up), and add an entry `L8` to `docs/differences.md`: Python
   allows the move, Rust does not.
4. Re-run `test_data/r226_bac120` and `scripts/compare_outputs.py`, and check that `f__Sporomusaceae` is no longer
   moved and that no other taxon is moved above a taxon of a higher rank.

## Measuring Jobs 11 to 14 against the curated releases

Jobs 11 to 13 came from testing against the curated releases (`docs/testing.md`). Each Python behaviour was switched
on, one at a time, in a scratch copy of the Rust code, and every run was scored against the curated truth. The
expected numbers below come from those runs, on `dev` at `0672bda` (before Job 10). Job 10 may shift them slightly.

**Caveat:** the curators started from the Python output, so the truth favours Python's behaviour. A score gain is
not proof that the current rule is wrong (see "Notes for agents" in `AGENTS.md`). Jobs 11 to 13 each need the user's
approval before they are started.

To measure a job, build with `cargo build --release`, then run each dataset (`stem` is
`gtdb_<dataset>.decorated.scaled.red` for r226 and `gtdb_<dataset>.decorated.red_decorated` for r232):

```shell
d=r232_bac120; stem=gtdb_$d.decorated.red_decorated; out=<scratch dir>/$d
./target/release/curatomatic --min-bs 95.0 \
    --taxonomy test_data/$d/rs_input/gtdb_$d.decorated_taxonomy.tsv \
    --novel-name-file test_data/$d/rs_input/gtdb_${d}_metadata_slim.tsv --novel-name-column ncbi_wgs_formatted \
    test_data/$d/rs_input/$stem.nwk test_data/$d/rs_input/gtdb_$d.median_red_rank_distribution.tsv $out
python3 scripts/compare_truth.py --name $d --truth <truth file, see docs/testing.md> \
    --baseline test_data/$d/rs_input/gtdb_$d.decorated_taxonomy.tsv --rs-out $out --out <scratch dir>/cmp/$d
```

"Matches" is the sum of `correct` over ranks p, c, o, f and g for version `rust` in `scores.tsv`. On `0672bda` they
are:

| | r226_ar53 | r226_bac120 | r232_ar53 | r232_bac120 | Total |
| --- | ---: | ---: | ---: | ---: | ---: |
| Genome-ranks (genomes x 5) | 34,840 | 683,230 | 50,610 | 949,005 | 1,717,685 |
| Input taxonomy (no curatomatic) | 31,527 | 623,982 | 48,000 | 890,093 | 1,593,602 |
| Python (stored output) | 31,825 | 624,791 | 49,041 | 898,390 | 1,604,047 |
| Rust at `0672bda` | 31,813 | 628,903 | 48,770 | 897,912 | 1,607,398 |

Not part of these jobs: the taxonomy file never overrides the tree's labels (see `AGENTS.md`). That rule already
holds, and is tested by `the_taxa_on_the_tree_are_never_changed_by_the_taxonomy_file` (`src/tree/tests.rs`). Python
kept a single-genome taxon from the file inside a tree taxon of the same rank. Rust correctly doesn't, so don't add
it back.

## Job 11: Species do not block a novel genus

**Status: done** (approved by the user on 2026-11-20). This reverses a rule of Job 9 ("valid name of a rank below
`R`, species included", and its test "a missing genus above a named species is not created").

### Problem

`fill_in_missing_ranks` (`src/tree.rs`) only creates a novel taxon where no valid name of a lower rank is at or below
its candidate node (Job 9, `Tree::lowest_valid_rank_below`). Every genome has a named species, so a missing genus is
almost never created, and those genomes reach the curator with `g__`. The rule was added for novel phyla above named
classes, orders and genera (the R237 `c__Chlorobiia` example). This job doesn't affect that case.

### Rule

In the Job 9 check, species do not count as a valid name below. Every other rank still counts, so no novel phylum,
class, order or family is created above a named lower rank. A novel genus is created where no genus (named or
WGS-style) is at or below its candidate node.

### Tasks

1. In `Tree::lowest_valid_rank_below`, ignore species (or add a separate "lowest valid rank above species" and use it
   in `fill_in_missing_ranks`). Keep `Taxon::is_valid_name` as it is.
2. Tests in `src/tree/tests.rs`:
   - Change the Job 9 test so that a missing genus above a named species **is** created.
   - A missing family above a named genus is still not created.
   - A missing phylum above a named class is still not created (R237 case).
3. Docs: `docs/refinement.md` "Novel taxa" (drop "species included"), `docs/differences.md` L7 (Python and Rust
   now agree at genus).
4. Measure (see above).

### Expected effect

| | r226_ar53 | r226_bac120 | r232_ar53 | r232_bac120 | Total |
| --- | ---: | ---: | ---: | ---: | ---: |
| Matches | 31,828 | 629,321 | 48,810 | 898,604 | 1,608,563 (+1,165) |

The new genera use the `g__G...` and WGS names the curators also use. With this job, 4 novel taxa in r226_bac120 and
1 in r232_bac120 become `UNRESOLVABLE` (see Job 14).

## Job 12: Keep filling lower ranks after a skipped rank

**Status: done** (approved by the user on 2026-11-20). This reverses a rule of Job 9 ("lower ranks are not filled
either", and its test "when a phylum is skipped, a missing family of the same genomes is not filled").

### Problem

When `fill_in_missing_ranks` skips a novel taxon (Job 9), it adds its leaves to `skipped_leaves`, and none of their
lower missing ranks are filled. For example, a skipped family leaves the genus empty, even when nothing named is below
the genus's candidate node.

### Rule

A skipped rank only skips that rank. Each lower missing rank is decided on its own, by the same Job 9 check on its own
candidate node.

### Tasks

1. Remove `skipped_leaves` from `fill_in_missing_ranks`, keeping the per-rank count of skipped taxa in the log.
2. Tests: change the Job 9 test so that the family **is** filled when nothing named is below its candidate node, and
   add one where the lower rank is still skipped because a valid name is below its own candidate.
3. Docs: `docs/refinement.md` "Novel taxa", `docs/differences.md` L7.
4. Measure (see above).

### Expected effect

| | r226_ar53 | r226_bac120 | r232_ar53 | r232_bac120 | Total |
| --- | ---: | ---: | ---: | ---: | ---: |
| Matches, this job alone | 31,813 | 628,953 | 48,771 | 898,170 | 1,607,707 (+309) |
| Matches, with Job 11 | 31,828 | 629,376 | 48,811 | 898,902 | 1,608,917 (+1,519) |

Counting only the genome-ranks Rust changed from the input, Jobs 11 and 12 together get 87.2% of the changes right.
Rust at `0672bda` gets 89.4%, and Python 45.2%.

## Job 13: A separate prior taxonomy for moving taxa down (`--prior-taxonomy`)

**Status: done** (approved by the user on 2026-11-20). Moving down moves labels that the decorated tree placed. Both versions are
designed to do this, but the tree's labels are authoritative (`AGENTS.md`), so the user must confirm that they want
it.

### Problem

Moving a taxon down (`docs/refinement.md` "Moving down") needs a prior taxonomy: the previous release. Rust reads it
from `--taxonomy` (`read_metadata`, `src/main.rs`). That is the decorated taxonomy already on the tree, so a move down
never happens. Python read the `gtdb_taxonomy` column of `--meta`, which is the previous release.

`docs/differences.md` L2 says that downward refinement didn't run in the Python reference outputs. That is true for
r226 only. Python ignored the r226 `metadata.tsv` because its column is named `taxonomy`. The r232 run was given
`curatomatic_metadata.tsv`, which has `gtdb_taxonomy`, so Python did move taxa down in r232.

### Rule

- New option `--prior-taxonomy <PATH>`, in the format of `--taxonomy` (read with `TaxonomyFile`, IDs made canonical).
- It is the only source of the prior taxonomy. `--taxonomy` no longer is: update its help text in `src/cli.rs`.
- A genome that isn't on the tree is ignored. A genome on the tree that isn't in the file has no prior. A placeholder
  rank (`p__`) has no prior at that rank. Log how many tree genomes have a prior (`taxonomy_coverage_messages`).
- Without `--prior-taxonomy`, downward refinement is disabled with the existing warning. The output must stay the same
  as at `0672bda` on all four datasets, because a move down never happens there today.
- The report's prior taxa (`ReportInput::prior_taxa`, the `displaced_from_prior` and `new_in_decoration` flags) come
  from `--prior-taxonomy` too. Today they come from `--taxonomy`, so `new_in_decoration` can never fire.

### Making prior files from the test data

Both domains of a release share one metadata file. Genomes that aren't on the tree are ignored.

- r232: `test_data/r232_*/curatomatic_metadata.tsv`, columns `accession` (1) and `gtdb_taxonomy` (3). It has 731,926
  named rows and 320,439 placeholder-only rows (genomes new in r232). Checked: for every genome named in both, it is
  identical to the curated r226 taxonomy, so it is the previous release.
  `awk -F'\t' 'NR>1{print $1"\t"$3}' test_data/r232_bac120/curatomatic_metadata.tsv > <scratch>/r232_prior.tsv`
- r226: `test_data/r226_*/metadata.tsv`, columns `accession` (1) and `taxonomy` (2). It has 596,626 named rows and
  258,061 `none`. It is probably r220 (not checked). `none` must be dropped:
  `awk -F'\t' 'NR>1 && $2!="none"{print $1"\t"$2}' test_data/r226_bac120/metadata.tsv > <scratch>/r226_prior.tsv`

### Tasks

1. Add `--prior-taxonomy` to `src/cli.rs`, and read it in `read_metadata` instead of `--taxonomy`.
2. Tests: a taxon is moved down with `--prior-taxonomy`, and isn't without it. A genome missing from the prior file
   has no prior.
3. Docs: `README.md` (inputs), `docs/refinement.md` "Moving down", `docs/differences.md` L2 (fix the r232 statement
   above).
4. Check that the output without `--prior-taxonomy` is unchanged on all four datasets.
5. Measure with `--prior-taxonomy` on all four datasets. r226 has never been measured with a prior, so report it
   separately.

### Expected effect (r232 only)

| | r232_ar53 | r232_bac120 |
| --- | ---: | ---: |
| Matches, this job alone | 48,790 (+20) | 897,919 (+7) |
| Matches, with Jobs 11 and 12 | 48,971 (+201) | 900,716 (+2,804) |

On its own, this job leaves the genomes that a taxon moved away from with an empty rank: 2,437 genome-ranks over
both. Job 12 fills them, mostly correctly. So this job should only ship together with Job 12. Jobs 11 to 13 together
give 1,610,891 matches in total, against 1,607,398 at `0672bda` and 1,604,047 for Python, but r226 includes no prior.

## Job 14: Report a namesake genome outside its taxon

**Status: approved.**

### Problem

`generate_unique_rank` (`src/taxon.rs`) names a novel taxon without a WGS name after its lowest genome ID whose name
isn't already a taxon. If every genome's name is taken, the taxon is `<rank>__UNRESOLVABLE`. This happens when the
input tree already has a taxon named after one of the genomes, but that genome isn't inside it. For example, with
Job 11, `G020350785` in r226_bac120 gets `g__UNRESOLVABLE`, because `g__G020350785` is already on the tree elsewhere
(the curators put the genome into `g__G020350785`). Refusing the name is right, but the report doesn't say why.

### Rule

When a novel taxon is `UNRESOLVABLE`, its report row gets a new flag, `namesake_outside_taxon`, and names the
existing taxa that are named after its genomes (e.g. `g__G020350785`). Nothing else changes.

### Tasks

1. Make `generate_unique_rank` (or its caller in `fill_in_missing_ranks`) return the taxa that blocked each name, and
   pass them to the report (`src/report.rs`, new `Flag` variant with a doc comment).
2. Tests in `src/taxon/tests.rs` and `src/report/tests.rs`.
3. Docs: `docs/refinement.md` "Novel taxa", the report columns in `README.md`.
4. With Job 11, check the r226_bac120 cases (`G020350785`, `G030406045`) and the r232_bac120 case (`G020350445`).
