# TODO


## Job 1: Write a review table for curators

### Problem
Some curator decisions can't be automated from RED and bootstrap, but curatomatic can point at them. The log is the
only record of what curatomatic changed, and it isn't machine readable.

### What to do
Write `<stem>.curatomatic.report.tsv` next to the output tree, one row per taxon that was changed or flagged:
`taxon, rank, action, flag, n_genomes_before, n_genomes_after, bootstrap_before, bootstrap_after, red_before,
red_after`. Actions: `moved_up`, `moved_down`, `absorbed` (Job 5), `novel`. Flags, all computed from data already in
memory:
- `nested_in_same_rank`: the taxon sits inside another taxon of the same rank (Job 2 skips these). Curators merged 8
  of the 13 found in bac120. See Job 9 for the pairs that are already nested in the input.
- `novel_low_red`: a novel taxon with two or more genomes whose node has `q < 0.5`. Curators split 38% of these
  genera in R226 bac120 (50/130) and 57% in R232 bac120 (130/227), against 2% at `q >= 0.5` in both.
- `sister_join_past_child_thresh`: two taxa of rank `r` whose common ancestor has RED above `child_thresh(r)`.
  For genera with a supported join, curators merged one into the other in 27% (R226 bac120) and 16% (R232 bac120)
  of cases, against about 1% for joins below `q = 0.5`. For phylum to order in R232 bac120 it was 56% (23/41).
- `unsupported_label`: the label stays on a node with bootstrap below `--min-bs`.
- `large_gain`: a move that multiplies the number of genomes by 5 or more. Curators kept 0 of 7 (R226 bac120) and 1
  of 12 (R232 bac120).

Flags added after the R232 review:
- `split_candidate`: an unmoved taxon on a supported node with `q < 0.5`, whose members from the previous release
  are all under one child while the other children hold only new genomes. Curators split 36–48% of these genera in
  R232 bac120 (Job 6, step 2).
- `displaced_from_prior`: a novel taxon whose members already had a taxon of that rank in the previous release
  (Job 8). Include the previous name in an extra column.
- `new_in_decoration`: a taxon that is on the input tree but not in any `gtdb_taxonomy` value in the metadata
  (names added before curatomatic runs, mostly Latin names from NCBI). Only 33–38% of these were kept unchanged in
  bac120 (42 of 129 in R232, 44 of 117 in R226), against 92–93% for names carried over from the previous release.

The report must not change the tree. Python and Rust must write identical files.


## Job 2: Name fallback novel taxa `<rank>__<accession>`, not `<rank>__NOVEL-<accession>`

### Problem
When no unused WGS prefix is available, `generate_unique_rank` names the taxon `g__NOVEL-G963674675`. Curators use the
bare accession: `g__G963674675`. In bac120, 298 novel taxa got a `NOVEL-` name; for 187 of them the curated name is
exactly the name with `NOVEL-` removed, on the same genomes (the other 111 were merged into a neighbour or given a
Latin name). ar53 has only 4 `NOVEL-` names, none kept.

### R232 check
Still valid. R232 bac120 has 365 `NOVEL-` names: 206 were curated to the bare accession on the same genomes, 72 were
merged into a neighbour and 87 got another name. ar53: 2 of 9. Many of the "other name" cases are WGS-style names
(`g__NOVEL-G020072925` became `g__JAIQGA01`), which means the genome has a WGS prefix that the metadata file didn't
supply. 13,146 of the 189,801 R232 bac120 genomes have an empty `ncbi_wgs_formatted`.

### Where
- Rust: `src/taxon.rs` → `generate_unique_rank` (last line).
- Python: `curatomatic/util/taxon.py` (the same function).

### What to do
1. Emit `{prefix}__{gid}`. If that name is already in `known_taxa`, fall back to the next lowest genome ID in the
   group, and only then to a suffixed name. Keep the result deterministic.
2. Update `src/taxon/tests.rs`.

### Done when
No `NOVEL-` appears in any output tree, and the bac120 exact count rises by about 187.
