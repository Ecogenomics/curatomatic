# TODO

Jobs 1–2 are bugs and Jobs 3–11 (under "Algorithmic improvements") are rule changes, all found by comparing
curatomatic's pre-curation trees with the manually curated taxonomies of two releases:
- R226: `test_data/r226_ar53/ar53_taxonomy_r226_curated.tsv`, `test_data/r226_bac120/bac120_taxonomy_r226_curated.tsv`
- R232: `test_data/r232_ar53/ar53_taxonomy_r232.tsv`, `test_data/r232_bac120/bac120_taxonomy_r232.tsv`

Each job stands alone and can be picked up in a fresh context window. Job 1's code fix is done, but its
`test_data` update is still pending (see Job 1). Do that first: Job 2's numbers assume Job 1's data fix.

Jobs 1–7 were written from R226 only and were re-checked against R232 on 2026-11-11. Each has an "R232 check" note.
Summary of that review: Jobs 1, 2, 3 and 7 stand as written; Job 4 gains a step; Job 5 loses rule (a) and gets
per-rank thresholds; Job 6 keeps step 1 (now measured on real downward moves) and step 2 becomes a report flag.
Jobs 8–11 are new.

Two things to know about the R232 data:
- The R232 runs are the production runs (v0.0.2, default `--min-bs 95`, a metadata file that has `gtdb_taxonomy`), so
  they are the only data with real downward moves: 1,451 in bac120 and 146 in ar53.
- `test_data/r232_*` holds only the output tree, the log and the curated taxonomy. There is no input tree, RED
  dictionary or metadata file, so R232 can't be used for parity yet (Job 11). RED was recomputed as the distance from
  the root (this matches the RED labels in the R226 input trees to within 5e-7, and identifies the origin node of
  every R232 "Moving" log line uniquely). The R232 median RED per rank was estimated as the median over internal nodes
  where the rank is the most specific label on the node (within 0.004 of the real dictionary on R226):
  ar53 p 0.286, c 0.435, o 0.570, f 0.764, g 0.925; bac120 p 0.350, c 0.482, o 0.655, f 0.800, g 0.942.

Ground rules (from `AGENTS.md`): Python is the source of truth for parity, so fix **both** implementations
(`curatomatic/` Python and `src/` Rust) and keep their outputs identical. Tests go in `foo/tests.rs`, never inline.
Use newtypes where possible. Don't run clippy or other lints.

# https://claude.ai/artifact/VkStq8kuS59BV5N9wqXZ4s?sk=Z9PaGddtyfhmVQ2qFivhnA#0f9ea517-82d7 
---

## Job 1 (done): A metadata file without a `gtdb_taxonomy` column silently disables downward refinement

Code fixed in both implementations. After the metadata is subset to the tree, both log (identical wording):
- WARNING "The metadata file has no gtdb_taxonomy column, downward refinement is disabled." if the column is missing.
- INFO "N of M tree genomes have a gtdb_taxonomy."
- WARNING "No tree genomes have a gtdb_taxonomy, downward refinement is disabled." if the column exists but N is 0.

Rust tests: `src/metadata/tests.rs` (no column parses, every taxonomy is `None`) and `src/main/tests.rs` (warning
messages). Verified on a copy of the metadata with the header fixed: Python and Rust trees are byte-identical for all
six R226 runs, and logs match apart from lines Python prints to stdout (`UNTESTED CODE`, `Shared name`). `Moving`
counts at bs45 / bs80 / bs95: ar53 41 / 48 / 53, bac120 940 / 966 / 1,080. `UNTESTED CODE` fires 113 times on
bac120 bs95. No Rust parity fixes were needed.

**Still pending (data only, `test_data` is gitignored):** rename the second header column of both
`test_data/r226_*/metadata.tsv` files from `taxonomy` to `gtdb_taxonomy` and regenerate the six reference outputs
`test_data/r226_{ar53,bac120}/curatomatic_bs{45.0,80.0,95.0}/` with Python, using the command on line 1 of each
existing `curatomatic.log` (move the old directory aside first; the output directory must not exist). Until this is
done the checked-in reference outputs still don't exercise the downward code path.

---

## Job 2: Upward refinement lets a nested taxon climb onto or past an ancestor of the same rank

### Problem
When a taxon is moved up, the conflict check (`would_cause_conflicts` in Rust, the `desc_node` loop in Python)
only looks at nodes **below** the candidate node. It never checks the candidate node itself or the nodes above it.
A taxon nested inside another taxon of the same rank can therefore:
- **Climb onto the ancestor's own node.** That node then carries two labels of the same rank, e.g.
  `'100.0:g__Mycobacterium; g__Mycolicibacter'`, and every genome below has two genera.
- **Climb part of the way.** The genomes on the way up, which belonged to the ancestor taxon, get relabelled to the
  nested taxon.

The legacy Python code avoided this by skipping nested taxa (`get_nested_ranks` / `find_taxa_for_refinement`, now
dead code below `sys.exit(0)` in `curatomatic/__main__.py`). The current `refine_existing_labels` lost that guard.
The input decorated trees contain no same-rank double labels; curatomatic creates them.

### Evidence (test_data, `--min-bs 95`)
- **bac120: 13 moves took genomes from another taxon of the same rank, about 4,950 genome-rank assignments in
  total.** 8 of them left two same-rank labels on one node:

  | Moved taxon | Ancestor taxon | Genomes taken |
  |---|---|---|
  | c__Bacilli_A | c__Bacilli | 3,528 |
  | g__Mycolicibacter | g__Mycobacterium (same node) | 593 |
  | f__RUG14156 | f__Gastranaerophilaceae (same node) | 488 |
  | g__Bombiscardovia | g__Bifidobacterium (same node) | 170 |
  | f__CALMFT01 | f__Absconditicoccaceae (same node) | 79 |
  | g__Saccharothrix | g__Actinosynnema (same node) | 33 |
  | g__Nitrospira | g__UBA8639 (same node) | 29 |
  | g__CANLOV01 | g__Cellulophaga (same node) | 19 |
  | g__Thalassotalea | g__Thalassotalea_A | 9 |
  | g__CAMLQK01 | g__JARIHD01 | 2 |
  | g__RYN-133 | g__JALZAI01 | 2 |
  | g__Rosettibacter | g__JARYLX01 (same node) | 1 |
  | g__40CM-2-70-7 | g__JACOVF01 | 1 |

  Curators removed 8 of these 13 names, merging the nested taxon into its ancestor.
- **ar53: 1 move.** g__Methanorbis (nested inside g__Methanocorpusculum) climbed and took 9 Methanocorpusculum
  genomes. Curators merged Methanorbis into Methanocorpusculum.
- Find the bac120 cases in the log with
  `grep -E "Moving (g__Mycolicibacter|c__Bacilli_A|f__RUG14156) " test_data/r226_bac120/curatomatic_bs95.0/curatomatic.log`.

### R232 check
Still valid, and the bug is in the production output. The R232 bac120 tree has 17 nodes with two labels of the same
rank, all from upward moves (e.g. `g__UBA817` onto the node of `g__JAYNJH01`, `f__JAAFJU01` onto
`f__Leptolyngbyaceae`). Curators removed the climbing name in 12 of the 17, kept it at the new node in 3 (the other
name was dropped instead) and did something else in 2. ar53 has none. Two side effects to expect when this is fixed:
- 12 upward moves in R232 bac120 grew a taxon 5 times or more; 10 of them are these nested cases and only 1 of the
  12 was kept (R226 bac120: 0 of 7). No separate size guard is needed once nested taxa are skipped.
- Species under a climbing genus keep the old genus in their name. This accounts for most species whose name
  disagrees with their genus (see Job 4, step 4).

The invariant in step 2 must stay "same node". A same-rank label on an **ancestor** is normal input: the R232 bac120
tree has 354 such pairs that curatomatic did not create (Job 9).

### Where
- Python: `curatomatic/model/tree.py` → `Tree.refine_existing_labels`, the "Check up the tree" `while` loop.
- Rust: `src/tree.rs` → `Tree::refine_existing_labels` (the upward `while let Some(cur_node)` loop) and
  `Tree::would_cause_conflicts`.

### What to do
1. Before walking up for `taxon` (rank `r`), check whether any **ancestor** of the taxon's current node carries a
   named taxon of rank `r` (names longer than the bare `x__` prefix) that differs from `taxon`. If one does, the taxon
   is nested: skip upward refinement for it. Downward refinement can stay as it is. This also covers the case where
   the candidate is the ancestor's node.
2. Add a post-condition after each rank's `process_refinement` and after `fill_in_missing_ranks`: no node may carry
   two named taxa of the same rank. Make a violation a hard error (`CuratomaticError::Invariant` in Rust,
   `logger.error` plus `sys.exit(1)` in Python), so this class of bug can't pass silently again.
3. Optional (needs agreement, because it changes output): log the nested taxa that were skipped at INFO, e.g.
   "Skipping nested taxon g__Mycolicibacter (inside g__Mycobacterium)". Curators merged most of these, so the line is
   a useful review hint.
4. Rust unit tests: add a small newick fixture where `g__B` is nested inside `g__A`, and a parent node of `g__B`
   would otherwise qualify. Assert that `g__B` isn't moved and that no node ends up with two genus labels. Put
   tree-level tests in a new `src/tree/tests.rs`, declared with `#[cfg(test)] mod tests;` in `src/tree.rs`.
5. Apply the same change to Python, regenerate the six `test_data` reference outputs (see Job 1, "Still pending"), and
   confirm Rust parity.

### Done when
- The regenerated bac120 bs95 log has no `Moving` line for any of the 13 taxa above, and the ar53 bs95 log has none
  for g__Methanorbis.
- No node in any curatomatic output tree has two labels of the same rank. A quick check: `python3 -c` over the tree
  text, splitting each quoted label on `;` and counting rank prefixes.
- Rust and Python outputs are identical.

---

# Algorithmic improvements

Jobs 3–11 change what curatomatic outputs, on purpose. They come from the same comparison against the curated R226
and R232 taxonomies. Do them after Jobs 1 and 2, in both implementations, and regenerate the `test_data` reference
outputs after each one. Every threshold below is a starting value taken from R226 and R232; expose each as a CLI
option so it can be tuned.

Terms used below, for a rank `r` with median RED `med(r)`:
- `parent_thresh(r)` = midpoint of `med(parent rank)` and `med(r)`; `child_thresh(r)` = midpoint of `med(r)` and
  `med(child rank)`. These are the existing `parent_red_thresh` / `child_red_thresh` in `refine_existing_labels`.
- `q(node, r)` = `(RED(node) - parent_thresh(r)) / (med(r) - parent_thresh(r))`: 0 at the parent threshold, 1 at the
  median. Using `q` makes thresholds comparable between domains and releases, whose medians differ.
- `accept(r)` = the RED at `q = accept_q(r)`. Curators keep a well-supported node at or above `accept(r)` as rank `r`
  and tend to split one below it. Defaults: **`accept_q` = 0.625 for genus and 0.25 for family and above** (one CLI
  option per rank group). The earlier R226-only value was 0.5 for every rank; R232 shows the same curve for genus,
  with the crossover between 0.5 and 0.625, and a much lower bar for families (evidence in Jobs 5 and 6).
- "Supported" = bootstrap >= `--min-bs`.

How to measure any job: for every curated taxon (ranks p to g, representative genomes on the tree), count it as
exact if the curatomatic tree has the same name on the same set of genomes. Baseline at `--min-bs 95`:

| | curated taxa | exact (as checked in) | after Job 1 |
|---|---|---|---|
| R226 ar53 | 2,940 | 2,694 (91.6%) | 2,771 |
| R226 bac120 | 35,403 | 31,804 (89.8%) | 33,135 |
| R232 ar53 | 3,643 | 3,370 (92.5%) | n/a (already has downward moves) |
| R232 bac120 | 43,800 | 40,249 (91.9%) | n/a |

A job is worth keeping if that count goes up in both domains and both releases. R232 can only be re-scored once its
inputs are checked in (Job 11); until then, use the R232 numbers quoted in each job as the expected effect.

## Job 3: Name fallback novel taxa `<rank>__<accession>`, not `<rank>__NOVEL-<accession>`

### Problem
When no unused WGS prefix is available, `generate_unique_rank` names the taxon `g__NOVEL-G963674675`. Curators use the
bare accession: `g__G963674675`. In bac120, 298 novel taxa got a `NOVEL-` name; for 187 of them the curated name is
exactly the name with `NOVEL-` removed, on the same genomes (the other 111 were merged into a neighbour or given a
Latin name). ar53 has only 4 `NOVEL-` names, none kept.

### R232 check
Still valid. R232 bac120 has 365 `NOVEL-` names: 206 were curated to the bare accession on the same genomes, 72 were
merged into a neighbour and 87 got another name. ar53: 2 of 9. Many of the "other name" cases are WGS-style names
(`g__NOVEL-G020072925` became `g__JAIQGA01`), which means the genome has a WGS prefix that the metadata file didn't
supply. 13,146 of the 189,801 R232 bac120 genomes have an empty `ncbi_wgs_formatted` (Job 10).

### Where
- Rust: `src/taxon.rs` → `generate_unique_rank` (last line).
- Python: `curatomatic/util/taxon.py` (the same function).

### What to do
1. Emit `{prefix}__{gid}`. If that name is already in `known_taxa`, fall back to the next lowest genome ID in the
   group, and only then to a suffixed name. Keep the result deterministic.
2. Update `src/taxon/tests.rs`.

### Done when
No `NOVEL-` appears in any output tree, and the bac120 exact count rises by about 187.

## Job 4: Name unresolved species `s__<genus> sp<accession digits>`

### Problem
`Tree::from_file` rewrites an `s__` placeholder to `s__{unresolved} <gid>` and nothing renames it later. Curators name
these `s__<genus> sp<digits of the accession>` (e.g. G016187405 in g__JACPJS01 becomes `s__JACPJS01 sp016187405`):
285 of 287 in ar53 and 5,267 of 5,276 in bac120. Using the genus that curatomatic itself assigns would already match
the curated species for 262 (ar53) and 4,548 (bac120).

### R232 check
Still valid: the predicted name matches the curated species for 7,199 of 8,074 unresolved species in bac120 and 675
of 750 in ar53. The misses are genera that curators changed. Step 2 was too broad and is replaced by steps 2 and 4:
the "fix 334, break 424" result on R226 came almost entirely from genera that Job 2's bug had moved up. Split by how
the genus got there, renaming a placeholder species (`s__<old genus> sp<digits>`) to follow its genus is:

| Genus of the species is | R226 bac120 right / wrong | R232 bac120 right / wrong | R232 ar53 right / wrong |
|---|---|---|---|
| novel (created by curatomatic) | 308 / 11 | 2,201 / 351 | 166 / 37 |
| moved up | 13 / 406 | 14 / 298 | 1 / 0 |
| untouched | 13 / 7 | 35 / 1 | 1 / 0 |

("wrong" = curators kept the original species name; cases where neither name is right are left out.)

### Where
- Rust: `src/tree.rs` → `Tree::from_file` (the "hack" block) and a new final step in `src/main.rs` → `curate()`, after
  `fill_in_missing_ranks` and before `remove_red_labels`.
- Python: the same two places in `curatomatic/model/tree.py` and `curatomatic/__main__.py`.

### What to do
1. After novel ranks are created, for every leaf whose species label is `s__{unresolved} <gid>`, read the genus on the
   path to the root and set the label to `s__<genus name> sp<gid without the leading letter>`.
2. Don't touch species with a Latin epithet, and don't touch any species whose genus was moved by refinement.
3. The check in `curate()` that the set of taxa only grew by the novel taxa must allow for the renamed species.
4. New: for a species named `s__<X> sp<digits>` whose genus is a **novel** genus `g__<Y>` with `Y != X`, rename it to
   `s__<Y> sp<digits>`. Do this after Job 5, so that a genome absorbed into an existing genus is handled by that
   genus instead. Log each rename at DEBUG.

### Done when
No `{unresolved}` label is left in any output tree, and no placeholder species under a novel genus carries another
genus name.

## Job 5: Absorb an orphan lineage into its only sister taxon instead of creating a novel taxon

### Problem
A genome with no label for rank `r` gets a novel taxon from `fill_in_missing_ranks`. Curators often don't want a new
taxon: they move the neighbouring taxon's label up one node so that it takes the orphan in. Curators did this for 408
of 4,437 novel taxa in bac120 (261 genera, 84 families, 61 orders, 2 phyla) and 15 of 248 in ar53. Upward refinement
doesn't do it today, because the move takes the label further from the median RED (`red_delta_worse`).

### R232 check
The pattern holds and is larger in R232: curators merged 939 of 9,343 novel taxa in bac120 into another taxon (558
of them into a taxon below `P`) and 113 of 873 in ar53 (75 below `P`). Two corrections to the R226 write-up:
- **Rule (a) is withdrawn.** The R226 figures for `RED(P) > child_thresh` (25 of 31 genera, 43 of 54 families, 60 of
  97 orders) counted every novel name that curators removed. Counting only "merged into a taxon that sits below `P`"
  gives 11 of 31, 0 of 54 and 4 of 97 in R226, and 24 of 140, 1 of 46 and 6 of 141 in R232. Most of the rest are
  genomes that already had a taxon of that rank in the previous release, and curators put them back into it (often
  with a new suffix: in R226 the novel orders and families around the insect symbionts all went to
  `o__Enterobacterales_A` / `f__Enterobacteriaceae_A`). That is Job 8, not an absorb.
- **The bar for rule (b) depends on the rank.** See the table: use `accept_q` (0.625 for genus, 0.25 for family).

### Evidence (`--min-bs 95`)
Let `P` be the first ancestor of the novel taxon's node that has another taxon of rank `r` below it. Novel taxa with
`P` supported, exactly one other taxon of rank `r` below `P`, and `RED(P) <= child_thresh(r)`; share that curators
merged into that sister, by `q(P, r)`:

| `q(P)` | R226 bac120 genus | R232 bac120 genus | R232 ar53 genus | R226 bac120 family | R232 bac120 family |
|---|---|---|---|---|---|
| below 0 | 3/610 | 10/994 | 4/111 | 5/48 | 12/114 |
| 0–0.25 | 15/150 (10%) | 20/288 (7%) | 2/18 | 0/5 | 14/31 (45%) |
| 0.25–0.5 | 33/150 (22%) | 60/246 (24%) | 13/22 | 4/8 | 13/15 (87%) |
| 0.5–0.625 | 30/64 (47%) | 60/140 (43%) | 11/12 | 7/7 | 8/9 |
| 0.625–0.75 | 33/50 (66%) | 87/139 (63%) | 8/12 | 2/2 | 2/2 |
| 0.75–1 | 32/40 (80%) | 74/92 (80%) | 3/5 | 1/1 | 2/2 |

- `P` not supported, or two or more other taxa below `P`: rarely absorbed (under 6% in every dataset), leave these alone.
- Orders, classes and phyla have too few cases (R232 bac120 orders: 7 of 9 at `q >= 0`); use the family value.

### Where
- Rust: `src/tree.rs` → `fill_in_missing_ranks` (before a novel taxon is created for a candidate node), reusing
  `transfer_taxon`. Python: `curatomatic/model/tree.py`, the matching function.

### What to do
1. For each rank, before creating a novel taxon at candidate node `c`: find `P` as above. Let `X` be the taxa of rank
   `r` below `P`.
2. Absorb if `X` holds exactly one taxon, `P` is supported, and `accept(r) <= RED(P) <= child_thresh(r)`. To absorb,
   move the label of that taxon to `P` and create no novel taxon. Log it:
   "Absorbing <n> genomes into <taxon> (moved to <P label>)". Don't absorb when `RED(P) > child_thresh(r)`.
3. Otherwise create the novel taxon as now.
4. Process ranks from highest to lowest as now, and re-check after each absorb, because several orphans can share `P`.
5. Unit tests in `src/tree/tests.rs`: one fixture each for absorb, and for no absorb (two sister taxa; unsupported
   `P`; RED below `accept`; RED above `child_thresh`).

### Done when
The exact count rises in both domains. Expected at the default `accept_q`, right / wrong: genus 68 / 26 (R226
bac120), 163 / 70 (R232 bac120), 11 / 6 (R232 ar53); family 14 / 4 (R226 bac120), 26 / 3 (R232 bac120). Each correct
absorb fixes two curated taxa (the novel name disappears and the sister gets the right genomes). Sweep `accept_q`
for genus over 0.5–0.75 before fixing the default.

## Job 6: Use the acceptance band when moving a label down

### Problem
Downward refinement (active after Job 1) ranks candidate nodes by distance to the median RED. Two patterns disagree
with curators (bac120, genus, original node supported):
- **Moved down, curators moved it back.** Original RED 0.88–0.90: 25 kept, 44 reverted. 0.90–0.92: 6 kept, 21
  reverted. Below 0.86: 114 kept, 5 reverted. Curators accept a supported node at or above `accept(r)` and let the
  taxon take in the new basal genomes.
- **Not moved, curators moved it down.** The taxon sits on a supported node below `accept(r)`, the members from the
  previous release are all under one child, and the other children hold only new genomes. Curators split 55% of these
  at RED 0.84–0.86 and 33% at 0.86–0.88, against 3% at 0.90–0.92. Curatomatic keeps the label, usually because the
  only candidate is a leaf (RED 1.0), which is further from the median than the original node.

Unsupported original nodes are fine as they are (94–100% of those moves were kept).

### R232 check
R232 is the first data with real downward moves, and it confirms step 1 and changes step 2.
- **Step 1 confirmed.** Downward moves from a supported node, share that curators reverted, by `q(original node)`:

  | `q` | R232 bac120 genus | R232 ar53 genus | R232 bac120 family |
  |---|---|---|---|
  | below 0 | 7/140 (5%) | 0/9 | 3/13 |
  | 0–0.25 | 10/129 (8%) | 1/6 | 4/7 |
  | 0.25–0.5 | 23/79 (29%) | 2/6 | 10/11 |
  | 0.5–0.625 | 16/48 (33%) | 4/5 | 4/4 |
  | 0.625–0.75 | 42/60 (70%) | 3/5 | 1/1 |
  | 0.75–1 | 27/32 (84%) | 1/2 | – |

  The genus crossover is at `q` of about 0.625, and for families curators revert nearly everything at `q >= 0`. Both
  match `accept_q`. Moves that end on a single genome are almost always kept (bac120 genus: 160 kept, 9 reverted),
  and they only occur below `q = 0.5`, so step 1 doesn't touch them. Moves from an unsupported node were kept 91–96%
  of the time when the node was below the median, and 75–79% above it (most of the rest were merged away).
- **Step 2 becomes a report flag (Job 7, `split_candidate`).** The R226 rates (55% split at RED 0.84–0.86) were
  measured on trees where no label could move down at all. In R232, after curatomatic's own downward moves, the
  taxa still matching the step 2 pattern (supported node, previous members under one child, other children all new
  genomes) were split by curators 48% of the time at `q` 0–0.25 (30/62), 36% at 0.25–0.5 (71/200), 16% at 0.5–0.625
  and under 6% above. That isn't a rule. The R226 "as checked in" rates for the same pattern (93% below `q = 0`,
  67% at 0–0.25) show what the existing downward code already handles once Job 1 is fixed.

### Where
- Rust: `src/tree.rs` → `refine_existing_labels` (the block after "A well supported taxon is not moved down") and the
  sort in `process_refinement`. Python: `Tree.refine_existing_labels` / `process_refinement` in
  `curatomatic/model/tree.py`.

### What to do
1. If the original node is supported and its RED is at or above `accept(r)`, skip downward refinement for the taxon.
2. Don't add the "prefer a descendant" rule from the R226 write-up. Report those taxa instead (Job 7).
3. Leave the unsupported case as it is.
4. Unit tests for each branch.

### Done when
The exact count rises in both domains. Expected from step 1 on R232 bac120, fixed / broken: genus 69 / 23, family 15 / 1
(ar53: genus 4 / 3). On the R226 re-run (after Job 1) the earlier estimate was 65 / 31 at `q = 0.5`.

## Job 7: Write a review table for curators

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

## Job 8: Don't mint a WGS name for genomes displaced from a taxon they already had

### Problem
`fill_in_missing_ranks` creates a novel taxon for any genome with no label at rank `r`, even when the metadata says
the genome had a taxon `T` of that rank in the previous release and `T` is still on the tree somewhere else. This
happens when `T` is no longer monophyletic (usually an unsupported backbone, e.g. the insect symbionts around
Enterobacteriaceae): the decoration step puts `T` on one clade and leaves the other members bare. Curators rarely
accept a new placeholder there. They put the genomes back in `T`, or in `T` with a new suffix.

### Evidence
Novel taxa, split by the previous-release taxonomy of their member genomes (`--min-bs 95`):

| | R226 bac120 | R232 bac120 | R232 ar53 |
|---|---|---|---|
| No member had a taxon at this rank: novel name kept | 3,642 / 4,214 (86%) | 7,608 / 8,802 (86%) | 733 / 855 (86%) |
| Members share one previous taxon `T`: novel name kept | 81 / 223 (36%) | 231 / 541 (43%) | 10 / 18 |
| ... curated back to `T` or a suffix variant of `T` | 111 / 223 (50%) | 224 / 541 (41%) | 5 / 18 |

In R232 bac120 the 224 are 120 "exactly `T`" and 104 "`T` with another suffix". The effect is strongest above genus
when the nearest other taxon of the rank joins past `child_thresh(r)`: orders 49 of 65 (R226) and 37 of 60 (R232),
phyla 11 of 14 (R232). No tree feature separates "kept" from "put back" well (support on the path to `T` and at the
common ancestor each move the rate by under 15 points), so this is a naming and reporting change, not a placement
rule.

### Where
- Rust: `src/tree.rs` → `fill_in_missing_ranks` and `src/taxon.rs` → `generate_unique_rank`. Python: the matching
  functions in `curatomatic/model/tree.py` and `curatomatic/util/taxon.py`.

### What to do
1. When naming a novel taxon of rank `r`, collect the previous-release taxa of rank `r` of its member genomes. If
   there is exactly one (`T`), name the novel taxon from `T` instead of from a WGS prefix: the first unused
   `<T without suffix>_<A, B, ...>` (the suffix rule GTDB uses for polyphyletic names). Keep it deterministic.
2. Log each at INFO: "Novel <taxon> holds <n> genomes previously in <T>". Add the `displaced_from_prior` flag (Job 7).
3. Needs agreement before it ships, because the chosen suffix will often differ from the curators' choice. If that
   isn't wanted, do step 2 only and keep the WGS name.
4. Unit tests in `src/taxon/tests.rs` for the suffix choice (no suffix in use, `_A` in use, `T` itself has a suffix).

### Done when
Every novel taxon with a single previous taxon is reported, and (if step 1 is agreed) carries a name derived from it.

## Job 9: Resolve same-rank taxa that are nested in the input tree

### Problem
The input tree often has a taxon `X` of rank `r` below another taxon `Y` of the same rank (label on an ancestor).
This is how the decoration step shows that `Y` is paraphyletic. Curatomatic leaves both in place, so every genome
under `X` reads as `X` and the rest of `Y` as `Y`. Curators always resolve it one way or the other.

### Evidence
Pairs not created by curatomatic: 228 in R226 bac120, 354 in R232 bac120, 13 and 4 in ar53.

| Outcome | R226 bac120 | R232 bac120 |
|---|---|---|
| `X` merged into `Y` | 36% | 48% |
| `X` kept, `Y` moved down or split so that it no longer contains `X` | 33% | 31% |
| Other (`Y` merged into `X`, renames) | 31% | 21% |

The best predictor is where `Y` sits (R232 bac120): `Y` unsupported and below `parent_thresh` → `Y` contracted 72%,
merged 12% (n=74); `Y` at or above the median → merged 58–73%, contracted 0–19% (n=111). In R226 the same split is
weaker (50% / 36%, n=14, for the first group). When `X` is not a name from the previous release, `X` almost never
survives (kept 5% in R232, 2% in R226).

### What to do
1. Report every pair with the `nested_in_same_rank` flag (Job 7), with both names, both REDs and the support of `Y`.
2. Don't add an automatic rule yet. The best single split found is 72% in one release and 50% in the other. Once
   the report exists and the R232 inputs are checked in (Job 11), sweep two candidates: move `Y` down when it is
   unsupported and below `parent_thresh(r)`; drop `X` when `Y` is supported and at or above the median.
3. Unit tests in `src/tree/tests.rs` for the report rows.

### Done when
The report lists every nested pair. A rule is added only if it raises the exact count in both releases.

## Job 10: Count genomes without a WGS prefix, and prefer complete metadata

See the R232 check in Job 3. Small, but it explains about a quarter of the `NOVEL-` names that curators replaced:
log "N of M tree genomes have no ncbi_wgs_formatted" at INFO in both implementations, and WARN when the column is
missing altogether (same pattern as Job 1).

## Job 11: Check in the R232 inputs and a scoring script

### Problem
R232 can't be used for parity or to score Jobs 3–9, because `test_data/r232_*` has no input tree, RED dictionary or
metadata. Every R232 number in this file was computed from the output tree and log with reconstructed RED values.
There is also no script in the repository that produces the "exact" count that every job is judged by.

### What to do
1. Add to `test_data/r232_ar53` and `test_data/r232_bac120`: the input `*.decorated.red_decorated.tree`, the
   `*.decorated.dict`, and a metadata file cut down to the three columns curatomatic reads (`accession` or
   `formatted_accession`, `gtdb_taxonomy`, `ncbi_wgs_formatted`). The paths are on line 1 of each `curatomatic.log`.
2. Confirm the Rust binary reproduces both R232 output trees byte for byte. This is the first parity test of the
   downward branch on production data (see Job 1, step 4).
3. Add `scripts/compare_to_curated.py`: given a curatomatic tree and a curated taxonomy, print the exact count per
   rank, and the same count split by how each taxon got there (unmoved, moved up, moved down, novel). Use it for the
   "Done when" of Jobs 3–9.
