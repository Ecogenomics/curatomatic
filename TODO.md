# TODO

Finished jobs are removed from this file: their rules are in `docs/refinement.md` and `docs/differences.md`, and
their measurements in `docs/testing.md`. Jobs that need the user's approval say so: don't start them before it is
given.

## Job 14: Report an unresolvable novel taxon, and never write it

**Status: approved** (re-worded by the user on 2026-11-23: an `UNRESOLVABLE` taxon is no longer written).

### Problem

`generate_unique_rank` (`src/taxon.rs`) names a novel taxon without a WGS name after its lowest genome ID whose name
isn't already a taxon. If every genome's name is taken, the taxon is the sentinel `<rank>__UNRESOLVABLE`, which
`fill_in_missing_ranks` (`src/tree.rs`) places on the tree like any other novel taxon. This happens when the input
tree already has a taxon named after one of the genomes, but that genome isn't inside it. For example, `G020350785`
in r226_bac120 gets `g__UNRESOLVABLE`, because `g__G020350785` is already on the tree elsewhere (the curators put the
genome into `g__G020350785`). In r237 bac120, G965247515 gets `g__UNRESOLVABLE` the same way.

Refusing the name is right, but:

- the sentinel is written to the curated tree and the taxonomy output files, as if it were a taxon. The user no
  longer wants it there.
- the report doesn't say why the name was refused.

### Rule

1. **Never written.** When `generate_unique_rank` finds no name, no taxon is placed on the candidate node. The
   genomes keep the placeholder for that rank (e.g. `g__`) in `<stem>.curatomatic.taxonomy.tsv`, exactly as for a
   novel taxon skipped because a valid name is below it, and nothing is added to the curated tree or the leaf
   taxonomy file. Only that rank is affected: the candidates and the `lowest_valid_rank_below` check of
   every rank are worked out before any novel taxon is created, so the other novel taxa are placed as before.
2. **Still reported.** Each such case is a row of the report: taxon `<rank>__UNRESOLVABLE`, a new action
   `unresolvable` (nothing was created, so not `novel`), the `*_after` columns those of the candidate node it would
   have been placed on, and a new flag `namesake_outside_taxon`. A new column, `namesakes`, lists the existing taxa
   named after its genomes (e.g. `g__G020350785`), which also says which genomes it is. The run logs the number of
   unresolvable taxa for each rank, beside the number skipped.
3. **The sentinel is no longer a taxon.** As it is never on the tree, it loses its exception to "every taxon is on
   exactly one node": remove the `is_unresolvable` exceptions in `Tree::get_taxon_labels` and
   `check_taxa_unchanged`, and the comments on `TaxonCounts`. An input tree with `g__UNRESOLVABLE` on two nodes is
   then a `NonUniqueTaxon` error like any other. `new_taxa_created` never holds it, so the taxa check and
   `check_every_change_has_a_row` need no special case; the unresolvable cases reach the report separately.

### Tasks

1. Make `generate_unique_rank` return either a name or the unresolvable case with the taxa that blocked each name
   (a new enum, not the sentinel `Taxon`), and make `fill_in_missing_ranks` return the unresolvable cases (rank,
   candidate node, blocking taxa) beside the taxa it created, without placing anything for them.
2. Report (`src/report.rs`): new `Action` and `Flag` variants with doc comments, the `namesakes` column, and the
   rows for the unresolvable cases. Every unresolvable case must have a row (extend
   `check_every_change_has_a_row`).
3. Remove the one-node exception (rule 3) and update its tests (`an_unresolvable_taxon_may_be_on_more_than_one_node`
   in `src/tree/tests.rs`, `taxa_check_accepts_an_unresolvable_novel_taxon_on_more_than_one_node` and the
   `g__UNRESOLVABLE` cases in `src/main/tests.rs`, the counts in `src/taxon/tests.rs`).
4. Tests: in `src/tree/tests.rs`, an unresolvable genus leaves its genomes with `g__` in the taxonomy output and
   nothing on the tree, and the other ranks of those genomes are still filled; in `src/taxon/tests.rs`, the blocking
   taxa are returned; in `src/report/tests.rs`, replace `every_node_given_the_unresolvable_sentinel_has_a_row` with
   one row per unresolvable case, with its action, flag and namesakes.
5. Docs: `docs/refinement.md` ("Novel taxa" and "Invariants"), `README.md` (the novel name section and the report
   columns), `docs/differences.md` (L4 and the "Every genome ID already taken" row: Python writes
   `<rank>__NOVEL-<gid>`, Rust now writes nothing and reports it). `docs/architecture.*` are snapshots: leave them.
6. Check the r226_bac120 cases (`G020350785`, `G030406045`), the r232_bac120 case (`G020350445`) and the r237 bac120
   case (`G965247515`): no `UNRESOLVABLE` in the tree or taxonomy outputs, one report row each with the right
   namesakes. Measure the curator workload on the four datasets (see `docs/testing.md`).

## Job 16: Don't create a novel taxon that breaks rank order

**Status: needs the user's approval** (proposed as P6 on 2026-11-20).

### Problem

`fill_in_missing_ranks` (`src/tree.rs`) can create a novel taxon that is out of rank order with the named taxa around
it. There are two forms:

- A novel taxon **below** a named taxon of a lower rank on its path. In r237 bac120, 10 novel orders are created
  inside `g__MED-G20`, and `o__G000829155` inside `g__Karelsulcia` (RED 0.891).
- A novel taxon **above** a named taxon of a higher rank. In r237 bac120, `g__CAZDDW01` is created above
  `f__DNKR01` and the novel `f__G965231015`.

With Jobs 10 to 13 there are 7 such taxa in r237 bac120 (Job 15 left that unchanged), and 13 / 15 in r226_bac120 /
r232_bac120 (6 / 10 at `32d3350`, before Jobs 11 and 12; not measured after Job 15). Every one checked is placed
elsewhere by the truth (`o__Flavobacteriales_D`, `o__Bacillales`, `f__CAZDDW01;g__G965231015`).

### Rule (proposed)

A novel taxon of rank `R` is skipped if its candidate node is below a named taxon of a rank lower than `R` on its
path, or has a named taxon of a rank higher than `R` below it. Novel taxa of higher ranks made earlier in the same
run count. Skipping only skips that rank (as for the existing check, `docs/refinement.md` "Novel taxa").

### Tasks

1. Add the check to `fill_in_missing_ranks`, next to the `lowest_valid_rank_below` check.
2. Tests in `src/tree/tests.rs`, one for each form, and one where a novel taxon of a higher rank made earlier in the
   run blocks a lower one.
3. Docs: `docs/refinement.md` "Novel taxa", `docs/differences.md` (new entry: Python has no such check).
4. Measure the curator workload on the four datasets (see `docs/testing.md`), and count the rank inversions left in
   r237 bac120.

## Job 17: Keep RED and bootstrap values as they were read when a label is rewritten

**Status: needs the user's decision.**

### Problem

`NodeLabel::to_label` (`src/node_label.rs`) writes the RED with 3 decimal places and the bootstrap with 1
(`100.0`), so a rewritten label differs from the input in format and precision (1,574 nodes in the r237 bac120
output with Jobs 10 to 13). The later refinement passes read the rounded RED back. This copies the Python behaviour,
but `AGENTS.md` says curatomatic never changes RED or bootstrap values.

### Question for the user

Should a rewritten label keep the RED and bootstrap exactly as read (no rounding), and should the later passes use
the unrounded RED? Either answer changes some outputs, so measure it on the four datasets.

## Job 18: Benchmark with a prior taxonomy and score the curator workload

**Status: not started. Confirm with the user before starting.**

### Problem

- `scripts/benchmark.sh` runs curatomatic without `--prior-taxonomy`, so it measures the behaviour without downward
  refinement (Jobs 13 and 15 are not measured).
- The curator workload, the user's metric of success, is scored by `wip/scripts/workload.py` and
  `wip/scripts/summary.py`, which are not committed (`wip/` is git-ignored).

### Tasks

1. Make the prior files in `benchmark.sh` from the Python metadata (recipes in `docs/testing.md`, "Prior taxonomy")
   and pass them with `--prior-taxonomy`.
2. Move the workload scorer into `scripts/`, run it in `benchmark.sh`, and document its metrics in
   `docs/testing.md`.
3. Check that `benchmark.sh` gives the numbers in `docs/testing.md`, "Curator workload".

## Job 25: Release a Linux binary on every push to main

**Status: approved** (decisions from the user on 2026-11-23 are in "Decisions").

### Problem

There are no releases: a user has to clone the repository and build curatomatic themselves (`README.md`, "Building"),
and the version in `Cargo.toml` (`0.1.0`, logged on every run as `curatomatic v{CARGO_PKG_VERSION}`) is never
changed. There is no `.github/` directory.

### Rule

On a push to `main`, a GitHub Actions workflow (`.github/workflows/release.yml`) releases a new version if the
Conventional Commits since the last tag call for one, and does nothing otherwise. The steps, in order:

1. **Check `main`.** Run `cargo fmt --check`, `cargo clippy --all-targets` and `cargo nextest run` (the same
   checks as "Development" in `AGENTS.md`). If any fails, the workflow stops and nothing is released.
2. **Work out the version.** semantic-release, in a dry run, reads the commits since the last `v*` tag:
   `fix`/`perf` give a patch, `feat` a minor, and `!` or a `BREAKING CHANGE:` footer a major. Other types (`docs`,
   `chore`, `refactor`, `test`, `ci`, ...) give no release, and the workflow stops. With no tag yet, the first
   release is `v1.0.0` (semantic-release's default, so no tag is made by hand first).
3. **Set the version, write the changelog and tag it.** semantic-release writes the version to `Cargo.toml` and
   `Cargo.lock` (`@semantic-release/exec`), adds the release notes to the top of `CHANGELOG.md`
   (`@semantic-release/changelog`), commits the three files to `main` as `chore(release): <version> [skip ci]`
   (`@semantic-release/git`), and tags that commit `v<version>`. The binary must be built after this, from the tag, so
   that the version it logs is the released one.
4. **Create the GitHub release** `v<version>` (`@semantic-release/github`), its description the release notes made
   from the commits (`conventional-changelog-conventionalcommits` preset: features, bug fixes, performance, breaking
   changes), the same notes as in `CHANGELOG.md`.
5. **Build the Linux binary and attach it.** Check out the tag, build a static `x86_64-unknown-linux-musl` release
   binary (no glibc dependency, so it runs on old cluster nodes), package it as
   `curatomatic-v<version>-x86_64-unknown-linux-musl.tar.gz` (binary, `README.md`, `CHANGELOG.md`, licence) with a
   `.sha256` beside it, and upload both with `gh release upload`. Only the musl binary is built.
6. **Bring `dev` up to date.** Merge the tag into `dev` as `chore(release): Merge v<version> into dev` and push, so
   `dev` has the version commit and the changelog. A conflict fails this job and pushes nothing. It runs beside the
   build (both after the release), so a failed build still leaves `dev` in sync.

The build runs in the same workflow (a `needs:` job after the release), because a release created with
`GITHUB_TOKEN` does not start other workflows.

Use the workflow pasted by the user from another Rust project as the reference (summarised here, as it is not in this
repository):

- `on: push: branches: [main]` and `workflow_dispatch`; `concurrency: {group: release, cancel-in-progress: false}`;
  `permissions: contents: read` at the top, `contents: write` only on the jobs that need it.
- semantic-release run with `npx --yes -p semantic-release@25 -p @semantic-release/changelog@7 -p
  @semantic-release/git@11 -p @semantic-release/exec@7 -p conventional-changelog-conventionalcommits@9
  semantic-release` (pinned majors), on Node 24 (`actions/setup-node`), with `actions/checkout` at `fetch-depth: 0`
  so it sees the tags.
- Jobs: `plan` (runs the checks of step 1, then `--dry-run`, outputs the next version, empty if none) → `release`
  (if a version: writes the changelog, commits, tags, creates the release, outputs the version) → `build` (checks out
  the tag, builds, uploads with `GH_TOKEN: ${{ github.token }}`) and `sync-dev` (the merge of step 6), both after
  `release`.
- cargo-nextest is installed from its prebuilt binary
  (`curl -LsSf https://get.nexte.st/latest/linux | tar zxf - -C "${CARGO_HOME:-$HOME/.cargo}/bin"`), not built from
  source.

The semantic-release config goes in `.releaserc.json` (`branches: ["main"]`, `tagFormat: "v${version}"`). Setting the
version needs no new crate or cargo tool: `sed` the `version` line of the `[package]` table, then
`cargo update --workspace --offline` to update `Cargo.lock`. This is not a crate dependency, so it doesn't need the
approval in "Rules", but no other new tool should be added without asking.

### Tasks

1. `.releaserc.json` and `.github/workflows/release.yml` as above.
2. Check the release config locally with `npx semantic-release --dry-run --no-ci` on a branch (it reports the
   version and the notes it would make without publishing), and check the musl build locally or in a test workflow
   run (`cargo build --release --locked --target x86_64-unknown-linux-musl`, then run the binary on the `example/`
   data).
3. Docs: `README.md` (a new "Installation" section: download the binary from the releases page, check the
   `.sha256`; keep "Building"), and `AGENTS.md` "Development" (pushes to `main` are checked and released; the commit
   type decides the version, so a breaking change to the CLI or the output files must be marked with `!`;
   `CHANGELOG.md` is written by the release workflow and never edited by hand).
4. Tell the user to check that `GITHUB_TOKEN` (`github-actions[bot]`) may push to `main` and `dev`, and create tags,
   if either branch is protected.

### Decisions

Made by the user on 2026-11-23:

1. **Checks before releasing:** yes. `main` is checked (step 1) before anything is released.
2. **The version commit:** yes. The version is committed to `main` and merged into `dev` (`sync-dev`, step 6).
3. **musl or glibc:** musl only.
4. **Starting version:** `v1.0.0`.
5. **Changelog file:** yes. `CHANGELOG.md` is written and committed with the version.

## Open questions (ask the user before starting)

- Moving down (Job 13) adds work as well as saving it: the operations added against the input went from 764 at
  `32d3350` to 1,897 with Jobs 10 to 13 (Python adds 1,814), and the taxon specificity in r226_bac120 dropped from
  99.95 to 99.56. Should the added operations be looked at?
- The effect of Job 10 was only measured without a prior (+3 operations, r226_bac120 only). It was not measured
  separately with Job 13.
