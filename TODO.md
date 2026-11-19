# TODO

## Job 7: Move the CLI definition into its own file, `cli.rs`

### Problem
The `Cli` struct (the clap argument definition) lives in `src/main.rs`, mixed in with `run`, `curate` and the other
pipeline steps. It should be in its own module.

### Where
- `src/cli.rs` (new): the `Cli` struct and its doc comments, unchanged, with the imports it needs (`clap::Parser`,
  `BootstrapValue`, `LogLevel`, `PathBuf`). Make `Cli` and its fields `pub(crate)`.
- `src/main.rs`: add `mod cli;`, import `crate::cli::Cli`, and remove the struct and any imports only it used.
- Tests: move the CLI parsing tests (`parse_cli` and the tests that use it) from `src/main/tests.rs` to
  `src/cli/tests.rs`, declared with `#[cfg(test)] mod tests;` in `src/cli.rs`.

### Done when
The command line (arguments, defaults, help text and `--version`) is unchanged, `src/main.rs` no longer defines `Cli`,
and `cargo fmt`, `cargo clippy --all-targets` and `cargo nextest run` pass.

## Job 8: Add `--novel-accession` to name the genome ID column of the novel name file

### Problem
`Metadata::from_reader` (`src/metadata.rs`) always looks for the genome IDs in a column named `accession`
(`COL_ACCESSION`). The example novel name file, `example/gtdb_r237_bac120_metadata_slim.tsv`, calls that column
`db_name`:
```
db_name	ncbi_wgs_formatted
G000008885	none
```
so the run stops with `Missing required column accession in: example/gtdb_r237_bac120_metadata_slim.tsv`.

### Rules
- A new option, `--novel-accession NAME`: the column of `--novel-name-file` that holds the genome IDs. It defaults to
  `accession`, so the GTDB metadata file and every current command line work unchanged.
- It requires `--novel-name-file`, like `--novel-name-column` does. Given on its own, it is a clap
  `MissingRequiredArgument` error.
- The column is found by header name, as now. A missing column is the existing `MissingColumn` error, naming the
  column that was asked for (e.g. `db_name`), not `accession`.
- The values are read as now: trimmed, then made canonical with `gid_to_canonical` (`G000008885` is kept as is,
  `GB_GCA_000008885.1` becomes `G000008885`).
- `--novel-accession` and `--novel-name-column` may name the same column only if that is intended, so don't reject it.

### Where
- `src/main.rs` (or `src/cli.rs` if Job 7 is done first): the `--novel-accession` argument, `value_name = "NAME"`,
  `default_value = "accession"`, `requires = "novel_name_file"`. Update the doc comment of `--novel-name-file`, which
  says it reads the `accession` column. Pass the column to `Metadata::from_file` in `read_metadata`, and log it with
  the other `debug!` lines in `run`.
- `src/metadata.rs`: `from_file` and `from_reader` take the accession column name as a parameter. `COL_ACCESSION`
  becomes the default only, used by the CLI. Update the doc comments.
- Tests:
  - `src/metadata/tests.rs`: a file with `db_name` as the accession column is read with `--novel-accession db_name`;
    a missing accession column names the requested column in the error. Existing tests pass `accession`.
  - `src/main/tests.rs` (or `src/cli/tests.rs`): the default is `accession`; `--novel-accession db_name` is parsed;
    `--novel-accession` without `--novel-name-file` is an error.
- Docs: `README.md` (the options list, and the "Novel name file" section, which says only `accession` is read).

### Done when
This runs to the end and names novel taxa from `ncbi_wgs_formatted`:
```shell
./target/release/curatomatic --min-bs 95.0 \
    --taxonomy example/gtdb_r237_decorated_taxonomy.tsv \
    --novel-name-file example/gtdb_r237_bac120_metadata_slim.tsv \
    --novel-accession db_name --novel-name-column ncbi_wgs_formatted \
    example/gtdb_r237_median_red_marked.nwk example/gtdb_r237_median_red_rank_distribution.tsv output
```
The log says `matched 229,307 IDs to the tree`. Commands without `--novel-accession` behave exactly as before, and
`cargo fmt`, `cargo clippy --all-targets` and `cargo nextest run` pass.
