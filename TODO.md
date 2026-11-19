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
