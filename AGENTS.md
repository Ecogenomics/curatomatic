# AGENTS.md


## About
Curatomatic was originally written as a Python library to perform pre-curation of a domain-specific prokaryotic phylogenetic tree.
The logic of the program follows several cases that have been identified as manual steps that curators of the tree perform. The aim is to remove as much manual work as possisble.

The goal of this project is to re-write this project using the Rust framework to ensure type safety and improve performance for large trees. 
The tree uses information based on Relative Evolutionary Divergence (RED) and a pre-decoration of the tree using a different software package.

## Goal
Ensure there is a 1:1 parity between the original Python library and the Rust implementation. The output of the Rust program should be identical to the output of the Python library.

## Test data
Test data has been generated from the Python library and exists under `./test_data`. You can find the command used to run curatomatic in the first line of each log file, e.g.: `test_data/r226_bac120/curatomatic_bs95.0/curatomatic.log`.

Each domain (ar53 - archaea, bac120 - bacteria) has been run three times using different minimum bootstrap values. All curatomatic output is under each test directory with the name `curatomatic_bsVALUE`. Input data is in the root of each directory - both `r226_ar53` and `r226_bac120`.

You are to compare the output of the Rust program to ensure that the decorated trees remain the same. The structure of the tree should always remain the same, bootstrap values should never change, the only changes to the final tree are the labels at each node (which form the taxonomy).

Ensure that there is 1:1 parity for each of the input/output data using the curomatic output as a source of truth.


## Rules
- Tests live in a separate file, never inline: declare `#[cfg(test)] mod tests;` in `foo.rs` and put the tests in `foo/tests.rs` (see `src/taxonomy.rs`).
- When possible, use newtypes.

## Rules (for development only)
- During this active stage of development, do not run lints (including clippy).
