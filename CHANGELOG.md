# Changelog

## 1.0.0 (2026-11-23)

### ⚠ BREAKING CHANGES

* **tree:** --taxonomy is now required, as leaves can only get
their species from it.

Also removes the unused log::warn import left in taxonomy.rs.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>

### Features

* **cli:** add --debug to write the loaded inputs with dummy nodes ([641518a](https://github.com/Ecogenomics/curatomatic/commit/641518a31eb59b7624bf40c36b3c5c05cf09ec1a))
* **cli:** add --prior-taxonomy as the only prior for moving taxa down (TODO Job 13) ([2622c2f](https://github.com/Ecogenomics/curatomatic/commit/2622c2f73a0a3814b5159231b82fcab5e05730c7))
* **input:** read gzipped input files (TODO Job 27) ([a82eab2](https://github.com/Ecogenomics/curatomatic/commit/a82eab2d9ed8ebd74742a05b60923fad683a895e))
* **main:** check that taxa are only ever moved or created ([a3aee1e](https://github.com/Ecogenomics/curatomatic/commit/a3aee1e6e1e7dc864d965e1686e77196cd403e19))
* **metadata:** Add --novel-accession to name the genome ID column of the novel name file (TODO Job 8) ([9ac6735](https://github.com/Ecogenomics/curatomatic/commit/9ac67353c7e2e74ca3d42b885f03b9ab2093c3fb))
* **meta:** Update metadata to use gtdb_taxonomy over taxonomy ([e72d923](https://github.com/Ecogenomics/curatomatic/commit/e72d923335eef14cfeefe4856d034432ec65f8a7))
* **novel:** report an unresolvable novel taxon, and never write it ([3dc483a](https://github.com/Ecogenomics/curatomatic/commit/3dc483a2f385d034b0e996e6d159a6ef70a82b89))
* **output:** write the taxa held by each leaf to a leaf taxonomy file ([825c317](https://github.com/Ecogenomics/curatomatic/commit/825c317b80802b69fff4f79327084b4ce2b3b446))
* **rank:** require every taxon to be `<prefix>__<name>` ([b1a4675](https://github.com/Ecogenomics/curatomatic/commit/b1a46758b036ca7f541f3d2c1f6a3eaa9f2b3bec))
* **report:** only list taxa that were moved or created ([537f048](https://github.com/Ecogenomics/curatomatic/commit/537f0483fe9cc1ad89e2ceee3cd1765940b6e6fc))
* **rust:** Re-write 1:1 parity in Rust. ([c1350e9](https://github.com/Ecogenomics/curatomatic/commit/c1350e9241ddb0df9633574514379044e25d0912))
* **scripts:** add a script to compare the Python and Rust outputs ([d4e0c75](https://github.com/Ecogenomics/curatomatic/commit/d4e0c75328be8b92340a6037eda434d2157868df))
* **scripts:** add a script to convert the Python inputs to the Rust format ([8622967](https://github.com/Ecogenomics/curatomatic/commit/862296750e44459eb461028570b6c04653c10bf8))
* **scripts:** add a script to score the outputs against a curated taxonomy ([241ec64](https://github.com/Ecogenomics/curatomatic/commit/241ec644cfef7f89f012f97668ca34c10651a1de))
* **scripts:** Add parsing script to create input file. ([6c2f396](https://github.com/Ecogenomics/curatomatic/commit/6c2f396f1fef28e7dcdb48e3ebce0f4531d9d1f0))
* **scripts:** score taxa on nodes with several taxa of a rank as ambiguous ([0672bda](https://github.com/Ecogenomics/curatomatic/commit/0672bda87b1a682aa46da7f4240ccbe215db00b7))
* **taxonomy_file:** require a full seven rank taxonomy in every prior taxonomy row ([2d562cc](https://github.com/Ecogenomics/curatomatic/commit/2d562cce4f9041ad5f95def6f267dd60cb16fc88))
* **taxonomy_file:** require rows in rank order, each taxon once, one domain and one species ([da1a5bf](https://github.com/Ecogenomics/curatomatic/commit/da1a5bfe5c794b6dc0cdd6acd77020b4c6039103))
* **taxonomy:** Add --taxonomy to name the ranks the tree can't hold (TODO Job 6) ([c1d66c2](https://github.com/Ecogenomics/curatomatic/commit/c1d66c296aa3a2a3294a0f9216f5029c17ec5ab7))
* **taxonomy:** ban placeholders on the tree and strip them from taxonomy files ([40e038d](https://github.com/Ecogenomics/curatomatic/commit/40e038dab692a88ad34c7a8443071dde89320633))
* **tree:** a genus is not moved down away from new genomes its species names ([f439325](https://github.com/Ecogenomics/curatomatic/commit/f43932599a25938cb266f6bccedfaa1696a44bf1))
* **tree:** a skipped novel taxon only skips its own rank (TODO Job 12) ([7513c7d](https://github.com/Ecogenomics/curatomatic/commit/7513c7dc0dcd16c289d0cde7c9590b07d8e5708d))
* **tree:** check each genome's taxa on the tree equal its row, and have a domain ([69e5607](https://github.com/Ecogenomics/curatomatic/commit/69e56077ce79108e48b113fcf4139f355b6f8b39))
* **tree:** don't create a novel taxon out of rank order (TODO Job 16) ([3be6fb4](https://github.com/Ecogenomics/curatomatic/commit/3be6fb44f5fe77ed0e3775bd8bc232fac40bd704))
* **tree:** keep a non-Latin genus on a supported node with the new genomes its species name ([5e27cd8](https://github.com/Ecogenomics/curatomatic/commit/5e27cd84485371543caffb6386a1df391b9f295d))
* **tree:** never move a taxon up above a taxon of a higher rank (TODO Job 10) ([9a6125a](https://github.com/Ecogenomics/curatomatic/commit/9a6125ad618d7bfd47c70fefbcecec02ccd60877))
* **tree:** only create novel taxa where no valid name is below (TODO Job 9) ([b71efa4](https://github.com/Ecogenomics/curatomatic/commit/b71efa4c9233dccb3d593f462afc8fd42e4f90bf))
* **tree:** only Latin genera inside the RED window keep new genomes their species name ([e073a16](https://github.com/Ecogenomics/curatomatic/commit/e073a164379dec89454210ebe0e1dc33a53b7b4d))
* **tree:** require every leaf to hold its named species ([7066cd2](https://github.com/Ecogenomics/curatomatic/commit/7066cd2bcf6691f46d3aea24b6c78960dc4e6fbc))
* **tree:** species no longer block a novel genus (TODO Job 11) ([fd77429](https://github.com/Ecogenomics/curatomatic/commit/fd7742986a3731bb74ec93cbf3f4aac9f5b396d6))

### Bug Fixes

* **meta:** Warn when the metadata has no gtdb_taxonomy column (TODO Job 1) ([796f56a](https://github.com/Ecogenomics/curatomatic/commit/796f56aa1ef660e6e306d7132c5b8d39cb85f769))
* **output:** label leaves of the output tree with their genome ID only ([266f503](https://github.com/Ecogenomics/curatomatic/commit/266f503d2a5b19f8d0bb3f3a11e9e8d4f05a7c39))
* **report:** write the bootstrap and RED columns exactly as read ([bea2984](https://github.com/Ecogenomics/curatomatic/commit/bea298460ffc96e68d14cfcbaba1bfef8bc8af05))
* satisfy rustfmt and clippy lints from Rust 1.100 ([7e4f9fb](https://github.com/Ecogenomics/curatomatic/commit/7e4f9fb3046e743c6ce5fa6f3b07ad3aa951024b))
* **tree:** keep the RED and bootstrap of a rewritten label exactly as read ([f13aed8](https://github.com/Ecogenomics/curatomatic/commit/f13aed8f7615c669a970459c3925054e2662ad32))
* **tree:** restore the dummy node taxa from the taxonomy file by taxon, not by rank ([cd09ece](https://github.com/Ecogenomics/curatomatic/commit/cd09ecef6475b8ae3930fcd65da1b0f1af7cfaff))
* **tree:** Skip upward refinement for nested taxa (TODO Job 2) ([e4968fa](https://github.com/Ecogenomics/curatomatic/commit/e4968fad5d31a2b1000ab5dc8c3c6712183baa1a))
