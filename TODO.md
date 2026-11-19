# TODO

## Job 10: Never move a taxon up above a taxon of a higher rank

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
