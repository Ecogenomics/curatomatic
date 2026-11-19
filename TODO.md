# TODO

## Job 9: Only create novel taxa where nothing below is named

### Problem

`fill_in_missing_ranks` (`src/tree.rs`) creates a novel taxon for every rank a genome is missing, and names
it after a WGS prefix (or genome ID), even when valid names are below it. In the R237 example, the
c__Chlorobiia clade had `p__` in the decorated taxonomy. Low bootstrap support (53) above it split the clade into six
novel phyla (`p__AASE01`, `p__AQXG01`, `p__ATTH01`, `p__CALRWV02`, `p__JBBSVG01`, `p__BAAFUN01`), all sitting
above named classes, orders, genera and species. A curator would never place a novel rank above existing names.

### Rule

A novel taxon may only be created (and given a WGS or genome ID name) when **no valid name is below it**. Otherwise
no taxon is created, and no other name is tried.

- For a missing rank `R` placed on candidate node `N` (the highest candidate, as today), the novel taxon is
  created only if no genome at or below `N` has a valid name at any rank below `R`, species included.
- A **valid name** is any taxon that is not a placeholder:
  - Placeholders (`c__`, `o__`, `f__`, `g__`) are not valid names.
  - `s__{unresolved} NAME` species are not valid names.
  - Every other taxon is a valid name, including named species (`s__Chlorobium limicola`), species named after a
    placeholder genus (`s__CAYYKD01 sp964471825`) and WGS-style taxa (`g__CAYYKD01`).
- The names are checked against the taxonomy as it is before novel creation (after refinement), i.e. the taxa on
  the nodes at or below `N` and in the `Taxonomy` of each `Leaf`. Taxa of lower ranks on `N` itself count (e.g.
  `c__Chlorobiia` on the same node as the would-be phylum).
- **Skipped genomes keep the placeholder** (`p__`) in the output taxonomy and tree. No report row or flag is added.
- **Lower ranks are not filled either**: once a rank is skipped for a genome, none of its lower missing ranks are
  filled. Ranks are still processed from highest to lowest.
- Unchanged: candidate placement (bootstrap, RED midpoint, stopping at an existing taxon of the rank), naming of
  the taxa that are created (`generate_unique_rank`), and `<rank>__UNRESOLVABLE`.

### Expected effect (R237 example)

A rough count from the current output, checking each novel taxon on its own and ignoring the "lower ranks are not
filled" cascade, so the real number kept will be lower:

| Rank | Kept | Skipped |
|---|---|---|
| p | 10 | 34 |
| c | 22 | 9 |
| o | 86 | 118 |
| f | 379 | 152 |
| g | 3,311 | 460 |

### Tasks

1. In `fill_in_missing_ranks`, after choosing the highest candidate node for a rank, skip it if any valid name of a
   lower rank is at or below it. Record the skipped leaves so their lower missing ranks are skipped too.
2. Add a predicate on `Taxon` for "valid name" (named, and not `s__{unresolved}`), and use it rather than
   repeating the check.
3. Log the number of novel taxa skipped per rank (INFO).
4. Tests in `src/tree/tests.rs`:
   - A missing phylum above a named class is not created, and the genomes keep `p__`.
   - A missing genus above only `s__{unresolved}` species is created.
   - A missing genus above a named species is not created.
   - A WGS-style name below (e.g. `g__CAYYKD01`) blocks creation.
   - A lower rank of the same name on the candidate node itself blocks creation.
   - When a phylum is skipped, a missing family of the same genomes is not filled, even when nothing named is
     below the family's candidate node.
5. Docs:
   - `docs/refinement.md` "Novel taxa": add the rule.
   - `docs/differences.md`: add a row (Python creates novel taxa above named taxa, Rust doesn't).
6. Check that the end-of-run check (taxa after minus novel taxa equals taxa before) still holds.
