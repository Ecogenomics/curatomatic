use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::error::{CuratomaticError, CuratomaticResult};
use crate::genome::{GenomeId, gid_to_canonical};
use crate::rank::Rank;
use crate::taxon::Taxon;
use crate::taxonomy::Taxonomy;

/// Column holding the genome accession.
const COL_ACCESSION: &str = "accession";
/// Values (compared after trimming, ignoring case) that mean a genome has no novel name.
const NO_NAME_VALUES: [&str; 4] = ["", "none", "n/a", "null"];

/// The metadata associated with a single genome.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct MetadataRow {
    pub(crate) accession: GenomeId,
    pub(crate) taxonomy: Option<Taxonomy>,
    /// The name used for novel taxa that contain this genome (e.g. its WGS prefix).
    pub(crate) wgs: Option<String>,
}

/// Metadata for a collection of genomes, keyed by canonical genome ID.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct Metadata {
    pub(crate) rows: BTreeMap<GenomeId, MetadataRow>,
    /// Whether a prior taxonomy was given. Without it no genome has a prior taxonomy, so no label
    /// can move down the tree.
    pub(crate) has_prior_taxonomy: bool,
}

impl Metadata {
    /// Read the novel names from a tab separated file with a header (e.g. the GTDB metadata file).
    ///
    /// Only the `accession` column and `column` are read, both found by header name; any other
    /// column is ignored. A value of `none`, `n/a`, `null` or empty (see [`NO_NAME_VALUES`]) means
    /// the genome has no name. No genome is given a prior taxonomy.
    pub(crate) fn from_file(path: &Path, column: &str) -> CuratomaticResult<Self> {
        let file = File::open(path)
            .map_err(|e| CuratomaticError::Io(format!("{}: {e}", path.display())))?;
        Self::from_reader(BufReader::new(file), path, column)
    }

    /// Read the novel names from a reader, see [`Metadata::from_file`]. `path` is only used in
    /// error messages.
    pub(crate) fn from_reader(
        reader: impl BufRead,
        path: &Path,
        column: &str,
    ) -> CuratomaticResult<Self> {
        let mut lines = reader.lines();
        let header_line = lines
            .next()
            .transpose()
            .map_err(|e| CuratomaticError::Io(format!("{}: {e}", path.display())))?
            .unwrap_or_default();
        let header: HashMap<&str, usize> = split_line(&header_line)
            .enumerate()
            .map(|(idx, name)| (name.trim(), idx))
            .collect();
        let find_column = |name: &str| {
            header
                .get(name)
                .copied()
                .ok_or_else(|| CuratomaticError::MissingColumn {
                    path: path.display().to_string(),
                    column: name.to_owned(),
                })
        };
        let accession_idx = find_column(COL_ACCESSION)?;
        let name_idx = find_column(column)?;

        let mut rows = BTreeMap::new();
        for (line_idx, line) in lines.enumerate() {
            let line =
                line.map_err(|e| CuratomaticError::Io(format!("{}: {e}", path.display())))?;
            if line.trim().is_empty() {
                continue;
            }
            // The header is line 1.
            let line_no = line_idx.saturating_add(2);
            // Values are not trimmed as a line, so empty trailing values keep their columns.
            let cols: Vec<&str> = split_line(&line).collect();
            let value = |idx: usize, name: &str| {
                cols.get(idx).map(|value| value.trim()).ok_or_else(|| {
                    CuratomaticError::MissingValue {
                        line: line_no,
                        column: name.to_owned(),
                    }
                })
            };

            let accession = gid_to_canonical(value(accession_idx, COL_ACCESSION)?);
            let wgs = named(value(name_idx, column)?).map(str::to_owned);
            rows.insert(
                accession.clone(),
                MetadataRow {
                    accession,
                    taxonomy: None,
                    wgs,
                },
            );
        }
        Ok(Metadata {
            rows,
            has_prior_taxonomy: false,
        })
    }

    /// Give `gid` a prior taxonomy, adding a row for it if there is none.
    pub(crate) fn set_taxonomy(&mut self, gid: GenomeId, taxonomy: Taxonomy) {
        self.has_prior_taxonomy = true;
        self.rows
            .entry(gid.clone())
            .or_insert_with(|| MetadataRow {
                accession: gid,
                taxonomy: None,
                wgs: None,
            })
            .taxonomy = Some(taxonomy);
    }

    /// Retain only the rows whose genome ID is in `ids`.
    pub(crate) fn subset_to_ids(&mut self, ids: &BTreeSet<GenomeId>) {
        self.rows.retain(|gid, _| ids.contains(gid));
    }

    /// Every named taxon in any prior taxonomy (placeholders such as `g__` are left out).
    pub(crate) fn all_taxa(&self) -> BTreeSet<Taxon> {
        self.rows
            .values()
            .filter_map(|row| row.taxonomy.as_ref())
            .flat_map(Taxonomy::iter_ranks)
            .filter(|taxon| taxon.is_named())
            .cloned()
            .collect()
    }

    /// The named taxa of `rank` that the metadata assigns to `gid`.
    pub(crate) fn taxa_of_rank(&self, gid: &GenomeId, rank: Rank) -> impl Iterator<Item = &Taxon> {
        self.rows
            .get(gid)
            .and_then(|row| row.taxonomy.as_ref())
            .map_or(&[][..], |taxonomy| taxonomy.get_rank(rank))
            .iter()
            .filter(|taxon| taxon.is_named())
    }

    /// True if the metadata gives `gid` a taxonomy, i.e. the genome was in the previous release.
    pub(crate) fn has_taxonomy(&self, gid: &GenomeId) -> bool {
        self.rows.get(gid).is_some_and(|row| row.taxonomy.is_some())
    }

    pub(crate) fn n_with_taxonomy(&self) -> usize {
        self.rows
            .values()
            .filter(|row| row.taxonomy.is_some())
            .count()
    }

    /// The taxonomy of each genome in `gids` that has one.
    pub(crate) fn get_taxonomy_from_gids(
        &self,
        gids: &BTreeSet<GenomeId>,
    ) -> BTreeMap<&GenomeId, &Taxonomy> {
        gids.iter()
            .filter_map(|gid| {
                let (gid, row) = self.rows.get_key_value(gid)?;
                Some((gid, row.taxonomy.as_ref()?))
            })
            .collect()
    }
}

/// The fields of a tab separated line, without its line ending.
fn split_line(line: &str) -> impl Iterator<Item = &str> {
    line.trim_end_matches(['\r', '\n']).split('\t')
}

/// `None` if the (trimmed) value is one of [`NO_NAME_VALUES`].
fn named(value: &str) -> Option<&str> {
    let is_placeholder = NO_NAME_VALUES
        .iter()
        .any(|placeholder| value.eq_ignore_ascii_case(placeholder));
    (!is_placeholder).then_some(value)
}

#[cfg(test)]
mod tests;
