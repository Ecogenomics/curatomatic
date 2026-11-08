use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::error::{CuratomaticError, CuratomaticResult};
use crate::genome::{GenomeId, gid_to_canonical};
use crate::taxonomy::Taxonomy;

/// Column holding the genome accession.
const COL_ACCESSION: &str = "accession";
/// Column holding the NCBI WGS prefix.
const COL_WGS: &str = "ncbi_wgs_formatted";
/// Column holding the GTDB taxonomy string.
const COL_TAXONOMY: &str = "gtdb_taxonomy";
/// Value used in the metadata file to indicate that there is no value.
const NONE_VALUE: &str = "none";

/// The metadata associated with a single genome.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct MetadataRow {
    pub(crate) accession: GenomeId,
    pub(crate) taxonomy: Option<Taxonomy>,
    pub(crate) wgs: Option<String>,
}

/// Metadata for a collection of genomes, keyed by canonical genome ID.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct Metadata {
    pub(crate) rows: BTreeMap<GenomeId, MetadataRow>,
}

impl Metadata {
    /// Read a tab separated metadata file.
    ///
    /// The `accession` column is required. The `ncbi_wgs_formatted` and `gtdb_taxonomy` columns
    /// are optional, and a value of `none` in either is treated as missing.
    pub(crate) fn from_file(path: &Path) -> CuratomaticResult<Self> {
        let file = File::open(path)
            .map_err(|e| CuratomaticError::Io(format!("{}: {e}", path.display())))?;
        Self::from_reader(BufReader::new(file))
    }

    /// Read tab separated metadata from a reader, see [`Metadata::from_file`].
    pub(crate) fn from_reader(reader: impl BufRead) -> CuratomaticResult<Self> {
        let mut lines = reader.lines();
        let header_line = lines
            .next()
            .transpose()
            .map_err(|e| CuratomaticError::Io(e.to_string()))?
            .unwrap_or_default();
        let header: HashMap<&str, usize> = header_line
            .trim()
            .split('\t')
            .enumerate()
            .map(|(idx, name)| (name, idx))
            .collect();

        let accession_idx = *header
            .get(COL_ACCESSION)
            .ok_or_else(|| CuratomaticError::MissingColumn(COL_ACCESSION.to_owned()))?;
        let wgs_idx = header.get(COL_WGS).copied();
        let taxonomy_idx = header.get(COL_TAXONOMY).copied();

        let mut rows = BTreeMap::new();
        for (line_idx, line) in lines.enumerate() {
            let line = line.map_err(|e| CuratomaticError::Io(e.to_string()))?;
            if line.trim().is_empty() {
                continue;
            }
            // The header is line 1.
            let line_no = line_idx.saturating_add(2);
            let cols: Vec<&str> = line.trim().split('\t').collect();
            let column = |idx: usize, name: &str| {
                cols.get(idx)
                    .copied()
                    .ok_or_else(|| CuratomaticError::MissingValue {
                        line: line_no,
                        column: name.to_owned(),
                    })
            };

            let accession = gid_to_canonical(column(accession_idx, COL_ACCESSION)?);

            let mut taxonomy = None;
            let mut wgs = None;
            if cols.len() > 1 {
                if let Some(idx) = wgs_idx {
                    wgs = not_none(column(idx, COL_WGS)?).map(str::to_owned);
                }
                if let Some(idx) = taxonomy_idx
                    && let Some(tax_string) = not_none(column(idx, COL_TAXONOMY)?)
                {
                    taxonomy = Some(Taxonomy::from_string(tax_string)?);
                }
            }

            rows.insert(
                accession.clone(),
                MetadataRow {
                    accession,
                    taxonomy,
                    wgs,
                },
            );
        }
        Ok(Metadata { rows })
    }

    /// Retain only the rows whose genome ID is in `ids`.
    pub(crate) fn subset_to_ids(&mut self, ids: &BTreeSet<GenomeId>) {
        self.rows.retain(|gid, _| ids.contains(gid));
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

/// `None` if the value is the `none` placeholder.
fn not_none(value: &str) -> Option<&str> {
    (value != NONE_VALUE).then_some(value)
}

#[cfg(test)]
mod tests;
