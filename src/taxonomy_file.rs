//! The taxonomy file given by `--taxonomy`: the full taxonomy of each genome, e.g. from gDecorate.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::error::{CuratomaticError, CuratomaticResult};
use crate::genome::{GenomeId, gid_to_canonical};
use crate::taxon::Taxon;
use crate::taxonomy::Taxonomy;

/// Separates the taxa of a taxonomy string, as read and written.
const TAXON_SEPARATOR: &str = ";";

/// The taxonomy of each genome, keyed by canonical genome ID.
///
/// Read from a tab separated file with no header and two columns, the genome ID and its taxonomy
/// string, e.g. `G003697165<TAB>d__Bacteria;p__Pseudomonadota;...;s__Escherichia coli`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct TaxonomyFile {
    rows: BTreeMap<GenomeId, Taxonomy>,
}

impl TaxonomyFile {
    /// Read a taxonomy file, see [`TaxonomyFile`].
    pub(crate) fn from_file(path: &Path) -> CuratomaticResult<Self> {
        let file = File::open(path)
            .map_err(|e| CuratomaticError::Io(format!("{}: {e}", path.display())))?;
        Self::from_reader(BufReader::new(file), path)
    }

    /// Read a taxonomy file from a reader. `path` is only used in error messages.
    ///
    /// Genome IDs are made canonical (see [`gid_to_canonical`]). Blank lines are skipped. A line
    /// without exactly two columns, an invalid taxonomy and a genome given twice are errors.
    /// Placeholders (e.g. `g__`) are not names: they are stripped, so each row holds named taxa
    /// only. A rank that holds both a placeholder and a name (`g__;g__X`) is an error.
    pub(crate) fn from_reader(reader: impl BufRead, path: &Path) -> CuratomaticResult<Self> {
        let error = |line: usize, message: String| CuratomaticError::TaxonomyFile {
            path: path.display().to_string(),
            line,
            message,
        };
        let mut rows = BTreeMap::new();
        for (line_idx, line) in reader.lines().enumerate() {
            let line_no = line_idx.saturating_add(1);
            let line =
                line.map_err(|e| CuratomaticError::Io(format!("{}: {e}", path.display())))?;
            let line = line.trim_end_matches(['\r', '\n']);
            if line.trim().is_empty() {
                continue;
            }

            let mut cols = line.split('\t');
            let (Some(gid), Some(taxonomy), None) = (cols.next(), cols.next(), cols.next()) else {
                return Err(error(
                    line_no,
                    "expected two tab separated columns".to_owned(),
                ));
            };
            let gid = gid_to_canonical(gid.trim());
            let taxonomy = parse_row(taxonomy).map_err(|message| error(line_no, message))?;
            if rows.insert(gid.clone(), taxonomy).is_some() {
                return Err(error(line_no, format!("{gid} is given more than once")));
            }
        }
        Ok(TaxonomyFile { rows })
    }

    /// The taxonomy of a genome.
    pub(crate) fn get(&self, gid: &GenomeId) -> Option<&Taxonomy> {
        self.rows.get(gid)
    }

    /// Every genome and its taxonomy, ordered by genome ID.
    pub(crate) fn iter(&self) -> impl Iterator<Item = (&GenomeId, &Taxonomy)> {
        self.rows.iter()
    }

    /// The genomes in `gids` that are not in the file, and the genomes in the file that are not in
    /// `gids`.
    pub(crate) fn compare(
        &self,
        gids: &BTreeSet<GenomeId>,
    ) -> (BTreeSet<GenomeId>, BTreeSet<GenomeId>) {
        let missing = gids
            .iter()
            .filter(|gid| !self.rows.contains_key(*gid))
            .cloned()
            .collect();
        let extra = self
            .rows
            .keys()
            .filter(|gid| !gids.contains(*gid))
            .cloned()
            .collect();
        (missing, extra)
    }
}

/// Parse the taxonomy string of a row, stripping its placeholders (see
/// [`TaxonomyFile::from_reader`]). The error is the message for the line.
fn parse_row(tax_string: &str) -> Result<Taxonomy, String> {
    let mut taxa = Vec::new();
    for taxon in tax_string.split(TAXON_SEPARATOR) {
        let taxon = Taxon::new(taxon.trim());
        let rank = taxon.rank().map_err(|e| format!("invalid taxonomy: {e}"))?;
        taxa.push((rank, taxon));
    }

    // Checked before stripping, which would hide the placeholder.
    for (rank, placeholder) in taxa.iter().filter(|(_, taxon)| !taxon.is_named()) {
        if let Some((_, name)) = taxa
            .iter()
            .find(|(other, taxon)| other == rank && taxon.is_named())
        {
            return Err(format!(
                "the {} rank holds both the placeholder {placeholder} and the name {name}",
                rank.name().to_lowercase()
            ));
        }
    }

    let mut out = Taxonomy::default();
    for (_, taxon) in taxa.into_iter().filter(|(_, taxon)| taxon.is_named()) {
        out.add_taxon(taxon)
            .map_err(|e| format!("invalid taxonomy: {e}"))?;
    }
    Ok(out)
}

/// Format rows in the format of the taxonomy file: the genome ID, a tab, and the taxa joined with
/// `;`, one row per line.
pub(crate) fn to_tsv<'a>(rows: impl IntoIterator<Item = (&'a GenomeId, &'a [Taxon])>) -> String {
    let mut out = String::new();
    for (gid, taxa) in rows {
        let taxa: Vec<&str> = taxa.iter().map(Taxon::as_str).collect();
        out.push_str(gid.as_str());
        out.push('\t');
        out.push_str(&taxa.join(TAXON_SEPARATOR));
        out.push('\n');
    }
    out
}

/// Write rows to `path`, see [`to_tsv`].
pub(crate) fn write<'a>(
    path: &Path,
    rows: impl IntoIterator<Item = (&'a GenomeId, &'a [Taxon])>,
) -> CuratomaticResult<()> {
    fs::write(path, to_tsv(rows))
        .map_err(|e| CuratomaticError::Io(format!("{}: {e}", path.display())))
}

#[cfg(test)]
mod tests;
