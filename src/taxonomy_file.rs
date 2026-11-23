//! The taxonomy files: the full taxonomy of each genome, given by `--taxonomy` (e.g. from
//! gDecorate), and the taxonomy of a previous release, given by `--prior-taxonomy`.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::Path;

use strum::IntoEnumIterator;

use crate::error::{CuratomaticError, CuratomaticResult};
use crate::genome::{GenomeId, gid_to_canonical};
use crate::rank::Rank;
use crate::taxon::Taxon;
use crate::taxonomy::Taxonomy;

/// Separates the taxa of a taxonomy string, as read and written.
const TAXON_SEPARATOR: &str = ";";

/// Which rules a row must follow, beyond the checks shared by every taxonomy file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RowFormat {
    /// `--taxonomy`: placeholders are stripped, and a rank can be missing or hold several taxa.
    Decorated,
    /// `--prior-taxonomy`: exactly one named taxon per rank, see [`check_prior_row`].
    Prior,
}

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
    /// only. A rank that holds both a placeholder and a name (`g__;g__X`) is an error, as is a row
    /// whose taxa are not in rank order, that gives a named taxon twice, or that has more than one
    /// named domain or species.
    pub(crate) fn from_reader(reader: impl BufRead, path: &Path) -> CuratomaticResult<Self> {
        Self::read_rows(reader, path, RowFormat::Decorated)
    }

    /// Read the rows of a taxonomy file, each following the rules of `format`.
    fn read_rows(reader: impl BufRead, path: &Path, format: RowFormat) -> CuratomaticResult<Self> {
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
            let taxonomy =
                parse_row(taxonomy, format).map_err(|message| error(line_no, message))?;
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

/// The prior taxonomy given by `--prior-taxonomy`: the taxonomy of each genome in a previous
/// release, keyed by canonical genome ID.
///
/// It has the format of [`TaxonomyFile`], but a previous release classifies every genome fully, so
/// each row must be exactly one named taxon per rank, see [`PriorTaxonomyFile::from_reader`]. A
/// newtype, so that the prior can't be read with the lenient reader of `--taxonomy`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct PriorTaxonomyFile(TaxonomyFile);

impl PriorTaxonomyFile {
    /// Read a prior taxonomy file, see [`PriorTaxonomyFile::from_reader`].
    pub(crate) fn from_file(path: &Path) -> CuratomaticResult<Self> {
        let file = File::open(path)
            .map_err(|e| CuratomaticError::Io(format!("{}: {e}", path.display())))?;
        Self::from_reader(BufReader::new(file), path)
    }

    /// Read a prior taxonomy file from a reader. `path` is only used in error messages.
    ///
    /// Every row is checked as in [`TaxonomyFile::from_reader`], and must also be exactly seven
    /// taxa, one per rank in rank order, each with a name, and a species that is a binomial name
    /// in the row's genus (`g__Escherichia;s__Escherichia coli`). A placeholder, a missing rank,
    /// two taxa of a rank, an unresolved species (`s__{unresolved} NAME`) and a species of another
    /// genus are errors. This applies to every row, including genomes that are not on the tree.
    pub(crate) fn from_reader(reader: impl BufRead, path: &Path) -> CuratomaticResult<Self> {
        TaxonomyFile::read_rows(reader, path, RowFormat::Prior).map(PriorTaxonomyFile)
    }

    /// The taxonomy of a genome.
    #[cfg(test)]
    pub(crate) fn get(&self, gid: &GenomeId) -> Option<&Taxonomy> {
        self.0.get(gid)
    }

    /// Every genome and its taxonomy, ordered by genome ID.
    pub(crate) fn iter(&self) -> impl Iterator<Item = (&GenomeId, &Taxonomy)> {
        self.0.iter()
    }
}

/// Check that a row of a prior taxonomy, as read (placeholders included), is exactly one named
/// taxon per rank, and that its species is a binomial name in its genus (see
/// [`PriorTaxonomyFile::from_reader`]). The rank order is checked before. The error is the message
/// for the line.
fn check_prior_row(taxa: &[(Rank, Taxon)]) -> Result<(), String> {
    for rank in Rank::iter() {
        let name = rank.name().to_lowercase();
        let of_rank: Vec<&Taxon> = taxa
            .iter()
            .filter(|(other, _)| *other == rank)
            .map(|(_, taxon)| taxon)
            .collect();
        if of_rank.is_empty() {
            return Err(format!(
                "the {name} rank is missing, but a prior taxonomy names every rank"
            ));
        }
        if let Some(placeholder) = of_rank.iter().find(|taxon| !taxon.is_named()) {
            return Err(format!(
                "the {name} rank is the placeholder {placeholder}, but a prior taxonomy names \
                 every rank"
            ));
        }
        if of_rank.len() > 1 {
            let joined: Vec<&str> = of_rank.iter().map(|taxon| taxon.as_str()).collect();
            return Err(format!(
                "the {name} rank has {}, but a prior taxonomy has one taxon per rank",
                joined.join(" and ")
            ));
        }
    }

    // Exactly one named taxon per rank, so the genus and species are the last two.
    let [.., (_, genus), (_, species)] = taxa else {
        return Err("the row has fewer than two taxa".to_owned());
    };
    let is_binomial = !species.is_unresolved_species()
        && species
            .as_str()
            .strip_prefix("s__")
            .and_then(|name| name.split_once(' '))
            .is_some_and(|(first, second)| {
                [first, second]
                    .iter()
                    .all(|word| !word.is_empty() && !word.contains(char::is_whitespace))
            });
    if !is_binomial {
        return Err(format!(
            "the species {species} is not a binomial name, but a prior taxonomy names every species"
        ));
    }
    if species.genus_of_species().as_ref() != Some(genus) {
        return Err(format!("the species {species} is not in the genus {genus}"));
    }
    Ok(())
}

/// Parse the taxonomy string of a row, stripping its placeholders (see
/// [`TaxonomyFile::from_reader`]), and check it follows the rules of `format`. The error is the
/// message for the line.
fn parse_row(tax_string: &str, format: RowFormat) -> Result<Taxonomy, String> {
    let mut taxa = Vec::new();
    for taxon in tax_string.split(TAXON_SEPARATOR) {
        let taxon = Taxon::new(taxon.trim());
        let rank = taxon.rank().map_err(|e| format!("invalid taxonomy: {e}"))?;
        taxa.push((rank, taxon));
    }

    // The row must be in rank order (placeholders included), checked before anything is kept.
    if let Some(((_, before), (_, after))) = taxa
        .iter()
        .zip(taxa.iter().skip(1))
        .find(|((before, _), (after, _))| before > after)
    {
        return Err(format!(
            "the taxa are not in rank order: {after} is after {before}"
        ));
    }

    let mut seen = BTreeSet::new();
    if let Some((_, taxon)) = taxa
        .iter()
        .find(|(_, taxon)| taxon.is_named() && !seen.insert(taxon))
    {
        return Err(format!("{taxon} is given more than once"));
    }

    for rank in [Rank::Domain, Rank::Species] {
        let named: Vec<&str> = taxa
            .iter()
            .filter(|(other, taxon)| *other == rank && taxon.is_named())
            .map(|(_, taxon)| taxon.as_str())
            .collect();
        if named.len() > 1 {
            return Err(format!(
                "a genome has one {} only, but the row has {}",
                rank.name().to_lowercase(),
                named.join(" and ")
            ));
        }
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

    if format == RowFormat::Prior {
        check_prior_row(&taxa)?;
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
