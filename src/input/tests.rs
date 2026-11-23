use std::fmt::Debug;
use std::fs;
use std::io::Write;
use std::path::PathBuf;

use flate2::Compression;
use flate2::write::GzEncoder;

use super::*;
use crate::metadata::{COL_ACCESSION, Metadata};
use crate::newick::NewickTree;
use crate::red::Red;
use crate::taxonomy_file::{PriorTaxonomyFile, TaxonomyFile};

/// A new, empty directory for the files of test `name`.
fn test_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("curatomatic-input-{name}-{}", std::process::id()));
    if dir.exists() {
        fs::remove_dir_all(&dir).unwrap();
    }
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// `bytes` as a single gzip member.
fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

/// Write `contents` to `dir/name` and gzipped to `dir/name.gz`, and return both paths.
fn write_plain_and_gzipped(dir: &Path, name: &str, contents: &str) -> (PathBuf, PathBuf) {
    let plain = dir.join(name);
    fs::write(&plain, contents).unwrap();
    let gzipped = dir.join(format!("{name}.gz"));
    fs::write(&gzipped, gzip(contents.as_bytes())).unwrap();
    (plain, gzipped)
}

/// Read `contents` plain and gzipped with `read`, check both give the same, and return it.
fn read_both<T: Debug + PartialEq>(
    test: &str,
    name: &str,
    contents: &str,
    read: impl Fn(&Path) -> CuratomaticResult<T>,
) -> T {
    let dir = test_dir(test);
    let (plain, gzipped) = write_plain_and_gzipped(&dir, name, contents);
    let from_plain = read(&plain).unwrap();
    let from_gzipped = read(&gzipped).unwrap();
    fs::remove_dir_all(&dir).unwrap();
    assert_eq!(from_plain, from_gzipped);
    from_plain
}

const NEWICK: &str = "((G1:0.1,G2:0.1)'100:g__A|red=0.95':0.1,G3:0.1)'d__D;p__P|red=0.0000';\n";

const RED: &str = "rank\tmedian\nPhylum\t0.3\nClass\t0.4\nOrder\t0.5\nFamily\t0.7\nGenus\t0.9\n";

const TAXONOMY: &str = "G1\td__D;p__P;c__C;o__O;f__F;g__A;s__A one\n\
                        G2\td__D;p__P;c__C;o__O;f__F;g__A;g__N;s__N two\n";

const PRIOR: &str = "G1\td__D;p__P;c__C;o__O;f__F;g__A;s__A one\n\
                     G2\td__D;p__P;c__C;o__O;f__F;g__N;s__N two\n";

const METADATA: &str = "accession\tncbi_wgs_formatted\tother\nG1\tWGS1\tx\nG2\tnone\ty\n";

// The rule: an input whose name ends in `.gz` is read decompressed, and gives exactly what the
// same content read plain gives.

#[test]
fn a_gzipped_tree_reads_like_the_plain_tree() {
    let tree = read_both("tree", "x.nwk", NEWICK, NewickTree::from_file);
    assert_eq!(tree, NewickTree::parse(NEWICK).unwrap());
}

#[test]
fn a_gzipped_median_red_file_reads_like_the_plain_file() {
    let red = read_both("red", "red.tsv", RED, Red::from_file);
    assert_eq!(red, Red::from_tsv(RED).unwrap());
}

#[test]
fn a_gzipped_taxonomy_file_reads_like_the_plain_file() {
    let taxonomy = read_both(
        "taxonomy",
        "taxonomy.tsv",
        TAXONOMY,
        TaxonomyFile::from_file,
    );
    assert_eq!(taxonomy.iter().count(), 2);
}

#[test]
fn a_gzipped_prior_taxonomy_file_reads_like_the_plain_file() {
    read_both("prior", "prior.tsv", PRIOR, PriorTaxonomyFile::from_file);
}

#[test]
fn a_gzipped_novel_name_file_reads_like_the_plain_file() {
    let metadata = read_both("metadata", "metadata.tsv", METADATA, |path| {
        Metadata::from_file(path, COL_ACCESSION, "ncbi_wgs_formatted")
    });
    assert_eq!(metadata.rows.len(), 2);
}

#[test]
fn any_name_ending_in_gz_is_gzipped() {
    // Not only `.tsv.gz`: the extension before `.gz` doesn't matter.
    let dir = test_dir("any-gz");
    for name in ["x.gz", "x.tree.gz", "x.txt.gz"] {
        let path = dir.join(name);
        fs::write(&path, gzip(b"text\n")).unwrap();
        assert_eq!(read_input_to_string(&path).unwrap(), "text\n", "{name}");
    }
    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn gzip_is_decided_by_the_name_never_by_the_content() {
    // Gzipped content in a file whose name doesn't end in `.gz` is read as it is: here, not text.
    let dir = test_dir("by-name");
    let path = dir.join("taxonomy.tsv");
    fs::write(&path, gzip(TAXONOMY.as_bytes())).unwrap();
    let result = TaxonomyFile::from_file(&path);
    fs::remove_dir_all(&dir).unwrap();
    assert!(matches!(result, Err(CuratomaticError::Io(_))), "{result:?}");

    // And `.GZ` or `gz` without the dot is not `.gz`.
    assert!(!is_gzipped(Path::new("x.GZ")));
    assert!(!is_gzipped(Path::new("xgz")));
    assert!(is_gzipped(Path::new("dir.gz/x.gz")));
    assert!(!is_gzipped(Path::new("dir.gz/x")));
}

#[test]
fn every_member_of_a_multi_member_gzip_is_read() {
    // `GzDecoder` would stop after the first member and silently drop the rest.
    let dir = test_dir("multi-member");
    let path = dir.join("taxonomy.tsv.gz");
    let (first, second) = TAXONOMY.split_at(TAXONOMY.find('\n').unwrap() + 1);
    let mut bytes = gzip(first.as_bytes());
    bytes.extend(gzip(second.as_bytes()));
    fs::write(&path, bytes).unwrap();
    let text = read_input_to_string(&path).unwrap();
    let taxonomy = TaxonomyFile::from_file(&path).unwrap();
    fs::remove_dir_all(&dir).unwrap();
    assert_eq!(text, TAXONOMY);
    assert_eq!(taxonomy.iter().count(), 2);
}

#[test]
fn a_bgzip_file_with_its_empty_end_member_is_read_in_full() {
    // bgzip ends a file with an empty member.
    let dir = test_dir("bgzip");
    let path = dir.join("x.nwk.gz");
    let mut bytes = gzip(NEWICK.as_bytes());
    bytes.extend(gzip(b""));
    fs::write(&path, bytes).unwrap();
    let text = read_input_to_string(&path).unwrap();
    fs::remove_dir_all(&dir).unwrap();
    assert_eq!(text, NEWICK);
}

/// The message of the [`CuratomaticError::Io`] of `result`.
fn io_message<T: Debug>(result: CuratomaticResult<T>) -> String {
    match result {
        Err(CuratomaticError::Io(message)) => message,
        other => panic!("expected an I/O error, got {other:?}"),
    }
}

#[test]
fn a_gz_file_that_is_not_gzip_is_an_io_error_naming_the_path() {
    let dir = test_dir("not-gzip");
    let path = dir.join("taxonomy.tsv.gz");
    fs::write(&path, TAXONOMY).unwrap();
    let path_text = path.display().to_string();
    let messages = [
        io_message(TaxonomyFile::from_file(&path)),
        io_message(PriorTaxonomyFile::from_file(&path)),
        io_message(Metadata::from_file(&path, COL_ACCESSION, "x")),
        io_message(Red::from_file(&path)),
        io_message(NewickTree::from_file(&path)),
    ];
    fs::remove_dir_all(&dir).unwrap();
    for message in messages {
        assert!(message.contains(&path_text), "{message}");
    }
}

#[test]
fn a_gzip_file_cut_short_is_an_io_error_naming_the_path() {
    // Never a partly read input.
    let dir = test_dir("cut-short");
    let path = dir.join("taxonomy.tsv.gz");
    let bytes = gzip(TAXONOMY.as_bytes());
    fs::write(&path, &bytes[..bytes.len() - 4]).unwrap();
    let message = io_message(TaxonomyFile::from_file(&path));
    let text = io_message(read_input_to_string(&path));
    fs::remove_dir_all(&dir).unwrap();
    assert!(message.contains(&path.display().to_string()), "{message}");
    assert!(text.contains(&path.display().to_string()), "{text}");
}

#[test]
fn an_empty_gz_file_is_an_io_error_naming_the_path() {
    let dir = test_dir("empty");
    let path = dir.join("taxonomy.tsv.gz");
    fs::write(&path, "").unwrap();
    let message = io_message(TaxonomyFile::from_file(&path));
    fs::remove_dir_all(&dir).unwrap();
    assert!(message.contains(&path.display().to_string()), "{message}");
}

#[test]
fn a_missing_file_is_an_io_error_naming_the_path() {
    for name in ["missing.tsv", "missing.tsv.gz"] {
        let path = Path::new("does/not/exist").join(name);
        let message = io_message(open_input(&path).map(|_| ()));
        assert!(message.contains(&path.display().to_string()), "{message}");
    }
}

#[test]
fn an_error_line_in_a_gzipped_taxonomy_file_is_the_line_of_the_plain_file() {
    let contents = format!("{TAXONOMY}\nG3\td__D;p__P;g__X;g__X;s__X three\n");
    let dir = test_dir("error-line");
    let (plain, gzipped) = write_plain_and_gzipped(&dir, "taxonomy.tsv", &contents);
    let line = |path: &Path| match TaxonomyFile::from_file(path) {
        Err(CuratomaticError::TaxonomyFile { line, .. }) => line,
        other => panic!("expected a taxonomy file error, got {other:?}"),
    };
    let (line_plain, line_gzipped) = (line(&plain), line(&gzipped));
    fs::remove_dir_all(&dir).unwrap();
    assert_eq!(line_plain, 4);
    assert_eq!(line_gzipped, line_plain);
}

#[test]
fn the_stem_drops_gz_and_then_the_last_extension() {
    let cases = [
        ("x.nwk.gz", "x"),
        ("x.nwk", "x"),
        ("x.tree.gz", "x"),
        ("dir/x.tree.gz", "x"),
        ("x.gz", "x"),
        ("x", "x"),
        ("x.y.nwk.gz", "x.y"),
        ("x.nwk.GZ", "x.nwk"),
    ];
    for (path, stem) in cases {
        assert_eq!(input_stem(Path::new(path)), stem, "{path}");
    }
}
