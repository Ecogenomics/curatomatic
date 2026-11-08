use super::*;

#[test]
fn strips_prefix_and_version() {
    assert_eq!(gid_to_canonical("GCA_947480885.1").as_str(), "G947480885");
    assert_eq!(gid_to_canonical("GCF_947480885.2").as_str(), "G947480885");
    assert_eq!(gid_to_canonical("GB_GCA_947480885.1").as_str(), "G947480885");
    assert_eq!(gid_to_canonical("RS_GCF_947480885.1").as_str(), "G947480885");
}

#[test]
fn leaves_other_ids_unchanged() {
    assert_eq!(gid_to_canonical("G947480885").as_str(), "G947480885");
    assert_eq!(gid_to_canonical("GB_GCF_947480885.1").as_str(), "GB_GCF_947480885.1");
    assert_eq!(gid_to_canonical("").as_str(), "");
}

#[test]
fn short_ids_after_prefix() {
    assert_eq!(gid_to_canonical("GCA_").as_str(), "G");
    assert_eq!(gid_to_canonical("GCA_1").as_str(), "G");
    assert_eq!(gid_to_canonical("GCA_12").as_str(), "G");
    assert_eq!(gid_to_canonical("GCA_123").as_str(), "G1");
}

// The expected values below come from the Python implementation
// (`curatomatic/util/genome.py`, `gid_to_canonical`), which slices by character.
#[test]
fn only_the_listed_prefixes_are_stripped() {
    assert_eq!(gid_to_canonical("RS_GCA_1.1").as_str(), "RS_GCA_1.1");
    assert_eq!(gid_to_canonical("GB_GCF_1.1").as_str(), "GB_GCF_1.1");
    assert_eq!(gid_to_canonical("gca_1.1").as_str(), "gca_1.1");
}

#[test]
fn prefix_only_ids() {
    assert_eq!(gid_to_canonical("GB_GCA_").as_str(), "G");
    assert_eq!(gid_to_canonical("RS_GCF_1").as_str(), "G");
}

#[test]
fn non_ascii_ids_are_sliced_by_character() {
    assert_eq!(gid_to_canonical("GCA_éé1").as_str(), "Gé");
    assert_eq!(gid_to_canonical("GCA_ab.é").as_str(), "Gab");
    assert_eq!(gid_to_canonical("GCA_é").as_str(), "G");
    assert_eq!(gid_to_canonical("GCA_ééé").as_str(), "Gé");
    assert_eq!(gid_to_canonical("GB_GCA_ééé.1").as_str(), "Gééé");
}

#[test]
fn display_and_as_str_give_the_id() {
    let gid = GenomeId::new("G002502135");
    assert_eq!(gid.as_str(), "G002502135");
    assert_eq!(gid.to_string(), "G002502135");
}

#[test]
fn zombies_start_with_d_dash() {
    assert!(GenomeId::new("D-G002502135").is_zombie());
    assert!(GenomeId::new("D-").is_zombie());
    assert!(!GenomeId::new("G002502135").is_zombie());
    assert!(!GenomeId::new("d-G1").is_zombie());
    assert!(!GenomeId::new("XD-G1").is_zombie());
    assert!(!GenomeId::new("").is_zombie());
}

// Python sorts genome IDs as plain strings (`util/taxon.py:26`).
#[test]
fn order_is_plain_string_order() {
    let mut gids: Vec<GenomeId> = ["G9", "G10", "D-G1", "G1", "g1", ""]
        .into_iter()
        .map(GenomeId::new)
        .collect();
    gids.sort();
    let sorted: Vec<&str> = gids.iter().map(GenomeId::as_str).collect();
    assert_eq!(sorted, ["", "D-G1", "G1", "G10", "G9", "g1"]);
}
