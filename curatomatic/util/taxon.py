from curatomatic.model.metadata import Metadata
from curatomatic.model.rank import RANKS

BAD_PAIRS = [("I", "l"), ("0", "O"), ("1", "I"), ("5", "S"), ("2", "Z")]


def score_candidate_taxon(taxon):
    for a, b in BAD_PAIRS:
        if a in taxon and b in taxon:
            return 1
    return 0


def generate_unique_rank(rank, leaf_node_gids, metadata: Metadata, known_taxa):
    wgs_available = set()
    for gid in leaf_node_gids:
        wgs_row = metadata.rows.get(gid)
        if wgs_row is not None:
            wgs = wgs_row.wgs
            if wgs is not None:
                wgs_rank = f'{rank}__{wgs}'
                if wgs_rank not in known_taxa:
                    wgs_available.add(wgs_rank)

    if len(wgs_available) == 0:
        return f'{rank}__NOVEL-{sorted(leaf_node_gids)[0]}'

    else:
        # Otherwise, score the WGS sequences by how readable they are
        taxa_sorted = sorted(wgs_available, key=lambda x: (score_candidate_taxon(x), x))
        return taxa_sorted[0]


def add_taxon_to_string(tax_string, taxon):
    if tax_string is None:
        return taxon
    tax_string_split = tax_string.split('; ')
    tax_string_split.append(taxon)
    tax_string_sorted = sorted(tax_string_split, key=lambda x: RANKS.index(x[0]))
    return '; '.join(tax_string_sorted)


def create_node_label(red, bs, taxon):
    if bs is not None:
        if taxon:
            return f'{bs}:{taxon}|RED={red}'
        else:
            return f'{bs}|RED={red}'
    else:
        if taxon:
            return f'{taxon}|RED={red}'
        else:
            # This will cause problems if no taxon is in there, or no BS values
            if red >= 1.0:
                return f'100|RED={red}'
