
PATH_METADATA = '/Users/aaron/gtdb/curatomatic/test_data/r237/gtdb_r237_bac120_metadata.tsv'
PATH_TAXONOMY = '/Users/aaron/gtdb/curatomatic/test_data/r237/final_taxonomy.tsv'
PATH_OUT = '/Users/aaron/gtdb/curatomatic/test_data/r237/curatomatic_input.tsv'

def read_metadata(path):
    out = dict()
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            cols = line.strip().split('\t')
            wgs = cols[header.index('ncbi_wgs_formatted')]
            gid = cols[header.index('db_name')]
            assert gid not in out
            out[gid] = wgs
    return out

def read_taxonomy(path):
    out = dict()
    with open(path) as f:
        for line in f:
            cols = line.strip().split('\t')
            gid = cols[0]
            taxonomy = cols[1]
            assert gid not in out
            out[gid] = taxonomy
    return out


def create_output(metadata, taxonomy, out):
    all_gids = set(metadata.keys()).union(set(taxonomy.keys()))
    with open(out, 'w') as f:
        f.write('accession\ttaxonomy\tncbi_wgs_formatted\n')
        for gid in sorted(all_gids):
            wgs = metadata[gid]
            tax = taxonomy[gid]
            f.write(f'{gid}\t{tax}\t{wgs}\n')
    print(f'Wrote to: {out}')

def main():
    metadata = read_metadata(PATH_METADATA)
    taxonomy = read_taxonomy(PATH_TAXONOMY)
    create_output(metadata, taxonomy, PATH_OUT)
    return

if __name__ == "__main__":
    main()
