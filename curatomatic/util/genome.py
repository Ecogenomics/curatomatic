def gid_to_canonical(gid: str) -> str:
    """
    Converts a genome ID to canonical format, e.g.
    GCA_947480885.1 -> G947480885
    """
    if gid.startswith('GB_GCA_') or gid.startswith('RS_GCF_'):
        return f'G{gid[7:-2]}'
    elif gid.startswith('GCA_') or gid.startswith('GCF_'):
        return f'G{gid[4:-2]}'
    return gid
