def read_marker_file(path):
    out = dict()
    with open(path) as f:
        header = {x: i for i, x in enumerate(f.readline().strip().split('\t'))}
        print(header)
        gid_idx = header['name']
        unq_idx = header['number_unique_genes']
        mul_idx = header['number_multiple_genes']
        muq_idx = header['number_multiple_unique_genes']
        mis_idx = header['number_missing_genes']
        for line in f.readlines():
            cols = line.strip().split('\t')
            out[cols[gid_idx]] = {
                'unq': int(cols[unq_idx]),
                'mul': int(cols[mul_idx]),
                'muq': int(cols[muq_idx]),
                'mis': int(cols[mis_idx]),
            }
    return out
