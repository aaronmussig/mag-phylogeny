def read_tax_file(path):
    out = dict()
    with open(path) as f:
        for line in f.readlines():
            gid, tax = line.strip().split('\t')
            out[gid] = tax
    return out
