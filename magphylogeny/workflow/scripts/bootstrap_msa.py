import gzip
import os
import random
import re
import sys


def bootstrap_alignment(msa, output_file, frac=1.0):
    alignment_len = len(msa[list(msa.keys())[0]])
    sample_len = int(alignment_len * frac)
    cols = [random.randint(0, alignment_len - 1) for _ in range(sample_len)]

    fout = open(output_file, 'w')
    for seq_id, seq in msa.items():
        fout.write('>' + seq_id + '\n')
        for col in cols:
            fout.write(seq[col])
        fout.write('\n')
    fout.close()


def read_msa_from_content(f):
    out = dict()
    hits = re.findall(r'>(.+)\n(.+)', f)
    for gid, seq in hits:
        out[gid] = seq.strip()
    return out


def read_msa(path):
    if path.endswith('gz'):
        with gzip.open(path) as f:
            return read_msa_from_content(f.read().decode('utf-8'))
    with open(path) as f:
        return read_msa_from_content(f.read())


def main():
    input_path, output_path = sys.argv[1], sys.argv[2]

    output_directory = os.path.dirname(output_path)
    os.makedirs(output_directory, exist_ok=True)

    d_gid_to_seq = read_msa(input_path)
    bootstrap_alignment(d_gid_to_seq, output_path)


if __name__ == '__main__':
    main()
