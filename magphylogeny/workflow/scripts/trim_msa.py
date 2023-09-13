import re
import sys
from collections import defaultdict


def get_hits(path):
    out = list()
    with open(path) as f:
        hits = re.findall(r'>(.+)\n(.+)', f.read())
        for gid, seq in hits:
            out.append((gid, seq.strip()))
    return out


def main():
    input_path, output_path = sys.argv[1], sys.argv[2]

    hits = get_hits(input_path)

    col_values = defaultdict(set)
    only_gaps = set()
    for gid, seq in hits:
        if set(seq) == {'-'}:
            only_gaps.add(gid)
        else:
            for col_idx, aa in enumerate(seq):
                col_values[col_idx].add(aa)

    cols_to_exclude = set()
    for col_idx, values in col_values.items():
        if len(values) == 1:
            cols_to_exclude.add(col_idx)

    # There were columns that contained only gaps, don't use them
    if len(cols_to_exclude) > 0:
        with open(output_path, 'w') as f:
            for gid, seq in hits:
                if gid not in only_gaps:
                    new_seq = ''.join([aa for idx, aa in enumerate(seq) if idx not in cols_to_exclude])
                    if len(new_seq) > 0:
                        f.write(f'>{gid}\n{new_seq}\n')

    # Only exclude those that contained gaps since all columns have information
    else:
        with open(output_path, 'w') as f:
            for gid, seq in hits:
                if gid not in only_gaps:
                    f.write(f'>{gid}\n{seq}\n')


if __name__ == '__main__':
    main()
