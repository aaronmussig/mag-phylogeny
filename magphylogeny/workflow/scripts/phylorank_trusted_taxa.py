import sys
from collections import defaultdict

RANKS = ('d', 'p', 'c', 'o', 'f', 'g', 's')


def get_num_ranks(tax):
    taxa = tax.split(';')
    seen = defaultdict(lambda: 0)
    for taxon in taxa:
        taxon = taxon.strip()
        seen[taxon[0]] += 1
    return seen


def main():
    input_path, output_path = sys.argv[1], sys.argv[2]

    trusted_taxa = set()
    untrusted = list()
    total_taxa = 0

    phyla = set()

    with open(input_path) as f:
        for line in f.readlines():
            gid, tax = line.strip().split('\t')
            total_taxa += 1

            phyla.add(tax.split(';')[1])

            n_seen = get_num_ranks(tax)

            n_seen_total = sum([n_seen[x] for x in RANKS])
            if n_seen_total != 7:
                untrusted.append((gid, tax))
            else:
                for taxon in tax.split(';'):
                    taxon = taxon.strip()
                    if len(taxon) > 3:
                        trusted_taxa.add(taxon)

    with open(output_path, 'w') as f:
        for gid in sorted(trusted_taxa):
            f.write(f'{gid}\n')

    return


if __name__ == '__main__':
    main()
