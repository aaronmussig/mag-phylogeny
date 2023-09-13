import sys
from pathlib import Path

import matplotlib.pyplot as plt

sys.path.insert(0, Path(__file__).absolute().parent.parent.parent.parent.as_posix())

from magphylogeny.methods.marker_summary_file import read_marker_file
from magphylogeny.methods.taxonomy_file import read_tax_file


def gen_plot(d_gid_to_marker, d_gid_to_tax, path_plot):
    gids = list()
    unq_markers = list()
    mul_markers = list()
    muq_markers = list()
    mis_markers = list()
    total_markers = list()
    for gid, tax in d_gid_to_tax.items():
        if gid.startswith('MFD'):
            marker = d_gid_to_marker[gid]
            gids.append(gid)
            unq_markers.append(marker['unq'])
            mul_markers.append(marker['mul'])
            muq_markers.append(marker['muq'])
            mis_markers.append(marker['mis'])
            total_markers.append(marker['unq'] + marker['muq'])

    fig, ax = plt.subplots()
    weights = [1 / len(total_markers) * 100 for _ in total_markers]
    ax.hist(total_markers, bins=30, weights=weights)
    ax.set_ylabel('Percentage of Samples (%)')
    ax.set_xlabel('Number of GTDB Markers')
    fig.suptitle(f'Distribution of marker genes used in tree inference')
    plt.savefig(path_plot)
    return


def main():
    path_markers, path_taxonomy, path_plot = sys.argv[1], sys.argv[2], sys.argv[3]

    d_gid_to_marker = read_marker_file(path_markers)
    d_gid_to_tax = read_tax_file(path_taxonomy)

    gen_plot(d_gid_to_marker, d_gid_to_tax, path_plot)

    return


if __name__ == '__main__':
    main()
