import sys
from pathlib import Path

import dendropy
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

sys.path.insert(0, Path(__file__).absolute().parent.parent.parent.parent.as_posix())

from magphylogeny.methods.phylorank_red_html import get_taxon_red_values
from magphylogeny.util.red import node_to_red, get_prop_lineages_at_threshold
from magphylogeny.util.tree import calculate_node_depths, d_node_to_leaves, calc_node_to_descendants


def generate_plot(d_node_to_red, d_node_to_descs, d_rank_to_percentiles, path_out, steps, d_rank_to_qty):
    x = list()
    y_mix = list()
    y_mfd = list()
    y_gtdb = list()

    red_colors = 'k'
    mfd_colour = '#1f77b4'
    mix_colour = '#cee9fb'
    gtdb_colour = '#63afe4'

    number_of_bins = steps
    n_taxa = list()
    print('Calculating proportions at RED thresholds')
    for thresh in tqdm(np.linspace(0, 1, number_of_bins)):
        props = get_prop_lineages_at_threshold(thresh, d_node_to_red, d_node_to_descs)
        x.append(thresh)
        y_mix.append(props['MIX'])
        y_mfd.append(props['MFD'])
        y_gtdb.append(props['GTDB'])
        n_taxa.append(props['TOTAL'])

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.set_xlim((0, 1.0))
    ax.set_ylim((0, 100))

    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

    ax.stackplot(
        x, y_mfd, y_gtdb, y_mix,
        labels=['MFD', 'GTDB', 'Mix'],
        colors=[mfd_colour, gtdb_colour, mix_colour]
    )

    # Top (2nd x axis) legend
    ax_x2 = ax.twiny()
    ax.get_shared_x_axes().join(ax, ax_x2)

    x2_labels = [get_prop_lineages_at_threshold(round(tick, 1), d_node_to_red, d_node_to_descs)['TOTAL'] for tick in
                 ax.get_xticks()]

    ax_x2.set_xticks(ax.get_xticks())
    ax_x2.set_xbound(ax.get_xbound())
    ax_x2.set_xlabel('Number of lineages')
    ax_x2.set_xticklabels(x2_labels)

    # Add the secondary y axis
    ax_y2 = ax.twinx()
    ax.get_shared_y_axes().join(ax, ax_y2)
    ax_y2.set_ybound(ax.get_ybound())

    phy_y = 10
    class_y = 30
    order_y = 50
    family_y = 70
    genus_y = 90
    y_delta = 2

    d_rank_to_axis_pos = {
        'Phylum': 10,
        'Class': 30,
        'Order': 50,
        'Family': 70,
        'Genus': 90,
    }

    ax_y2.set_yticks((phy_y, class_y, order_y, family_y, genus_y))
    ax_y2.grid(False)

    new_tick_labels = list()
    for rank in d_rank_to_axis_pos:
        cur_qty = d_rank_to_qty[rank]
        new_tick_labels.append(f'{rank} ({cur_qty:,})')

    ax_y2.set_yticklabels(new_tick_labels)

    for cur_rank, cur_qty in d_rank_to_qty.items():
        cur_key = cur_rank[0].lower()
        axis_pos = d_rank_to_axis_pos[cur_rank]
        cur_percentile = d_rank_to_percentiles[cur_key]
        ax.hlines(axis_pos, cur_percentile[0], cur_percentile[2], colors=red_colors)
        ax.vlines(cur_percentile[1], axis_pos - y_delta, axis_pos + y_delta, colors=red_colors)

    fig.legend(loc='upper right')
    ax.set_xlabel('Relative Evolutionary Divergence (RED)')
    ax.set_ylabel('Percentage of lineages (%)')
    plt.title('Percentage of lineages comprised entirely of a given source at increasing RED thresholds')
    ax.grid(True)
    fig.tight_layout()
    fig.subplots_adjust(top=0.80)

    plt.savefig(path_out)


def main():
    # path_tree = '/Users/aaron/obsidian/phd/PhD/Lab Stays/Aalborg/data/phylorank_denovo_arc_trusted/gtdbtk.ar53.decorated.scaled.tree'
    # path_red = '/Users/aaron/obsidian/phd/PhD/Lab Stays/Aalborg/data/phylorank_denovo_arc_trusted/gtdbtk.ar53.decorated.node_rd.tsv'
    # path_red_taxa = '/Users/aaron/obsidian/phd/PhD/Lab Stays/Aalborg/data/phylorank_denovo_arc_trusted/gtdbtk.ar53.decorated.html'
    # path_out = '/tmp/output.png'
    # steps = 10000

    path_tree = sys.argv[1]
    path_red = sys.argv[2]
    path_red_taxa = sys.argv[3]
    path_out = sys.argv[4]
    steps = int(sys.argv[5])

    tree = dendropy.Tree.get(path=path_tree, schema='newick', preserve_underscores=True)

    d_depth_to_node, d_node_to_depth = calculate_node_depths(tree)
    d_node_to_descs, d_leaf_taxon_to_node = d_node_to_leaves(d_depth_to_node)
    d_node_to_red = node_to_red(d_node_to_descs, d_node_to_depth, d_leaf_taxon_to_node, path_red)
    d_node_to_descs = calc_node_to_descendants(d_node_to_red)
    d_rank_to_percentiles, d_rank_to_qty = get_taxon_red_values(path_red_taxa)

    generate_plot(d_node_to_red, d_node_to_descs, d_rank_to_percentiles, path_out, steps, d_rank_to_qty)

    return


if __name__ == '__main__':
    main()
