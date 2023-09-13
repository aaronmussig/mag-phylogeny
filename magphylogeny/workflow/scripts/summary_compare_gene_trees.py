import os
import sys
from collections import Counter
from pathlib import Path

import dendropy
import numpy as np
import seaborn as sns
from tqdm import tqdm

sys.path.insert(0, Path(__file__).absolute().parent.parent.parent.parent.as_posix())

from magphylogeny.methods.taxonomy_file import read_tax_file
from magphylogeny.util.config import read_config_file
from magphylogeny.util.pdm import get_pdm
from magphylogeny.util.stats import weighted_jaccard
from magphylogeny.util.tree import RANKS


def get_marker_paths(marker_dir, domain):
    config = read_config_file()
    markers = set(config['gtdbtk']['markers'][domain])
    out = dict()
    for file in os.listdir(marker_dir):
        if file.endswith('.tree'):
            # This is horribly inefficient but it's compatible with my dev environment
            for marker in markers:
                if marker in file:
                    out[marker] = os.path.join(marker_dir, file)
                    break
    return out


def process_marker(target_gid, dm, labels, n, d_gtdb_tax, target_taxon):
    # Get the index of the taxon
    try:
        label_idx = labels.index(target_gid)
    except ValueError:
        return dict()

    rank_idx = RANKS.index(target_taxon)

    label_row = dm[label_idx]
    label_row_argsort = np.argsort(label_row)

    taxa_of_closest = list()
    for cur_closest_idx in label_row_argsort[1:]:
        cur_gid = labels[cur_closest_idx]
        if cur_gid.startswith('MFD'):
            continue
        taxa_of_closest.append(d_gtdb_tax[cur_gid].split(';')[rank_idx])
        if len(taxa_of_closest) >= n:
            break

    counts_of_closest = Counter(taxa_of_closest)
    return counts_of_closest


def process_all_markers(n, gid_of_interest, d_marker_to_path, path_out, ref_dm, ref_labels, d_gtdb_tax, target_taxon):
    d_marker_to_res = dict()
    ref_result = process_marker(gid_of_interest, ref_dm, ref_labels, n, d_gtdb_tax, target_taxon)
    d_marker_to_res['ref'] = ref_result
    all_taxa = set()
    all_taxa.update(ref_result.keys())

    for marker, path_marker_tree in d_marker_to_path.items():
        tree = dendropy.Tree.get(path=path_marker_tree, schema='newick', preserve_underscores=True)
        dm, labels = get_pdm(tree)
        result = process_marker(gid_of_interest, dm, labels, n, d_gtdb_tax, target_taxon)
        d_marker_to_res[marker] = result
        all_taxa.update(result.keys())

    vectors = dict()
    for marker, result in d_marker_to_res.items():
        vector = list()
        for taxon in all_taxa:
            vector.append(result.get(taxon, 0))
        vectors[marker] = vector

    heatmap_labels = sorted(vectors.keys())
    heatmap_arr = np.zeros((len(heatmap_labels), len(heatmap_labels)))
    print('Processing each marker')
    for i, marker_i in tqdm(enumerate(heatmap_labels), total=len(heatmap_labels)):
        for j, marker_j in enumerate(heatmap_labels):
            vector_i = vectors[marker_i]
            vector_j = vectors[marker_j]
            wj = weighted_jaccard(vector_i, vector_j) * 100
            heatmap_arr[i, j] = wj

    # fig, ax = plt.subplots(figsize=(10, 10))
    d_taxon_prefix_to_rank = {'d': 'domain', 'p': 'phylum', 'c': 'class', 'o': 'order', 'f': 'family', 'g': 'genus',
                              's': 'species'}
    target_rank = d_taxon_prefix_to_rank[target_taxon]
    g = sns.clustermap(heatmap_arr, xticklabels=heatmap_labels, yticklabels=heatmap_labels, annot=True, fmt=".0f")
    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle(
        f'Weighted jaccard similarity index of {gid_of_interest} to the {n:,} closest neighbours \n(by patristic distance) classified at the {target_rank} level'
    )
    g.savefig(path_out, dpi=300)

    return


def main():
    path_ref_tree = '/Users/aaron/obsidian/phd/PhD/Lab Stays/Aalborg/mfd/data/gtdbtk.ar53.decorated.bootstrap.tree'
    path_marker = '/Users/aaron/obsidian/phd/PhD/Lab Stays/Aalborg/mfd/data/ar53_marker'
    path_gtdb_tax = '/Users/aaron/obsidian/phd/PhD/Lab Stays/Aalborg/mfd/data/ar53_taxonomy_r214.tsv'
    path_out = '/tmp/clustermap.png'
    gid_of_interest = 'MFD04408.bin.c.2'
    n_neighbours = 50
    target_taxon = 'c'
    domain = 'ar53'

    d_gtdb_tax = read_tax_file(path_gtdb_tax)
    d_marker_to_path = get_marker_paths(path_marker, domain)

    ref_tree = dendropy.Tree.get(path=path_ref_tree, schema='newick', preserve_underscores=True)
    ref_dm, ref_labels = get_pdm(ref_tree)

    process_all_markers(
        n_neighbours,
        gid_of_interest,
        d_marker_to_path,
        path_out,
        ref_dm,
        ref_labels,
        d_gtdb_tax,
        target_taxon
    )

    return


if __name__ == '__main__':
    main()
