from collections import Counter

from tqdm import tqdm

from magphylogeny.util.tree import mrca_faster


def get_prop_lineages_at_threshold(thresh_max, d_node_to_red, d_node_to_descs):
    out = list()

    # 1. Filter to a subset of nodes that satifsy the RED criteria
    candidate_nodes = set()
    for node, red in d_node_to_red.items():
        if red >= thresh_max:
            candidate_nodes.add(node)

    # Filter through to make sure these are not duplicated
    for node in candidate_nodes:
        if node.parent_node not in candidate_nodes:
            out.append(d_node_to_descs[node])

    counts = Counter(out)
    percentages = dict()
    total = sum(counts.values())
    for key, value in counts.items():
        percentages[key] = (value / total) * 100
    counts_out = {
        'MIX': percentages.get('MIX', 0),
        'MFD': percentages.get('MFD', 0),
        'GTDB': percentages.get('GTDB', 0),
        'MIX_count': counts.get('MIX', 0),
        'MFD_count': counts.get('MFD', 0),
        'GTDB_count': counts.get('GTDB', 0),
        'TOTAL': total
    }
    return counts_out


def node_to_red(d_node_to_descs, d_node_to_depth, d_leaf_taxon_to_node, path_red):
    out = dict()
    print('Calculating RED of each node')
    with open(path_red) as f:
        for line in tqdm(f.readlines()):
            mrca, red = line.strip().split('\t')
            mrca_split = mrca.split('|')

            if len(mrca_split) > 1:
                # Find the node in the tree using MRCA
                # node = TREE.mrca(taxon_labels=mrca.split('|'))
                node = mrca_faster(mrca_split, d_node_to_descs, d_node_to_depth)
                # assert(node==node2)
            else:
                # node = TREE.find_node_with_taxon_label(mrca)
                node = d_leaf_taxon_to_node[mrca]
            out[node] = float(red)
    return out
