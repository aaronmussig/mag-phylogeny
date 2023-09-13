from collections import defaultdict, deque

from tqdm import tqdm

from magphylogeny.util.values import is_float

RANKS = ('d', 'p', 'c', 'o', 'f', 'g', 's')


def parse_label(label):
    """Parse a Newick label which may contain a support value, taxon, and/or auxiliary information.

    Parameters
    ----------
    label : str
        Internal label in a Newick tree.

    Returns
    -------
    float
        Support value specified by label, or None
    str
        Taxon specified by label, or None
    str
        Auxiliary information, on None
    """

    support = None
    taxon = None
    auxiliary_info = None

    if label:
        label = label.strip()
        if '|' in label:
            label, auxiliary_info = label.split('|')

        if ':' in label:
            support, taxon = label.split(':')
            support = float(support)
        else:
            if is_float(label):
                support = float(label)
            elif label != '':
                taxon = label

    return support, taxon, auxiliary_info


def calculate_node_depths(tree):
    out = defaultdict(list)
    out_rev = dict()
    queue = deque([(tree.seed_node, 0)])
    while queue:
        node, cur_depth = queue.popleft()
        for child in node.child_node_iter():
            queue.append((child, cur_depth + 1))
        out[cur_depth].append(node)
        out_rev[node] = cur_depth
    return out, out_rev


def d_node_to_leaves(d_depth_to_node):
    out = defaultdict(set)
    out2 = dict()
    for cur_depth in sorted(d_depth_to_node.keys(), reverse=True):
        depth_nodes = d_depth_to_node[cur_depth]
        for cur_node in depth_nodes:
            if cur_node.is_leaf():
                out[cur_node] = {cur_node.taxon.label}
                out2[cur_node.taxon.label] = cur_node
            for child in cur_node.child_node_iter():
                out[cur_node].update(out[child])
    return out, out2


def mrca_faster(mrca_split, d_node_to_descs, d_node_to_depth):
    # 47 seconds time to beat
    mrca_set = set(mrca_split)

    mrca_node = None
    mrca_depth = None
    for node, node_descs in d_node_to_descs.items():
        if mrca_set.issubset(node_descs):
            cur_depth = d_node_to_depth[node]
            if mrca_node is None:
                mrca_node = node
                mrca_depth = cur_depth
            else:
                if mrca_depth < cur_depth:
                    mrca_node = node
                    mrca_depth = cur_depth

    # Get the deepest node
    return mrca_node


def calc_node_to_descendants(d_node_to_red):
    out = dict()
    print('Calculating descendants of each node')
    for node in tqdm(d_node_to_red):
        leaf_nodes = set()
        for leaf_node in node.leaf_nodes():
            label = leaf_node.taxon.label
            if label.startswith('MFD'):
                leaf_nodes.add('MFD')
            else:
                leaf_nodes.add('GTDB')
        if len(leaf_nodes) == 1:
            out[node] = leaf_nodes.pop()
        else:
            out[node] = 'MIX'
    return out


def d_node_to_desc_labels(tree):
    out = defaultdict(set)
    for leaf_node in tqdm(tree.leaf_nodes()):
        current_node = leaf_node
        desc_taxa = set()
        while current_node is not None:
            _, taxon, _ = parse_label(current_node.label)
            if taxon:
                desc_taxa.add(taxon)
            out[current_node].update(desc_taxa)
            current_node = current_node.parent_node
    return out


def get_parent_taxon(start_node):
    # Go up the tree until we find the highest rank we're contained in
    current_node = start_node
    parent_taxon = None
    while current_node is not None or parent_taxon is None:
        _, taxon, _ = parse_label(current_node.label)

        # Found a taxon label
        if taxon is not None:
            parent_taxon = taxon.split(';')[-1].strip()
            break

        # Keep going up
        current_node = current_node.parent_node

    return current_node, parent_taxon


def get_contained_taxa(start_node, target_rank, d_node_to_red):
    assert (not target_rank.startswith('g__'))

    queue = deque([start_node])
    d_found_node_to_taxonomy = dict()
    while len(queue) > 0:
        node = queue.popleft()

        _, taxon, _ = parse_label(node.label)

        # Found a taxon label, check it
        if taxon and taxon.split(';')[-1].strip().startswith(f'{target_rank}'):
            d_found_node_to_taxonomy[node] = taxon

        # Otherwise, keep going down
        else:
            for child in node.child_node_iter():
                queue.append(child)

    all_red = list()
    for node, taxonomy in d_found_node_to_taxonomy.items():
        current_red = d_node_to_red[node]
        all_red.append(current_red)
    return all_red
