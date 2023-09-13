import sys
from collections import defaultdict
from pathlib import Path

import dendropy
from tqdm import tqdm

sys.path.insert(0, Path(__file__).absolute().parent.parent.parent.parent.as_posix())

from magphylogeny.methods.phylorank_red_html import get_taxon_red_values
from magphylogeny.util.red import node_to_red
from magphylogeny.util.tree import calculate_node_depths, d_node_to_leaves, d_node_to_desc_labels, get_parent_taxon, \
    RANKS


def find_candidate_placement_nodes(parent_node, start_node, target_taxon, iqr_min, d_node_desc_labels, d_node_to_descs,
                                   d_node_to_red):
    current_node = start_node
    previous_node = None
    candidate_node = None
    while current_node is not None or current_node == parent_node:

        taxon_descs = d_node_desc_labels[current_node]
        gids_descs = d_node_to_descs[current_node]
        node_red = d_node_to_red[current_node]

        # There is nothing left for us to choose from
        if node_red < iqr_min:
            break

        # print(gids_descs)
        # print(taxon_descs)
        # print(node_red)
        # print()

        # No candidate nodes are allowed outside of the IRQ range
        within_red_range = iqr_min <= node_red
        has_taxon_conflict = any(x for x in taxon_descs if f'{target_taxon}__' in x)
        # has_gid_conflict = any(x for x in gids_descs if not x.startswith('MFD'))

        """
        Taxon conflict is not so bad, if there is only one taxon in the set. If
        that is the case, then we could move that node higher.
        EDIT: Changing this as it's too complicated otherwise
        """
        # n_desc_taxa_at_rank = 0
        # for cur_taxon in taxon_descs:
        #     if f'{target_taxon}__' in cur_taxon:
        #         n_desc_taxa_at_rank += 1
        # has_taxon_conflict = n_desc_taxa_at_rank > 1

        if has_taxon_conflict:
            # print('Taxon conflict')
            break

        # Check if we would cause any taxon conflicts
        # if has_gid_conflict or has_taxon_conflict:
        #     candidate_nodes[depth].add(current_node)
        #     candidate_nodes[depth].add(previous_node)
        #     break

        if within_red_range:
            candidate_node = current_node

        # Keep searching
        previous_node = current_node
        current_node = current_node.parent_node
    return candidate_node


def get_proposed_rank_for_gid(gid, d_leaf_taxon_to_node, d_novelty_to_range, d_node_to_descs, d_node_desc_labels,
                              d_node_to_red):
    node = d_leaf_taxon_to_node[gid]

    # Get the parent taxon
    parent_node, parent_taxon = get_parent_taxon(node)

    target_taxon = RANKS[RANKS.index(parent_taxon[0]) + 1]

    if target_taxon == 's':
        return None, None

    # Get find all named nodes that contain the rank below the parent taxon
    # e.g. if p__A, find all p__A; c__X nodes
    # contained_red = get_contained_taxa(parent_node, target_taxon)

    # q1, q3 = np.percentile(contained_red, [25, 75])
    # iqr = (q3 - q1) * 1.5
    # iqr_min = q1 - iqr
    # iqr_max = q3 + iqr
    # min_red = min(contained_red) if contained_red else 0
    # max_red = max(contained_red) if contained_red else 1
    min_red = d_novelty_to_range[target_taxon][0]
    max_red = d_novelty_to_range[target_taxon][2]

    """
    Find all eligble nodes in the tree that could host a label for this.
    Note: The node may still be unique within it (e.g. novel class), but there is no
    place to put the node. This would be if the RED is off.
    """
    # parent_node, start_node, target_taxon, iqr_min, d_node_desc_labels, d_node_to_descs,
    #                                    d_node_to_red
    candidate_node = find_candidate_placement_nodes(parent_node, node, target_taxon, min_red, d_node_desc_labels,
                                                    d_node_to_descs, d_node_to_red)
    candidate_desc_gids = d_node_to_descs[candidate_node]
    # print('----')
    # print(candidate_node)
    # print(candidate_desc_gids)
    # print('---')
    has_gtdb_genomes = any(x for x in candidate_desc_gids if not x.startswith('MFD'))
    novelty = None

    # Easy case, if there was no candidate node within range, this is simply novel
    if candidate_node is None:
        print('????', gid)
    else:
        if has_gtdb_genomes:
            novelty = f'mix_{target_taxon}'
        else:
            novelty = f'mfd_{target_taxon}'

    return candidate_node, novelty


# get_proposed_rank_for_gid('MFD01689.bin.1.29')
# get_proposed_rank_for_gid('MFD01476.bin.2.32')
# get_proposed_rank_for_gid('MFD03726.bin.1.205')  # novel order
# get_proposed_rank_for_gid('MFD04408.bin.2.191')  # novel family
# get_proposed_rank_for_gid('MFD04408.bin.2.183')  # this sucks, should be higher
# get_proposed_rank_for_gid('MFD02416.bin.c.18')
# get_proposed_rank_for_gid('')

def get_novelty(tree, d_leaf_taxon_to_node, d_novelty_to_range, d_node_to_descs, d_node_desc_labels, d_node_to_red):
    out = dict()
    non_novel_gids = set()
    d_candidate_node_to_gids = defaultdict(list)
    for leaf_node in tqdm(tree.leaf_nodes()):
        gid = leaf_node.taxon.label
        if gid.startswith('MFD'):
            candidate_node, novelty = get_proposed_rank_for_gid(gid, d_leaf_taxon_to_node, d_novelty_to_range,
                                                                d_node_to_descs, d_node_desc_labels, d_node_to_red)

            if candidate_node is None and novelty is None:
                non_novel_gids.add(gid)
            else:
                # print(candidate_node, novelty, gid)
                out[gid] = (candidate_node, novelty)
                d_candidate_node_to_gids[candidate_node].append((gid, novelty))
    return out, non_novel_gids, d_candidate_node_to_gids


from collections import Counter


def report_on_results(d_candidate_nodes_to_gids, non_novel_gids, d_gid_to_novelty, path_out):
    rows = list()
    gid_to_cluster = dict()
    d_rank_to_novel_gids = defaultdict(set)
    d_cluster_id_to_gids = defaultdict(set)

    # Index the clusters
    gids_seen = set()
    for i, lst_gid_novelty in enumerate(d_candidate_nodes_to_gids.values()):
        for gid, novelty in lst_gid_novelty:
            assert (gid not in gids_seen)
            grouping, rank = novelty.split('_')
            is_novel = grouping == 'mfd'
            gid_to_cluster[gid] = i
            gids_seen.add(gid)

            if is_novel:
                d_rank_to_novel_gids[rank].add(i)
                d_cluster_id_to_gids[i].add(gid)

    for rank, novel_set in sorted(d_rank_to_novel_gids.items(), key=lambda x: RANKS.index(x[0])):
        print(f'Novel {rank} contains {len(novel_set):,} clades')
        gids_in_each_clade = sorted([len(d_cluster_id_to_gids[x]) for x in novel_set])
        counts = Counter(gids_in_each_clade)
        for k, v in counts.items():
            print(f'{k}(x{v})') if v > 1 else print(f'{k}')
        print()
        # gids_in_each_clade = ', '.join(map(str, gids_in_each_clade))
        # print(f'Clades contain {gids_in_each_clade}')

    for gid in non_novel_gids:
        rows.append(('s', False, 'N/A', gid))

    for gid, (candidate_node, novelty) in d_gid_to_novelty.items():
        grouping, rank = novelty.split('_')
        is_novel = grouping == 'mfd'
        rows.append((rank, is_novel, gid_to_cluster[gid], gid))

    with open(path_out, 'w') as f:
        rows = sorted(rows, key=lambda x: (RANKS.index(x[0]), x[2]))
        rows.insert(0, ('rank', 'is_novel', 'cluster_id', 'gids'))
        for row in rows:
            f.write('\t'.join(map(str, row)) + '\n')

    return


# print(D_CANDIDATE_NODE_TO_GIDS)


def main():
    # path_tree = '/Users/aaron/obsidian/phd/PhD/Lab Stays/Aalborg/data/phylorank_denovo_arc_trusted/gtdbtk.ar53.decorated.scaled.tree'
    # path_red = '/Users/aaron/obsidian/phd/PhD/Lab Stays/Aalborg/data/phylorank_denovo_arc_trusted/gtdbtk.ar53.decorated.node_rd.tsv'
    # path_red_taxa = '/Users/aaron/obsidian/phd/PhD/Lab Stays/Aalborg/data/phylorank_denovo_arc_trusted/gtdbtk.ar53.decorated.html'
    # path_out = '/tmp/output.png'

    # path_tree = '/tmp/mfd_tmp/gtdbtk.bac120.rooted.decorated.bootstrap.scaled.tree'
    # path_red = '/tmp/mfd_tmp/gtdbtk.bac120.rooted.decorated.bootstrap.node_rd.tsv'
    # path_red_taxa = '/tmp/mfd_tmp/gtdbtk.bac120.rooted.decorated.bootstrap.html'
    # path_out = '/tmp/output.png'
    #
    path_tree = sys.argv[1]
    path_red = sys.argv[2]
    path_red_taxa = sys.argv[3]
    path_out = sys.argv[4]

    tree = dendropy.Tree.get(path=path_tree, schema='newick', preserve_underscores=True)

    d_depth_to_node, d_node_to_depth = calculate_node_depths(tree)
    d_node_to_descs, d_leaf_taxon_to_node = d_node_to_leaves(d_depth_to_node)
    node_to_desc_labels = d_node_to_desc_labels(tree)
    d_node_to_red = node_to_red(d_node_to_descs, d_node_to_depth, d_leaf_taxon_to_node, path_red)
    d_novelty_to_range, _ = get_taxon_red_values(path_red_taxa)

    d_gid_to_novelty, non_novel_gids, d_candidate_nodes_to_gids = get_novelty(
        tree,
        d_leaf_taxon_to_node,
        d_novelty_to_range,
        d_node_to_descs,
        node_to_desc_labels,
        d_node_to_red
    )

    report_on_results(
        d_candidate_nodes_to_gids,
        non_novel_gids,
        d_gid_to_novelty,
        path_out
    )

    return


if __name__ == '__main__':
    main()
