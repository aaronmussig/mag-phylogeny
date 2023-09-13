import json
import re


def get_from_group(name, line):
    if name == 'data05':
        hit = re.findall(r'"data05": \[\[([\d.]+), ([\d.]+), ([\d.]+)\]', line)
    else:
        hit = re.findall(r'"' + name + '": \[\[([\d.]+), ([\d.]+), ([\d.]+), ([\d.]+)\], \[.+?\]\]', line)
    print(hit)
    hit = hit[0]
    return float(hit[0]), float(hit[2]), float(hit[3])


def get_taxon_red_values(path_red_taxa):
    out = dict()
    out_d_rank_to_qty = dict()
    with open(path_red_taxa) as f:
        for line in f.readlines():
            if line.strip().startswith('mpld3.draw_figure('):
                content = line.strip()
                break

    json_str = re.findall(r'({.+})', content)[0]
    json_obj = json.loads(json_str)

    # Find the content of interest
    group_to_rank = {'data01': 'p', 'data02': 'c', 'data03': 'o', 'data04': 'f', 'data05': 'g'}
    for group, rank in group_to_rank.items():
        cur_json = json_obj['data'][group][0]
        median = cur_json[0]
        minimum = cur_json[2]
        if group == 'data05':
            maximum = 1.0
        else:
            maximum = cur_json[3]
        out[rank] = (minimum, median, maximum)

    re_rank_count = re.compile(r'(\w+) \(([\d,]+)\)')
    for rank_count in json_obj['axes'][0]['axes'][1]['tickformat']:
        hit = re_rank_count.match(rank_count)
        out_d_rank_to_qty[hit.group(1)] = int(hit.group(2).replace(',', ''))
    return out, out_d_rank_to_qty
