import sys
from pathlib import Path

sys.path.insert(0,Path(workflow.basedir).parent.parent.as_posix())

from magphylogeny.model.marker_genes.rep_markers import MarkerGenesRepDownload
from magphylogeny.model.gtdbtk.align import GTDBTkAlign
from magphylogeny.model.marker_genes.concat_marker import MarkerGenesConcatMarker, MarkerGenesConcatMarkerTrimmed, \
    MarkerGenesConcatGeneTree

MARKER_REP_DL = MarkerGenesRepDownload(config)
GTDBTK_ALIGN = GTDBTkAlign(config)
CONCAT_MARKER = MarkerGenesConcatMarker(config)
CONCAT_MARKER_TRIMMED = MarkerGenesConcatMarkerTrimmed(config)
CONCAT_MARKER_TREE = MarkerGenesConcatGeneTree(config)

rule download_rep_markers_ar53:
    threads: 1
    resources:
        mem_mb=1,
        runtime=1,
        queue='quick'
    input:
        archive=MARKER_REP_DL.archive('ar53')
    output:
        root=directory(MARKER_REP_DL.dir('ar53')),
        markers=MARKER_REP_DL.markers('ar53')
    benchmark:
        "benchmarks/marker_genes/download_rep_markers_ar53.tsv"
    shell:
        "tar -xzf {input.archive} -C {output.root}"

rule download_rep_markers_bac120:
    threads: 1
    resources:
        mem_mb=1,
        runtime=1,
        queue='quick'
    input:
        archive=MARKER_REP_DL.archive('bac120')
    output:
        root=directory(MARKER_REP_DL.dir('bac120')),
        markers=MARKER_REP_DL.markers('bac120')
    benchmark:
        "benchmarks/marker_genes/download_rep_markers_bac120.tsv"
    shell:
        "tar -xzf {input.archive} -C {output.root}"

rule concat_marker:
    threads: 1
    resources:
        mem_mb=1,
        runtime=1,
        queue='quick'
    input:
        rep=MARKER_REP_DL.marker('{domain}','{marker}'),
        mag=GTDBTK_ALIGN.marker('{domain}','{marker}')
    output:
        faa=CONCAT_MARKER.marker('{domain}','{marker}')
    benchmark:
        "benchmarks/marker_genes/concat_marker_{domain}_{marker}.tsv"
    shell:
        "cat {input.mag} {input.rep} > {output.faa}"

rule remove_blank_sequences:
    threads: 1
    resources:
        mem_mb=10,
        runtime=1,
        queue='quick'
    input:
        faa=CONCAT_MARKER.marker('{domain}','{marker}')
    output:
        faa=CONCAT_MARKER_TRIMMED.marker('{domain}','{marker}')
    benchmark:
        "benchmarks/marker_genes/remove_blank_sequences_{domain}_{marker}.tsv"
    shell:
        """
        python "{workflow.basedir}/scripts/trim_msa.py" "{input.faa}" "{output.faa}"
        """

rule create_marker_gene_tree:
    conda: "fasttree-2.1.11"
    threads: config["marker_genes"]["fasttree"]["cpus"]
    resources:
        mem_mb=lambda wildcards: CONCAT_MARKER_TREE.mem_mb(wildcards.domain),
        runtime=lambda wildcards: CONCAT_MARKER_TREE.runtime(wildcards.domain),
        queue='memory'
    input:
        faa=CONCAT_MARKER_TRIMMED.marker('{domain}','{marker}')
    output:
        log=CONCAT_MARKER_TREE.log('{domain}','{marker}'),
        tree=CONCAT_MARKER_TREE.tree('{domain}','{marker}'),
    benchmark:
        "benchmarks/marker_genes/create_marker_gene_tree_{domain}_{marker}.tsv"
    shell:
        """
        OMP_NUM_THREADS={threads} FastTreeMP -wag -log {output.log} {input.faa} > {output.tree}
        """
