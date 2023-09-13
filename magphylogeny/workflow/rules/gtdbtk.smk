import sys
from pathlib import Path

sys.path.insert(0,Path(workflow.basedir).parent.parent.as_posix())

from magphylogeny.model.gtdbtk.identify import GTDBTkIdentify
from magphylogeny.model.gtdbtk.align import GTDBTkAlign
from magphylogeny.model.gtdbtk.infer import GTDBTkInfer
from magphylogeny.model.gtdbtk.root import GTDBTkRoot
from magphylogeny.model.gtdbtk.decorate import GTDBTkDecorate
from magphylogeny.model.gtdbtk.bootstrap import GTDBTkBootstrapMsaReplicate, GTDBTkBootstrapTreeReplicate, \
    GTDBTkBootstrapMergeTrees
from magphylogeny.model.gtdbtk.phylorank import GTDBTkPhylorankTrustedTaxa, GTDBTkPhylorankOutliers

GTDBTK_IDENTIFY = GTDBTkIdentify(config)
GTDBTK_ALIGN = GTDBTkAlign(config)
GTDBTK_INFER = GTDBTkInfer(config)
GTDBTK_ROOT = GTDBTkRoot(config)
GTDBTK_DECORATE = GTDBTkDecorate(config)
GTDBTK_BS_MSA_REP = GTDBTkBootstrapMsaReplicate(config)
GTDBTK_BS_TREE_REP = GTDBTkBootstrapTreeReplicate(config)
GTDBTK_BS_MERGE = GTDBTkBootstrapMergeTrees(config)
GTDBTK_PHYLORANK_TRUSTED_TAXA = GTDBTkPhylorankTrustedTaxa(config)
GTDBTK_PHYLORANK_OUTLIERS = GTDBTkPhylorankOutliers(config)

rule gtdbtk_identify:
    conda: "gtdbtk-2.3.0"
    threads: config["gtdbtk"]["identify"]["cpus"]
    resources:
        queue='slow'
    benchmark:
        "benchmarks/gtdbtk/identify.tsv"
    output:
        root=directory(GTDBTK_IDENTIFY.dir),
        ar53=GTDBTK_IDENTIFY.markers_summary('ar53'),
        bac120=GTDBTK_IDENTIFY.markers_summary('bac120')
    shell:
        "gtdbtk identify --genome_dir {config[genome][dir]} --extension {config[genome][ext]} --cpus {threads} --out_dir={output.root}"

rule gtdbtk_align:
    conda: "gtdbtk-2.3.0"
    threads: config["gtdbtk"]["align"]["cpus"]
    resources:
        queue='slow'
    input:
        identify_dir=GTDBTK_IDENTIFY.dir
    output:
        root=directory(GTDBTK_ALIGN.dir),
        ar53_msa=GTDBTK_ALIGN.msa('ar53'),
        bac120_msa=GTDBTK_ALIGN.msa('bac120'),
        bac120_markers=GTDBTK_ALIGN.markers('bac120'),
        ar53_markers=GTDBTK_ALIGN.markers('ar53'),
    shell:
        "gtdbtk align --identify_dir {input.identify_dir} --out_dir {output.root} --cpus {config[gtdbtk][align][cpus]} --debug"

rule gtdbtk_infer:
    conda: "gtdbtk-2.3.0"
    threads: config["gtdbtk"]["infer"]["cpus"]
    resources:
        queue='slow'
    input:
        msa=GTDBTK_ALIGN.msa('{domain}')
    output:
        root=directory(GTDBTK_INFER.dir('{domain}')),
        tree=GTDBTK_INFER.tree('{domain}')
    shell:
        "gtdbtk infer --msa_file {input.msa} --out_dir {output.root} --cpus {config[gtdbtk][infer][cpus]}"

rule gtdbtk_root:
    conda: "gtdbtk-2.3.0"
    threads: 1
    resources:
        queue='quick'
    input:
        tree=GTDBTK_INFER.tree('{domain}')
    output:
        tree=GTDBTK_ROOT.tree('{domain}')
    params:
        outgroup=lambda wildcards: GTDBTK_ROOT.outgroup(wildcards.domain)
    shell:
        "gtdbtk root --input_tree {input.tree} --outgroup_taxon {params.outgroup} --output_tree {output.tree}"

rule gtdbtk_decorate:
    conda: "gtdbtk-2.3.0"
    threads: 1
    resources:
        queue='quick'
    input:
        tree=GTDBTK_ROOT.tree('{domain}')
    output:
        tree=GTDBTK_DECORATE.tree('{domain}'),
        tree_table=GTDBTK_DECORATE.tree_table('{domain}'),
        taxonomy=GTDBTK_DECORATE.taxonomy('{domain}')
    shell:
        "gtdbtk decorate --input_tree {input.tree} --output_tree {output.tree}"

rule bootstrap_create_msa:
    threads: 1
    resources:
        queue='quick'
    input:
        msa=GTDBTK_ALIGN.msa('{domain}')
    output:
        faa=GTDBTK_BS_MSA_REP.faa('{domain}','{replicate}')
    shell:
        """
        python "{workflow.basedir}/scripts/bootstrap_msa.py" "{input.msa}" "{output.faa}"
        """

rule bootstrap_create_tree:
    conda: "fasttree-2.1.11"
    threads: config["marker_genes"]["fasttree"]["cpus"]
    resources:
        queue='memory'
    input:
        faa=GTDBTK_BS_MSA_REP.faa('{domain}','{replicate}')
    output:
        log=GTDBTK_BS_TREE_REP.log('{domain}','{replicate}'),
        tree=GTDBTK_BS_TREE_REP.tree('{domain}','{replicate}'),
    shell:
        """
        OMP_NUM_THREADS={threads} FastTreeMP -wag -log {output.log} {input.faa} > {output.tree}
        """

rule bootstrap_merge_trees:
    conda: "genometreetk-0.1.8"
    threads: 1
    resources:
        queue='quick'
    input:
        domain_tree=GTDBTK_DECORATE.tree('{domain}'),
        replicate_trees=GTDBTK_BS_TREE_REP.trees('{{domain}}')
    output:
        root=directory(GTDBTK_BS_MERGE.dir('{domain}')),
        tree=GTDBTK_BS_MERGE.tree('{domain}')
    params:
        boot_dir=lambda wildcards: GTDBTK_BS_TREE_REP.dir(wildcards.domain),
        n_replicates=config["gtdbtk"]["bootstrap"]["n_boot"]
    shell:
        """
        genometreetk bootstrap --boot_dir {params.boot_dir} -r {params.n_replicates} {input.domain_tree} NONE {output.root}
        """

rule create_phylorank_trusted_taxa:
    threads: 1
    resources:
        queue='quick'
    input:
        tax=GTDBTK_DECORATE.taxonomy('{domain}')
    output:
        txt=GTDBTK_PHYLORANK_TRUSTED_TAXA.trusted_taxa('{domain}')
    shell:
        """
        python "{workflow.basedir}/scripts/phylorank_trusted_taxa.py" "{input.tax}" "{output.txt}"
        """

rule phylorank_scale_bootstrap:
    conda: "phylorank-0.1.12"
    resources:
        queue='quick'
    threads: 1
    input:
        tree=GTDBTK_BS_MERGE.tree('{domain}'),
        trusted_taxa=GTDBTK_PHYLORANK_TRUSTED_TAXA.trusted_taxa('{domain}'),
        tax=GTDBTK_DECORATE.taxonomy('{domain}')
    output:
        root=directory(GTDBTK_PHYLORANK_OUTLIERS.dir('{domain}')),
        node_rd=GTDBTK_PHYLORANK_OUTLIERS.node_rd('{domain}'),
        red_scaled_tree=GTDBTK_PHYLORANK_OUTLIERS.red_scaled_tree('{domain}'),
        red_scaled_decorated_tree=GTDBTK_PHYLORANK_OUTLIERS.red_scaled_decorated_tree('{domain}'),
        red_html=GTDBTK_PHYLORANK_OUTLIERS.red_html_plot('{domain}'),
    shell:
        """
        phylorank outliers -t {input.trusted_taxa} {input.tree} {input.tax} {output.root}
        """
