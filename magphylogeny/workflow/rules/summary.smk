import sys
from pathlib import Path


sys.path.insert(0,Path(workflow.basedir).parent.parent.as_posix())

from magphylogeny.model.gtdbtk.identify import GTDBTkIdentify
from magphylogeny.model.gtdbtk.decorate import GTDBTkDecorate
from magphylogeny.model.plots.marker_summary import PlotMarkerSummary
from magphylogeny.model.gtdbtk.phylorank import GTDBTkPhylorankOutliers
from magphylogeny.model.summary.novelty_of_genomes import SummaryNoveltyOfGenomes


GTDBTK_PHYLORANK_OUTLIERS = GTDBTkPhylorankOutliers(config)
GTDBTK_IDENTIFY = GTDBTkIdentify(config)
GTDBTK_DECORATE = GTDBTkDecorate(config)
PLOT_MARKER_SUMMARY = PlotMarkerSummary(config)
SUMMARY_NOVELTY_GENOMES = SummaryNoveltyOfGenomes(config)

rule summary_novelty_of_genomes:
    input:
        tree=GTDBTK_PHYLORANK_OUTLIERS.red_scaled_tree('{domain}'),
        red=GTDBTK_PHYLORANK_OUTLIERS.node_rd('{domain}'),
        red_html=GTDBTK_PHYLORANK_OUTLIERS.red_html_plot('{domain}')
    resources:
        queue='quick'
    output:
        tsv=SUMMARY_NOVELTY_GENOMES.path('{domain}'),
    shell:
        """
        python "{workflow.basedir}/scripts/summary_novelty_of_genomes.py" "{input.tree}" "{input.red}" "{input.red_html}" "{output.tsv}"
        """
